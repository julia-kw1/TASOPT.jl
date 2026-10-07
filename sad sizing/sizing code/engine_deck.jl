module GenerateE175CFM56EngineDeck

using TASOPT
using CSV
using Printf

using ..WorkflowPaths

include(__TASOPTindices__)

const ENGINE_MODEL = "CFM56-7B27"
const RATED_THRUST_LBF = 27_300.0
const OUTPUT_FILE = WorkflowPaths.output_path("e175_cfm56_7b27_engine_deck.csv")

# Deck limits stay fixed while the runner controls grid spacing.
const MAX_ALTITUDE_FT = 45_000.0
const MAX_MACH = 0.85
const DEFAULT_ALTITUDE_STEP_FT = 5_000.0
const DEFAULT_MACH_STEP = 0.10
const DEFAULT_THROTTLE_STEP = 0.05

# The lowest converged state above this Tt4 approximates idle.
const IDLE_PROXY_TT4_K = 900.0
const TT4_SEARCH_STEP_K = 25.0
const CORRECTED_SPEED_TOLERANCE = 1.0e-5
const ROOT_ITERATIONS = 16

const IMISSION = 1
const ENGINE_POINT = iptakeoff # tfwrap! uses fixed-Tt4 off-design mode here.
const KG_TO_LB = 2.2046226218487757

function grid_to_limit(limit::Float64, step::Float64)
    step > 0.0 || error("Grid spacing must be positive.")
    values = collect(0.0:step:limit)
    if values[end] < limit - 10eps(limit)
        push!(values, limit)
    end
    return values
end

function quiet(f::Function)
    redirect_stdout(devnull) do
        redirect_stderr(devnull) do
            f()
        end
    end
end

function load_base_aircraft()
    ac = quickload_aircraft(WorkflowPaths.require_model(:hybrid))
    ac.parm[imDeltaTatm, IMISSION] = 0.0

    # Engine decks use clear-flow conditions without aircraft offtakes.
    ac.parg[igfBLIw] = 0.0
    ac.parg[igfBLIf] = 0.0
    ac.parg[igmofWpay] = 0.0
    ac.parg[igmofWMTO] = 0.0
    ac.parg[igPofWpay] = 0.0
    ac.parg[igPofWMTO] = 0.0
    return ac
end

function condition_aircraft(base, altitude_ft::Float64, mach::Float64)
    ac = deepcopy(base)
    ac.para[iaalt, ENGINE_POINT, IMISSION] = altitude_ft * ft_to_m
    ac.para[iaMach, ENGINE_POINT, IMISSION] = mach
    TASOPT.set_ambient_conditions!(ac, ENGINE_POINT, mach; im = IMISSION)
    ac.pare[ieConvFail, ENGINE_POINT, IMISSION] = 0.0
    return ac
end

function evaluate_turbine_temperature(condition, tt4_K::Float64)
    ac = deepcopy(condition)
    ac.pare[ieTt4, ENGINE_POINT, IMISSION] = tt4_K
    ac.pare[ieConvFail, ENGINE_POINT, IMISSION] = 0.0

    try
        quiet() do
            ac.engine.model.enginecalc!(
                ac, "off_design", IMISSION, ENGINE_POINT, true
            )
        end
    catch exception
        exception isa DomainError || rethrow()
        return nothing
    end

    outputs = (
        ac.pare[ieFe, ENGINE_POINT, IMISSION],
        ac.pare[ieNbf, ENGINE_POINT, IMISSION],
        ac.pare[ieTSFC, ENGINE_POINT, IMISSION],
    )
    converged = ac.pare[ieConvFail, ENGINE_POINT, IMISSION] == 0.0
    return converged && all(isfinite, outputs) ? ac : nothing
end

corrected_fan_speed(ac) = ac.pare[ieNbf, ENGINE_POINT, IMISSION]

function first_converged_endpoint(condition, temperatures)
    for tt4_K in temperatures
        ac = evaluate_turbine_temperature(condition, tt4_K)
        if !isnothing(ac)
            return ac, tt4_K
        end
    end
    error("No converged engine point was found in the Tt4 range.")
end

function converged_endpoints(condition, maximum_tt4_K::Float64)
    low_ac, low_tt4 = first_converged_endpoint(
        condition,
        IDLE_PROXY_TT4_K:TT4_SEARCH_STEP_K:maximum_tt4_K,
    )
    high_ac, high_tt4 = first_converged_endpoint(
        condition,
        maximum_tt4_K:-TT4_SEARCH_STEP_K:IDLE_PROXY_TT4_K,
    )

    low_nbf = corrected_fan_speed(low_ac)
    high_nbf = corrected_fan_speed(high_ac)
    high_tt4 > low_tt4 || error("The converged Tt4 interval is empty.")
    high_nbf > low_nbf || error("Corrected fan speed is not monotonic with Tt4.")
    return (
        low_ac = low_ac,
        high_ac = high_ac,
        low_tt4 = low_tt4,
        high_tt4 = high_tt4,
        low_nbf = low_nbf,
        high_nbf = high_nbf,
    )
end

function solve_throttle(condition, endpoints, throttle::Float64)
    throttle <= 0.0 && return endpoints.low_ac, endpoints.low_tt4
    throttle >= 1.0 && return endpoints.high_ac, endpoints.high_tt4

    target = endpoints.low_nbf +
             throttle * (endpoints.high_nbf - endpoints.low_nbf)
    lo_t, hi_t = endpoints.low_tt4, endpoints.high_tt4
    lo_n, hi_n = endpoints.low_nbf, endpoints.high_nbf
    best_ac = nothing
    best_t = NaN
    best_error = Inf

    for _ in 1:ROOT_ITERATIONS
        fraction = clamp((target - lo_n) / (hi_n - lo_n), 0.05, 0.95)
        tt4 = lo_t + fraction * (hi_t - lo_t)
        trial = evaluate_turbine_temperature(condition, tt4)
        isnothing(trial) && error("The throttle solve reached an unconverged engine point.")

        nbf = corrected_fan_speed(trial)
        residual = abs(nbf - target)
        if residual < best_error
            best_ac, best_t, best_error = trial, tt4, residual
        end
        residual <= CORRECTED_SPEED_TOLERANCE && return trial, tt4

        if nbf < target
            lo_t, lo_n = tt4, nbf
        else
            hi_t, hi_n = tt4, nbf
        end
    end

    best_error <= 10.0 * CORRECTED_SPEED_TOLERANCE ||
        error("The throttle solve did not reach the corrected-speed tolerance.")
    return best_ac, best_t
end

function output_row(ac, altitude_ft, mach, throttle, endpoints, tt4)
    pare = ac.pare
    nbf = pare[ieNbf, ENGINE_POINT, IMISSION]
    achieved = (nbf - endpoints.low_nbf) /
               (endpoints.high_nbf - endpoints.low_nbf)
    neng = ac.parg[igneng]
    fuel_kg_s = pare[iemfuel, ENGINE_POINT, IMISSION] / neng
    fan_pr = pare[iepif, ENGINE_POINT, IMISSION]
    lpc_pr = pare[iepilc, ENGINE_POINT, IMISSION]
    hpc_pr = pare[iepihc, ENGINE_POINT, IMISSION]
    thrust_N = pare[ieFe, ENGINE_POINT, IMISSION]

    return (
        engine_model = ENGINE_MODEL,
        rated_thrust_lbf = RATED_THRUST_LBF,
        altitude_ft = altitude_ft,
        mach = mach,
        isa_deviation_K = 0.0,
        ambient_temperature_K = pare[ieT0, ENGINE_POINT, IMISSION],
        ambient_pressure_Pa = pare[iep0, ENGINE_POINT, IMISSION],
        ambient_density_kg_m3 = pare[ierho0, ENGINE_POINT, IMISSION],
        flight_speed_m_s = pare[ieu0, ENGINE_POINT, IMISSION],
        throttle_command = throttle,
        throttle_achieved = achieved,
        corrected_fan_speed = nbf,
        corrected_fan_speed_design_fraction =
            nbf / pare[ieNbfD, ENGINE_POINT, IMISSION],
        idle_corrected_fan_speed = endpoints.low_nbf,
        maximum_corrected_fan_speed = endpoints.high_nbf,
        idle_tt4_K = endpoints.low_tt4,
        maximum_tt4_K = endpoints.high_tt4,
        tt4_K = tt4,
        net_thrust_N_per_engine = thrust_N,
        net_thrust_lbf_per_engine = thrust_N / lbf_to_N,
        fuel_flow_kg_s_per_engine = fuel_kg_s,
        fuel_flow_lb_hr_per_engine = fuel_kg_s * KG_TO_LB * 3600.0,
        tsfc_per_hr = pare[ieTSFC, ENGINE_POINT, IMISSION] * 3600.0,
        bypass_ratio = pare[ieBPR, ENGINE_POINT, IMISSION],
        fan_pressure_ratio = fan_pr,
        lpc_pressure_ratio = lpc_pr,
        hpc_pressure_ratio = hpc_pr,
        overall_pressure_ratio = lpc_pr * hpc_pr,
        fan_speed = pare[ieNf, ENGINE_POINT, IMISSION],
        low_spool_speed = pare[ieN1, ENGINE_POINT, IMISSION],
        high_spool_speed = pare[ieN2, ENGINE_POINT, IMISSION],
        core_mass_flow_kg_s = pare[iemcore, ENGINE_POINT, IMISSION],
        converged = true,
        status = "ok",
    )
end

function generate_deck(
    output_file::AbstractString = OUTPUT_FILE;
    altitude_step_ft::Float64 = DEFAULT_ALTITUDE_STEP_FT,
    mach_step::Float64 = DEFAULT_MACH_STEP,
    throttle_step::Float64 = DEFAULT_THROTTLE_STEP,
)
    altitudes = grid_to_limit(MAX_ALTITUDE_FT, altitude_step_ft)
    machs = grid_to_limit(MAX_MACH, mach_step)
    throttles = grid_to_limit(1.0, throttle_step)
    rows = NamedTuple[]
    base = load_base_aircraft()
    maximum_tt4_K = base.pare[ieTt4, ipstatic, IMISSION]

    total_conditions = length(altitudes) * length(machs)
    condition_number = 0
    for altitude_ft in altitudes, mach in machs
        condition_number += 1
        @printf("[%3d/%3d] altitude %7.0f ft, Mach %.2f\n",
                condition_number, total_conditions, altitude_ft, mach)
        condition = condition_aircraft(base, altitude_ft, mach)
        endpoints = converged_endpoints(condition, maximum_tt4_K)

        for throttle in throttles
            ac, tt4 = solve_throttle(condition, endpoints, throttle)
            push!(rows, output_row(
                ac, altitude_ft, mach, throttle, endpoints, tt4
            ))
        end
    end

    mkpath(dirname(output_file))
    CSV.write(output_file, rows)
    println("Saved $(length(rows)) engine-deck points to: $output_file")
    println("All engine-deck points converged.")
    return rows
end

end

