module CFM56Installation

using TASOPT
using Printf

using ..AirframeChanges
using ..WorkflowPaths

include(TASOPT.__TASOPTindices__)

const CFM56 = (
    model = "CFM56-7B27",
    takeoff_thrust_lbf = 27_300.0,
    bpr = 5.1,
    opr = 32.7,
    fan_diameter_in = 61.0,
    dry_weight_lb = 5_359.0,
)
const CRUISE_ALTITUDE_FT = 35_000.0
const CRUISE_MACH = 0.78
const RANGE_UPPER_NMI = 4_000.0
const RANGE_TOLERANCE_NMI = 5.0
const RATED_TT4_BOUNDS_K = (1_600.0, 2_000.0)
const RATED_THRUST_TOLERANCE_LBF = 0.05
const RATED_THRUST_ITERATIONS = 24

const HARDWARE_PARAMETERS = (
    igHTRf, igHTRlc, igHTRhc,
    igrSnace, igrVnace,
    igfeadd, igfpylon,
    igGearf, igTmetal,
)
const MAP_FLOW_PARAMETERS = (iembfD, iemblcD, iembhcD, iembhtD, iembltD)
const MAP_AREA_PARAMETERS = (ieA25, ieA5, ieA7)

struct FixedAirframe
    mtow_N::Float64
    fuel_capacity_N::Float64
    nonengine_bow_N::Float64
    payload_max_N::Float64
end

lbf(weight_N) = weight_N / lbf_to_N

function quiet_output(f::Function)
    redirect_stdout(devnull) do
        redirect_stderr(devnull) do
            f()
        end
    end
end

function fixed_airframe(ac)
    return FixedAirframe(
        ac.parg[igWMTO],
        ac.parg[igWfmax],
        ac.parg[igWMTO] - ac.parg[igWfuel] - ac.parg[igWpay] - ac.parg[igWeng],
        ac.parg[igWpaymax],
    )
end

function set_cruise_condition!(ac)
    for mission in axes(ac.parm, 2)
        ac.para[iaalt, ipclimbn:ipcruise1, mission] .= CRUISE_ALTITUDE_FT * ft_to_m
        ac.para[iaMach, ipclimbn:ipdescent1, mission] .= CRUISE_MACH
        TASOPT.set_ambient_conditions!(ac, ipclimbn; im = mission)
        TASOPT.set_ambient_conditions!(ac, ipcruise1; im = mission)
    end
    return ac
end

function scale_donor_map!(ac, donor)
    fan_diameter_m = CFM56.fan_diameter_in * in_to_m
    fan_area_m2 = π * fan_diameter_m^2 * (1.0 - ac.parg[igHTRf]^2) / 4.0
    area_scale = fan_area_m2 / donor.pare[ieA2, ipcruise1, 1]

    ac.pare[ieA2, :, :] .= fan_area_m2
    for index in MAP_FLOW_PARAMETERS
        ac.pare[index, :, :] .*= area_scale
    end
    for index in MAP_AREA_PARAMETERS
        ac.pare[index, :, :] .*= area_scale
    end

    length_scale = sqrt(area_scale)
    ac.parg[igdfan] = fan_diameter_m
    ac.parg[igdlcomp] = donor.parg[igdlcomp] * length_scale
    ac.parg[igdhcomp] = donor.parg[igdhcomp] * length_scale
    ac.parg[igA5] = donor.parg[igA5] * area_scale
    ac.parg[igA7] = donor.parg[igA7] * area_scale
    return ac
end

function set_design_opr!(ac, donor)
    donor_opr = donor.pare[iepilcD, ipcruise1, 1] * donor.pare[iepihcD, ipcruise1, 1]
    ac.pare[iepihcD, :, :] .*= CFM56.opr / donor_opr
    return ac
end

function install_donor_map!(ac, donor)
    ac.engine = deepcopy(donor.engine)
    for mission in axes(ac.pare, 3)
        ac.pare[:, :, mission] .= donor.pare[:, :, 1]
    end
    for index in HARDWARE_PARAMETERS
        ac.parg[index] = donor.parg[index]
    end

    scale_donor_map!(ac, donor)
    set_design_opr!(ac, donor)
    TASOPT.engine.tfweightwrap!(ac)
    target_bare_N = CFM56.dry_weight_lb * lbf_to_N * ac.parg[igneng]
    ac.parg[igWebare] = target_bare_N
    ac.parg[igWeng] =
        (target_bare_N * (1.0 + ac.parg[igfeadd]) + ac.parg[igWnace]) *
        (1.0 + ac.parg[igfpylon])
    return ac
end

function rated_static_state(ac, tt4_K::Float64)
    trial = deepcopy(ac)
    mission = 1
    trial.parm[imDeltaTatm, mission] = 0.0
    trial.para[iaalt, ipstatic, mission] = 0.0
    trial.para[iaMach, ipstatic, mission] = 0.0
    trial.parg[igfBLIw] = 0.0
    trial.parg[igfBLIf] = 0.0
    trial.parg[igmofWpay] = 0.0
    trial.parg[igmofWMTO] = 0.0
    trial.parg[igPofWpay] = 0.0
    trial.parg[igPofWMTO] = 0.0
    TASOPT.set_ambient_conditions!(trial, ipstatic, 0.0; im = mission)
    trial.pare[ieTt4, ipstatic, mission] = tt4_K
    trial.pare[ieConvFail, ipstatic, mission] = 0.0
    quiet_output() do
        trial.engine.model.enginecalc!(trial, "off_design", mission, ipstatic, true)
    end
    trial.pare[ieConvFail, ipstatic, mission] == 0.0 ||
        error("The rated sea-level-static engine point did not converge.")
    return (
        thrust_lbf = trial.pare[ieFe, ipstatic, mission] / lbf_to_N,
        bpr = trial.pare[ieBPR, ipstatic, mission],
        opr = trial.pare[iepilc, ipstatic, mission] * trial.pare[iepihc, ipstatic, mission],
    )
end

function calibrate_rated_turbine_temperature!(ac)
    lower_K, upper_K = RATED_TT4_BOUNDS_K
    lower_thrust = rated_static_state(ac, lower_K).thrust_lbf
    upper_thrust = rated_static_state(ac, upper_K).thrust_lbf
    target = CFM56.takeoff_thrust_lbf
    lower_thrust <= target <= upper_thrust || error(
        "CFM56 rated thrust is not bracketed by the configured Tt4 bounds."
    )

    rated_tt4_K = 0.5 * (lower_K + upper_K)
    for _ in 1:RATED_THRUST_ITERATIONS
        rated_tt4_K = 0.5 * (lower_K + upper_K)
        thrust_lbf = rated_static_state(ac, rated_tt4_K).thrust_lbf
        abs(thrust_lbf - target) <= RATED_THRUST_TOLERANCE_LBF && break
        if thrust_lbf < target
            lower_K = rated_tt4_K
        else
            upper_K = rated_tt4_K
        end
    end

    ac.pare[ieTt4, ipstatic:ipclimbn, :] .= rated_tt4_K
    return rated_tt4_K
end

function set_reference_missions!(ac, airframe, knee_range_nmi, knee_payload_N, ferry_range_nmi)
    ac.parg[igWMTO] = airframe.mtow_N
    ac.parg[igWfmax] = airframe.fuel_capacity_N
    ac.parm[imRange, 1] = knee_range_nmi * nmi_to_m
    ac.parm[imWpay, 1] = knee_payload_N
    ac.parm[imRange, 2] = ferry_range_nmi * nmi_to_m
    ac.parm[imWpay, 2] = 0.0
    ac.parg[igRange] = ac.parm[imRange, 1]
    ac.parg[igWpay] = ac.parm[imWpay, 1]
    return ac
end

function full_fuel_range(ac, airframe, mission::Int, payload_N::Float64, itermax::Int)
    low_nmi = 0.0
    high_nmi = RANGE_UPPER_NMI
    while high_nmi - low_nmi > RANGE_TOLERANCE_NMI
        range_nmi = 0.5 * (low_nmi + high_nmi)
        trial = deepcopy(ac)
        trial.parm[imRange, mission] = range_nmi * nmi_to_m
        trial.parm[imWpay, mission] = payload_N
        quiet_output() do
            TASOPT.fly_mission!(
                trial,
                mission;
                itermax,
                initializes_engine = true,
                opt_prescribed_cruise_parameter = "altitude",
            )
        end
        feasible =
            trial.parm[imWfuel, mission] <= airframe.fuel_capacity_N &&
            trial.para[iagamV, ipclimbn, mission] >= 0.0 &&
            sum(trial.pare[ieConvFail, :, mission]) == 0.0
        if feasible
            low_nmi = range_nmi
        else
            high_nmi = range_nmi
        end
    end
    return low_nmi
end

function fly_reference_missions!(ac, itermax::Int)
    for mission in 1:2
        quiet_output() do
            TASOPT.fly_mission!(
                ac,
                mission;
                itermax,
                initializes_engine = true,
                opt_prescribed_cruise_parameter = "altitude",
            )
        end
        sum(ac.pare[ieConvFail, :, mission]) == 0.0 ||
            error("Reference mission $mission has an unconverged engine point.")
    end
    ac.parg[igWfuel] = ac.parm[imWfuel, 1]
    ac.parg[igWpay] = ac.parm[imWpay, 1]
    return ac
end

function prepare_airframe(; composite::Bool, no_interior::Bool, density_factor::Float64)
    ac = TASOPT.quickload_aircraft(WorkflowPaths.require_model(:e175))
    composite && AirframeChanges.apply_composites!(ac, density_factor)
    no_interior && AirframeChanges.remove_interior!(ac)
    return ac
end

function run(
    output_model::Symbol;
    composite::Bool = false,
    no_interior::Bool = false,
    density_factor::Float64 = 0.85,
    itermax::Int = 100,
)
    ac = prepare_airframe(; composite, no_interior, density_factor)
    size(ac.parm, 2) >= 2 || error("The E175 model must contain two missions.")
    airframe = fixed_airframe(ac)

    # Install last so airframe sizing cannot overwrite the CFM56 data.
    donor = TASOPT.quickload_aircraft(WorkflowPaths.require_model(:b737))
    install_donor_map!(ac, donor)
    set_cruise_condition!(ac)
    rated_tt4_K = calibrate_rated_turbine_temperature!(ac)

    knee_payload_N = clamp(
        airframe.mtow_N - airframe.nonengine_bow_N - ac.parg[igWeng] - airframe.fuel_capacity_N,
        0.0,
        airframe.payload_max_N,
    )
    set_reference_missions!(ac, airframe, 1_000.0, knee_payload_N, 1_000.0)
    knee_range_nmi = full_fuel_range(ac, airframe, 1, knee_payload_N, itermax)
    ferry_range_nmi = full_fuel_range(ac, airframe, 2, 0.0, itermax)
    set_reference_missions!(ac, airframe, knee_range_nmi, knee_payload_N, ferry_range_nmi)
    fly_reference_missions!(ac, itermax)

    output_file = WorkflowPaths.model_path(output_model)
    TASOPT.quicksave_aircraft(ac, output_file)
    rated = rated_static_state(ac, rated_tt4_K)
    println("Saved $(WorkflowPaths.model_title(output_model)) to: $output_file")
    @printf("  Rated thrust: %.1f lbf per engine\n", rated.thrust_lbf)
    @printf("  Fan diameter: %.1f in\n", ac.parg[igdfan] / in_to_m)
    return ac
end

end
