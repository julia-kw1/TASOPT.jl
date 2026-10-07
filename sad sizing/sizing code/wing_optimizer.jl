module WingOptimizer

using TASOPT
using NLopt
using Logging
using CSV

using ..CFM56Installation
using ..OptimizationCommon
using ..PayloadMachSweep
using ..SO2Payload
using ..WorkflowPaths

include(TASOPT.__TASOPTindices__)

# x = [
#     mach,
#     AR,
#     Sref,
#     sweep,
#     altitude_ft,
#     range_nmi,
#     ]

LOWER_BOUNDS = [0.65, 7.5, 55.0, 10.0, 45_000.0, 1.0]
UPPER_BOUNDS = [0.89, 18.0, 175.0, 45.0, 70_000.0, 4_500.0]
INITIAL_STEP = [0.01, 0.5, 5.0, 1.0, 1_000.0, 100.0]

BASELINE_MODEL = :hybrid_reduced_radius_composite_no_interior
INITIAL_CL_LIMIT = 0.80

function reference_ceiling_point()
    path = WorkflowPaths.output_path(WorkflowPaths.payload_mach_file(BASELINE_MODEL))
    isfile(path) || error("Missing input-model ceiling sweep: $path")
    best = nothing
    for row in CSV.File(path)
        if row.ceiling_ft isa Number && isfinite(row.ceiling_ft) &&
           (best === nothing || row.ceiling_ft > best.ceiling_ft)
            best = (
                ceiling_ft = row.ceiling_ft,
                payload_lb = row.so2_payload_lb,
                mach = row.mach,
                range_nmi = row.range_nmi,
                fuel_reserve_fraction = row.fuel_reserve_fraction,
                minimum_roc_fpm = row.service_ceiling_roc_fpm,
            )
        end
    end
    best === nothing && error("No feasible input-model ceiling point in $path")
    return best
end

struct FixedInstalledCFM56Map end

function (calculator::FixedInstalledCFM56Map)(
    ac,
    case::String,
    imission::Int64,
    ip::Int64,
    initializes_engine::Bool,
    iterw::Int64 = 0,
)
    fixed_case = case == "design" ? "off_design" : case
    TASOPT.engine.tfwrap!(ac, fixed_case, imission, ip, initializes_engine, iterw)
    return nothing
end

function install_fixed_cfm56!(ac)
    old_model = ac.engine.model
    installed = (
        engine_weight = ac.parg[igWeng],
        bare_engine_weight = ac.parg[igWebare],
        nacelle_weight = ac.parg[igWnace],
        nacelle_area_fraction = ac.parg[igfSnace],
        nacelle_length = ac.parg[iglnace],
    )
    fixed_weight! = function (trial)
        trial.parg[igWeng] = installed.engine_weight
        trial.parg[igWebare] = installed.bare_engine_weight
        trial.parg[igWnace] = installed.nacelle_weight
        trial.parg[igfSnace] = installed.nacelle_area_fraction
        trial.parg[iglnace] = installed.nacelle_length
        return nothing
    end
    ac.engine = TASOPT.engine.Engine(
        TASOPT.engine.TurbofanModel(
            old_model.model_name,
            FixedInstalledCFM56Map(),
            old_model.weight_model_name,
            fixed_weight!,
            old_model.has_BLI_cores,
        ),
        ac.engine.data,
        ac.engine.heat_exchangers,
    )
    return ac
end

function prepare_baseline(tank_spec, payload_lb)
    ac = TASOPT.quickload_aircraft(WorkflowPaths.require_model(BASELINE_MODEL))
    target_fan_diameter_m = CFM56Installation.CFM56.fan_diameter_in * in_to_m
    if !isapprox(ac.parg[igdfan], target_fan_diameter_m; rtol = 1.0e-3)
        donor = TASOPT.quickload_aircraft(WorkflowPaths.require_model(:b737))
        CFM56Installation.install_donor_map!(ac, donor)
        CFM56Installation.calibrate_rated_turbine_temperature!(ac)
    end

    usable_diameter_m = 2.0 * (ac.fuselage.layout.radius - ac.fuselage.skin.thickness)
    tank = SO2Payload.size_tank(payload_lb, tank_spec;
        usable_inner_diameter_m = usable_diameter_m)
    ac.parm[imWpay, 1] = tank.total_payload_lb * lbf_to_N
    ac.parg[igWpay] = ac.parm[imWpay, 1]
    installed_engine = deepcopy(ac.engine)
    install_fixed_cfm56!(ac)
    return ac, installed_engine
end

function set_wing_planform!(ac, area_m2, aspect_ratio, sweep_deg)
    wing = ac.wing
    wing.layout.S = area_m2
    wing.layout.AR = aspect_ratio
    wing.layout.sweep = sweep_deg
    wing.layout.span = sqrt(area_m2 * aspect_ratio)
    wing.layout.ηs = max(wing.layout.ηs, wing.layout.ηo)

    ηo = wing.layout.ηo
    ηs = wing.layout.ηs
    λs = wing.inboard.λ
    λt = wing.outboard.λ
    Kc = ηo + 0.5 * (1.0 + λs) * (ηs - ηo) +
        0.5 * (λs + λt) * (1.0 - ηs)
    wing.layout.root_chord = area_m2 / (Kc * wing.layout.span)
    wing.inboard.co = wing.layout.root_chord
    wing.outboard.co = wing.inboard.co * λs
    TASOPT.structures.calculate_centroid_offset!(wing; calc_cma = true)
    return ac
end

function set_design_cl!(ac)
    TASOPT.set_ambient_conditions!(ac, ipcruise1)
    weight_N = ac.parg[igWMTO] * ac.para[iafracW, ipcruise1, 1]
    dynamic_pressure_Pa = 0.5 * ac.pare[ierho0, ipcruise1, 1] *
        ac.pare[ieu0, ipcruise1, 1]^2
    cruise_cl = weight_N / (dynamic_pressure_Pa * ac.wing.layout.S)
    ac.para[iaCL, ipclimb1+1:ipdescentn-1, 1] .= cruise_cl
    return cruise_cl
end

function size_candidate!(ac; iter = 40, printiter = false)
    set_design_cl!(ac)
    with_logger(NullLogger()) do
        TASOPT.size_aircraft!(
            ac;
            iter,
            initwgt = true,
            printiter = false,
            fixed_geometry = true,
        )
    end
    set_design_cl!(ac)
    with_logger(NullLogger()) do
        TASOPT.size_aircraft!(
            ac;
            iter = max(50, iter ÷ 2),
            initwgt = true,
            printiter,
            fixed_geometry = true,
        )
    end
    if sum(ac.pare[ieConvFail, :, 1]) != 0.0
        set_design_cl!(ac)
        with_logger(NullLogger()) do
            TASOPT.size_aircraft!(
                ac;
                iter = 50,
                initwgt = true,
                printiter = false,
                fixed_geometry = true,
            )
        end
    end
    return sum(ac.pare[ieConvFail, :, 1]) == 0.0
end

function apply_design!(ac, x)
    ac.para[iaMach, ipclimbn:ipdescent1, 1] .= x[1]
    set_wing_planform!(ac, x[3], x[2], x[4])
    ac.para[iaalt, ipclimbn:ipcruise1, 1] .= x[5] * ft_to_m
    ac.parm[imRange, 1] = x[6] * nmi_to_m
    ac.parg[igRange] = ac.parm[imRange, 1]
    return ac
end

function initial_design(ac, reference_altitude_ft)
    trial = deepcopy(ac)
    im = 2
    initial_altitude_ft = clamp(ac.para[iaalt, ipcruise1, 1] / ft_to_m,
        LOWER_BOUNDS[5], UPPER_BOUNDS[5])
    trial.para[iaMach, ipcruise1, im] = trial.para[iaMach, ipcruise1, 1]
    trial.para[iaalt, ipcruise1, im] = reference_altitude_ft * ft_to_m
    TASOPT.set_ambient_conditions!(trial, ipcruise1; im)
    weight_N = trial.parg[igWMTO] * trial.para[iafracW, ipcruise1, 1]
    dynamic_pressure_Pa = 0.5 * trial.pare[ierho0, ipcruise1, im] *
        trial.pare[ieu0, ipcruise1, im]^2
    required_area_m2 = 1.08 * weight_N /
        (dynamic_pressure_Pa * INITIAL_CL_LIMIT)

    initial_area_m2 = sqrt(ac.wing.layout.S * max(ac.wing.layout.S, required_area_m2))
    initial = [
        ac.para[iaMach, ipcruise1, 1],
        ac.wing.layout.AR * initial_area_m2 / ac.wing.layout.S,
        initial_area_m2,
        ac.wing.layout.sweep,
        initial_altitude_ft,
        ac.parm[imRange, 1] / nmi_to_m,
    ]
    return OptimizationCommon.clamp_initial!(initial, LOWER_BOUNDS, UPPER_BOUNDS)
end

function reference_mission(ac, reference, altitude_ft, settings)
    return PayloadMachSweep.evaluate_mission(
        ac, reference.payload_lb, reference.mach, altitude_ft; settings,
    )
end

function ceiling_above_reference(ac, reference, settings;
    lower_ft = reference.ceiling_ft,
    upper_ft = settings.altitude_search_limit_ft,
    tolerance_ft = settings.altitude_tolerance_ft,
)
    low = lower_ft
    high = upper_ft
    while high - low > tolerance_ft
        altitude_ft = (low + high) / 2.0
        feasible = try
            reference_mission(ac, reference, altitude_ft, settings).feasible
        catch err
            (err isa DomainError ||
                (err isa ErrorException && occursin("NaN", err.msg))) || rethrow()
            false
        end
        if feasible
            low = altitude_ft
        else
            high = altitude_ft
        end
    end
    return low
end

function run(;
    maxeval::Int = 500,
    tolerance::Float64 = 1.0e-5,
    tank::SO2Payload.TankSpec = SO2Payload.TankSpec(),
)
    WorkflowPaths.ensure_directories()
    reference = reference_ceiling_point()
    base, installed_engine = prepare_baseline(tank, reference.payload_lb)
    mission_settings_for_range(range_nmi) = PayloadMachSweep.Settings(
        range_nmi = range_nmi,
        fuel_reserve_fraction = reference.fuel_reserve_fraction,
        minimum_roc_fpm = reference.minimum_roc_fpm,
        tank = tank,
    )
    mission_settings = mission_settings_for_range(reference.range_nmi)
    println("Input-model ceiling to beat: $(reference.ceiling_ft) ft at " *
        "$(reference.payload_lb) lb SO2, M$(reference.mach), $(reference.range_nmi) nmi")
    iteration = Ref(0)
    phase = Ref(:altitude)
    best_aircraft = Ref{Any}(nothing)
    best_design = Ref{Any}(nothing)
    best_altitude_ft = Ref(-Inf)
    best_range_nmi = Ref(-Inf)
    ceiling_finalists = NamedTuple[]
    coarse_tolerance_ft = max(500.0, mission_settings.altitude_tolerance_ft)

    function objective(x, _)
        trial = apply_design!(deepcopy(base), x)
        try
        converged = size_candidate!(trial)
        if !converged
            iteration[] += 1
            diagnostics = (
                "Altitude" => x[5],
                "Range" => x[6],
                "ROC" => 0.0,
                "Fuel" => 0.0,
                "Sref" => trial.wing.layout.S / 0.092903,
                "AR" => trial.wing.layout.AR,
                "Mach" => trial.para[iaMach, ipcruise1, 1],
            )
            violations = [OptimizationCommon.Constraint(
                "Sizing failed", 0.0, 1.0, 1.0e10
            )]
            OptimizationCommon.print_progress(iteration[], diagnostics, violations)
            return -1.0e9
        end

        check_altitude_ft = phase[] == :altitude ?
            max(reference.ceiling_ft, best_altitude_ft[]) : best_altitude_ft[]
        check_settings = phase[] == :altitude ?
            mission_settings : mission_settings_for_range(x[6])
        mission = reference_mission(
            trial, reference, check_altitude_ft, check_settings,
        )
        roc_fpm = mission.toc_roc_fpm
        violations = OptimizationCommon.Constraint[]

        if !mission.feasible
            push!(violations, OptimizationCommon.Constraint(
                "Reference mission", 0.0, 1.0, 1.0e10
            ))
        end
        if trial.parg[igWfuel] > trial.parg[igWfmax]
            push!(violations, OptimizationCommon.Constraint(
                "Fuel volume", trial.parg[igWfuel], trial.parg[igWfmax], 1.0
            ))
        end
        if trial.parm[imlBF] > 1_720.0
            push!(violations, OptimizationCommon.Constraint(
                "Field length", trial.parm[imlBF], 1_720.0, 1.0
            ))
        end
        if trial.wing.layout.span > 51.9
            push!(violations, OptimizationCommon.Constraint(
                "Wing span", trial.wing.layout.span, 51.9, 1.0
            ))
        end

        ceiling_ft = if isempty(violations) && phase[] == :altitude
            ceiling_above_reference(trial, reference, mission_settings;
                lower_ft = check_altitude_ft,
                tolerance_ft = coarse_tolerance_ft,
            )
        else
            NaN
        end
        if isfinite(ceiling_ft) &&
           (length(ceiling_finalists) < 3 || ceiling_ft >= ceiling_finalists[end].ceiling_ft)
            push!(ceiling_finalists, (
                ceiling_ft = ceiling_ft,
                design = copy(x),
                aircraft = deepcopy(trial),
            ))
            sort!(ceiling_finalists; by = candidate -> candidate.ceiling_ft, rev = true)
            length(ceiling_finalists) > 3 && pop!(ceiling_finalists)
        end
        score = if isempty(violations)
            phase[] == :altitude ? ceiling_ft : x[6]
        elseif mission.mission_completed && mission.engine_converged
            minimum((
                roc_fpm - reference.minimum_roc_fpm,
                mission.cruise_range_m / nmi_to_m,
                mission.mtow_margin_lb,
                (trial.parg[igWfmax] - mission.fuel_N) / lbf_to_N,
            ))
        else
            -1.0e9
        end
        if isempty(violations) && (
            phase[] == :altitude ? ceiling_ft > best_altitude_ft[] : x[6] > best_range_nmi[]
        )
            phase[] == :altitude && (best_altitude_ft[] = ceiling_ft)
            best_range_nmi[] = phase[] == :altitude ? reference.range_nmi : x[6]
            best_design[] = copy(x)
            best_aircraft[] = deepcopy(trial)
        end
        iteration[] += 1
        diagnostics = (
            "DesignAlt" => x[5],
            "Ceiling" => ceiling_ft,
            "Range" => x[6],
            "ROC" => roc_fpm,
            "Fuel" => trial.parg[igWfuel] / lbf_to_N,
            "Sref" => trial.wing.layout.S / 0.092903,
            "AR" => trial.wing.layout.AR,
            "Mach" => trial.para[iaMach, ipcruise1, 1],
        )
        OptimizationCommon.print_progress(iteration[], diagnostics, violations)
        return score
        catch
            iteration[] += 1
            diagnostics = (
                "Altitude" => x[5],
                "Range" => x[6],
                "ROC" => 0.0,
                "Fuel" => 0.0,
                "Sref" => x[3] / 0.092903,
                "AR" => x[2],
                "Mach" => x[1],
            )
            violations = [OptimizationCommon.Constraint(
                "Engine map or sizing failure", 0.0, 1.0, 1.0e10
            )]
            OptimizationCommon.print_progress(iteration[], diagnostics, violations)
            return -1.0e9
        end
    end

    initial = initial_design(base, reference.ceiling_ft)
    results = Symbol[]
    for current_phase in (:altitude, :range)
        if current_phase == :range
            best_aircraft[] === nothing && error("No feasible altitude and range found.")
            initial = copy(best_design[])
            initial[6] = reference.range_nmi
        end
        phase[] = current_phase
        lower = copy(LOWER_BOUNDS)
        upper = copy(UPPER_BOUNDS)
        if current_phase == :range
            lower[5] = upper[5] = best_design[][5]
        end
        optimizer = NLopt.Opt(:LN_COBYLA, length(initial))
        optimizer.lower_bounds = lower
        optimizer.upper_bounds = upper
        optimizer.max_objective = objective
        optimizer.initial_step = INITIAL_STEP
        optimizer.ftol_rel = tolerance
        optimizer.maxeval = maxeval
        println("Optimizing $(current_phase)")
        OptimizationCommon.print_setup(
            optimizer, initial, lower, upper, tolerance, maxeval
        )
        _, _, result = NLopt.optimize(optimizer, initial)
        push!(results, result)
        if current_phase == :altitude
            isempty(ceiling_finalists) && error("No feasible altitude and range found.")
            best_altitude_ft[] = -Inf
            for candidate in ceiling_finalists
                refined_ft = ceiling_above_reference(
                    candidate.aircraft, reference, mission_settings;
                    lower_ft = candidate.ceiling_ft,
                    upper_ft = min(candidate.ceiling_ft + coarse_tolerance_ft,
                        mission_settings.altitude_search_limit_ft),
                )
                if refined_ft > best_altitude_ft[]
                    best_altitude_ft[] = refined_ft
                    best_range_nmi[] = reference.range_nmi
                    best_design[] = candidate.design
                    best_aircraft[] = candidate.aircraft
                end
            end
        end
    end

    output = deepcopy(best_aircraft[])
    output.engine = installed_engine
    TASOPT.weight_buildup(output)
    TASOPT.geometry(output)
    best_altitude_ft[] > reference.ceiling_ft || error(
        "Optimized ceiling $(best_altitude_ft[]) ft did not exceed input-model ceiling $(reference.ceiling_ft) ft."
    )
    verified = PayloadMachSweep.evaluate_mission(
        output, reference.payload_lb, reference.mach, best_altitude_ft[];
        settings = mission_settings_for_range(best_range_nmi[]),
    )
    verified.feasible || error("Saved aircraft does not meet the optimized range at the optimized altitude.")
    output_file = WorkflowPaths.model_path(:wing_optimized)
    TASOPT.quicksave_aircraft(output, output_file)
    println("Saved optimized wing model to: $output_file")
    println("Best altitude: $(best_altitude_ft[]) ft")
    println("Best range at that altitude: $(best_range_nmi[]) nmi")
    println("Reference-mission rate of climb: $(verified.toc_roc_fpm) ft/min")
    println("NLopt results: $(results)")
    return output
end

end
