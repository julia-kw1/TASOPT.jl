module E175Sizing

using TASOPT
using NLopt

using ..OptimizationCommon
using ..WorkflowPaths

include(TASOPT.__TASOPTindices__)

const LOWER_BOUNDS = [
    8.75, 0.42, 21.0, 10_134.6, 0.58, 0.15, 0.10, 0.090,
    1.00, 0.60, 1_300.0, 1_500.0, 15.5, 1.30, 4.9,
]
const UPPER_BOUNDS = [
    8.9, 0.65, 23.5, 11_201.4, 0.605, 0.29, 0.125, 0.11,
    1.35, 0.95, 1_550.0, 1_700.0, 22.0, 1.70, 5.3,
]
const INITIAL_STEP = [
    0.05, 0.05, 0.1, 150.0, 0.01, 0.01, 0.01, 0.01,
    0.01, 0.01, 50.0, 50.0, 0.5, 0.05, 1.0,
]

function apply_design!(ac, x; final = false)
    ac.wing.layout.AR = x[1]
    ac.wing.layout.sweep = x[3]
    ac.wing.inboard.λ = x[5]
    ac.wing.outboard.λ = x[6]
    ac.wing.inboard.cross_section.thickness_to_chord = x[7]
    ac.wing.outboard.cross_section.thickness_to_chord = x[8]

    cl_start = final ? ipcruise1 : ipclimb1 + 1
    ac.para[iaCL, cl_start:ipdescentn-1, 1] .= x[2]
    ac.para[iaalt, ipclimbn:ipcruise1, 1] .= x[4]
    ac.para[iarcls, ipclimb2:ipdescent4, 1] .= x[9]
    ac.para[iarclt, ipclimb2:ipdescent4, 1] .= x[10]
    ac.pare[ieTt4, ipcruise1:ipcruise2, 1] .= x[11]
    ac.pare[ieTt4, iptakeoff:ipclimbn, 1] .= x[12]
    ac.pare[iepihc, ipcruise1:ipcruise2, 1] .= x[13]
    ac.pare[iepif, ipcruise1:ipcruise2, 1] .= x[14]
    ac.pare[ieBPR, ipcruise1:ipcruise2, 1] .= x[15]
    return ac
end

function initial_design(ac)
    initial = [
        ac.wing.layout.AR,
        ac.para[iaCL, ipcruise1, 1],
        ac.wing.layout.sweep,
        ac.para[iaalt, ipcruise1, 1],
        ac.wing.inboard.λ,
        ac.wing.outboard.λ,
        ac.wing.inboard.cross_section.thickness_to_chord,
        ac.wing.outboard.cross_section.thickness_to_chord,
        ac.para[iarcls, ipcruise1, 1],
        ac.para[iarclt, ipcruise1, 1],
        ac.pare[ieTt4, ipcruise1, 1],
        ac.pare[ieTt4, iptakeoff, 1],
        ac.pare[iepihc, ipcruise1, 1],
        ac.pare[iepif, ipcruise1, 1],
        5.0,
    ]
    return OptimizationCommon.clamp_initial!(initial, LOWER_BOUNDS, UPPER_BOUNDS)
end

function penalty(ac)
    violations = OptimizationCommon.Constraint[]
    total = 0.0

    span_limit = 85.33 * ft_to_m
    if ac.wing.span > span_limit
        error = ac.wing.span / span_limit - 1.0
        value = ac.parg[igWpay] * error^2
        total += value
        push!(violations, OptimizationCommon.Constraint("Wing span", ac.wing.span, span_limit, value))
    end

    climb_limit = ac.parg[iggtocmin]
    climb_gradient = ac.para[iagamV, ipclimbn, 1]
    if climb_gradient < climb_limit
        error = 1.0 - climb_gradient / climb_limit
        value = 100.0 * ac.parg[igWpay] * error^2
        total += value
        push!(violations, OptimizationCommon.Constraint("Climb gradient", climb_gradient, climb_limit, value))
    end

    if ac.parg[igWfuel] > ac.parg[igWfmax]
        error = ac.parg[igWfuel] / ac.parg[igWfmax] - 1.0
        value = 10.0 * ac.parg[igWpay] * error^2
        total += value
        push!(violations, OptimizationCommon.Constraint("Fuel volume", ac.parg[igWfuel], ac.parg[igWfmax], value))
    end

    vane_limit_K = 1_280.33
    vane_temperature_K = maximum(ac.pare[ieTmet1, :, 1])
    if vane_temperature_K > vane_limit_K
        error = vane_temperature_K / vane_limit_K - 1.0
        total += 5.0 * ac.parg[igWpay] * error^2
    end

    if ac.parm[imWTO, 1] > ac.parg[igWMTO]
        error = ac.parm[imWTO, 1] / ac.parg[igWMTO] - 1.0
        total += 10.0 * ac.parg[igWpay] * error^2
    end

    fan_limit_m = 1.25
    if ac.parg[igdfan] > fan_limit_m
        error = ac.parg[igdfan] / fan_limit_m - 1.0
        value = 50.0 * ac.parg[igWpay] * error^2
        total += value
        push!(violations, OptimizationCommon.Constraint("Fan diameter", ac.parg[igdfan], fan_limit_m, value))
    end

    field_length_limit_m = 1_720.0
    if ac.parm[imlBF] > field_length_limit_m
        error = ac.parm[imlBF] / field_length_limit_m - 1.0
        value = 10.0 * ac.parg[igWpay] * error^2
        total += value
        push!(violations, OptimizationCommon.Constraint("Field length", ac.parm[imlBF], field_length_limit_m, value))
    end

    area_target_m2 = 783.0 * 0.092903
    area_error = abs(ac.wing.layout.S / area_target_m2 - 1.0)
    if area_error > 0.01
        value = 50.0 * ac.parg[igWpay] * area_error^2
        total += value
        push!(violations, OptimizationCommon.Constraint("Wing area", ac.wing.layout.S, area_target_m2, value))
    end

    fuel_target_N = 20_785.0 * lbf_to_N
    fuel_error = abs(ac.parg[igWfmax] / fuel_target_N - 1.0)
    if fuel_error > 0.02
        value = 25.0 * ac.parg[igWpay] * fuel_error^2
        total += value
        push!(violations, OptimizationCommon.Constraint("Target fuel", ac.parg[igWfmax], fuel_target_N, value))
    end

    mtow_target_N = 82_673.0 * lbf_to_N
    mtow_error = abs(ac.parg[igWMTO] / mtow_target_N - 1.0)
    if mtow_error > 0.01
        value = 50.0 * ac.parg[igWpay] * mtow_error^2
        total += value
        push!(violations, OptimizationCommon.Constraint("Target MTOW", ac.parg[igWMTO], mtow_target_N, value))
    end

    return total, violations
end

function run(; maxeval::Int = 500, tolerance::Float64 = 1.0e-5)
    input_file = WorkflowPaths.require_model(:e175_input)
    source = TASOPT.read_aircraft_model(input_file; templatefile = input_file)
    base = deepcopy(source)
    iteration = Ref(0)

    function objective(x, _)
        trial = apply_design!(deepcopy(base), x)
        try
            TASOPT.size_aircraft!(trial; iter = 50, printiter = false)
        catch error
            println("Aircraft sizing failed: $error")
            return 1.0e10
        end

        pfei = trial.parm[imPFEI]
        penalty_value, violations = penalty(trial)
        iteration[] += 1
        diagnostics = (
            "PFEI" => pfei,
            "FTO" => trial.pare[ieFe, iptakeoff, 1] / lbf_to_N,
            "TSFC" => trial.pare[ieTSFC, ipcruise1, 1] * 3_600.0,
            "OPR" => trial.pare[iept3, ipcruise1, 1] / trial.pare[iept2, ipcruise1, 1],
            "Wfuel" => trial.parg[igWfmax] / lbf_to_N,
            "CL" => trial.para[iaCL, ipcruise1, 1],
            "MTOW" => trial.parg[igWMTO] / lbf_to_N,
            "Sref" => trial.wing.layout.S / 0.092903,
            "Wwing" => trial.wing.weight / lbf_to_N,
        )
        OptimizationCommon.print_progress(iteration[], diagnostics, violations)
        return pfei + penalty_value
    end

    initial = initial_design(source)
    optimizer = NLopt.Opt(:LN_NELDERMEAD, length(initial))
    optimizer.lower_bounds = LOWER_BOUNDS
    optimizer.upper_bounds = UPPER_BOUNDS
    optimizer.min_objective = objective
    optimizer.initial_step = INITIAL_STEP
    optimizer.ftol_rel = tolerance
    optimizer.maxeval = maxeval
    OptimizationCommon.print_setup(
        optimizer, initial, LOWER_BOUNDS, UPPER_BOUNDS, tolerance, maxeval
    )

    _, optimum, result = NLopt.optimize(optimizer, initial)
    aircraft = apply_design!(deepcopy(source), optimum; final = true)
    TASOPT.size_aircraft!(aircraft; iter = 100, printiter = true)
    TASOPT.weight_buildup(aircraft)
    TASOPT.geometry(aircraft)
    output_file = WorkflowPaths.model_path(:e175)
    TASOPT.quicksave_aircraft(aircraft, output_file)
    println("Saved E175 model to: $output_file")
    println("NLopt result: $result")
    return aircraft
end

end
