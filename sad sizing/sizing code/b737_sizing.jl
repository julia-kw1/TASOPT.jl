module B737Sizing

using TASOPT
using NLopt

using ..OptimizationCommon
using ..WorkflowPaths

include(TASOPT.__TASOPTindices__)

const LOWER_BOUNDS = [
    8.0, 0.45, 25.0, 10_000.0, 0.55, 0.2, 0.115,
    0.09, 0.9, 0.7, 1_400.0, 9.0, 1.25, 4.5,
]
const UPPER_BOUNDS = [
    15.0, 0.75, 30.0, 11_000.0, 0.85, 0.4, 0.125,
    0.115, 1.3, 1.0, 1_650.0, 11.0, 2.0, 5.5,
]
const INITIAL_STEP = [
    0.5, 0.05, 0.1, 10.0, 0.01, 0.01, 0.01,
    0.01, 0.01, 0.01, 100.0, 0.5, 0.05, 1.0,
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
    ac.pare[iepihc, ipcruise1:ipcruise2, 1] .= x[12]
    ac.pare[iepif, ipcruise1:ipcruise2, 1] .= x[13]
    ac.pare[iepilc, ipcruise1, 1] = 3.0
    ac.pare[ieBPR, ipcruise1:ipcruise2, 1] .= final ? x[14] : 5.142
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
        ac.pare[iepihc, ipcruise1, 1],
        ac.pare[iepif, ipcruise1, 1],
        5.10,
    ]
    return OptimizationCommon.clamp_initial!(initial, LOWER_BOUNDS, UPPER_BOUNDS)
end

function penalty(ac)
    violations = OptimizationCommon.Constraint[]
    total = 0.0

    span_limit_m = 34.3
    if ac.wing.span > span_limit_m
        difference = ac.wing.span / span_limit_m - 1.0
        value = 25.0 * ac.parg[igWpay] * difference^2
        total += value
        push!(violations, OptimizationCommon.Constraint("Wing span", ac.wing.span, span_limit_m, value))
    end

    climb_limit = ac.parg[iggtocmin]
    climb_gradient = ac.para[iagamV, ipclimbn, 1]
    if climb_gradient < climb_limit
        difference = 1.0 - climb_gradient / climb_limit
        value = ac.parg[igWpay] * difference^2
        total += value
        push!(violations, OptimizationCommon.Constraint("Climb gradient", climb_gradient, climb_limit, value))
    end

    compressor_limit_K = 1_050.0
    compressor_temperature_K = maximum(ac.pare[ieTt3, :, 1])
    if compressor_temperature_K > compressor_limit_K
        difference = compressor_temperature_K / compressor_limit_K - 1.0
        value = 5.0 * ac.parg[igWpay] * difference^2
        total += value
        push!(violations, OptimizationCommon.Constraint("Turbine Tt3", compressor_temperature_K, compressor_limit_K, value))
    end

    if ac.parg[igWfuel] > ac.parg[igWfmax]
        difference = ac.parg[igWfuel] / ac.parg[igWfmax] - 1.0
        value = 10.0 * ac.parg[igWpay] * difference^2
        total += value
        push!(violations, OptimizationCommon.Constraint("Fuel volume", ac.parg[igWfuel], ac.parg[igWfmax], value))
    end

    if ac.parm[imWTO, 1] > ac.parg[igWMTO]
        difference = ac.parm[imWTO, 1] / ac.parg[igWMTO] - 1.0
        total += 10.0 * ac.parg[igWpay] * difference^2
    end

    fan_limit_m = 2.0
    if ac.parg[igdfan] > fan_limit_m
        difference = ac.parg[igdfan] / fan_limit_m - 1.0
        value = ac.parg[igWpay] * difference^2
        total += value
        push!(violations, OptimizationCommon.Constraint("Fan diameter", ac.parg[igdfan], fan_limit_m, value))
    end

    field_length_limit_m = 2_300.0
    if ac.parm[imlBF] > field_length_limit_m
        difference = ac.parm[imlBF] / field_length_limit_m - 1.0
        value = ac.parg[igWpay] * difference^2
        total += value
        push!(violations, OptimizationCommon.Constraint("Field length", ac.parm[imlBF], field_length_limit_m, value))
    end

    thrust_target_N = 27_300.4 * lbf_to_N
    thrust_error = abs(ac.pare[ieFe, iptakeoff, 1] / thrust_target_N - 1.0)
    if thrust_error > 0.02
        value = 50.0 * ac.parg[igWpay] * thrust_error^2
        total += value
        push!(violations, OptimizationCommon.Constraint("Target FTO", ac.pare[ieFe, iptakeoff, 1], thrust_target_N, value))
    end

    mtow_target_N = 155_500.0 * lbf_to_N
    mtow_error = abs(ac.parg[igWMTO] / mtow_target_N - 1.0)
    if mtow_error > 0.01
        value = 50.0 * ac.parg[igWpay] * mtow_error^2
        total += value
        push!(violations, OptimizationCommon.Constraint("Target MTOW", ac.parg[igWMTO], mtow_target_N, value))
    end

    return total, violations
end

function run(; maxeval::Int = 500, tolerance::Float64 = 1.0e-5)
    input_file = WorkflowPaths.require_model(:b737_input)
    aircraft = TASOPT.read_aircraft_model(input_file; templatefile = input_file)
    iteration = Ref(0)

    function objective(x, _)
        apply_design!(aircraft, x)
        try
            TASOPT.size_aircraft!(aircraft; iter = 50, printiter = false)
        catch error
            println("Aircraft sizing failed: $error")
            return 1.0e10
        end

        pfei = aircraft.parm[imPFEI]
        penalty_value, violations = penalty(aircraft)
        iteration[] += 1
        diagnostics = (
            "PFEI" => pfei,
            "FTO" => aircraft.pare[ieFe, iptakeoff, 1] / lbf_to_N,
            "TSFC" => aircraft.pare[ieTSFC, ipcruise1, 1] * 3_600.0,
            "OPR" => aircraft.pare[iept3, ipcruise1, 1] / aircraft.pare[iept2, ipcruise1, 1],
            "Wfuel" => aircraft.parg[igWfuel] / lbf_to_N,
            "CL" => aircraft.para[iaCL, ipcruise1, 1],
            "MTOW" => aircraft.parg[igWMTO] / lbf_to_N,
        )
        OptimizationCommon.print_progress(iteration[], diagnostics, violations)
        return pfei + penalty_value
    end

    initial = initial_design(aircraft)
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
    output = apply_design!(deepcopy(aircraft), optimum; final = true)
    TASOPT.size_aircraft!(output; iter = 100, printiter = true)
    TASOPT.weight_buildup(output)
    TASOPT.geometry(output)
    output_file = WorkflowPaths.model_path(:b737)
    TASOPT.quicksave_aircraft(output, output_file)
    println("Saved B737 donor model to: $output_file")
    println("NLopt result: $result")
    return output
end

end
