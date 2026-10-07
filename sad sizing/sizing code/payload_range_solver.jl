function ensure_off_design_mission!(ac::TASOPT.aircraft)
    if size(ac.parm, 2) >= 2
        return ac
    end

    ac.parm = cat(ac.parm[:, 1], ac.parm[:, 1], dims = 2)
    ac.para = cat(ac.para[:, :, 1], ac.para[:, :, 1], dims = 3)
    ac.pare = cat(ac.pare[:, :, 1], ac.pare[:, :, 1], dims = 3)

    for HX in ac.engine.heat_exchangers
        HX.HXgas_mission = cat(HX.HXgas_mission[:, 1], HX.HXgas_mission[:, 1], dims = 2)
    end

    return ac
end

function configure_off_design_mission!(
    ac::TASOPT.aircraft,
    payload_N::Float64,
    range_nmi::Float64,
    altitude_ft::Float64,
    mach::Float64;
    imission::Int = 2,
)
    altitude_m = altitude_ft * ft_to_m
    range_m = range_nmi * nmi_to_m

    ac.parm[imWpay, imission] = payload_N
    ac.parm[imRange, imission] = range_m

    ac.para[iaMach, ipclimbn:ipdescent1, imission] .= mach
    ac.para[iaalt, ipclimbn:ipcruise1, imission] .= altitude_m
    ac.para[iaCL, :, imission] .= ac.para[iaCL, :, 1]

    TASOPT.set_ambient_conditions!(ac, ipclimbn; im = imission)
    TASOPT.set_ambient_conditions!(ac, ipcruise1; im = imission)

    return ac
end

function run_off_design_mission(
    ac_base::TASOPT.aircraft,
    payload_N::Float64,
    range_nmi::Float64,
    settings::Settings;
    imission::Int = 2,
)
    ac = deepcopy(ac_base)
    configure_off_design_mission!(ac, payload_N, range_nmi, settings.altitude_ft, settings.mach; imission)

    WBOW = empty_weight(ac)
    WMTO = ac.parg[igWMTO]
    Wusable_fuel = ac.parg[igWfmax] * ac.parg[igrWfmax]

    try
        redirect_stdout(devnull) do
            redirect_stderr(devnull) do
                TASOPT.fly_mission!(ac, imission; itermax = 50, initializes_engine = true,
                    opt_prescribed_cruise_parameter = "altitude")
            end
        end
    catch err
        return MissionResult(
            feasible = false,
            range_nmi = range_nmi,
            payload_N = payload_N,
            fuel_N = NaN,
            WTO_N = NaN,
            fuel_margin_N = -Inf,
            mtow_margin_N = -Inf,
            pfei = NaN,
            message = sprint(showerror, err),
        )
    end

    fuel_N = ac.parm[imWfuel, imission]
    WTO_N = WBOW + payload_N + fuel_N
    fuel_margin_N = Wusable_fuel - fuel_N
    mtow_margin_N = WMTO - WTO_N

    feasible = isfinite(fuel_N) &&
               isfinite(WTO_N) &&
               fuel_N >= -settings.weight_tolerance_N &&
               WTO_N >= -settings.weight_tolerance_N &&
               fuel_margin_N >= -settings.weight_tolerance_N &&
               mtow_margin_N >= -settings.weight_tolerance_N

    return MissionResult(
        feasible = feasible,
        range_nmi = range_nmi,
        payload_N = payload_N,
        fuel_N = fuel_N,
        WTO_N = WTO_N,
        fuel_margin_N = fuel_margin_N,
        mtow_margin_N = mtow_margin_N,
        pfei = ac.parm[imPFEI, imission],
        message = feasible ? "" : "Fuel or takeoff weight limit exceeded.",
    )
end

function find_first_feasible_range(ac::TASOPT.aircraft, payload_N::Float64, settings::Settings)
    first_trial = run_off_design_mission(ac, payload_N, 0.0, settings)
    first_trial.feasible && return first_trial

    lower_nmi = 0.0
    upper_nmi = min(100.0, settings.range_search_limit_nmi)
    while true
        trial = run_off_design_mission(ac, payload_N, upper_nmi, settings)
        if trial.feasible
            best = trial
            while upper_nmi - lower_nmi > settings.range_tolerance_nmi
                range_nmi = 0.5 * (lower_nmi + upper_nmi)
                candidate = run_off_design_mission(ac, payload_N, range_nmi, settings)
                if candidate.feasible
                    upper_nmi, best = range_nmi, candidate
                else
                    lower_nmi = range_nmi
                end
            end
            return best
        end
        upper_nmi >= settings.range_search_limit_nmi && break
        lower_nmi = upper_nmi
        upper_nmi = min(settings.range_search_limit_nmi, max(2.0 * upper_nmi, upper_nmi + 100.0))
    end
    error("No feasible mission was found for payload $(payload_N / lbf_to_N) lb.")
end

function find_max_range(ac::TASOPT.aircraft, payload_N::Float64, settings::Settings)
    best = find_first_feasible_range(ac, payload_N, settings)
    low = best.range_nmi
    low >= settings.range_search_limit_nmi && return best

    seed_range_nmi = max(ac.parm[imRange, 1] / nmi_to_m, 100.0)
    high = min(settings.range_search_limit_nmi, max(seed_range_nmi, low + 100.0))
    trial = run_off_design_mission(ac, payload_N, high, settings)

    while trial.feasible
        low = high
        best = trial
        high >= settings.range_search_limit_nmi && return best
        high = min(settings.range_search_limit_nmi, max(high * 1.5, high + 100.0))
        trial = run_off_design_mission(ac, payload_N, high, settings)
    end

    while (high - low) > settings.range_tolerance_nmi
        mid = 0.5 * (low + high)
        trial = run_off_design_mission(ac, payload_N, mid, settings)
        if trial.feasible
            low = mid
            best = trial
        else
            high = mid
        end
    end

    return best
end

function find_range_for_target_fuel(
    ac::TASOPT.aircraft,
    payload_N::Float64,
    target_fuel_N::Float64,
    settings::Settings,
)
    low = find_first_feasible_range(ac, payload_N, settings)
    low.fuel_N <= target_fuel_N + settings.fuel_target_tolerance_N ||
        error("Target fuel is below the first feasible mission requirement.")

    best_above = low.fuel_N >= target_fuel_N - settings.fuel_target_tolerance_N ? low : nothing
    seed_range_nmi = max(ac.parm[imRange, 1] / nmi_to_m, 100.0)
    high_range = min(settings.range_search_limit_nmi, max(seed_range_nmi, low.range_nmi + 100.0))
    high = run_off_design_mission(ac, payload_N, high_range, settings)

    while high.feasible && high.fuel_N < target_fuel_N - settings.fuel_target_tolerance_N &&
          high_range < settings.range_search_limit_nmi
        low = high
        high_range = min(settings.range_search_limit_nmi, max(high_range * 1.5, high_range + 100.0))
        high = run_off_design_mission(ac, payload_N, high_range, settings)
    end

    if high.feasible && high.fuel_N >= target_fuel_N - settings.fuel_target_tolerance_N
        best_above = high
    end

    search_high_range = high_range
    while (search_high_range - low.range_nmi) > settings.range_tolerance_nmi
        mid_range = 0.5 * (low.range_nmi + search_high_range)
        mid = run_off_design_mission(ac, payload_N, mid_range, settings)

        if !mid.feasible
            search_high_range = mid_range
            continue
        end

        if mid.fuel_N < target_fuel_N - settings.fuel_target_tolerance_N
            low = mid
        else
            best_above = mid
            search_high_range = mid_range
        end
    end

    if best_above === nothing
        abs(low.fuel_N - target_fuel_N) <= settings.fuel_target_tolerance_N ||
            error("Could not reach the target fuel at FL$(round(Int, settings.altitude_ft / 100.0)) and M$(settings.mach).")
        return low
    end

    return abs(low.fuel_N - target_fuel_N) <= abs(best_above.fuel_N - target_fuel_N) ? low : best_above
end

function point_from_result(label::String, result::MissionResult; note::String = "")
    return PayloadRangePoint(
        label,
        result.range_nmi,
        result.payload_N,
        result.fuel_N,
        result.WTO_N,
        result.fuel_margin_N,
        result.mtow_margin_N,
        result.pfei,
        note,
    )
end

function solve_payload_points(
    ac::TASOPT.aircraft,
    payload_N::Float64,
    usable_fuel_payload_limit_N::Float64,
    WBOW::Float64,
    maximum_fuel_N::Float64,
    settings::Settings,
)
    # Maximum fuel respects both tank capacity and MTOW.
    usable_fuel_payload_N = min(payload_N, usable_fuel_payload_limit_N)
    max_payload_fuel_N = min(maximum_fuel_N, ac.parg[igWMTO] - WBOW - payload_N)
    max_payload_fuel_N >= 0.0 || return nothing

    max_payload_result = try
        find_range_for_target_fuel(ac, payload_N, max_payload_fuel_N, settings)
    catch
        return nothing
    end
    max_fuel_result = try
        find_range_for_target_fuel(ac, usable_fuel_payload_N, maximum_fuel_N, settings)
    catch
        return nothing
    end

    return max_payload_result, max_fuel_result
end

function build_payload_range_points(ac::TASOPT.aircraft, settings::Settings)
    WMTO = ac.parg[igWMTO]
    WBOW = empty_weight(ac)
    Wusable_fuel = ac.parg[igWfmax] * ac.parg[igrWfmax]
    maximum_fuel_N = min(Wusable_fuel, max(0.0, WMTO - WBOW))
    Wpaymax = ac.parg[igWpaymax]
    MRW = WMTO - maximum_fuel_N
    usable_fuel_payload_requested_N = clamp(MRW - WBOW, 0.0, Wpaymax)
    payload_step_N = settings.payload_reduction_step_lb * lbf_to_N
    payload_step_N > 0.0 || error("Payload reduction step must be positive.")

    ferry_result = find_max_range(ac, 0.0, settings)

    requested_payload_N = Wpaymax
    solved_points = solve_payload_points(
        ac, requested_payload_N, usable_fuel_payload_requested_N, WBOW, maximum_fuel_N, settings,
    )

    if solved_points === nothing
        settings.enable_payload_recovery || error(
            "No feasible Origin, Max Payload, and Usable Fuel points at the requested payload of " *
            "$(fmt_number(lbf(Wpaymax))) lb. Enable payload recovery to reduce payload."
        )

        lower_payload_N = 0.0
        solved_points = solve_payload_points(
            ac, lower_payload_N, usable_fuel_payload_requested_N, WBOW, maximum_fuel_N, settings,
        )
        solved_points === nothing && error(
            "Could not recover feasible Origin, Max Payload, and Usable Fuel points by reducing payload to zero."
        )

        upper_payload_N = Wpaymax
        while (upper_payload_N - lower_payload_N) > payload_step_N
            trial_payload_N = 0.5 * (lower_payload_N + upper_payload_N)
            trial_points = solve_payload_points(
                ac, trial_payload_N, usable_fuel_payload_requested_N, WBOW, maximum_fuel_N, settings,
            )
            if trial_points === nothing
                upper_payload_N = trial_payload_N
            else
                lower_payload_N = trial_payload_N
                solved_points = trial_points
            end
        end
        requested_payload_N = lower_payload_N
    end

    max_payload_result, max_fuel_result = solved_points

    payload_reduction_N = Wpaymax - requested_payload_N
    recovery_note = payload_reduction_N > settings.weight_tolerance_N ?
        "Payload reduced by $(fmt_number(lbf(payload_reduction_N))) lb to recover feasibility." : ""
    usable_fuel_payload_N = min(requested_payload_N, usable_fuel_payload_requested_N)
    usable_fuel_payload_reduction_N = usable_fuel_payload_requested_N - usable_fuel_payload_N
    usable_fuel_note = usable_fuel_payload_reduction_N > settings.weight_tolerance_N ?
        "Payload reduced by $(fmt_number(lbf(usable_fuel_payload_reduction_N))) lb to recover feasibility." : ""
    if Wusable_fuel - maximum_fuel_N > settings.weight_tolerance_N
        usable_fuel_note *= " Fuel limited to $(fmt_number(lbf(maximum_fuel_N))) lb by MTOW."
        usable_fuel_note = String(strip(usable_fuel_note))
    end

    return [
        PayloadRangePoint("Origin", 0.0, requested_payload_N, 0.0, WBOW + requested_payload_N, maximum_fuel_N, WMTO - (WBOW + requested_payload_N), NaN, recovery_note),
        point_from_result("Max Payload", max_payload_result; note = recovery_note),
        point_from_result("Usable Fuel", max_fuel_result; note = usable_fuel_note),
        point_from_result("Ferry", ferry_result),
    ], Wpaymax, usable_fuel_payload_requested_N
end
