"""Payload/Mach service-ceiling sweep for selectable E175 quicksaves."""

module PayloadMachSweep

using TASOPT
using Printf
using CSV

using ..SO2Payload
using ..WorkflowPaths

include(TASOPT.__TASOPTindices__)

Base.@kwdef struct Settings
    payload_min_lb::Float64 = 0.0
    payload_max_lb::Float64 = 10_000.0
    mach_min::Float64 = 0.70
    mach_max::Float64 = 0.80
    range_nmi::Float64 = 500.0
    fuel_reserve_fraction::Float64 = 0.05
    minimum_roc_fpm::Float64 = 100.0
    ps_altitude_ft::Float64 = 35_000.0
    altitude_lower_ft::Float64 = 35_000.0
    altitude_initial_upper_ft::Float64 = 65_000.0
    altitude_search_limit_ft::Float64 = 70_000.0
    altitude_tolerance_ft::Float64 = 50.0
    weight_tolerance_N::Float64 = 10.0
    mission_itermax::Int = 75
    ceiling_requires_mission_feasibility::Bool = true
    tank::SO2Payload.TankSpec = SO2Payload.TankSpec()
end

weight_lb(weight_N) = weight_N / lbf_to_N

function payload_package(base, so2_payload_lb, settings::Settings)
    spec = settings.tank
    usable_inner_diameter_m = max(
        0.0,
        2.0 * (base.fuselage.layout.radius - base.fuselage.skin.thickness),
    )
    tank = SO2Payload.size_tank(so2_payload_lb, spec; usable_inner_diameter_m)
    fit = SO2Payload.fuselage_fit(base, tank, spec)
    structural_payload_limit_lb = weight_lb(base.parg[igWpaymax])
    structural_fit = tank.total_payload_lb <= structural_payload_limit_lb
    feasible = fit.feasible && structural_fit

    reasons = String[]
    fit.radial_fit || push!(reasons, "tank diameter plus radial clearance exceeds cabin diameter")
    fit.longitudinal_fit || push!(reasons, "tank length exceeds cylindrical cabin length")
    structural_fit || push!(reasons, "SO2 plus tank exceeds structural payload capacity")
    return (;
        feasible,
        reason = join(reasons, "; "),
        tank,
        fit,
        structural_fit,
        structural_payload_limit_lb,
    )
end

function active_configuration(configuration::Symbol)
    spec = WorkflowPaths.model_spec(configuration)
    isnothing(spec.result_stem) && error("Model $configuration does not support payload analysis.")
    return (
        model = configuration,
        title = spec.title,
        data_output = WorkflowPaths.payload_mach_file(configuration),
        show_radius = spec.show_radius,
        input_file = WorkflowPaths.require_model(configuration),
    )
end

"""Replace the radius-title marker with the radius stored in the quicksave."""
function resolve_configuration_title(config, ac::TASOPT.aircraft)
    config.show_radius || return config
    radius_in = ac.fuselage.layout.radius / in_to_m
    title = @sprintf("%s (radius %.0f in.)", config.title, radius_in)
    return merge(config, (; title))
end

function save_sweep_data_csv(
    configuration::Symbol,
    config,
    payload_lb,
    mach,
    ceiling_ft,
    cruise_cl,
    specific_excess_power_fpm,
    ceiling_fuel_lb,
    ceiling_mtow_margin_lb,
    ceiling_toc_roc_fpm,
    payload_packages;
    structural_payload_limit_lb::Float64,
    settings::Settings,
    payload_grid_mode::String,
    mach_grid_mode::String,
    data_file::String = WorkflowPaths.output_path(config.data_output),
)
    rows = NamedTuple[]
    for (j, current_mach) in pairs(mach), (i, payload) in pairs(payload_lb)
        package = payload_packages[i]
        tank = package.tank
        fit = package.fit
        push!(rows, (
            schema_version = 4,
            configuration = String(configuration),
            input_file = basename(config.input_file),
            title = String(config.title),
            payload_definition = "SO2 mass; tank mass is added to mission payload",
            payload_grid_mode,
            mach_grid_mode,
            configured_payload_min_lb = settings.payload_min_lb,
            configured_payload_max_lb = settings.payload_max_lb,
            configured_mach_min = settings.mach_min,
            configured_mach_max = settings.mach_max,
            range_nmi = settings.range_nmi,
            fuel_reserve_fraction = settings.fuel_reserve_fraction,
            service_ceiling_roc_fpm = settings.minimum_roc_fpm,
            ps_reference_altitude_ft = settings.ps_altitude_ft,
            so2_temperature_K = settings.tank.temperature_K,
            so2_liquid_density_kg_m3 = settings.tank.liquid_density_kg_m3,
            tank_ullage_fraction = settings.tank.ullage_fraction,
            tank_count = 1,
            tank_minimum_aspect_ratio = settings.tank.minimum_aspect_ratio,
            tank_pressure_margin_Pa = settings.tank.pressure_margin_Pa,
            tank_material_allowable_stress_Pa =
                settings.tank.material_allowable_stress_Pa,
            tank_weld_efficiency = settings.tank.weld_efficiency,
            tank_material_density_kg_m3 = settings.tank.material_density_kg_m3,
            tank_radial_clearance_m = settings.tank.radial_clearance_m,
            payload_lb = payload,
            so2_payload_lb = payload,
            tank_mass_lb = tank.tank_mass_lb,
            total_payload_lb = tank.total_payload_lb,
            tank_saturation_pressure_Pa = tank.saturation_pressure_Pa,
            tank_design_pressure_Pa = tank.design_pressure_Pa,
            tank_aspect_ratio = tank.aspect_ratio,
            tank_cylinder_wall_thickness_m = tank.cylinder_wall_thickness_m,
            tank_head_wall_thickness_m = tank.head_wall_thickness_m,
            tank_installed_diameter_m = tank.installed_diameter_m,
            tank_installed_length_m = tank.installed_length_m,
            tank_radial_fit = fit.radial_fit,
            tank_longitudinal_fit = fit.longitudinal_fit,
            structural_payload_fit = package.structural_fit,
            payload_package_feasible = package.feasible,
            structural_payload_limit_lb,
            mach = current_mach,
            ceiling_ft = ceiling_ft[j, i],
            cruise_cl = cruise_cl[j, i],
            specific_excess_power_fpm = specific_excess_power_fpm[j, i],
            fuel_lb_at_ceiling = ceiling_fuel_lb[j, i],
            mtow_margin_lb_at_ceiling = ceiling_mtow_margin_lb[j, i],
            toc_roc_fpm_at_ceiling = ceiling_toc_roc_fpm[j, i],
        ))
    end
    mkpath(dirname(data_file))
    CSV.write(data_file, rows)
    return data_file
end

"""Add mission 2 when the quicksave supplies only its sized design mission."""
function ensure_off_design_mission!(ac::TASOPT.aircraft)
    size(ac.parm, 2) >= 2 && return ac
    ac.parm = cat(ac.parm[:, 1], ac.parm[:, 1], dims = 2)
    ac.para = cat(ac.para[:, :, 1], ac.para[:, :, 1], dims = 3)
    ac.pare = cat(ac.pare[:, :, 1], ac.pare[:, :, 1], dims = 3)
    for hx in ac.engine.heat_exchangers
        hx.HXgas_mission = cat(hx.HXgas_mission[:, 1], hx.HXgas_mission[:, 1], dims = 2)
    end
    return ac
end

"""Configure mission 2 for one payload/Mach/altitude point."""
function configure_mission!(ac, total_payload_lb, mach, altitude_ft, settings::Settings)
    im = 2
    ac.parm[imWpay, im] = total_payload_lb * lbf_to_N
    ac.parm[imRange, im] = settings.range_nmi * nmi_to_m
    ac.para[iaMach, ipclimbn:ipdescent1, im] .= mach
    ac.para[iaalt, ipclimbn:ipcruise1, im] .= altitude_ft * ft_to_m
    ac.para[iaCL, :, im] .= ac.para[iaCL, :, 1]
    TASOPT.set_ambient_conditions!(ac, ipclimbn; im)
    TASOPT.set_ambient_conditions!(ac, ipcruise1; im)
    return ac
end

"""Return the total available thrust and aerodynamic drag at top of climb."""
function ceiling_performance(ac::TASOPT.aircraft, im::Int = 2)
    ip = ipclimbn
    total_thrust_available_lbf =
        ac.pare[ieFe, ip, im] * ac.parg[igneng] / lbf_to_N
    dynamic_pressure_Pa =
        0.5 * ac.pare[ierho0, ip, im] * ac.pare[ieu0, ip, im]^2
    drag_lbf =
        dynamic_pressure_Pa * ac.wing.layout.S * ac.para[iaCD, ip, im] / lbf_to_N
    return (; total_thrust_available_lbf, drag_lbf)
end

"""Evaluate all service-ceiling feasibility constraints at one altitude."""
function evaluate_mission(
    base::TASOPT.aircraft,
    so2_payload_lb::Float64,
    mach::Float64,
    altitude_ft::Float64;
    settings::Settings,
    minimum_roc_fpm::Float64 = settings.minimum_roc_fpm,
)
    package = payload_package(base, so2_payload_lb, settings)
    if !package.feasible
        return (; feasible = false, payload_feasible = false, reason = package.reason,
            fuel_lb = NaN, toc_roc_fpm = NaN, cruise_cl = NaN,
            mtow_margin_lb = -Inf, total_thrust_available_lbf = NaN,
            drag_lbf = NaN, cruise_range_m = 0.0, fuel_N = Inf,
            takeoff_weight_N = Inf, reached_altitude_ft = 0.0,
            mission_completed = false, engine_converged = false, package)
    end

    ac = deepcopy(base)
    ac.parg[igfreserve] = settings.fuel_reserve_fraction
    ac.pare[:, :, 2] .= ac.pare[:, :, 1]

    reached_altitude_ft = 0.0
    mission_completed = false
    engine_converged = true
    reason = ""
    bow_N = base.parg[igWMTO] - base.parg[igWfuel] - base.parg[igWpay]
    fuel_N = Inf
    wto_N = Inf
    fuel_lb = Inf
    toc_roc_fpm = 0.0
    cruise_cl = 0.0
    mtow_margin_lb = -Inf
    cruise_range_m = 0.0
    total_thrust_available_lbf = NaN
    drag_lbf = NaN
    altitude_steps_ft = sort!(unique(filter(
        altitude -> altitude <= altitude_ft,
        [settings.altitude_lower_ft, 45_000.0, 55_000.0, 60_000.0,
            62_000.0, 64_000.0, altitude_ft],
    )))

    for step_altitude_ft in altitude_steps_ft
        configure_mission!(
            ac,
            package.tank.total_payload_lb,
            mach,
            step_altitude_ft,
            settings,
        )
        weight_N = base.parg[igWMTO] * ac.para[iafracW, ipcruise1, 2]
        dynamic_pressure_Pa = 0.5 * ac.pare[ierho0, ipcruise1, 2] *
            ac.pare[ieu0, ipcruise1, 2]^2
        ac.para[iaCL, ipclimb1+1:ipdescentn-1, 2] .=
            weight_N / (dynamic_pressure_Pa * ac.wing.layout.S)

        try
            redirect_stdout(devnull) do
                redirect_stderr(devnull) do
                    TASOPT.fly_mission!(ac, 2;
                        itermax = settings.mission_itermax,
                        initializes_engine = false,
                        opt_prescribed_cruise_parameter = "altitude")
                end
            end
        catch err
            if err isa DomainError ||
               (err isa ErrorException && occursin("NaN", err.msg))
                engine_converged = false
                reason = "mission calculation outside valid domain"
                break
            end
            rethrow()
        end

        engine_converged = sum(ac.pare[ieConvFail, :, 2]) == 0.0
        if !engine_converged
            reason = "engine convergence constraint"
            break
        end
        reached_altitude_ft = step_altitude_ft
        mission_completed = step_altitude_ft == altitude_ft
        fuel_N = ac.parm[imWfuel, 2]
        wto_N = bow_N + package.tank.total_payload_lb * lbf_to_N + fuel_N
        fuel_lb = weight_lb(fuel_N)
        toc_roc_fpm = ac.para[iaROC, ipclimbn, 2]
        cruise_cl = ac.para[iaCL, ipcruise1, 2]
        mtow_margin_lb = weight_lb(base.parg[igWMTO] - wto_N)
        cruise_range_m = ac.para[iaRange, ipcruisen, 2] -
            ac.para[iaRange, ipcruise1, 2]
        performance = ceiling_performance(ac)
        total_thrust_available_lbf = performance.total_thrust_available_lbf
        drag_lbf = performance.drag_lbf
    end

    feasible = mission_completed && engine_converged &&
        all(isfinite, (fuel_lb, toc_roc_fpm, cruise_cl, mtow_margin_lb, cruise_range_m)) &&
        cruise_range_m > 0.0 &&
        fuel_lb <= weight_lb(base.parg[igWfmax]) + weight_lb(settings.weight_tolerance_N) &&
        mtow_margin_lb >= -weight_lb(settings.weight_tolerance_N) &&
        toc_roc_fpm >= minimum_roc_fpm

    return (; feasible, payload_feasible = true,
        reason = feasible ? "" : isempty(reason) ?
            "fuel, MTOW, ROC, or engine convergence constraint" : reason,
        fuel_lb, toc_roc_fpm, cruise_cl, mtow_margin_lb,
        total_thrust_available_lbf, drag_lbf,
        cruise_range_m, fuel_N, takeoff_weight_N = wto_N,
        reached_altitude_ft, mission_completed, engine_converged, package)
end

"""Evaluate specific excess power at the common fixed-altitude reference condition."""
function evaluate_specific_excess_power(base, payload_lb, mach, settings::Settings)
    result = evaluate_mission(
        base,
        payload_lb,
        mach,
        settings.ps_altitude_ft;
        settings,
        minimum_roc_fpm = -Inf,
    )
    # TASOPT reports specific excess power here as climb rate in ft/min.
    return (; feasible = result.feasible, ps_fpm = result.toc_roc_fpm, result)
end

"""Bisection-search the highest feasible cruise altitude for one grid point."""
function find_service_ceiling(base, payload_lb, mach, settings::Settings)
    ceiling_feasible(result) = settings.ceiling_requires_mission_feasibility ?
        result.feasible :
        result.payload_feasible && result.mission_completed &&
        result.engine_converged && isfinite(result.toc_roc_fpm) &&
        result.toc_roc_fpm >= settings.minimum_roc_fpm

    low = settings.altitude_lower_ft
    best = evaluate_mission(base, payload_lb, mach, low; settings)
    while !ceiling_feasible(best) && low < settings.altitude_initial_upper_ft
        low = min(low + 5_000.0, settings.altitude_initial_upper_ft)
        best = evaluate_mission(base, payload_lb, mach, low; settings)
    end
    ceiling_feasible(best) || return (; altitude_ft = NaN, result = best)

    high = max(settings.altitude_initial_upper_ft, low + 5_000.0)
    high = min(high, settings.altitude_search_limit_ft)
    high == low && return (; altitude_ft = low, result = best)
    while true
        trial = evaluate_mission(base, payload_lb, mach, high; settings)
        !ceiling_feasible(trial) && break
        low, best = high, trial
        high >= settings.altitude_search_limit_ft && return (; altitude_ft = low, result = best)
        high = min(settings.altitude_search_limit_ft, high + 5_000.0)
    end

    while high - low > settings.altitude_tolerance_ft
        altitude_ft = 0.5 * (low + high)
        trial = evaluate_mission(base, payload_lb, mach, altitude_ft; settings)
        if ceiling_feasible(trial)
            low, best = altitude_ft, trial
        else
            high = altitude_ft
        end
    end
    return (; altitude_ft = low, result = best)
end

function run(
    configuration::Symbol;
    payload_values_lb = nothing,
    mach_values = nothing,
    settings::Settings = Settings(),
)
    isnothing(payload_values_lb) && error("Provide a payload grid in run_sad_sizing.jl.")
    isnothing(mach_values) && error("Provide a Mach grid in run_sad_sizing.jl.")
    config = active_configuration(configuration)
    base = ensure_off_design_mission!(TASOPT.quickload_aircraft(config.input_file))
    config = resolve_configuration_title(config, base)
    base.is_sized[1] ||
        error("$(config.input_file) is not sized; this sweep intentionally performs no sizing.")
    base.parg[igfreserve] = settings.fuel_reserve_fraction

    structural_payload_limit_lb = weight_lb(base.parg[igWpaymax])
    payload_lb = Float64.(collect(payload_values_lb))
    mach = Float64.(collect(mach_values))
    all((settings.payload_min_lb .<= payload_lb) .& (payload_lb .<= settings.payload_max_lb)) ||
        error("payload_values_lb lies outside the configured payload limits.")
    all((settings.mach_min .<= mach) .& (mach .<= settings.mach_max)) ||
        error("mach_values lies outside the configured Mach limits.")
    payload_packages = [payload_package(base, payload, settings) for payload in payload_lb]
    tank_mass_lb = [package.tank.tank_mass_lb for package in payload_packages]
    total_payload_lb = [package.tank.total_payload_lb for package in payload_packages]
    tank_installed_diameter_m =
        [package.tank.installed_diameter_m for package in payload_packages]
    tank_installed_length_m =
        [package.tank.installed_length_m for package in payload_packages]
    payload_package_feasible = [package.feasible for package in payload_packages]

    ceiling_ft = fill(NaN, length(mach), length(payload_lb))
    cruise_cl = fill(NaN, length(mach), length(payload_lb))
    specific_excess_power_fpm = fill(NaN, length(mach), length(payload_lb))
    ceiling_fuel_lb = fill(NaN, length(mach), length(payload_lb))
    ceiling_mtow_margin_lb = fill(NaN, length(mach), length(payload_lb))
    ceiling_toc_roc_fpm = fill(NaN, length(mach), length(payload_lb))
    println("$(config.title) SO2 payload / Mach service-ceiling sweep")
    @printf("  SO2 payload: %.0f to %.0f lb\n", first(payload_lb), last(payload_lb))
    @printf("  Mach: M%.2f to M%.2f\n", first(mach), last(mach))
    @printf("  Mission: %.0f nmi, %.0f%% fuel reserve, ROC >= %.0f ft/min\n",
        settings.range_nmi, 100 * settings.fuel_reserve_fraction, settings.minimum_roc_fpm)
    @printf("  Specific-excess-power reference altitude: %.0f ft\n",
        settings.ps_altitude_ft)

    for (j, current_mach) in pairs(mach), (i, payload) in pairs(payload_lb)
        package = payload_packages[i]
        if !package.feasible
            @printf("  M%.2f, SO2 %6.0f lb -> rejected: %s\n",
                current_mach, payload, package.reason)
            continue
        end

        outcome = find_service_ceiling(base, payload, current_mach, settings)
        ceiling_ft[j, i] = outcome.altitude_ft
        if isfinite(outcome.altitude_ft)
            cruise_cl[j, i] = outcome.result.cruise_cl
            ceiling_fuel_lb[j, i] = outcome.result.fuel_lb
            ceiling_mtow_margin_lb[j, i] = outcome.result.mtow_margin_lb
            ceiling_toc_roc_fpm[j, i] = outcome.result.toc_roc_fpm
            @printf("  M%.2f, SO2 %6.0f lb + tank %6.1f lb -> %6.0f ft, CL %.3f\n",
                current_mach, payload, package.tank.tank_mass_lb,
                outcome.altitude_ft, outcome.result.cruise_cl)
        else
            @printf("  M%.2f, SO2 %6.0f lb -> infeasible at %.0f ft\n",
                current_mach, payload, settings.altitude_lower_ft)
        end

        ps_outcome = evaluate_specific_excess_power(base, payload, current_mach, settings)
        if ps_outcome.feasible
            specific_excess_power_fpm[j, i] = ps_outcome.ps_fpm
            @printf("      Ps at %.0f ft: %7.1f ft/min\n",
                settings.ps_altitude_ft, ps_outcome.ps_fpm)
        else
            @printf("      Ps at %.0f ft: infeasible\n", settings.ps_altitude_ft)
        end
    end

    data_path = save_sweep_data_csv(
        configuration,
        config,
        payload_lb,
        mach,
        ceiling_ft,
        cruise_cl,
        specific_excess_power_fpm,
        ceiling_fuel_lb,
        ceiling_mtow_margin_lb,
        ceiling_toc_roc_fpm,
        payload_packages;
        structural_payload_limit_lb,
        settings,
        payload_grid_mode = "configured",
        mach_grid_mode = "configured",
    )
    println("Saved sweep data to: $data_path")
    return (; payload_lb, so2_payload_lb = payload_lb, mach, ceiling_ft, cruise_cl,
        specific_excess_power_fpm, ceiling_fuel_lb, ceiling_mtow_margin_lb,
        ceiling_toc_roc_fpm, tank_mass_lb, total_payload_lb,
        tank_installed_diameter_m, tank_installed_length_m,
        payload_package_feasible, data_path)
end

"""Run the common grid for every presentation configuration."""
function run_all(
    configurations;
    payload_values_lb = nothing,
    mach_values = nothing,
    settings::Settings = Settings(),
)
    results = Dict{Symbol,Any}()
    for configuration in configurations
        results[configuration] = run(
            configuration;
            payload_values_lb,
            mach_values,
            settings,
        )
    end
    return results
end

end
