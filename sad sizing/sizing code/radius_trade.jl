# Run the fixed-length E175/CFM56 fuselage-radius trade study.

module RadiusTrade

using TASOPT
using Printf
using CSV

using ..CFM56Installation
using ..SO2Payload
using ..WorkflowPaths

include(TASOPT.__TASOPTindices__)

const DATA_FILE = WorkflowPaths.output_path("e175_cfm56_radius_payload_ceiling_sweep.csv")
const REDUCED_RADIUS_MODEL = :hybrid_reduced_radius_no_interior
const REDUCED_RADIUS_MODEL_FILE = WorkflowPaths.model_path(REDUCED_RADIUS_MODEL)

Base.@kwdef struct Settings
    design_so2_payload_lb::Float64 = 10_000.0
    baseline_radius_in::Float64 = 59.0
    cruise_mach::Float64 = 0.78
    range_nmi::Float64 = 500.0
    fuel_reserve_fraction::Float64 = 0.05
    payload_reference_altitude_ft::Float64 = 35_000.0
    payload_capacity_tolerance_lb::Float64 = 25.0
    reference_payload_lb::Float64 = 2_000.0
    mtow_closure_altitude_ft::Float64 = 35_000.0
    mtow_closure_tolerance_lb::Float64 = 10.0
    mtow_closure_max_iterations::Int = 12
    mtow_closure_relaxation::Float64 = 1.0
    mtow_closure_mission_iterations::Int = 35
    reference_payload_tolerance_lb::Float64 = 10.0
    service_ceiling_roc_fpm::Float64 = 100.0
    ceiling_lower_ft::Float64 = 35_000.0
    ceiling_initial_upper_ft::Float64 = 60_000.0
    ceiling_search_limit_ft::Float64 = 70_000.0
    ceiling_tolerance_ft::Float64 = 50.0
    mission_iterations::Int = 75
    tank::SO2Payload.TankSpec = SO2Payload.TankSpec()
end

include(joinpath(@__DIR__, "radius_model.jl"))

function payload_package(variant, so2_payload_lb, settings::Settings)
    spec = settings.tank
    usable_inner_diameter_m = max(
        0.0,
        2.0 * (variant.ac.fuselage.layout.radius - variant.ac.fuselage.skin.thickness),
    )
    tank = SO2Payload.size_tank(
        so2_payload_lb,
        spec;
        usable_inner_diameter_m,
    )
    fit = SO2Payload.fuselage_fit(variant.ac, tank, spec)
    structural_payload_limit_lb = weight_lb(variant.structural_payload_N)
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

function configure_ceiling_mission!(ac, total_payload_lb, altitude_ft, settings::Settings)
    im = 2
    ac.parm[imRange, im] = settings.range_nmi * nmi_to_m
    ac.parm[imWpay, im] = total_payload_lb * lbf_to_N
    ac.para[iaMach, ipclimbn:ipdescent1, im] .= settings.cruise_mach
    ac.para[iaalt, ipclimbn:ipcruise1, im] .= altitude_ft * ft_to_m
    TASOPT.set_ambient_conditions!(ac, ipclimbn; im)
    TASOPT.set_ambient_conditions!(ac, ipcruise1; im)
    return ac
end

function ceiling_trial(variant, so2_payload_lb, altitude_ft, settings::Settings)
    package = payload_package(variant, so2_payload_lb, settings)
    if !package.feasible
        return (; feasible = false, payload_feasible = false, reason = package.reason,
            roc_fpm = NaN, specific_excess_thrust = NaN, fuel_lb = NaN,
            mtow_margin_lb = -Inf, package)
    end

    ac = deepcopy(variant.ac)
    configure_ceiling_mission!(
        ac,
        package.tank.total_payload_lb,
        altitude_ft,
        settings,
    )

    try
        redirect_stdout(devnull) do
            redirect_stderr(devnull) do
                TASOPT.fly_mission!(
                    ac,
                    2;
                    itermax = settings.mission_iterations,
                    initializes_engine = true,
                    opt_prescribed_cruise_parameter = "altitude",
                )
            end
        end
    catch err
        return (; feasible = false, payload_feasible = true,
            reason = sprint(showerror, err), roc_fpm = NaN,
            specific_excess_thrust = NaN, fuel_lb = NaN,
            mtow_margin_lb = -Inf, package)
    end

    fuel_N = ac.parm[imWfuel, 2]
    wto_N = variant.bow_N + package.tank.total_payload_lb * lbf_to_N + fuel_N
    roc_fpm = ac.para[iaROC, ipclimbn, 2]
    specific_excess_thrust = sin(ac.para[iagamV, ipclimbn, 2])
    fuel_lb = weight_lb(fuel_N)
    mtow_margin_lb = weight_lb(variant.mtow_N - wto_N)
    engine_converged = sum(ac.pare[ieConvFail, :, 2]) == 0.0

    feasible = engine_converged &&
        all(isfinite, (roc_fpm, specific_excess_thrust, fuel_lb, mtow_margin_lb)) &&
        roc_fpm >= settings.service_ceiling_roc_fpm &&
        fuel_lb <= weight_lb(variant.ac.parg[igWfmax]) &&
        mtow_margin_lb >= 0.0

    reason = feasible ? "" : "ROC, fuel reserve, MTOW, or engine-convergence constraint"
    return (; feasible, payload_feasible = true, reason, roc_fpm,
        specific_excess_thrust, fuel_lb, mtow_margin_lb, package)
end

function find_max_payload_capacity(variant, settings::Settings)
    low = 0.0
    high = weight_lb(variant.structural_payload_N)

    function payload_fits(so2_payload_lb)
        trial = ceiling_trial(
            variant, so2_payload_lb, settings.payload_reference_altitude_ft, settings
        )
        return trial.payload_feasible && isfinite(trial.fuel_lb) &&
            trial.fuel_lb <= weight_lb(variant.ac.parg[igWfmax]) &&
            trial.mtow_margin_lb >= 0.0
    end

    payload_fits(low) || return NaN
    while high - low > settings.payload_capacity_tolerance_lb
        so2_payload_lb = 0.5 * (low + high)
        if payload_fits(so2_payload_lb)
            low = so2_payload_lb
        else
            high = so2_payload_lb
        end
    end
    return low
end

function find_service_ceiling(variant, so2_payload_lb, settings::Settings)
    low = settings.payload_reference_altitude_ft
    best = ceiling_trial(variant, so2_payload_lb, low, settings)
    if !best.feasible
        for altitude_ft in (low - 2_500.0):-2_500.0:settings.ceiling_lower_ft
            trial = ceiling_trial(variant, so2_payload_lb, altitude_ft, settings)
            if trial.feasible
                low, best = altitude_ft, trial
                break
            end
        end
    end
    best.feasible || return (; altitude_ft = NaN, result = best)

    high = max(settings.ceiling_initial_upper_ft, low + 5_000.0)
    while true
        trial = ceiling_trial(variant, so2_payload_lb, high, settings)
        !trial.feasible && break
        low, best = high, trial
        high >= settings.ceiling_search_limit_ft && return (; altitude_ft = low, result = best)
        high = min(settings.ceiling_search_limit_ft, high + 5_000.0)
    end

    while high - low > settings.ceiling_tolerance_ft
        altitude_ft = 0.5 * (low + high)
        trial = ceiling_trial(variant, so2_payload_lb, altitude_ft, settings)
        if trial.feasible
            low, best = altitude_ft, trial
        else
            high = altitude_ft
        end
    end
    return (; altitude_ft = low, result = best)
end

function save_sweep_data_csv(
    radius_in,
    so2_payload_lb,
    ceiling_ft,
    mtow_lb,
    so2_payload_capacity_lb,
    structural_payload_limit_lb,
    sizing_payload_lb,
    sizing_fuel_lb,
    tank_sizing,
    tank_fit,
    structural_payload_fit,
    payload_package_feasible,
    design_payload_feasible,
    design_service_ceiling_ft,
    selected_radius,
    fixed_gear_tailstrike_angle_deg,
    gear_was_redesigned,
    required_tailstrike_angle_deg;
    settings::Settings,
    radius_grid_mode::String,
    payload_grid_mode::String,
    data_file::String = DATA_FILE,
)
    rows = NamedTuple[]
    for (i, radius) in pairs(radius_in), (j, payload) in pairs(so2_payload_lb)
        tank = tank_sizing[j, i]
        fit = tank_fit[j, i]
        push!(rows, (
            schema_version = 5,
            study = "e175_cfm56_radius_payload_ceiling",
            radius_grid_mode,
            payload_grid_mode,
            payload_definition = "SO2 mass; tank mass is added to mission payload",
            design_so2_payload_lb = settings.design_so2_payload_lb,
            baseline_radius_in = settings.baseline_radius_in,
            cruise_mach = settings.cruise_mach,
            range_nmi = settings.range_nmi,
            fuel_reserve_fraction = settings.fuel_reserve_fraction,
            payload_capacity_reference_altitude_ft = settings.payload_reference_altitude_ft,
            service_ceiling_roc_fpm = settings.service_ceiling_roc_fpm,
            excess_thrust_reference_so2_payload_lb = settings.reference_payload_lb,
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
            radius_in = radius,
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
            usable_inner_fuselage_diameter_m = fit.usable_inner_diameter_m,
            required_diameter_with_clearance_m = fit.required_diameter_m,
            cabin_cylinder_length_m = fit.cabin_length_m,
            tank_radial_fit = fit.radial_fit,
            tank_longitudinal_fit = fit.longitudinal_fit,
            structural_payload_fit = structural_payload_fit[j, i],
            payload_package_feasible = payload_package_feasible[j, i],
            design_payload_feasible = design_payload_feasible[i],
            design_service_ceiling_ft = design_service_ceiling_ft[i],
            selected_radius = selected_radius[i],
            ceiling_ft = ceiling_ft[j, i],
            mtow_lb = mtow_lb[i],
            payload_capacity_lb = so2_payload_capacity_lb[i],
            so2_payload_capacity_lb = so2_payload_capacity_lb[i],
            structural_payload_limit_lb = structural_payload_limit_lb[i],
            mtow_closure_payload_lb = sizing_payload_lb[i],
            mtow_closure_fuel_lb = sizing_fuel_lb[i],
            mtow_closure_altitude_ft = settings.mtow_closure_altitude_ft,
            fixed_gear_tailstrike_angle_deg = fixed_gear_tailstrike_angle_deg[i],
            required_tailstrike_angle_deg,
            gear_was_redesigned = gear_was_redesigned[i],
        ))
    end
    mkpath(dirname(data_file))
    CSV.write(data_file, rows)
    return data_file
end

function build_at_selected_radius(
    source_model::Symbol,
    output_model::Symbol;
    settings::Settings = Settings(),
)
    selected = TASOPT.quickload_aircraft(WorkflowPaths.require_model(REDUCED_RADIUS_MODEL))
    radius_in = selected.fuselage.layout.radius / in_to_m
    base_ac = TASOPT.quickload_aircraft(WorkflowPaths.require_model(source_model))
    airframe = CFM56Installation.fixed_airframe(base_ac)
    reference_payload_N = find_reference_sizing_payload(base_ac, airframe, settings)
    variant = size_radius_variant(base_ac, airframe, radius_in, reference_payload_N, settings)
    design_package = payload_package(variant, settings.design_so2_payload_lb, settings)
    design_package.feasible || error(
        "Selected radius $(round(radius_in, digits = 1)) in cannot carry " *
        "$(round(settings.design_so2_payload_lb, digits = 0)) lb of SO2: " *
        design_package.reason,
    )
    output_file = WorkflowPaths.model_path(output_model)
    quicksave_radius_variant(variant, base_ac.engine, output_file, settings)
    @printf("Saved %s: radius %.1f in, MTOW %.0f lb -> %s\n",
        WorkflowPaths.model_title(output_model), radius_in, weight_lb(variant.mtow_N), output_file)
    return variant
end

function run(;
    radius_values_in = nothing,
    payload_values_lb = nothing,
    settings::Settings = Settings(),
)
    isnothing(radius_values_in) && error("Provide a radius grid in run_sad_sizing.jl.")
    isnothing(payload_values_lb) && error("Provide a payload grid in run_sad_sizing.jl.")
    radius_in_values = Float64.(collect(radius_values_in))
    payload_lb_values = Float64.(collect(payload_values_lb))
    isempty(radius_in_values) && error("Radius grid cannot be empty.")
    isempty(payload_lb_values) && error("Payload grid cannot be empty.")
    all(payload_lb_values .>= 0.0) || error("SO2 payload grid cannot contain negative values.")
    settings.design_so2_payload_lb >= 0.0 ||
        error("Design SO2 payload cannot be negative.")
    base_ac = TASOPT.quickload_aircraft(
        WorkflowPaths.require_model(:hybrid_no_interior)
    )
    airframe = CFM56Installation.fixed_airframe(base_ac)
    reference_sizing_payload_N = find_reference_sizing_payload(base_ac, airframe, settings)
    ceiling_ft = fill(NaN, length(payload_lb_values), length(radius_in_values))
    mtow_lb = fill(NaN, length(radius_in_values))
    so2_payload_capacity_lb = fill(NaN, length(radius_in_values))
    structural_payload_limit_lb = fill(NaN, length(radius_in_values))
    sizing_payload_lb = fill(NaN, length(radius_in_values))
    sizing_fuel_lb = fill(NaN, length(radius_in_values))
    result_size = (length(payload_lb_values), length(radius_in_values))
    tank_sizing = Matrix{SO2Payload.TankSizing}(undef, result_size)
    tank_fit = Matrix{SO2Payload.FuselageFit}(undef, result_size)
    structural_payload_fit = fill(false, result_size)
    payload_package_feasible = fill(false, result_size)
    design_payload_feasible = fill(false, length(radius_in_values))
    design_service_ceiling_ft = fill(NaN, length(radius_in_values))
    selected_radius = fill(false, length(radius_in_values))
    fixed_gear_tailstrike_angle_deg = fill(NaN, length(radius_in_values))
    gear_was_redesigned = fill(false, length(radius_in_values))
    selected_variant = nothing
    selected_radius_in = NaN
    selected_design_ceiling_ft = -Inf

    println("E175 / CFM56 fixed-length radius-payload service-ceiling sweep")
    @printf("  Radius: %.1f to %.1f in\n", first(radius_in_values), last(radius_in_values))
    @printf("  SO2 payload: %.0f to %.0f lb; design requirement: %.0f lb\n",
        first(payload_lb_values), last(payload_lb_values), settings.design_so2_payload_lb)
    @printf("  Tank: %d hemispherical-ended aluminum pressure vessel, minimum L/D %.2f, %.1f%% ullage, %.2f bar margin, %.0f MPa allowable, weld efficiency %.2f, %.3f m radial clearance\n",
        1, settings.tank.minimum_aspect_ratio,
        100 * settings.tank.ullage_fraction,
        settings.tank.pressure_margin_Pa / 1.0e5,
        settings.tank.material_allowable_stress_Pa / 1.0e6,
        settings.tank.weld_efficiency,
        settings.tank.radial_clearance_m)
    @printf("  Ceiling mission: M%.2f, %.0f nmi, ROC >= %.0f ft/min\n",
        settings.cruise_mach, settings.range_nmi, settings.service_ceiling_roc_fpm)
    @printf("  59-in sizing payload closed from input MTOW: %.0f lb\n",
        weight_lb(reference_sizing_payload_N))

    for (i, radius_in) in pairs(radius_in_values)
        variant = size_radius_variant(
            base_ac,
            airframe,
            radius_in,
            reference_sizing_payload_N,
            settings,
        )
        mtow_lb[i] = weight_lb(variant.mtow_N)
        structural_payload_limit_lb[i] = weight_lb(variant.structural_payload_N)
        sizing_payload_lb[i] = weight_lb(variant.sizing_payload_N)
        sizing_fuel_lb[i] = weight_lb(variant.sizing_fuel_N)
        fixed_gear_tailstrike_angle_deg[i] =
            rad2deg(variant.fixed_gear_tailstrike_angle_rad)
        gear_was_redesigned[i] = variant.gear_was_redesigned
        so2_payload_capacity_lb[i] = find_max_payload_capacity(variant, settings)
        design_package = payload_package(variant, settings.design_so2_payload_lb, settings)
        design_outcome = design_package.feasible ?
            find_service_ceiling(variant, settings.design_so2_payload_lb, settings) : nothing
        if design_outcome !== nothing && isfinite(design_outcome.altitude_ft)
            design_payload_feasible[i] = true
            design_service_ceiling_ft[i] = design_outcome.altitude_ft
            ceiling_improvement_ft = design_outcome.altitude_ft - selected_design_ceiling_ft
            ceiling_tied = abs(ceiling_improvement_ft) <= settings.ceiling_tolerance_ft
            lower_mtow = selected_variant === nothing || variant.mtow_N < selected_variant.mtow_N
            if selected_variant === nothing ||
                    ceiling_improvement_ft > settings.ceiling_tolerance_ft ||
                    (ceiling_tied && lower_mtow)
                selected_variant = variant
                selected_radius_in = radius_in
                selected_design_ceiling_ft = design_outcome.altitude_ft
            end
        end
        @printf("\nRadius %.1f in: MTOW %.0f lb, max SO2 %.0f lb, structural payload %.0f lb, sizing payload %.0f lb, sizing fuel %.0f lb, fuselage %.0f lb, gear %.0f lb\n",
            radius_in, mtow_lb[i], so2_payload_capacity_lb[i],
            structural_payload_limit_lb[i],
            sizing_payload_lb[i], sizing_fuel_lb[i],
            weight_lb(variant.fuselage_weight_N), weight_lb(variant.landing_gear_weight_N))
        @printf("  Design SO2 %.0f lb + tank %.1f lb = %.1f lb: D %.3f m + clearance -> %.3f m available %.3f m; L %.3f m available %.3f m; %s\n",
            settings.design_so2_payload_lb,
            design_package.tank.tank_mass_lb,
            design_package.tank.total_payload_lb,
            design_package.tank.installed_diameter_m,
            design_package.fit.required_diameter_m,
            design_package.fit.usable_inner_diameter_m,
            design_package.tank.installed_length_m,
            design_package.fit.cabin_length_m,
            design_payload_feasible[i] ?
                "accepted; ceiling $(round(design_service_ceiling_ft[i], digits = 0)) ft" :
                "rejected ($(design_package.feasible ? design_outcome.result.reason : design_package.reason))")
        @printf("  Fixed-gear tail-strike angle: %.2f deg (required %.2f deg); gear redesign: %s\n",
            fixed_gear_tailstrike_angle_deg[i],
            rad2deg(variant.ac.landing_gear.tailstrike_angle),
            gear_was_redesigned[i] ? "yes" : "no")

        for (j, payload_lb) in pairs(payload_lb_values)
            package = payload_package(variant, payload_lb, settings)
            tank_sizing[j, i] = package.tank
            tank_fit[j, i] = package.fit
            structural_payload_fit[j, i] = package.structural_fit
            payload_package_feasible[j, i] = package.feasible

            if !package.feasible
                @printf("  SO2 %6.0f lb -> rejected: %s\n", payload_lb, package.reason)
                continue
            end
            payload_lb > so2_payload_capacity_lb[i] && continue
            outcome = payload_lb == settings.design_so2_payload_lb && design_outcome !== nothing ?
                design_outcome : find_service_ceiling(variant, payload_lb, settings)
            ceiling_ft[j, i] = outcome.altitude_ft
            if isfinite(outcome.altitude_ft)
                @printf("  SO2 %6.0f lb + tank %6.1f lb -> %6.0f ft\n",
                    payload_lb, package.tank.tank_mass_lb, outcome.altitude_ft)
            else
                @printf("  SO2 %6.0f lb -> infeasible at %.0f ft\n",
                    payload_lb, settings.ceiling_lower_ft)
            end
        end
    end

    selected_variant === nothing && error(
        "No radius in the configured trade can fly the service-ceiling mission with " *
        "$(round(settings.design_so2_payload_lb, digits = 0)) lb of SO2 plus its tank.",
    )
    selected_index = findfirst(==(selected_radius_in), radius_in_values)
    selected_radius[selected_index] = true

    data_path = save_sweep_data_csv(
        radius_in_values,
        payload_lb_values,
        ceiling_ft,
        mtow_lb,
        so2_payload_capacity_lb,
        structural_payload_limit_lb,
        sizing_payload_lb,
        sizing_fuel_lb,
        tank_sizing,
        tank_fit,
        structural_payload_fit,
        payload_package_feasible,
        design_payload_feasible,
        design_service_ceiling_ft,
        selected_radius,
        fixed_gear_tailstrike_angle_deg,
        gear_was_redesigned,
        rad2deg(base_ac.landing_gear.tailstrike_angle);
        settings,
        radius_grid_mode = "configured",
        payload_grid_mode = "configured",
    )
    quicksave_radius_variant(
        selected_variant,
        base_ac.engine,
        REDUCED_RADIUS_MODEL_FILE,
        settings,
    )

    println("Saved sweep data to: $data_path")
    @printf("Saved reduced-radius model selected by highest %.0f-lb-SO2 ceiling: radius %.1f in, ceiling %.0f ft, MTOW %.0f lb -> %s\n",
        settings.design_so2_payload_lb, selected_radius_in, selected_design_ceiling_ft,
        weight_lb(selected_variant.mtow_N), REDUCED_RADIUS_MODEL_FILE)
    tank_mass_lb = map(tank -> tank.tank_mass_lb, tank_sizing)
    total_payload_lb = map(tank -> tank.total_payload_lb, tank_sizing)
    tank_installed_diameter_m = map(tank -> tank.installed_diameter_m, tank_sizing)
    tank_installed_length_m = map(tank -> tank.installed_length_m, tank_sizing)
    tank_radial_fit = map(fit -> fit.radial_fit, tank_fit)
    tank_longitudinal_fit = map(fit -> fit.longitudinal_fit, tank_fit)

    return (; radius_in = radius_in_values, payload_lb = payload_lb_values,
        so2_payload_lb = payload_lb_values, ceiling_ft, mtow_lb,
        payload_capacity_lb = so2_payload_capacity_lb, so2_payload_capacity_lb,
        structural_payload_limit_lb, sizing_payload_lb, sizing_fuel_lb,
        tank_mass_lb, total_payload_lb, tank_installed_diameter_m,
        tank_installed_length_m,
        tank_radial_fit, tank_longitudinal_fit,
        structural_payload_fit, payload_package_feasible,
        design_payload_feasible, design_service_ceiling_ft, selected_radius,
        fixed_gear_tailstrike_angle_deg,
        gear_was_redesigned, data_path,
        selected = (radius_in = selected_radius_in,
            design_service_ceiling_ft = selected_design_ceiling_ft,
            mtow_lb = weight_lb(selected_variant.mtow_N),
            filepath = REDUCED_RADIUS_MODEL_FILE))
end

end
