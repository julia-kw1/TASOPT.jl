module ExportSummary

using TASOPT
using CSV

using ..WorkflowPaths
using ..PayloadMachSweep
using ..SO2Payload

include(TASOPT.__TASOPTindices__)

const OUTPUT_FILE = WorkflowPaths.output_path("e175_cfm56_configuration_summary.csv")

weight_lb(weight_N) = weight_N / lbf_to_N
finite_or_missing(value) = isfinite(value) ? value : missing

function finite_rows(rows, column; minimum_payload_lb = -Inf)
    return filter(rows) do row
        value = getproperty(row, column)
        !ismissing(value) && isfinite(value) && row.payload_lb >= minimum_payload_lb
    end
end

function maximum_row(rows, column; minimum_payload_lb = -Inf)
    candidates = finite_rows(rows, column; minimum_payload_lb)
    return candidates[argmax(getproperty.(candidates, column))]
end

function weight_breakdown(ac)
    mtow_lb = weight_lb(ac.parg[igWMTO])
    bow_lb = weight_lb(ac.parg[igWMTO] - ac.parg[igWfuel] - ac.parg[igWpay])

    fuselage_structure_lb = weight_lb(
        ac.fuselage.shell.weight.W +
        ac.fuselage.cone.weight.W +
        ac.fuselage.floor.weight.W +
        ac.fuselage.window.W +
        ac.fuselage.insulation.W +
        ac.fuselage.bendingmaterial_h.weight.W +
        ac.fuselage.bendingmaterial_v.weight.W,
    )

    total_engine_lb = weight_lb(ac.parg[igWeng])
    nacelle_lb = weight_lb(ac.parg[igWnace])
    bare_engine_lb =
        (total_engine_lb / (1.0 + ac.parg[igfpylon]) - nacelle_lb) /
        (1.0 + ac.parg[igfeadd])
    engine_accessories_lb = ac.parg[igfeadd] * bare_engine_lb
    pylon_lb = ac.parg[igfpylon] * (bare_engine_lb + engine_accessories_lb + nacelle_lb)

    landing_gear_lb = weight_lb(
        ac.landing_gear.nose_gear.weight.W + ac.landing_gear.main_gear.weight.W,
    )
    electrical_system_lb = ac.fuselage.HPE_sys.W * mtow_lb
    apu_lb = weight_lb(ac.fuselage.APU.W)
    systems_installations_lb = weight_lb(
        ac.fuselage.seat.W +
        ac.fuselage.fixed.W +
        ac.fuselage.added_payload.W +
        ac.parg[igWtesys] +
        ac.parg[igWftank] +
        ac.parg[igWinsftank],
    ) + engine_accessories_lb + pylon_lb

    groups = (
        wing_lb = weight_lb(ac.wing.weight),
        horizontal_tail_lb = weight_lb(ac.htail.weight),
        vertical_tail_lb = weight_lb(ac.vtail.weight),
        fuselage_structure_lb,
        bare_engine_lb,
        nacelle_lb,
        landing_gear_lb,
        electrical_system_lb,
        apu_lb,
        systems_installations_lb,
    )
    grouped_weight_lb = sum(groups)
    other_unallocated_weight_lb = bow_lb - grouped_weight_lb
    breakdown_sum_lb = grouped_weight_lb + other_unallocated_weight_lb

    return (; mtow_lb, bow_lb, groups..., other_unallocated_weight_lb, breakdown_sum_lb,
        breakdown_residual_lb = bow_lb - breakdown_sum_lb)
end

function engine_point_metrics(ac, ip)
    lpc_pressure_ratio = ac.pare[iepilc, ip, 1]
    hpc_pressure_ratio = ac.pare[iepihc, ip, 1]
    return (
        bypass_ratio = ac.pare[ieBPR, ip, 1],
        overall_pressure_ratio = lpc_pressure_ratio * hpc_pressure_ratio,
        thrust_lbf_per_engine = weight_lb(ac.pare[ieFe, ip, 1]),
        tsfc_lb_lbf_hr = ac.pare[ieTSFC, ip, 1] * 3600.0,
    )
end

function ceiling_point_metrics(base, row, tank)
    settings = PayloadMachSweep.Settings(
        range_nmi = row.range_nmi,
        fuel_reserve_fraction = row.fuel_reserve_fraction,
        minimum_roc_fpm = row.service_ceiling_roc_fpm,
        ps_altitude_ft = row.ps_reference_altitude_ft,
        tank = tank,
    )
    ac = PayloadMachSweep.ensure_off_design_mission!(deepcopy(base))
    ac.parg[igfreserve] = settings.fuel_reserve_fraction
    result = PayloadMachSweep.evaluate_mission(
        ac,
        Float64(row.payload_lb),
        Float64(row.mach),
        Float64(row.ceiling_ft);
        settings,
    )
    return (
        total_thrust_available_lbf =
            finite_or_missing(result.total_thrust_available_lbf),
        drag_lbf = finite_or_missing(result.drag_lbf),
        ps_fpm = finite_or_missing(row.toc_roc_fpm_at_ceiling),
    )
end

function summary_row(configuration::Symbol, composite_density_factor, tank)
    spec = WorkflowPaths.model_spec(configuration)
    sweep_file = WorkflowPaths.payload_mach_file(configuration)
    quicksave_path = WorkflowPaths.require_model(configuration)
    sweep_path = WorkflowPaths.output_path(sweep_file)
    ac = TASOPT.quickload_aircraft(quicksave_path)
    redirect_stdout(devnull) do
        TASOPT.weight_buildup(ac)
    end

    rows = collect(CSV.File(sweep_path))
    absolute_ceiling = maximum_row(rows, :ceiling_ft)
    loaded_ceiling = maximum_row(rows, :ceiling_ft; minimum_payload_lb = 2_000.0)
    max_ps = maximum_row(rows, :specific_excess_power_fpm)
    absolute_ceiling_performance = ceiling_point_metrics(ac, absolute_ceiling, tank)
    loaded_ceiling_performance = ceiling_point_metrics(ac, loaded_ceiling, tank)
    weights = weight_breakdown(ac)
    inherited_payload_limit_lb = weight_lb(ac.parg[igWpaymax])
    payload_at_zero_fuel_lb = min(
        inherited_payload_limit_lb,
        max(0.0, weights.mtow_lb - weights.bow_lb),
    )
    payload_at_design_fuel_lb = min(
        inherited_payload_limit_lb,
        max(0.0, weights.mtow_lb - weights.bow_lb - weight_lb(ac.parg[igWfuel])),
    )
    payload_at_full_fuel_capacity_lb = min(
        inherited_payload_limit_lb,
        max(0.0, weights.mtow_lb - weights.bow_lb - weight_lb(ac.parg[igWfmax])),
    )
    cruise_engine = engine_point_metrics(ac, ipcruise1)
    takeoff_engine = engine_point_metrics(ac, ipstatic)
    feasible_rows = finite_rows(rows, :ceiling_ft)

    return (
        configuration_id = spec.result_stem,
        configuration = spec.title,
        input_quicksave = basename(quicksave_path),
        sweep_source = sweep_file,
        composite_density_factor = spec.composite ? composite_density_factor : missing,
        mtow_lb = weights.mtow_lb,
        bow_lb = weights.bow_lb,
        wing_weight_lb = weights.wing_lb,
        horizontal_tail_weight_lb = weights.horizontal_tail_lb,
        vertical_tail_weight_lb = weights.vertical_tail_lb,
        fuselage_structure_weight_lb = weights.fuselage_structure_lb,
        total_bare_engine_weight_lb = weights.bare_engine_lb,
        nacelle_weight_lb = weights.nacelle_lb,
        landing_gear_weight_lb = weights.landing_gear_lb,
        electrical_system_weight_lb = weights.electrical_system_lb,
        apu_weight_lb = weights.apu_lb,
        systems_installations_weight_lb = weights.systems_installations_lb,
        other_unallocated_weight_lb = weights.other_unallocated_weight_lb,
        weight_breakdown_sum_lb = weights.breakdown_sum_lb,
        weight_breakdown_residual_lb = weights.breakdown_residual_lb,
        saved_reference_mission_payload_lb = weight_lb(ac.parg[igWpay]),
        inherited_airframe_payload_limit_lb = inherited_payload_limit_lb,
        mtow_limited_payload_at_zero_fuel_lb = payload_at_zero_fuel_lb,
        mtow_limited_payload_at_design_fuel_lb = payload_at_design_fuel_lb,
        mtow_limited_payload_at_full_fuel_capacity_lb =
            payload_at_full_fuel_capacity_lb,
        sampled_max_feasible_payload_lb = maximum(getproperty.(feasible_rows, :payload_lb)),
        usable_fuel_capacity_lb = weight_lb(ac.parg[igWfmax]),
        saved_reference_mission_fuel_lb = weight_lb(ac.parg[igWfuel]),
        saved_reference_mission_range_nmi = ac.parm[imRange, 1] / nmi_to_m,
        pfei = finite_or_missing(ac.parm[imPFEI, 1]),
        absolute_service_ceiling_ft = absolute_ceiling.ceiling_ft,
        absolute_ceiling_mach = absolute_ceiling.mach,
        absolute_ceiling_payload_lb = absolute_ceiling.payload_lb,
        absolute_ceiling_cruise_cl = absolute_ceiling.cruise_cl,
        absolute_ceiling_fuel_lb = absolute_ceiling.fuel_lb_at_ceiling,
        absolute_ceiling_total_thrust_available_lbf =
            absolute_ceiling_performance.total_thrust_available_lbf,
        absolute_ceiling_drag_lbf = absolute_ceiling_performance.drag_lbf,
        absolute_ceiling_ps_fpm = absolute_ceiling_performance.ps_fpm,
        max_loaded_service_ceiling_ft = loaded_ceiling.ceiling_ft,
        max_loaded_ceiling_mach = loaded_ceiling.mach,
        max_loaded_ceiling_payload_lb = loaded_ceiling.payload_lb,
        max_loaded_ceiling_cruise_cl = loaded_ceiling.cruise_cl,
        max_loaded_ceiling_fuel_lb = loaded_ceiling.fuel_lb_at_ceiling,
        max_loaded_ceiling_total_thrust_available_lbf =
            loaded_ceiling_performance.total_thrust_available_lbf,
        max_loaded_ceiling_drag_lbf = loaded_ceiling_performance.drag_lbf,
        max_loaded_ceiling_ps_fpm = loaded_ceiling_performance.ps_fpm,
        max_specific_excess_power_fpm = max_ps.specific_excess_power_fpm,
        max_ps_mach = max_ps.mach,
        max_ps_payload_lb = max_ps.payload_lb,
        ps_reference_altitude_ft = max_ps.ps_reference_altitude_ft,
        service_ceiling_roc_threshold_fpm = absolute_ceiling.service_ceiling_roc_fpm,
        sweep_range_nmi = absolute_ceiling.range_nmi,
        sweep_fuel_reserve_fraction = absolute_ceiling.fuel_reserve_fraction,
        fuselage_radius_in = ac.fuselage.layout.radius / in_to_m,
        fuselage_length_ft = ac.fuselage.layout.x_end / ft_to_m,
        wing_area_ft2 = ac.wing.layout.S / ft_to_m^2,
        wing_span_ft = ac.wing.layout.span / ft_to_m,
        wing_aspect_ratio = ac.wing.layout.span^2 / ac.wing.layout.S,
        engine_count = Int(round(ac.parg[igneng])),
        fan_diameter_in = ac.parg[igdfan] / in_to_m,
        cruise_bypass_ratio = cruise_engine.bypass_ratio,
        cruise_overall_pressure_ratio = cruise_engine.overall_pressure_ratio,
        cruise_thrust_lbf_per_engine = cruise_engine.thrust_lbf_per_engine,
        cruise_tsfc_lb_lbf_hr = cruise_engine.tsfc_lb_lbf_hr,
        takeoff_bypass_ratio = takeoff_engine.bypass_ratio,
        takeoff_overall_pressure_ratio = takeoff_engine.overall_pressure_ratio,
        takeoff_thrust_lbf_per_engine = takeoff_engine.thrust_lbf_per_engine,
        takeoff_tsfc_lb_lbf_hr = takeoff_engine.tsfc_lb_lbf_hr,
    )
end

function run(
    configurations;
    composite_density_factor::Float64,
    tank::SO2Payload.TankSpec = SO2Payload.TankSpec(),
)
    rows = [
        summary_row(configuration, composite_density_factor, tank)
        for configuration in configurations
    ]
    CSV.write(OUTPUT_FILE, rows)
    println("Saved E175/CFM56 configuration summary to: $OUTPUT_FILE")

    return rows
end

end
