# Configure and run the SAD aircraft-sizing workflow.

import Pkg
import TASOPT
import Plots
Pkg.activate(normpath(joinpath(@__DIR__, "..")); io = devnull)

BUILD_MODELS = Symbol[ # Listed models are always rebuilt.
    # :e175,
    # :b737,
    # :hybrid,
    # :hybrid_no_interior,
    # :hybrid_composite,
    # :hybrid_composite_no_interior,
    # :hybrid_reduced_radius_no_interior,
    # :hybrid_reduced_radius_composite_no_interior,
    :wing_optimized,
]

PAYLOAD_MACH_MODELS = Symbol[
    # :e175,
    # :hybrid,
    # :hybrid_no_interior,
    # :hybrid_composite,
    # :hybrid_composite_no_interior,
    # :hybrid_reduced_radius_no_interior,
    # :hybrid_reduced_radius_composite_no_interior,
    # :wing_optimized,
]

STICKFIG_MODELS = Symbol[
    # :e175,
    # :b737,
    # :hybrid,
    # :hybrid_no_interior,
    # :hybrid_composite,
    # :hybrid_composite_no_interior,
    # :hybrid_reduced_radius_no_interior,
    # :hybrid_reduced_radius_composite_no_interior,
    :wing_optimized,
]

PAYLOAD_RANGE_MODELS = copy(STICKFIG_MODELS)
SUMMARY_MODELS = Symbol[
    :e175,
    :hybrid,
    :hybrid_no_interior,
    :hybrid_composite,
    :hybrid_composite_no_interior,
    :hybrid_reduced_radius_no_interior,
    :hybrid_reduced_radius_composite_no_interior,
    :wing_optimized,
]
RUN_ENGINE_DECK = false
RUN_SUMMARY = true
RUN_STICKFIGURES = true

CONFIG = (
    composite_density_factor = 0.85,
    e175_maxeval = 500,
    b737_maxeval = 500,
    hybrid_mission_iterations = 100,
    wing_optimizer = (
        maxeval = 100,
        tolerance = 1.0e-5,
    ),
    payload_mach = (
        # payload_lb = collect(0.0:1_000.0:20_000.0),
        # mach = collect(0.72:0.01:0.82),
        payload_lb = collect(2_500.0:2_500.0:10_000.0),
        mach = collect(0.76:0.02:0.80),
        range_nmi = 500.0,
        fuel_reserve_fraction = 0.05,
        minimum_roc_fpm = 100.0,
        ps_altitude_ft = 35_000.0,
        altitude_lower_ft = 35_000.0,
        altitude_initial_upper_ft = 65_000.0,
        altitude_search_limit_ft = 70_000.0,
        altitude_tolerance_ft = 50.0,
        weight_tolerance_N = 10.0,
        mission_iterations = 75,
    ),
    payload_range = (
        altitude_ft = [35_000.0, 41_000.0],
        mach = [0.78],
        range_tolerance_nmi = 1.0,
        range_search_limit_nmi = 4_500.0,
        weight_tolerance_N = 5.0,
        fuel_target_tolerance_N = 250.0,
        enable_payload_recovery = true,
        payload_reduction_step_lb = 50.0,
        published_range_nmi = [0.0, 950.0, 1_980.0, 2_270.0],
        published_zfw_lb = [69_887.0, 69_887.0, 61_955.0, 47_399.0],
    ),
    radius_trade = (
        radius_in = collect(47.0:6.0:59.0),
        payload_lb = collect(0.0:10_000.0:20_000.0),
        design_so2_payload_lb = 10_000.0,
        baseline_radius_in = 59.0,
        cruise_mach = 0.78,
        range_nmi = 500.0,
        fuel_reserve_fraction = 0.05,
        payload_reference_altitude_ft = 35_000.0,
        payload_capacity_tolerance_lb = 25.0,
        reference_payload_lb = 2_000.0,
        mtow_closure_altitude_ft = 35_000.0,
        mtow_closure_tolerance_lb = 10.0,
        mtow_closure_max_iterations = 12,
        mtow_closure_relaxation = 1.0,
        mtow_closure_mission_iterations = 35,
        reference_payload_tolerance_lb = 10.0,
        service_ceiling_roc_fpm = 100.0,
        ceiling_lower_ft = 35_000.0,
        ceiling_initial_upper_ft = 60_000.0,
        ceiling_search_limit_ft = 70_000.0,
        ceiling_tolerance_ft = 50.0,
        mission_iterations = 75,
        tank = (
            temperature_K = 303.15,
            liquid_density_kg_m3 = 1_350.20,
            ullage_fraction = 0.04,
            minimum_aspect_ratio = 4.0,
            pressure_margin_Pa = 5.0e5,
            material_allowable_stress_Pa = 1.0e8,
            weld_efficiency = 0.85,
            material_density_kg_m3 = 2_700.0,
            radial_clearance_m = 0.1000,
        ),
    ),
    engine_deck = (
        altitude_step_ft = 5_000.0,
        mach_step = 0.10,
        throttle_step = 0.05,
    ),
)

const CODE_DIR = joinpath(@__DIR__, "sizing code")
include(joinpath(CODE_DIR, "workflow_paths.jl"))
include(joinpath(CODE_DIR, "so2_payload.jl"))
include(joinpath(CODE_DIR, "optimization_common.jl"))
include(joinpath(CODE_DIR, "e175_sizing.jl"))
include(joinpath(CODE_DIR, "b737_sizing.jl"))
include(joinpath(CODE_DIR, "airframe_changes.jl"))
include(joinpath(CODE_DIR, "cfm56_installation.jl"))
include(joinpath(CODE_DIR, "payload_mach_sweep.jl"))
include(joinpath(CODE_DIR, "wing_optimizer.jl"))
include(joinpath(CODE_DIR, "radius_trade.jl"))
include(joinpath(CODE_DIR, "payload_range.jl"))
include(joinpath(CODE_DIR, "engine_deck.jl"))
include(joinpath(CODE_DIR, "export_summary.jl"))

function main()
    WorkflowPaths.ensure_directories()

    tank = SO2Payload.TankSpec(; CONFIG.radius_trade.tank...)

    :e175 in BUILD_MODELS && E175Sizing.run(maxeval = CONFIG.e175_maxeval)
    :b737 in BUILD_MODELS && B737Sizing.run(maxeval = CONFIG.b737_maxeval)

    recipes = (
        (name = :hybrid, composite = false, no_interior = false),
        (name = :hybrid_no_interior, composite = false, no_interior = true),
        (name = :hybrid_composite, composite = true, no_interior = false),
        (name = :hybrid_composite_no_interior, composite = true, no_interior = true),
    )
    for recipe in recipes
        recipe.name in BUILD_MODELS || continue
        CFM56Installation.run(
            recipe.name;
            composite = recipe.composite,
            no_interior = recipe.no_interior,
            density_factor = CONFIG.composite_density_factor,
            itermax = CONFIG.hybrid_mission_iterations,
        )
    end

    reduced_radius_models = (
        :hybrid_reduced_radius_no_interior,
        :hybrid_reduced_radius_composite_no_interior,
    )
    if any(model -> model in BUILD_MODELS, reduced_radius_models)
        radius_settings = RadiusTrade.Settings(
            design_so2_payload_lb = CONFIG.radius_trade.design_so2_payload_lb,
            baseline_radius_in = CONFIG.radius_trade.baseline_radius_in,
            cruise_mach = CONFIG.radius_trade.cruise_mach,
            range_nmi = CONFIG.radius_trade.range_nmi,
            fuel_reserve_fraction = CONFIG.radius_trade.fuel_reserve_fraction,
            payload_reference_altitude_ft = CONFIG.radius_trade.payload_reference_altitude_ft,
            payload_capacity_tolerance_lb = CONFIG.radius_trade.payload_capacity_tolerance_lb,
            reference_payload_lb = CONFIG.radius_trade.reference_payload_lb,
            mtow_closure_altitude_ft = CONFIG.radius_trade.mtow_closure_altitude_ft,
            mtow_closure_tolerance_lb = CONFIG.radius_trade.mtow_closure_tolerance_lb,
            mtow_closure_max_iterations = CONFIG.radius_trade.mtow_closure_max_iterations,
            mtow_closure_relaxation = CONFIG.radius_trade.mtow_closure_relaxation,
            mtow_closure_mission_iterations = CONFIG.radius_trade.mtow_closure_mission_iterations,
            reference_payload_tolerance_lb = CONFIG.radius_trade.reference_payload_tolerance_lb,
            service_ceiling_roc_fpm = CONFIG.radius_trade.service_ceiling_roc_fpm,
            ceiling_lower_ft = CONFIG.radius_trade.ceiling_lower_ft,
            ceiling_initial_upper_ft = CONFIG.radius_trade.ceiling_initial_upper_ft,
            ceiling_search_limit_ft = CONFIG.radius_trade.ceiling_search_limit_ft,
            ceiling_tolerance_ft = CONFIG.radius_trade.ceiling_tolerance_ft,
            mission_iterations = CONFIG.radius_trade.mission_iterations,
            tank = tank,
        )
        if :hybrid_reduced_radius_no_interior in BUILD_MODELS
            RadiusTrade.run(
                radius_values_in = CONFIG.radius_trade.radius_in,
                payload_values_lb = CONFIG.radius_trade.payload_lb,
                settings = radius_settings,
            )
        end
        if :hybrid_reduced_radius_composite_no_interior in BUILD_MODELS
            RadiusTrade.build_at_selected_radius(
                :hybrid_composite_no_interior,
                :hybrid_reduced_radius_composite_no_interior;
                settings = radius_settings,
            )
        end
    end

    if :wing_optimized in BUILD_MODELS
        WingOptimizer.run(
            maxeval = CONFIG.wing_optimizer.maxeval,
            tolerance = CONFIG.wing_optimizer.tolerance,
            tank = tank,
        )
    end

    if RUN_STICKFIGURES
        for model in STICKFIG_MODELS
            aircraft = TASOPT.quickload_aircraft(WorkflowPaths.require_model(model))
            figure = TASOPT.stickfig(aircraft)
            model_stem = splitext(WorkflowPaths.model_spec(model).file)[1]
            image_file = WorkflowPaths.image_path("$(model_stem)_stickfig.png")
            Plots.savefig(figure, image_file)
            println("Saved stick figure to: $image_file")
        end
    end

    payload_values_lb = copy(CONFIG.payload_mach.payload_lb)
    mach_values = copy(CONFIG.payload_mach.mach)
    if :wing_optimized in PAYLOAD_MACH_MODELS
        reference = WingOptimizer.reference_ceiling_point()
        payload_values_lb = sort!(unique([payload_values_lb; reference.payload_lb]))
        mach_values = sort!(unique([mach_values; reference.mach]))
    end
    payload_mach_settings = PayloadMachSweep.Settings(
        payload_min_lb = minimum(payload_values_lb),
        payload_max_lb = maximum(payload_values_lb),
        mach_min = minimum(mach_values),
        mach_max = maximum(mach_values),
        range_nmi = CONFIG.payload_mach.range_nmi,
        fuel_reserve_fraction = CONFIG.payload_mach.fuel_reserve_fraction,
        minimum_roc_fpm = CONFIG.payload_mach.minimum_roc_fpm,
        ps_altitude_ft = CONFIG.payload_mach.ps_altitude_ft,
        altitude_lower_ft = CONFIG.payload_mach.altitude_lower_ft,
        altitude_initial_upper_ft = CONFIG.payload_mach.altitude_initial_upper_ft,
        altitude_search_limit_ft = CONFIG.payload_mach.altitude_search_limit_ft,
        altitude_tolerance_ft = CONFIG.payload_mach.altitude_tolerance_ft,
        weight_tolerance_N = CONFIG.payload_mach.weight_tolerance_N,
        mission_itermax = CONFIG.payload_mach.mission_iterations,
        tank = tank,
    )
    PayloadMachSweep.run_all(
        PAYLOAD_MACH_MODELS;
        payload_values_lb = payload_values_lb,
        mach_values = mach_values,
        settings = payload_mach_settings,
    )

    payload_range_settings = PayloadRangeStudy.StudySettings(
        altitude_ft = copy(CONFIG.payload_range.altitude_ft),
        mach = copy(CONFIG.payload_range.mach),
        range_tolerance_nmi = CONFIG.payload_range.range_tolerance_nmi,
        range_search_limit_nmi = CONFIG.payload_range.range_search_limit_nmi,
        weight_tolerance_N = CONFIG.payload_range.weight_tolerance_N,
        fuel_target_tolerance_N = CONFIG.payload_range.fuel_target_tolerance_N,
        enable_payload_recovery = CONFIG.payload_range.enable_payload_recovery,
        payload_reduction_step_lb = CONFIG.payload_range.payload_reduction_step_lb,
        pub_ranges_nmi = copy(CONFIG.payload_range.published_range_nmi),
        pub_zfw_lb = copy(CONFIG.payload_range.published_zfw_lb),
    )
    for model in PAYLOAD_RANGE_MODELS
        PayloadRangeStudy.run(model; study = payload_range_settings)
    end

    if RUN_ENGINE_DECK
        GenerateE175CFM56EngineDeck.generate_deck(
            GenerateE175CFM56EngineDeck.OUTPUT_FILE;
            altitude_step_ft = CONFIG.engine_deck.altitude_step_ft,
            mach_step = CONFIG.engine_deck.mach_step,
            throttle_step = CONFIG.engine_deck.throttle_step,
        )
    end

    RUN_SUMMARY && !isempty(SUMMARY_MODELS) && ExportSummary.run(
        SUMMARY_MODELS;
        composite_density_factor = CONFIG.composite_density_factor,
        tank = tank,
    )
    return nothing
end

main()
