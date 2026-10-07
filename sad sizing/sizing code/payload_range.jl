module PayloadRangeStudy

using TASOPT
using Printf
using CSV

using ..WorkflowPaths

include(TASOPT.__TASOPTindices__)

Base.@kwdef struct StudySettings
    altitude_ft::Vector{Float64} = [35_000.0, 41_000.0]
    mach::Vector{Float64} = [0.78]
    range_tolerance_nmi::Float64 = 1.0
    range_search_limit_nmi::Float64 = 3_500.0
    weight_tolerance_N::Float64 = 5.0
    fuel_target_tolerance_N::Float64 = 250.0
    enable_payload_recovery::Bool = true
    payload_reduction_step_lb::Float64 = 50.0
    pub_ranges_nmi::Vector{Float64} = [0.0, 950.0, 1_980.0, 2_270.0]
    pub_zfw_lb::Vector{Float64} = [69_887.0, 69_887.0, 61_955.0, 47_399.0]
end

Base.@kwdef struct Settings
    model_file::String
    title::String
    altitude_ft::Float64
    mach::Float64
    data_file::String
    range_tolerance_nmi::Float64 = 1.0
    range_search_limit_nmi::Float64 = 3_500.0
    weight_tolerance_N::Float64 = 5.0
    fuel_target_tolerance_N::Float64 = 250.0
    enable_payload_recovery::Bool = true
    payload_reduction_step_lb::Float64 = 50.0
    pub_ranges_nmi::Vector{Float64} = [0.0, 950.0, 1_980.0, 2_270.0]
    pub_zfw_lb::Vector{Float64} = [69_887.0, 69_887.0, 61_955.0, 47_399.0]
end

Base.@kwdef struct MissionResult
    feasible::Bool
    range_nmi::Float64
    payload_N::Float64
    fuel_N::Float64
    WTO_N::Float64
    fuel_margin_N::Float64
    mtow_margin_N::Float64
    pfei::Float64
    message::String = ""
end

struct PayloadRangePoint
    label::String
    range_nmi::Float64
    payload_N::Float64
    fuel_N::Float64
    WTO_N::Float64
    fuel_margin_N::Float64
    mtow_margin_N::Float64
    pfei::Float64
    note::String
end

lbf(weight_N) = weight_N / lbf_to_N
fmt_number(x) = replace(@sprintf("%.0f", x), r"(?<=[0-9])(?=(?:[0-9]{3})+$)" => ",")

function active_configuration(configuration::Symbol)
    spec = WorkflowPaths.model_spec(configuration)
    isnothing(spec.result_stem) && error("Model $configuration does not support payload analysis.")
    return (
        model = configuration,
        title = spec.title,
        data_output = WorkflowPaths.payload_range_file(configuration),
        show_radius = spec.show_radius,
        input_file = WorkflowPaths.require_model(configuration),
    )
end

function resolve_configuration_title(config, ac::TASOPT.aircraft)
    config.show_radius || return config
    radius_in = ac.fuselage.layout.radius / in_to_m
    return merge(config, (; title = @sprintf("%s (radius %.0f in.)", config.title, radius_in)))
end

function empty_weight(ac::TASOPT.aircraft)
    return ac.fuselage.weight + ac.wing.weight + ac.wing.strut.weight +
           ac.htail.weight + ac.vtail.weight + ac.parg[igWeng] +
           ac.parg[igWtesys] + ac.parg[igWftank] +
           ac.landing_gear.nose_gear.weight.W + ac.landing_gear.main_gear.weight.W +
           ac.parg[igWMTO] * ac.fuselage.HPE_sys.W
end

function build_settings(
    config,
    study::StudySettings,
)
    isempty(study.altitude_ft) && error("Altitude grid cannot be empty.")
    isempty(study.mach) && error("Mach grid cannot be empty.")

    return [
        Settings(
            ; model_file = config.input_file,
            title = String(config.title),
            altitude_ft,
            mach,
            data_file = WorkflowPaths.output_path(config.data_output),
            range_tolerance_nmi = study.range_tolerance_nmi,
            range_search_limit_nmi = study.range_search_limit_nmi,
            weight_tolerance_N = study.weight_tolerance_N,
            fuel_target_tolerance_N = study.fuel_target_tolerance_N,
            enable_payload_recovery = study.enable_payload_recovery,
            payload_reduction_step_lb = study.payload_reduction_step_lb,
            pub_ranges_nmi = copy(study.pub_ranges_nmi),
            pub_zfw_lb = copy(study.pub_zfw_lb),
        )
        for altitude_ft in study.altitude_ft, mach in study.mach
    ]
end

include(joinpath(@__DIR__, "payload_range_solver.jl"))

function save_payload_range_data_csv(payload_range_results, fallback_settings::Settings)
    rows = NamedTuple[]
    for result in payload_range_results
        points = result.points
        settings = result.settings
        oew_lb = lbf(points[1].WTO_N - points[1].payload_N)
        for point in points
            push!(rows, (
                series_type = "model",
                title = settings.title,
                altitude_ft = settings.altitude_ft,
                mach = settings.mach,
                point_label = point.label,
                range_nmi = point.range_nmi,
                oew_plus_payload_lb = oew_lb + lbf(point.payload_N),
            ))
        end
    end

    settings = isempty(payload_range_results) ? fallback_settings : first(payload_range_results).settings
    for (range_nmi, oew_plus_payload_lb) in zip(settings.pub_ranges_nmi, settings.pub_zfw_lb)
        push!(rows, (
            series_type = "published_reference",
            title = settings.title,
            altitude_ft = NaN,
            mach = NaN,
            point_label = "",
            range_nmi,
            oew_plus_payload_lb,
        ))
    end
    mkpath(dirname(settings.data_file))
    CSV.write(settings.data_file, rows)
    return settings.data_file
end

function print_summary(
    points::Vector{PayloadRangePoint},
    settings::Settings,
    max_payload_requested_N::Float64,
    usable_fuel_payload_requested_N::Float64,
)
    println("Model file: $(settings.model_file)")
    @printf("Cruise altitude: %.0f ft\n", settings.altitude_ft)
    @printf("Cruise Mach: %.2f\n\n", settings.mach)

    @printf("%-12s %12s %14s %14s %14s %14s %14s %10s\n",
        "Point", "Range [nmi]", "Requested [lb]", "Payload [lb]", "Reduction [lb]", "Fuel [lb]", "WTO [lb]", "PFEI")

    for point in points
        requested_payload_N = point.label == "Usable Fuel" ? usable_fuel_payload_requested_N :
                              point.label == "Ferry" ? 0.0 : max_payload_requested_N
        pfei_text = isfinite(point.pfei) ? @sprintf("%.4f", point.pfei) : "-"
        @printf("%-12s %12.1f %14.1f %14.1f %14.1f %14.1f %14.1f %10s\n",
            point.label,
            point.range_nmi,
            lbf(requested_payload_N),
            lbf(point.payload_N),
            lbf(requested_payload_N - point.payload_N),
            lbf(point.fuel_N),
            lbf(point.WTO_N),
            pfei_text,
        )
        isempty(point.note) || println("  note: $(point.note)")
    end

    println()
end

function run(
    configuration::Symbol;
    study::StudySettings = StudySettings(),
)
    config = active_configuration(configuration)
    ac = TASOPT.quickload_aircraft(config.input_file)
    config = resolve_configuration_title(config, ac)
    settings_list = build_settings(config, study)
    ensure_off_design_mission!(ac)

    if !ac.is_sized[1]
        TASOPT.size_aircraft!(ac; iter = 50, printiter = false)
    end

    payload_range_results = NamedTuple[]
    for settings in settings_list
        try
            points, max_payload_requested_N, usable_fuel_payload_requested_N =
                build_payload_range_points(ac, settings)
            print_summary(points, settings, max_payload_requested_N, usable_fuel_payload_requested_N)
            push!(payload_range_results, (; points, settings))
        catch err
            @warn "Skipping infeasible payload-range case" altitude_ft = settings.altitude_ft mach = settings.mach reason = sprint(showerror, err)
        end
    end

    isempty(payload_range_results) && @warn "No feasible payload-range cases were found; saving the published reference only" configuration
    data_file = save_payload_range_data_csv(payload_range_results, first(settings_list))
    println("Saved payload-range data to: $data_file")
    return payload_range_results
end

end
