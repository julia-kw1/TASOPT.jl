module WorkflowPaths

const ROOT = normpath(joinpath(@__DIR__, ".."))
const MODELS = joinpath(ROOT, "models")
const OUTPUTS = joinpath(ROOT, "results", "outputs")
const IMAGES = joinpath(ROOT, "results", "images")

const MODEL_SPECS = Dict(
    :e175_input => (
        file = "e175_input.toml", title = "E175 input", result_stem = nothing,
        composite = false, show_radius = false,
    ),
    :b737_input => (
        file = "b737_input.toml", title = "B737 input", result_stem = nothing,
        composite = false, show_radius = false,
    ),
    :e175 => (
        file = "e175.toml", title = "E175", result_stem = "e175",
        composite = false, show_radius = false,
    ),
    :b737 => (
        file = "b737_cfm56_donor.toml", title = "B737 CFM56 donor", result_stem = nothing,
        composite = false, show_radius = false,
    ),
    :hybrid => (
        file = "e175_cfm56.toml", title = "E175 with CFM56", result_stem = "e175_cfm56",
        composite = false, show_radius = false,
    ),
    :hybrid_no_interior => (
        file = "e175_cfm56_no_interior.toml", title = "E175 with CFM56, Without Interior",
        result_stem = "e175_cfm56_no_interior", composite = false, show_radius = false,
    ),
    :hybrid_composite => (
        file = "e175_cfm56_composite.toml", title = "E175 with CFM56, Composites",
        result_stem = "e175_cfm56_composite", composite = true, show_radius = false,
    ),
    :hybrid_composite_no_interior => (
        file = "e175_cfm56_composite_no_interior.toml",
        title = "E175 with CFM56, Composites, Without Interior",
        result_stem = "e175_cfm56_composite_no_interior", composite = true,
        show_radius = false,
    ),
    :hybrid_reduced_radius_no_interior => (
        file = "e175_cfm56_reduced_radius_no_interior.toml",
        title = "E175 with CFM56, Without Interior",
        result_stem = "e175_cfm56_reduced_radius_no_interior", composite = false,
        show_radius = true,
    ),
    :hybrid_reduced_radius_composite_no_interior => (
        file = "e175_cfm56_reduced_radius_composite_no_interior.toml",
        title = "E175 with CFM56, Composites, Without Interior",
        result_stem = "e175_cfm56_reduced_radius_composite_no_interior",
        composite = true, show_radius = true,
    ),
    :wing_optimized => (
        file = "e175_cfm56_reduced_radius_composite_no_interior_wing_optimized.toml",
        title = "E175 with CFM56, Optimized Wing and Tail",
        result_stem = "e175_cfm56_reduced_radius_composite_no_interior_wing_optimized",
        composite = true, show_radius = true,
    ),
)

function model_spec(name::Symbol)
    haskey(MODEL_SPECS, name) || error("Unknown model: $name")
    return MODEL_SPECS[name]
end

model_path(name::Symbol) = joinpath(MODELS, model_spec(name).file)
model_title(name::Symbol) = model_spec(name).title
payload_mach_file(name::Symbol) = "$(model_spec(name).result_stem)_payload_mach_sweep.csv"
payload_range_file(name::Symbol) = "$(model_spec(name).result_stem)_payload_range_comparison.csv"
output_path(filename::AbstractString) = joinpath(OUTPUTS, filename)
image_path(filename::AbstractString) = joinpath(IMAGES, filename)

function require_model(name::Symbol)
    path = model_path(name)
    isfile(path) || error("Missing model: $path")
    return path
end

function ensure_directories()
    mkpath(MODELS)
    mkpath(OUTPUTS)
    mkpath(IMAGES)
    return nothing
end

end
