"""SO2 payload-tank sizing and fuselage-fit utilities for the SAD workflow."""

module SO2Payload

const LB_TO_KG = 0.45359237
const KG_TO_LB = 1.0 / LB_TO_KG

Base.@kwdef struct TankSpec
    temperature_K::Float64 = 303.15
    liquid_density_kg_m3::Float64 = 1_350.20
    ullage_fraction::Float64 = 0.04
    minimum_aspect_ratio::Float64 = 4.0
    pressure_margin_Pa::Float64 = 5.0e5
    material_allowable_stress_Pa::Float64 = 1.0e8
    weld_efficiency::Float64 = 0.85
    material_density_kg_m3::Float64 = 2_700.0
    radial_clearance_m::Float64 = 0.1000
end

Base.@kwdef struct TankSizing
    so2_mass_lb::Float64
    saturation_pressure_Pa::Float64
    design_pressure_Pa::Float64
    aspect_ratio::Float64
    installed_diameter_m::Float64
    installed_length_m::Float64
    cylinder_wall_thickness_m::Float64
    head_wall_thickness_m::Float64
    tank_mass_lb::Float64
    total_payload_lb::Float64
end

struct FuselageFit
    feasible::Bool
    radial_fit::Bool
    longitudinal_fit::Bool
    usable_inner_diameter_m::Float64
    required_diameter_m::Float64
    cabin_length_m::Float64
end

"""
Size one longitudinal SO2 pressure vessel with hemispherical ends.

SO2 saturation pressure is evaluated with the NIST Antoine correlation at the
specified temperature, and the design pressure adds `pressure_margin_Pa`. The
cylindrical barrel and hemispherical heads are sized independently from membrane
stress. If `usable_inner_diameter_m` is finite, the aspect ratio is increased
from `minimum_aspect_ratio` just enough to satisfy radial clearance. Reported
mass is bare pressure-shell mass; no installation allowance is applied.
"""
function size_tank(
    so2_mass_lb::Real,
    spec::TankSpec = TankSpec();
    usable_inner_diameter_m::Real = Inf,
)
    so2_mass_lb = Float64(so2_mass_lb)
    so2_mass_kg = so2_mass_lb * LB_TO_KG
    liquid_volume_m3 = so2_mass_kg / spec.liquid_density_kg_m3
    tank_internal_volume_m3 = liquid_volume_m3 / (1.0 - spec.ullage_fraction)

    # NIST Chemistry WebBook Antoine equation: log10(P/bar) = A - B/(T + C).
    saturation_pressure_Pa =
        1.0e5 * 10.0^(4.37798 - 966.575 / (spec.temperature_K - 42.071))
    design_pressure_Pa = saturation_pressure_Pa + spec.pressure_margin_Pa
    cylinder_denominator_Pa =
        spec.material_allowable_stress_Pa * spec.weld_efficiency -
        0.6 * design_pressure_Pa
    head_denominator_Pa =
        2.0 * spec.material_allowable_stress_Pa * spec.weld_efficiency -
        0.2 * design_pressure_Pa
    aspect_ratio = spec.minimum_aspect_ratio
    if isfinite(usable_inner_diameter_m)
        maximum_outer_diameter_m =
            Float64(usable_inner_diameter_m) - 2.0 * spec.radial_clearance_m
        maximum_internal_diameter_m = maximum_outer_diameter_m /
            (1.0 + design_pressure_Pa / cylinder_denominator_Pa)
        aspect_ratio = max(
            aspect_ratio,
            (1.0 + 12.0 * tank_internal_volume_m3 /
                (pi * maximum_internal_diameter_m^3)) / 3.0,
        )
    end

    # V = pi*D^2*(L-D)/4 + pi*D^3/6, with L = aspect_ratio*D.
    internal_diameter_m = cbrt(
        12.0 * tank_internal_volume_m3 /
        (pi * (3.0 * aspect_ratio - 1.0)),
    )
    internal_radius_m = 0.5 * internal_diameter_m
    internal_length_m = aspect_ratio * internal_diameter_m
    barrel_length_m = internal_length_m - internal_diameter_m
    cylinder_wall_thickness_m =
        design_pressure_Pa * internal_radius_m / cylinder_denominator_Pa
    head_wall_thickness_m =
        design_pressure_Pa * internal_radius_m / head_denominator_Pa
    installed_diameter_m = internal_diameter_m + 2.0 * cylinder_wall_thickness_m
    installed_length_m = internal_length_m + 2.0 * head_wall_thickness_m
    shell_volume_m3 =
        pi * ((internal_radius_m + cylinder_wall_thickness_m)^2 - internal_radius_m^2) *
        barrel_length_m +
        (4.0 / 3.0) * pi *
        ((internal_radius_m + head_wall_thickness_m)^3 - internal_radius_m^3)
    tank_mass_kg = shell_volume_m3 * spec.material_density_kg_m3
    tank_mass_lb = tank_mass_kg * KG_TO_LB

    return TankSizing(;
        so2_mass_lb,
        saturation_pressure_Pa,
        design_pressure_Pa,
        aspect_ratio,
        installed_diameter_m,
        installed_length_m,
        cylinder_wall_thickness_m,
        head_wall_thickness_m,
        tank_mass_lb,
        total_payload_lb = so2_mass_lb + tank_mass_lb,
    )
end

"""Check one longitudinal tank envelope against the empty cylindrical cabin."""
function fuselage_fit(ac, tank::TankSizing, spec::TankSpec = TankSpec())
    layout = ac.fuselage.layout
    skin_thickness_m = max(0.0, ac.fuselage.skin.thickness)
    usable_inner_diameter_m = max(0.0, 2.0 * (layout.radius - skin_thickness_m))
    required_diameter_m = iszero(tank.installed_diameter_m) ? 0.0 :
        tank.installed_diameter_m + 2.0 * spec.radial_clearance_m
    cabin_length_m = layout.l_cabin_cylinder
    if cabin_length_m <= 0.0
        cabin_length_m = layout.x_end_cylinder - layout.x_start_cylinder
    end

    radial_fit = required_diameter_m <= usable_inner_diameter_m
    longitudinal_fit = tank.installed_length_m <= cabin_length_m
    return FuselageFit(
        radial_fit && longitudinal_fit,
        radial_fit,
        longitudinal_fit,
        usable_inner_diameter_m,
        required_diameter_m,
        cabin_length_m,
    )
end

end
