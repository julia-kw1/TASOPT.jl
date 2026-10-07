module AirframeChanges

using TASOPT

include(TASOPT.__TASOPTindices__)

const RETAINED_COCKPIT_AND_PILOTS_LB = 3_000.0
const CABIN_PRESSURE_ALTITUDE_FT = 8_000.0

function structural_materials(ac)
    return (
        ac.wing.inboard.caps.material,
        ac.wing.outboard.caps.material,
        ac.wing.inboard.webs.material,
        ac.wing.outboard.webs.material,
        ac.wing.strut.material,
        ac.htail.inboard.caps.material,
        ac.htail.outboard.caps.material,
        ac.htail.inboard.webs.material,
        ac.htail.outboard.webs.material,
        ac.htail.strut.material,
        ac.vtail.inboard.caps.material,
        ac.vtail.outboard.caps.material,
        ac.vtail.inboard.webs.material,
        ac.vtail.outboard.webs.material,
        ac.vtail.strut.material,
        ac.fuselage.skin.material,
        ac.fuselage.cone.material,
        ac.fuselage.bendingmaterial_h.material,
        ac.fuselage.bendingmaterial_v.material,
        ac.fuselage.floor.material,
    )
end

function apply_composites!(ac, density_factor::Float64)
    0.0 < density_factor <= 1.0 || error("Composite density factor must be in (0, 1].")
    scaled = IdDict()
    for material in structural_materials(ac)
        haskey(scaled, material) && continue
        material.ρ *= density_factor
        scaled[material] = nothing
    end

    TASOPT.size_aircraft!(ac; iter = 100, printiter = false, fixed_geometry = true)
    return ac
end

function cabin_centroid(ac)
    layout = ac.fuselage.layout
    tank = ac.fuse_tank
    tank_length = ac.parg[iglftank]
    tank.tank_count == 1 || return 0.5 * (
        layout.x_pressure_shell_fwd + layout.x_pressure_shell_aft
    )

    if tank.placement == "front"
        return 0.5 * (
            layout.x_pressure_shell_fwd + tank_length + 2.0 * ft_to_m +
            layout.x_pressure_shell_aft
        )
    end

    return 0.5 * (
        layout.x_pressure_shell_fwd + layout.x_pressure_shell_aft -
        tank_length - 2.0 * ft_to_m
    )
end

function remove_interior!(ac)
    ac.is_sized[1] || error("The E175 must be sized before its interior is removed.")
    fuse = ac.fuselage
    removed_seat_N = fuse.seat.W
    removed_floor_N = fuse.floor.weight.W
    retained_fixed_N = min(fuse.fixed.W, RETAINED_COCKPIT_AND_PILOTS_LB * lbf_to_N)
    removed_fixed_N = fuse.fixed.W - retained_fixed_N
    removed_total_N = removed_seat_N + removed_floor_N + removed_fixed_N

    removed_moment_Nm =
        removed_seat_N * cabin_centroid(ac) +
        removed_floor_N * 0.5 * (
            fuse.layout.x_pressure_shell_fwd + fuse.layout.x_pressure_shell_aft
        ) +
        removed_fixed_N * fuse.fixed.r.x

    fuse.seat.W = 0.0
    fuse.floor.weight.W = 0.0
    fuse.fixed.W = retained_fixed_N
    fuse.weight -= removed_total_N
    fuse.moment -= removed_moment_Nm
    # ac.parg[igWMTO] -= removed_total_N
    ac.parm[imWpay, :] .= 0.0

    _, cabin_pressure_Pa, _, _, _ =
        TASOPT.atmos(CABIN_PRESSURE_ALTITUDE_FT * ft_to_m / 1_000.0)
    ac.parg[igpcabin] = cabin_pressure_Pa
    ac.parg[igmofWpay] = 0.0
    ac.parg[igPofWpay] = 0.0
    return ac
end

end
