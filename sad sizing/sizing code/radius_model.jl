weight_lb(weight_N) = weight_N / lbf_to_N

struct RadiusVariant
    ac::TASOPT.aircraft
    bow_N::Float64
    mtow_N::Float64
    structural_payload_N::Float64
    fuselage_weight_N::Float64
    landing_gear_weight_N::Float64
    sizing_payload_N::Float64
    sizing_fuel_N::Float64
    fixed_gear_tailstrike_angle_rad::Float64
    gear_was_redesigned::Bool
end

struct BaselineGearGeometry
    nose_length_m::Float64
    main_length_m::Float64
end

struct FixedInstalledCFM56Map end

function (calculator::FixedInstalledCFM56Map)(
    ac,
    case::String,
    imission::Int64,
    ip::Int64,
    initializes_engine::Bool,
    iterw::Int64 = 0,
)
    fixed_case = case == "design" ? "off_design" : case
    TASOPT.engine.tfwrap!(ac, fixed_case, imission, ip, initializes_engine, iterw)
    return nothing
end

function install_fixed_cfm56!(ac::TASOPT.aircraft)
    old_model = ac.engine.model
    installed = (
        engine_weight = ac.parg[igWeng],
        bare_engine_weight = ac.parg[igWebare],
        nacelle_weight = ac.parg[igWnace],
        nacelle_area_fraction = ac.parg[igfSnace],
        nacelle_length = ac.parg[iglnace],
    )
    fixed_weight! = function (trial)
        trial.parg[igWeng] = installed.engine_weight
        trial.parg[igWebare] = installed.bare_engine_weight
        trial.parg[igWnace] = installed.nacelle_weight
        trial.parg[igfSnace] = installed.nacelle_area_fraction
        trial.parg[iglnace] = installed.nacelle_length
        return nothing
    end
    ac.engine = TASOPT.engine.Engine(
        TASOPT.engine.TurbofanModel(
            old_model.model_name,
            FixedInstalledCFM56Map(),
            old_model.weight_model_name,
            fixed_weight!,
            old_model.has_BLI_cores,
        ),
        ac.engine.data,
        ac.engine.heat_exchangers,
    )
    return ac
end

function set_fuselage_radius!(ac::TASOPT.aircraft, radius_in::Float64)
    scale = radius_in * in_to_m / ac.fuselage.layout.radius
    section = ac.fuselage.layout.cross_section
    section.radius *= scale
    section.bubble_lower_downward_shift *= scale
    if hasfield(typeof(section), :bubble_center_y_offset)
        setfield!(section, :bubble_center_y_offset,
            getfield(section, :bubble_center_y_offset) * scale)
    end
    return ac
end

function move_wing_inboard!(ac::TASOPT.aircraft, radius_reduction_m::Float64)
    wing = ac.wing
    span_reduction_m = 2.0 * radius_reduction_m
    break_span_m = wing.layout.break_span - span_reduction_m
    wing.layout.root_span -= span_reduction_m
    wing.layout.span -= span_reduction_m
    wing.layout.ηs = break_span_m / wing.layout.span
    wing.layout.AR = wing.layout.span^2 / wing.layout.S
    return ac
end

function update_fuselage_structure!(ac::TASOPT.aircraft)
    parg, parm, para, pare = ac.parg, ac.parm, ac.para, ac.pare
    fuse, fuse_tank = ac.fuselage, ac.fuse_tank
    wing, htail, vtail = ac.wing, ac.htail, ac.vtail

    TASOPT.set_ambient_conditions!(ac, ipcruise1; im = 1)
    TASOPT.fuselage_drag!(fuse, @view(parm[:, 1]), @view(para[:, :, 1]), ipcruise1)
    TASOPT.broadcast_fuselage_drag!(@view(para[:, :, 1]), ipcruise1)

    qne = 0.5 * TASOPT.ρSL * parg[igVne]^2
    Lhmax = qne * htail.layout.S * htail.CL_max
    Lvmax = qne * vtail.layout.S * vtail.CL_max / vtail.ntails
    wcd = para[iafracW, ipcruisen, 1] / para[iafracW, ipcruise1, 1]
    delta_p = parg[igpcabin] - pare[iep0, ipcruise1, 1] * wcd
    parg[igdeltap] = delta_p

    engine_location = lowercase(strip(ac.options.opt_engine_location))
    Wengtail = engine_location == "wing" ? 0.0 :
        engine_location == "fuselage" ? parg[igWeng] :
        error("Unknown engine location: $(ac.options.opt_engine_location)")

    nftanks = fuse_tank.tank_count
    Waftfuel = 0.0
    Wftank_single = 0.0
    ltank = 0.0
    xftank_fuse = 0.0
    tank_placement = ""
    if !ac.options.has_wing_fuel
        tank_placement = fuse_tank.placement
        Wftank_single = parg[igWftank] / nftanks
        ltank = parg[iglftank]
        if tank_placement == "rear"
            Waftfuel = parg[igWfuel]
            xftank_fuse = parg[igxftankaft]
        elseif tank_placement == "both"
            Waftfuel = parg[igWfuel] / 2.0
            xftank_fuse = parg[igxftankaft]
        elseif tank_placement == "front"
            Waftfuel = parg[igWfuel]
            xftank_fuse = fuse.layout.x_end - parg[igxftank]
        end
    end

    cbox = wing.layout.root_chord * wing.inboard.cross_section.width_to_chord
    parg[igcabVol] = TASOPT.fusew!(
        fuse, parg[igNland], parg[igWpaymax], Wengtail,
        nftanks, Waftfuel, Wftank_single, ltank, xftank_fuse, tank_placement,
        delta_p, htail.weight, vtail.weight, parg[igrMh], parg[igrMv],
        Lhmax, Lvmax, vtail.layout.span, vtail.outboard.λ, vtail.ntails,
        htail.layout.x, vtail.layout.x, wing.layout.x, wing.layout.box_x,
        cbox, parg[igxeng],
    )
    return ac
end

landing_gear_weight(ac::TASOPT.aircraft) =
    ac.landing_gear.nose_gear.weight.W + ac.landing_gear.main_gear.weight.W

function available_tailstrike_angle(
    ac::TASOPT.aircraft,
    equivalent_gear_length_m::Float64,
)
    main_gear = ac.landing_gear.main_gear
    x_lg = ac.parg[igxCGaft] + main_gear.distance_CG_to_landing_gear
    tail_arm = ac.fuselage.layout.x_end - x_lg
    tail_arm > 0.0 || error("Main gear must be forward of the fuselage end.")
    return atan(
        (equivalent_gear_length_m + 2.0 * ac.fuselage.layout.radius) /
        tail_arm,
    )
end

function update_fixed_length_landing_gear_weights!(
    ac::TASOPT.aircraft,
    baseline::BaselineGearGeometry,
)
    landing_gear = ac.landing_gear
    nose_gear = landing_gear.nose_gear
    main_gear = landing_gear.main_gear
    WMTO = ac.parg[igWMTO]

    if lowercase(landing_gear.model) != "historical_correlations"
        TASOPT.size_landing_gear!(ac)
        return ac
    end

    load_factor = 4.5
    Vstall = ac.pare[ieu0, iprotate, 1]
    main_nstruts = main_gear.number_struts
    main_nwheels = main_nstruts * main_gear.wheels_per_strut
    nose_nwheels = nose_gear.number_struts * nose_gear.wheels_per_strut

    Wmain = 0.0106 * (WMTO * lb_N)^0.888 * load_factor^0.25 *
        (baseline.main_length_m / in_to_m)^0.4 * main_nwheels^0.321 *
        main_nstruts^(-0.5) * (Vstall / kts_to_mps)^0.1 / lb_N
    Wnose = 0.032 * (WMTO * lb_N)^0.646 * load_factor^0.2 *
        (baseline.nose_length_m / in_to_m)^0.5 * nose_nwheels^0.45 / lb_N

    nose_gear.length = baseline.nose_length_m
    main_gear.length = baseline.main_length_m
    nose_gear.overall_mass_fraction = Wnose / WMTO
    main_gear.overall_mass_fraction = Wmain / WMTO
    nose_gear.weight = TASOPT.Weight(W = Wnose, x = nose_gear.weight.r[1])
    main_gear.weight = TASOPT.Weight(
        W = Wmain,
        x = ac.parg[igxCGaft] + main_gear.distance_CG_to_landing_gear,
        y = ac.wing.span / 2 * main_gear.y_offset_halfspan_fraction,
    )
    return ac
end

function update_landing_gear_for_radius!(
    ac::TASOPT.aircraft,
    baseline::BaselineGearGeometry,
)
    available_angle = available_tailstrike_angle(ac, baseline.nose_length_m)
    required_angle = ac.landing_gear.tailstrike_angle
    redesigned = available_angle + 1.0e-10 < required_angle
    if redesigned
        TASOPT.size_landing_gear!(ac)
    else
        update_fixed_length_landing_gear_weights!(ac, baseline)
    end
    return redesigned, available_angle
end

function quicksave_radius_variant(
    variant::RadiusVariant,
    installed_engine::TASOPT.engine.Engine,
    filepath::String,
    settings::Settings,
)
    ac = deepcopy(variant.ac)
    ac.engine = deepcopy(installed_engine)
    ac.parg[igWMTO] = variant.mtow_N

    payload_N = variant.sizing_payload_N
    fuel_N = variant.sizing_fuel_N
    ac.parg[igWpay] = payload_N
    ac.parg[igWfuel] = fuel_N
    ac.parm[imWTO, 1] = variant.mtow_N
    ac.parm[imWpay, 1] = payload_N
    ac.parm[imWfuel, 1] = fuel_N
    ac.parm[imRange, 1] = settings.range_nmi * nmi_to_m

    mkpath(dirname(filepath))
    TASOPT.quicksave_aircraft(ac, filepath)
    return filepath
end

function sizing_mission_fuel!(
    ac::TASOPT.aircraft,
    bow_N::Float64,
    payload_N::Float64,
    settings::Settings,
)
    fuel_seed_N = ac.parg[igWfmax]
    ac.parg[igWMTO] = bow_N + payload_N + fuel_seed_N
    ac.parg[igWpay] = payload_N
    ac.parg[igWfuel] = fuel_seed_N
    configure_ceiling_mission!(
        ac,
        weight_lb(payload_N),
        settings.mtow_closure_altitude_ft,
        settings,
    )

    redirect_stdout(devnull) do
        redirect_stderr(devnull) do
            TASOPT.fly_mission!(
                ac,
                2;
                itermax = settings.mtow_closure_mission_iterations,
                initializes_engine = true,
                opt_prescribed_cruise_parameter = "altitude",
            )
        end
    end

    sum(ac.pare[ieConvFail, :, 2]) == 0.0 ||
        error("The fixed-engine sizing mission did not converge.")
    fuel_N = ac.parm[imWfuel, 2]
    isfinite(fuel_N) && fuel_N > 0.0 ||
        error("The sizing mission returned an invalid fuel weight.")
    fuel_N <= ac.parg[igWfmax] || error(
        "Sizing mission requires $(round(weight_lb(fuel_N), digits = 1)) lb of fuel, " *
        "above the $(round(weight_lb(ac.parg[igWfmax]), digits = 1))-lb capacity.")
    return fuel_N
end

function find_reference_sizing_payload(
    base_ac::TASOPT.aircraft,
    airframe::CFM56Installation.FixedAirframe,
    settings::Settings,
)
    ac = deepcopy(base_ac)
    install_fixed_cfm56!(ac)
    ac.parg[igfreserve] = settings.fuel_reserve_fraction
    update_radius_dependent_weights!(ac, airframe.mtow_N)
    bow_N = airframe.nonengine_bow_N + ac.parg[igWeng]

    function payload_fits(payload_N)
        trial = deepcopy(ac)
        try
            fuel_N = sizing_mission_fuel!(
                trial,
                bow_N,
                payload_N,
                settings,
            )
            return bow_N + payload_N + fuel_N <= airframe.mtow_N
        catch
            return false
        end
    end

    reference_tank = SO2Payload.size_tank(
        settings.reference_payload_lb,
        settings.tank;
        usable_inner_diameter_m = max(
            0.0,
            2.0 * (base_ac.fuselage.layout.radius - base_ac.fuselage.skin.thickness),
        ),
    )
    low_N = reference_tank.total_payload_lb * lbf_to_N
    high_N = min(
        airframe.payload_max_N,
        airframe.mtow_N - bow_N - 1.0 * lbf_to_N,
    )
    payload_fits(low_N) || error(
        "The reference 500-nmi sizing mission is infeasible at " *
        "$(round(settings.reference_payload_lb, digits = 0)) lb SO2 plus " *
        "$(round(reference_tank.tank_mass_lb, digits = 1)) lb of tank mass.")
    if payload_fits(high_N)
        return high_N
    end

    tolerance_N = settings.reference_payload_tolerance_lb * lbf_to_N
    while high_N - low_N > tolerance_N
        payload_N = 0.5 * (low_N + high_N)
        if payload_fits(payload_N)
            low_N = payload_N
        else
            high_N = payload_N
        end
    end
    return low_N
end

function update_radius_dependent_weights!(
    ac::TASOPT.aircraft,
    mtow_N::Float64;
    baseline_gear::Union{Nothing,BaselineGearGeometry} = nothing,
)
    ac.parg[igWMTO] = mtow_N
    update_fuselage_structure!(ac)
    if isnothing(baseline_gear)
        TASOPT.size_landing_gear!(ac)
    else
        update_landing_gear_for_radius!(ac, baseline_gear)
    end
    return ac
end

function size_radius_variant(
    base_ac::TASOPT.aircraft,
    airframe::CFM56Installation.FixedAirframe,
    radius_in::Float64,
    reference_sizing_payload_N::Float64,
    settings::Settings,
)
    reference = deepcopy(base_ac)
    install_fixed_cfm56!(reference)
    reference.parg[igfreserve] = settings.fuel_reserve_fraction
    update_radius_dependent_weights!(reference, airframe.mtow_N)
    reference_fuse_N = reference.fuselage.weight
    reference_gear_N = landing_gear_weight(reference)
    baseline_gear = BaselineGearGeometry(
        reference.landing_gear.nose_gear.length,
        reference.landing_gear.main_gear.length,
    )

    ac = deepcopy(base_ac)
    install_fixed_cfm56!(ac)
    ac.parg[igfreserve] = settings.fuel_reserve_fraction
    baseline_radius_m = ac.fuselage.layout.radius
    set_fuselage_radius!(ac, radius_in)
    move_wing_inboard!(ac, baseline_radius_m - ac.fuselage.layout.radius)

    payload_scale = (radius_in / settings.baseline_radius_in)^2
    ac.parg[igWpaymax] = airframe.payload_max_N * payload_scale

    sizing_payload_N = min(
        ac.parg[igWpaymax],
        reference_sizing_payload_N * payload_scale,
    )
    sizing_payload_N > 0.0 || error("Radius $radius_in in has no structural payload capacity.")
    baseline_bow_N = airframe.nonengine_bow_N + ac.parg[igWeng]
    update_radius_dependent_weights!(
        ac,
        airframe.mtow_N;
        baseline_gear,
    )
    initial_bow_N = baseline_bow_N +
        (ac.fuselage.weight - reference_fuse_N) +
        (landing_gear_weight(ac) - reference_gear_N)
    mtow_N = max(
        airframe.mtow_N,
        initial_bow_N + sizing_payload_N + ac.parg[igWfmax],
    )
    sizing_fuel_N = NaN
    bow_N = NaN
    converged = false
    for _ in 1:settings.mtow_closure_max_iterations
        update_radius_dependent_weights!(ac, mtow_N; baseline_gear)
        delta_empty_N = (ac.fuselage.weight - reference_fuse_N) +
            (landing_gear_weight(ac) - reference_gear_N)
        bow_N = baseline_bow_N + delta_empty_N
        sizing_fuel_N = sizing_mission_fuel!(
            ac,
            bow_N,
            sizing_payload_N,
            settings,
        )
        new_mtow_N = bow_N + sizing_payload_N + sizing_fuel_N
        if abs(new_mtow_N - mtow_N) <= settings.mtow_closure_tolerance_lb * lbf_to_N
            mtow_N = new_mtow_N
            converged = true
            break
        end
        mtow_N = (1.0 - settings.mtow_closure_relaxation) * mtow_N +
            settings.mtow_closure_relaxation * new_mtow_N
    end
    converged || error("MTOW mission closure did not converge for radius $radius_in in.")

    update_radius_dependent_weights!(ac, mtow_N; baseline_gear)
    bow_N = baseline_bow_N + (ac.fuselage.weight - reference_fuse_N) +
        (landing_gear_weight(ac) - reference_gear_N)
    sizing_fuel_N = sizing_mission_fuel!(
        ac,
        bow_N,
        sizing_payload_N,
        settings,
    )
    mtow_N = bow_N + sizing_payload_N + sizing_fuel_N
    ac.parg[igWMTO] = mtow_N
    ac.parg[igWpay] = sizing_payload_N
    ac.parg[igWfuel] = sizing_fuel_N

    fixed_gear_tailstrike_angle_rad = available_tailstrike_angle(
        ac,
        baseline_gear.nose_length_m,
    )
    gear_was_redesigned = fixed_gear_tailstrike_angle_rad + 1.0e-10 <
        ac.landing_gear.tailstrike_angle

    return RadiusVariant(
        ac, bow_N, mtow_N, ac.parg[igWpaymax],
        ac.fuselage.weight, landing_gear_weight(ac),
        sizing_payload_N, sizing_fuel_N,
        fixed_gear_tailstrike_angle_rad, gear_was_redesigned,
    )
end
