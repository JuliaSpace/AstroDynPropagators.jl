function rotational!(δv, t, ep, r, v, params, propagator)
    ep = epoch(params.s0) + t * seconds
    b = propagator.center
    δα = right_ascension_rate(b, ep)
    δ = declination(b, ep)
    δδ = declination_rate(b, ep)
    ω = rotation_angle(b, ep)
    δω = rotation_rate(b, ep)
    χ = π/2 - δ

    ωr = [
        δα * sin(χ) * sin(ω) - δδ * cos(ω),
        δα * sin(χ) * cos(ω) + δδ * sin(ω),
        δα * cos(χ) + δω,
    ]
    δv .-= ωr × (ωr × r) .+ 2(ωr × v)
end
