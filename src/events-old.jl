ImpulsiveManeuver(;radial=0.0, along=0.0, cross=0.0) = ImpulsiveManeuver([radial, along, cross])
deltav(man::ImpulsiveManeuver) = norm(man.Δv)

function apply!(man::ImpulsiveManeuver, t, y, params, propagator)
    m = rotation_matrix(RAC, propagator.frame, y)
    y[4:6] += m*man.Δv
end
