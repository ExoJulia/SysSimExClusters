"""
    photoevap_boundary_Carrera2018(R, P)

Compute whether a planet is above (1) or below (0) the photoevaporation boundary defined in Carrera et al. 2018 (Eq. 5).

# Arguments:
- `R::Float64`: planet radius (Earth radii).
- `P::Float64`: period (days).

# Returns:
1 if the planet is above, or 0 if the planet is below, the photevaporation boundary.
"""
function photoevap_boundary_Carrera2018(R::Float64, P::Float64)
    @assert R > 0.
    @assert P > 0.
    Rtrans = 2.6*P^(-0.1467)
    if R >= Rtrans
        return 1
    elseif R < Rtrans
        return 0
    end
end



"""
    calc_transit_duration_central_circ_obs(P; Mstar, Rstar)

Calculate the expected transit duration assuming a circular orbit with impact parameter b=0. Analogous to the function `ExoplanetsSysSim.calc_transit_duration_central_circ_small_angle_approx`.

# Arguments:
- `P::Float64`: period (days).
- `Mstar::Float64`: stellar mass (solar masses).
- `Rstar::Float64`: stellar radius (solar radii).

# Returns:
- `duration::Float64`: transit duration (days) assuming a circular orbit with impact parameter b=0, using observed properties.
"""
function calc_transit_duration_central_circ_obs(P::Float64; Mstar::Float64, Rstar::Float64)
    duration = (ExoplanetsSysSim.rsol_in_au*Rstar*P)/(pi*semimajor_axis(P, Mstar))
    return duration
end



"""
    calc_transit_duration_central_circ_obs(t, pl)

Calculate the expected transit duration assuming a circular orbit with impact parameter b=0. Analogous to the function `ExoplanetsSysSim.calc_transit_duration_central_circ_small_angle_approx` except this function takes a `KeplerTargetObs` instead of a `KeplerTarget` (i.e., it uses the observed properties (stellar mass, radius, and orbital period) instead of the true values).

# Arguments:
- `t::KeplerTargetObs`: object containing the observed planets of a Kepler target.
- `pl::Integer`: index of the observed planet we wish to compute the duration for.

# Returns:
- `duration::Float64`: transit duration (days) assuming a circular orbit with impact parameter b=0, using observed properties.
"""
function calc_transit_duration_central_circ_obs(t::KeplerTargetObs, pl::Integer)
    return calc_transit_duration_central_circ_obs(t.obs[pl].period; Mstar=t.star.mass, Rstar=t.star.radius)
end
