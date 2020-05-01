if !@isdefined ExoplanetsSysSim
    using ExoplanetsSysSim
end

function generate_stable_cluster(star::StarT, sim_param::SimParam; n::Int64=1) where {StarT<:StarAbstract}

    @assert(n >= 1)

    # Load functions and model parameters:
    generate_sizes = get_function(sim_param, "generate_sizes")
    min_radius::Float64 = get_real(sim_param, "min_radius")
    max_radius::Float64 = get_real(sim_param, "max_radius")
    sigma_log_radius_in_cluster = get_real(sim_param, "sigma_log_radius_in_cluster")
    generate_planet_mass_from_radius = get_function(sim_param, "generate_planet_mass_from_radius")

    min_period::Float64 = get_real(sim_param, "min_period")
    max_period::Float64 = get_real(sim_param, "max_period")
    max_period_ratio = max_period/min_period
    sigma_logperiod_per_pl_in_cluster = get_real(sim_param, "sigma_logperiod_per_pl_in_cluster")

    # If single planet in cluster:
    if n==1
        R = generate_sizes(star, sim_param)
        mass = [generate_planet_mass_from_radius(R[1], sim_param)]
        P = [1.0] # generate_periods_power_law(star,sim_param)
        return (P, R, mass)
    end

    # If reach here, then at least 2 planets in cluster

    # Draw radii and masses:
    mean_R = generate_sizes(star,sim_param)[1]
    Rdist = Truncated(LogNormal(log(mean_R),sigma_log_radius_in_cluster), min_radius, max_radius) # clustered planet sizes
    R = rand(Rdist, n)
    mass = map(r -> generate_planet_mass_from_radius(r, sim_param), R)
    #println("# Rp = ", R)
    #println("# mass = ", mass)

    # Draw unscaled periods, checking for mutual-Hill stability (assuming circular orbits) of the entire cluster as a precondition:
    log_mean_P = 0.0
    Pdist = Truncated(LogNormal(log_mean_P,sigma_logperiod_per_pl_in_cluster*n), 1/sqrt(max_period_ratio), sqrt(max_period_ratio)) # truncated unscaled period distribution to ensure that the cluster can fit in the period range [min_period, max_period] after scaling by a period scale
    local P



    # Draw unscaled periods first, checking for mutual Hill separation stability assuming circular and coplanar orbits

    #= New sampling:
    P = zeros(n)
    for i in 1:n # Draw periods one at a time
        if any(isnan.(P))
            P[i:end] .= NaN
            #println("Cannot fit any more planets in cluster.")
            break
        end
        P[i] = draw_period_lognormal_allowed_regions_mutualHill(P[1:i-1], mass[1:i-1], mass[i], star.mass, sim_param; μ=log_mean_P, σ=n*sigma_logperiod_per_pl_in_cluster, x_min=1/sqrt(max_period_ratio), x_max=sqrt(max_period_ratio))
    end
    found_good_periods = all(isnan.(P)) ? false : true
    @assert(test_stability(P, mass, star.mass, sim_param)) # should always be true if our unscaled period draws are correct
    =#

    # Old rejection sampling:
    found_good_periods = false # will be true if entire cluster is likely to be stable assuming circular and coplanar orbits (given sizes/masses and periods)
    max_attempts = 100
    attempts_periods = 0
    while !found_good_periods && attempts_periods < max_attempts
        attempts_periods += 1
        P = rand(Pdist, n)
        if test_stability(P, mass, star.mass, sim_param)
            # If pass mutual Hill criteria, also check circular MMR overlap criteria:
            a = semimajor_axis.(P, mass .+star.mass)
            μ = mass ./star.mass
            if test_stability_circular_MMR_overlap(μ, a)
                found_good_periods = true
            else
                @info("Found set of periods passing mutual Hill criteria but not circular MMR overlap criteria.")
            end
        end
    end # while trying to draw periods
    #

    return (P, R, mass) # NOTE: can also return earlier if only one planet in cluster; also, the planets are NOT sorted at this point
end

function generate_num_clusters_poisson(s::Star, sim_param::SimParam)
    lambda::Float64 = exp(get_real(sim_param, "log_rate_clusters"))
    max_clusters_in_sys::Int64 = get_int(sim_param, "max_clusters_in_sys")
    return draw_truncated_poisson(lambda, min=0, max=max_clusters_in_sys, n=1)[1]
    #return ExoplanetsSysSim.generate_num_planets_poisson(lambda, max_clusters_in_sys) ##### Use this if setting max_clusters_in_sys > 20
end

function generate_num_clusters_ZTP(s::Star, sim_param::SimParam)
    lambda::Float64 = exp(get_real(sim_param, "log_rate_clusters"))
    max_clusters_in_sys::Int64 = get_int(sim_param, "max_clusters_in_sys")
    return draw_truncated_poisson(lambda, min=1, max=max_clusters_in_sys, n=1)[1]
end

function generate_num_planets_in_cluster_poisson(s::Star, sim_param::SimParam)
    lambda::Float64 = exp(get_real(sim_param, "log_rate_planets_per_cluster"))
    max_planets_in_cluster::Int64 = get_int(sim_param, "max_planets_in_cluster")
    return draw_truncated_poisson(lambda, min=0, max=max_planets_in_cluster, n=1)[1]
end

function generate_num_planets_in_cluster_ZTP(s::Star, sim_param::SimParam)
    lambda::Float64 = exp(get_real(sim_param, "log_rate_planets_per_cluster"))
    max_planets_in_cluster::Int64 = get_int(sim_param, "max_planets_in_cluster")
    return draw_truncated_poisson(lambda, min=1, max=max_planets_in_cluster, n=1)[1]
end

function generate_planetary_system_clustered(star::StarAbstract, sim_param::SimParam; verbose::Bool=false)

    # To include a dependence on stellar color for the fraction of stars with planets:
    #
    global stellar_catalog
    star_color = stellar_catalog[:bp_rp][star.id] - stellar_catalog[:e_bp_min_rp_interp][star.id]
    f_stars_with_planets_attempted_color_slope = get_real(sim_param, "f_stars_with_planets_attempted_color_slope")
    f_stars_with_planets_attempted_at_med_color = get_real(sim_param, "f_stars_with_planets_attempted_at_med_color")
    med_color = get_real(sim_param, "med_color")
    @assert(0 <= f_stars_with_planets_attempted_at_med_color <= 1)

    f_stars_with_planets_attempted = f_stars_with_planets_attempted_color_slope*(star_color - med_color) + f_stars_with_planets_attempted_at_med_color
    f_stars_with_planets_attempted = min(f_stars_with_planets_attempted, 1.)
    f_stars_with_planets_attempted = max(f_stars_with_planets_attempted, 0.)
    @assert(0 <= f_stars_with_planets_attempted <= 1)
    #

    # Load functions and model parameters for drawing planet properties:
    #=
    if haskey(sim_param, "f_stars_with_planets_attempted")
        f_stars_with_planets_attempted = get_real(sim_param, "f_stars_with_planets_attempted")
        @assert(0 <= f_stars_with_planets_attempted <= 1)
    else
        f_stars_with_planets_attempted = 1.
    end
    =#

    generate_num_clusters = get_function(sim_param, "generate_num_clusters")
    generate_num_planets_in_cluster = get_function(sim_param, "generate_num_planets_in_cluster")
    power_law_P = get_real(sim_param, "power_law_P")
    min_period = get_real(sim_param, "min_period")
    max_period = get_real(sim_param, "max_period")

    # Decide whether to assign a planetary system to the star at all:
    if rand() > f_stars_with_planets_attempted
        #println("Star not assigned a planetary system.")
        return PlanetarySystem(star)
    end

    # Assign a reference (invariant) plane for the system:
    # NOTE: aligning sky plane to x-y plane (z-axis is unit normal)
    vec_z = [0.,0.,1.]

    vec_sys = draw_random_normal_vector() # unit normal of reference plane
    incl_sys = calc_angle_between_vectors(vec_z, vec_sys) # inclination of reference plane (rad) relative to sky plane
    Ω_sys = calc_Ω_in_sky_plane(vec_sys) # argument of ascending node of reference plane (rad) relative to sky plane

    # Generate a set of periods, planet radii, and planet masses:
    local num_pl, clusteridlist, Plist, Rlist, masslist, ecclist, ωlist, Ωlist, meananomlist, inclmutlist, inclskylist, Ωskylist

    # First, generate number of clusters (to attempt) and planets (to attempt) in each cluster:
    num_clusters = generate_num_clusters(star, sim_param)::Int64
    num_pl_in_cluster = map(x -> generate_num_planets_in_cluster(star, sim_param)::Int64, 1:num_clusters)
    num_pl_in_cluster_true = zeros(Int64, num_clusters) # true numbers of planets per cluster, after subtracting the number of NaNs
    num_pl = sum(num_pl_in_cluster)
    #println("num_clusters: ", num_clusters, " ; num_pl_in_clusters", num_pl_in_cluster)

    if num_pl==0
        return PlanetarySystem(star)
    end

    clusteridlist::Array{Int64,1} = Array{Int64}(undef, num_pl)
    Plist::Array{Float64,1} = Array{Float64}(undef, num_pl)
    Rlist::Array{Float64,1} = Array{Float64}(undef, num_pl)
    masslist::Array{Float64,1} = Array{Float64}(undef, num_pl)

    @assert(num_pl_in_cluster[1] >= 1)
    pl_start = 1
    pl_stop = 0
    for c in 1:num_clusters
        n = num_pl_in_cluster[c]
        pl_stop += n

        # Draw a stable cluster (with unscaled periods):
        Plist_tmp::Array{Float64,1}, Rlist_tmp::Array{Float64,1}, masslist_tmp::Array{Float64,1} = generate_stable_cluster(star, sim_param, n=num_pl_in_cluster[c])

        clusteridlist[pl_start:pl_stop] = ones(Int64, num_pl_in_cluster[c])*c
        Rlist[pl_start:pl_stop], masslist[pl_start:pl_stop] = Rlist_tmp, masslist_tmp

        #= New sampling:
        idx = .!isnan.(Plist[1:pl_stop-n])
        idy = .!isnan.(Plist_tmp)
        if any(idy)
            min_period_scale::Float64 = min_period/minimum(Plist_tmp[idy])
            max_period_scale::Float64 = max_period/maximum(Plist_tmp[idy])

            if min_period_scale < max_period_scale
                period_scale = draw_periodscale_power_law_allowed_regions_mutualHill(num_pl_in_cluster_true[1:c-1], Plist[1:pl_stop-n][idx], masslist[1:pl_stop-n][idx], Plist_tmp[idy], masslist_tmp[idy], star.mass, sim_param; x0=min_period_scale, x1=max_period_scale, α=power_law_P)
            else # cluster cannot fit at all
                period_scale = NaN
            end
        else # void cluster; all NaNs
            period_scale = NaN
        end
        Plist[pl_start:pl_stop] = Plist_tmp .* period_scale
        if isnan(period_scale)
            Plist[pl_stop:end] .= NaN
            #println("Cannot fit cluster into system; returning clusters that did fit.")
            break
        end
        @assert(test_stability(view(Plist,1:pl_stop), view(masslist,1:pl_stop), star.mass, sim_param)) # should always be true if our period scale draws are correct
        =#

        # Old rejection sampling:
        valid_cluster = !any(isnan.(Plist_tmp)) # if the cluster has any nans, the whole cluster is discarded
        valid_period_scale = false
        max_attempts_period_scale = 100
        attempts_period_scale = 0
        min_period_scale::Float64 = min_period/minimum(Plist_tmp)
        max_period_scale::Float64 = max_period/maximum(Plist_tmp)
        while !valid_period_scale && attempts_period_scale<max_attempts_period_scale && valid_cluster
            attempts_period_scale += 1

            period_scale::Array{Float64,1} = draw_power_law(power_law_P, min_period_scale, max_period_scale, 1)

            Plist[pl_start:pl_stop] = Plist_tmp .* period_scale

            if test_stability(view(Plist,1:pl_stop), view(masslist,1:pl_stop), star.mass, sim_param)
                # If pass mutual Hill criteria, also check circular MMR overlap criteria:
                alist = semimajor_axis.(view(Plist,1:pl_stop), view(masslist,1:pl_stop) .+star.mass)
                μlist = view(masslist,1:pl_stop) ./star.mass
                if test_stability_circular_MMR_overlap(μlist, alist)
                    valid_period_scale = true
                else
                    @info("Found period scale passing mutual Hill criteria but not circular MMR overlap criteria.")
                end
            end
        end  # while !valid_period_scale...

        #if attempts_period_scale > 1
            #println("attempts_period_scale: ", attempts_period_scale)
        #end

        if !valid_period_scale
            Plist[pl_start:pl_stop] .= NaN
        end
        #

        num_pl_in_cluster_true[c] = sum(.!isnan.(Plist[pl_start:pl_stop]))
        pl_start += n
    end # for c in 1:num_clusters

    # Discard failed planets and sort the remaining planets by period:

    isnanPlist::Array{Bool,1} = isnan.(Plist::Array{Float64,1})
    if any(isnanPlist) # if any loop failed to generate valid planets, it should set a NaN in the period list
        keep::Array{Bool,1} = .!(isnanPlist) # currently, keeping clusters that could be fit, rather than throwing out entire systems and starting from scratch
        num_pl = sum(keep)

        clusteridlist = clusteridlist[keep]
        Plist = Plist[keep]
        Rlist = Rlist[keep]
        masslist = masslist[keep]
    end

    # Return the system (star) early if failed to draw any planets:
    if num_pl==0
        return PlanetarySystem(star)
    end

    idx = sortperm(Plist)
    Plist = Plist[idx]
    Rlist = Rlist[idx]
    masslist = masslist[idx]

    # Now compute the critical AMD for the system and distribute it between the planets (to draw eccentricities and inclinations):

    f_amd_crit = get_real(sim_param, "f_amd_crit")
    f_amd_crit = f_amd_crit + (1 - f_amd_crit)*rand() # uniform between f_amd_crit and 1
    @assert(0 <= f_amd_crit <= 1)

    ecclist::Array{Float64,1} = Array{Float64}(undef, num_pl)
    ωlist::Array{Float64,1} = Array{Float64}(undef, num_pl)
    Ωlist::Array{Float64,1} = Array{Float64}(undef, num_pl)
    meananomlist::Array{Float64,1} = Array{Float64}(undef, num_pl)
    inclmutlist::Array{Float64,1} = Array{Float64}(undef, num_pl)
    AMDlist::Array{Float64,1} = Array{Float64}(undef, num_pl)

    μlist = masslist ./ star.mass # mass ratios
    alist = map(P -> semimajor_axis(P, star.mass), Plist)
    if num_pl == 1
        generate_e_omega = get_function(sim_param, "generate_e_omega")
        sigma_ecc::Float64 = get_real(sim_param, "sigma_hk")
        ecc::Float64, ω::Float64 = generate_e_omega(sigma_ecc)
        AMDlist, ecclist, ωlist, inclmutlist = μlist .*sqrt.(alist) .*(1 - sqrt(1 - ecc^2)), [ecc], [ω], [0.]
    else
        AMDlist, ecclist, ωlist, inclmutlist = draw_ecc_incl_system_critical_AMD(μlist, alist; f_amd_crit=f_amd_crit, check_stability=false)
    end

    Ωlist, meananomlist = 2π .* rand(num_pl), 2π .* rand(num_pl) # relative to the reference plane

    # Calculate the sky orientations of each orbit:
    inclskylist::Array{Float64,1} = Array{Float64}(undef, num_pl)
    Ωskylist::Array{Float64,1} = Array{Float64}(undef, num_pl)

    inclskylist, Ωskylist = calc_sky_incl_Ω_orbits_given_system_vector(inclmutlist, Ωlist, vec_sys)

    #= For debugging:
    println("P (d): ", Plist)
    println("R (R⊕): ", Rlist ./ ExoplanetsSysSim.earth_radius)
    println("M (M⊕): ", masslist ./ ExoplanetsSysSim.earth_mass)
    println("e: ", ecclist)
    println("ω (deg): ", ωlist .* 180/π)
    println("i_m (deg): ", inclmutlist .* 180/π)
    println("i (deg): ", inclskylist .* 180/π)
    println("AMD: ", AMDlist)
    =#

    # Final checks and returning the system:

    # This final mutual-Hill stability test should be unnecessary if the period scales were drawn properly
    # NOTE: NOT including eccentricities in final mutual-Hill stability test since they were drawn by distributing AMD after the periods were set
    @assert(test_stability(Plist, masslist, star.mass, sim_param))

    pl = Array{Planet}(undef, num_pl)
    orbit = Array{Orbit}(undef, num_pl)
    for i in 1:num_pl
        pl[i] = Planet(Rlist[i], masslist[i], clusteridlist[i])
        orbit[i] = Orbit(Plist[i], ecclist[i], inclskylist[i], ωlist[i], Ωskylist[i], meananomlist[i])
    end
    sys_ref_plane = SystemPlane(incl_sys, Ω_sys)

    return PlanetarySystem(star, pl, orbit, sys_ref_plane)
end
