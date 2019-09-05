if !@isdefined ExoplanetsSysSim
    using ExoplanetsSysSim
end

function generate_stable_cluster(star::StarT, sim_param::SimParam, incl_sys::Float64, sigma_incl_use::Float64; n::Int64=1) where {StarT<:StarAbstract}

    @assert n >= 1
    @assert sigma_incl_use >= 0

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

    generate_e_omega = get_function(sim_param, "generate_e_omega")
    sigma_ecc::Float64 = haskey(sim_param, "sigma_hk") ? get_real(sim_param, "sigma_hk") : 0.0
    sigma_incl = deg2rad(get_real(sim_param, "sigma_incl"))
    sigma_incl_near_mmr = deg2rad(get_real(sim_param, "sigma_incl_near_mmr"))

    # If single planet in cluster:
    if n==1
        R = generate_sizes(star, sim_param)
        mass = [generate_planet_mass_from_radius(R[1], sim_param)]
        P = [1.0] # generate_periods_power_law(star,sim_param)
        eccentricity::Float64, omega_e::Float64 = generate_e_omega(sigma_ecc)::Tuple{Float64,Float64}
        ecc, omega = [eccentricity], [omega_e]
        asc_node, mean_anom, incl = [2pi*rand()], [2pi*rand()], [incl_sys]
        return (P, R, mass, ecc, omega, asc_node, mean_anom, incl)
    end

    # If reach here, then at least 2 planets in cluster

    # Draw radii and masses:
    mean_R = generate_sizes(star,sim_param)[1]
    Rdist = Truncated(LogNormal(log(mean_R),sigma_log_radius_in_cluster), min_radius, max_radius) # clustered planet sizes
    R = rand(Rdist, n)
    #R = generate_sizes(star,sim_param, num_pl=n) ##### if want non-clustered planet sizes instead
    mass = map(r -> generate_planet_mass_from_radius(r, sim_param), R)
    #println("# Rp = ", R)
    #println("# mass = ", mass)

    # Draw unscaled periods, eccentricities, and mutual inclinations, checking for stability of the entire cluster:
    log_mean_P = 0.0 # log(generate_periods_power_law(star,sim_param))
    Pdist = Truncated(LogNormal(log_mean_P,sigma_logperiod_per_pl_in_cluster*n), 1/sqrt(max_period_ratio), sqrt(max_period_ratio)) # truncated unscaled period distribution to ensure that the cluster can fit in the period range [min_period, max_period] after scaling by a period scale
    local P
    local ecc, omega, asc_nodes, mean_anoms, incls
    ecc = Array{Float64}(undef, n)
    omega = Array{Float64}(undef, n)
    asc_nodes = Array{Float64}(undef, n)
    mean_anoms = Array{Float64}(undef, n)
    incls_mut = Array{Float64}(undef, n)
    incls = Array{Float64}(undef, n)

    found_stable_cluster = false # will be true if entire cluster is stable (given sizes/masses, periods, eccentricities, and mutual inclinations)
    max_attempts = 100
    attempts_cluster = 0
    while !found_stable_cluster && attempts_cluster < max_attempts
        attempts_cluster += 1

        # Draw unscaled periods first, checking for mutual Hill separation stability assuming circular and coplanar orbits
        found_good_periods = false # will be true if entire cluster is likely to be stable assuming circular and coplanar orbits (given sizes/masses and periods)
        attempts_periods = 0
        while !found_good_periods && attempts_periods < max_attempts
            attempts_periods += 1
            P = rand(Pdist, n)
            if test_stability(P, mass, star.mass, sim_param)
                found_good_periods = true
            end
        end # while trying to draw periods

        # Draw eccentricities and mutual inclinations (configure orbits), checking for stability of the entire cluster:
        if found_good_periods
            idx = sortperm(P)
            is_near_resonance = calc_if_near_resonance(P[idx], sim_param) # NOTE: this only checks if each planet is near a resonance with another planet in the same cluster; it does not know about whether each planet is near a resonance with a planet in a different cluster (since period scales have not been assigned yet)
            for i in 1:n # NOTE: this 'i' is indexing the planets in this cluster sorted by their periods
                (ecc[idx[i]], omega[idx[i]]) = generate_e_omega(sigma_ecc)::Tuple{Float64,Float64}

                incl_mut = sigma_incl_use*sqrt(randn()^2+randn()^2) # rand(Distributions.Rayleigh(sigma_incl_use))
                #
                if is_near_resonance[i]
                    incl_mut = min(sigma_incl, sigma_incl_near_mmr)*sqrt(randn()^2+randn()^2)
                end
                #
                asc_node = 2pi*rand()
                mean_anom = 2pi*rand()
                incl = incl_mut != zero(incl_mut) ? acos(cos(incl_sys)*cos(incl_mut) + sin(incl_sys)*sin(incl_mut)*cos(asc_node)) : incl_sys

                asc_nodes[idx[i]], mean_anoms[idx[i]], incls_mut[idx[i]], incls[idx[i]] = asc_node, mean_anom, incl_mut, incl
            end

            μ = mass./star.mass
            a = map(P -> semimajor_axis(P, star.mass), P)
            stats, pairs, ratios = AMD_stability(μ[idx], a[idx], ecc[idx], incls_mut[idx])
            #if test_stability(P, mass, star.mass, sim_param; ecc=ecc)
            if all(stats .== :Stable)
                found_stable_cluster = true
            end
        end # if found_good_periods

    end # while trying to make the cluster stable
    println("attempts for cluster: ", attempts_cluster)

    if !found_stable_cluster
        #println("# Warning: Did not find a good set of sizes, masses, periods, eccentricities, and mutual inclinations for one cluster.")
        return (fill(NaN,n), R, mass, ecc, omega, asc_nodes, mean_anoms, incls)  # Return NaNs for periods to indicate failed
    end
    return (P, R, mass, ecc, omega, asc_nodes, mean_anoms, incls) # NOTE: can also return earlier if only one planet in cluster or if fail to generate a good set of values; also, the planets are NOT sorted at this point
end

function generate_num_planets_in_cluster_poisson(s::Star, sim_param::SimParam)
    lambda::Float64 = exp(get_real(sim_param, "log_rate_planets_per_cluster"))
    max_planets_in_cluster::Int64 = get_int(sim_param, "max_planets_in_cluster")
    return draw_truncated_poisson(lambda, min=1, max=max_planets_in_cluster, n=1)[1]
end

function generate_num_clusters_poisson(s::Star, sim_param::SimParam)
    lambda::Float64 = exp(get_real(sim_param, "log_rate_clusters"))
    max_clusters_in_sys::Int64 = get_int(sim_param, "max_clusters_in_sys")
    return draw_truncated_poisson(lambda, min=0, max=max_clusters_in_sys, n=1)[1]
    #return ExoplanetsSysSim.generate_num_planets_poisson(lambda, max_clusters_in_sys) ##### Use this if setting max_clusters_in_sys > 20
end

function generate_planetary_system_clustered(star::StarAbstract, sim_param::SimParam; verbose::Bool=false)

    # Load functions and model parameters for drawing planet properties:
    generate_num_clusters = get_function(sim_param, "generate_num_clusters")
    generate_num_planets_in_cluster = get_function(sim_param, "generate_num_planets_in_cluster")
    power_law_P = get_real(sim_param, "power_law_P")
    min_period = get_real(sim_param, "min_period")
    max_period = get_real(sim_param, "max_period")

    sigma_incl = deg2rad(get_real(sim_param, "sigma_incl"))
    sigma_incl_near_mmr = deg2rad(get_real(sim_param, "sigma_incl_near_mmr"))
    max_incl_sys = get_real(sim_param, "max_incl_sys")
    f_high_incl = get_real(sim_param, "f_high_incl")
    sigma_incl_use = rand() < f_high_incl ? max(sigma_incl, sigma_incl_near_mmr) : min(sigma_incl, sigma_incl_near_mmr)

    # Assign a reference plane for the system:
    incl_sys = acos(cos(max_incl_sys*pi/180)*rand()) # acos(rand()) for isotropic distribution of system inclinations; acos(cos(X*pi/180)*rand()) gives angles from X (deg) to 90 (deg)

    # Generate a set of periods, planet radii, and planet masses:
    attempts_system = 0
    max_attempts_system = 1 # NOTE: currently this should not matter; each system is always attempted just once
    local num_pl, clusteridlist, Plist, Rlist, masslist, ecclist, omegalist, ascnodelist, meananomlist, incllist
    valid_system = false
    while !valid_system && attempts_system < max_attempts_system

        # First, generate number of clusters (to attempt) and planets (to attempt) in each cluster:
        num_clusters = generate_num_clusters(star, sim_param)::Int64
        num_pl_in_cluster = map(x -> generate_num_planets_in_cluster(star, sim_param)::Int64, 1:num_clusters)
        num_pl = sum(num_pl_in_cluster)
        #println("num_clusters: ", num_clusters, " ; num_pl_in_clusters", num_pl_in_cluster)

        if num_pl==0
            return PlanetarySystem(star)
        end

        clusteridlist::Array{Int64,1} = Array{Int64}(undef, num_pl)
        Plist::Array{Float64,1} = Array{Float64}(undef, num_pl)
        Rlist::Array{Float64,1} = Array{Float64}(undef, num_pl)
        masslist::Array{Float64,1} = Array{Float64}(undef, num_pl)
        ecclist::Array{Float64,1} = Array{Float64}(undef, num_pl)
        omegalist::Array{Float64,1} = Array{Float64}(undef, num_pl)
        ascnodelist::Array{Float64,1} = Array{Float64}(undef, num_pl)
        meananomlist::Array{Float64,1} = Array{Float64}(undef, num_pl)
        incllist::Array{Float64,1} = Array{Float64}(undef, num_pl)

        @assert num_pl_in_cluster[1] >= 1
        pl_start = 1
        pl_stop = 0
        for c in 1:num_clusters
            pl_stop += num_pl_in_cluster[c]

            # Draw a stable cluster (with unscaled periods):
            Plist_tmp::Array{Float64,1}, Rlist_tmp::Array{Float64,1}, masslist_tmp::Array{Float64,1}, ecclist_tmp::Array{Float64,1}, omegalist_tmp::Array{Float64,1}, ascnodelist_tmp::Array{Float64,1}, meananomlist_tmp::Array{Float64,1}, incllist_tmp::Array{Float64,1} = generate_stable_cluster(star, sim_param, incl_sys, sigma_incl_use, n=num_pl_in_cluster[c])

            clusteridlist[pl_start:pl_stop] = ones(Int64, num_pl_in_cluster[c])*c
            Rlist[pl_start:pl_stop], masslist[pl_start:pl_stop] = Rlist_tmp, masslist_tmp
            ecclist[pl_start:pl_stop], omegalist[pl_start:pl_stop] = ecclist_tmp, omegalist_tmp
            ascnodelist[pl_start:pl_stop], meananomlist[pl_start:pl_stop], incllist[pl_start:pl_stop] = ascnodelist_tmp, meananomlist_tmp, incllist_tmp

            valid_cluster = !any(isnan.(Plist_tmp)) # if the cluster has any nans, the whole cluster is discarded
            valid_period_scale = false
            max_attempts_period_scale = 100
            attempts_period_scale = 0
            while !valid_period_scale && attempts_period_scale<max_attempts_period_scale && valid_cluster
                attempts_period_scale += 1

                period_scale::Array{Float64,1} = draw_power_law(power_law_P, min_period/minimum(Plist_tmp), max_period/maximum(Plist_tmp), 1)
                # NOTE: this ensures that the minimum and maximum periods will be in the range [min_period, max_period]
                # WARNING: not sure about the behaviour when min_period/minimum(Plist_tmp) > max_period/maximum(Plist_tmp) (i.e. when the cluster cannot fit in the given range)?
                # TODO OPT: could draw period_scale more efficiently by computing the allowed regions in [min_period, max_period] given the previous cluster draws (difficult)

                Plist[pl_start:pl_stop] = Plist_tmp .* period_scale

                if test_stability(view(Plist,1:pl_stop), view(masslist,1:pl_stop), star.mass, sim_param; ecc=view(ecclist,1:pl_stop))
                    valid_period_scale = true
                end
            end  # while !valid_period_scale...

            #if attempts_period_scale > 1
                #println("attempts_period_scale: ", attempts_period_scale)
            #end

            if !valid_period_scale
                Plist[pl_start:pl_stop] .= NaN
            end
            pl_start += num_pl_in_cluster[c]
        end # for c in 1:num_clusters

        isnanPlist::Array{Bool,1} = isnan.(Plist::Array{Float64,1})
        if any(isnanPlist)  # if any loop failed to generate valid planets, it should set a NaN in the period list
            keep::Array{Bool,1} = .!(isnanPlist)  # currently, keeping clusters that could be fit, rather than throwing out entire systems and starting from scratch
            num_pl = sum(keep)
            clusteridlist = clusteridlist[keep]
            Plist = Plist[keep]
            Rlist = Rlist[keep]
            masslist = masslist[keep]
            ecclist = ecclist[keep]
            omegalist = omegalist[keep]
            ascnodelist = ascnodelist[keep]
            meananomlist = meananomlist[keep]
            incllist = incllist[keep]
        end

        valid_system = false

        #println("P: ", Plist)
        #println("R: ", Rlist)
        #println("M: ", masslist)
        #println("ecc: ", ecclist)
        #println("omega: ", omegalist)

        # NOTE: this would be for drawing each cluster separately and then accepting or rejecting the whole lot. By testing for stability before adding each cluster, this last test should be unnecessary.
        if length(Plist) > 0
            if test_stability(Plist, masslist, star.mass, sim_param; ecc=ecclist)
                valid_system = true
            else
                println("Warning: re-attempting system because it fails stability test even though its clusters each pass the test.")
                # NOTE: this should never happen because we check for stability before adding each cluster, and unstable additions are set to NaN and then discarded
            end
        else
            valid_system = true # this else statement is to allow for systems with no planets to pass
        end

        attempts_system += 1
    end # while !valid_system...

    if attempts_system > 1
        println("attempts_system: ", attempts_system)
    end

    # To print out periods, radii, and masses (for troubleshooting):
    #=
    i_sort = sortperm(Plist)
    Plist_sorted = sort(Plist)
    if length(Plist) > 1
        ratio_list = Plist_sorted[2:end]./Plist_sorted[1:end-1]
        if minimum(ratio_list) < 1.1
            println("P: ", Plist_sorted)
            println(Rlist[i_sort])
            println(masslist[i_sort])
            println(ecclist[i_sort])
            println(omegalist[i_sort])
        end
    end
    =#

    pl = Array{Planet}(undef, num_pl)
    orbit = Array{Orbit}(undef, num_pl)
    idx = sortperm(Plist) # TODO OPT: Check to see if sorting is significant time sink. If so, could reduce redundant sortperm
    for i in 1:num_pl
        orbit[i] = Orbit(Plist[idx[i]], ecclist[idx[i]], incllist[idx[i]], omegalist[idx[i]], ascnodelist[idx[i]], meananomlist[idx[i]])
        pl[i] = Planet(Rlist[idx[i]], masslist[idx[i]], clusteridlist[idx[i]])
    end

    return PlanetarySystem(star, pl, orbit)
end

function generate_planetary_system_non_clustered(star::StarAbstract, sim_param::SimParam; verbose::Bool=false)
    # Load functions to use for drawing parameters:
    lambda::Float64 = exp(get_real(sim_param, "log_rate_clusters"))
    max_clusters_in_sys::Int64 = get_int(sim_param, "max_clusters_in_sys")

    generate_periods = get_function(sim_param, "generate_periods")
    generate_sizes = get_function(sim_param, "generate_sizes")
    generate_planet_mass_from_radius = get_function(sim_param, "generate_planet_mass_from_radius")
    generate_e_omega = get_function(sim_param, "generate_e_omega")
    sigma_ecc::Float64 = haskey(sim_param, "sigma_hk") ? get_real(sim_param, "sigma_hk") : 0.0

    # Generate a set of periods, planet radii, and planet masses:
    attempts_system = 0
    max_attempts_system = 100 # Note: Should this be a parameter?
    local num_pl, Plist, Rlist, masslist, ecclist, omegalist
    valid_system = false
    while !valid_system && attempts_system < max_attempts_system

        num_pl = ExoplanetsSysSim.generate_num_planets_poisson(lambda, max_clusters_in_sys)

        if num_pl==0
            valid_system = true
            return PlanetarySystem(star)
        end

        if num_pl==1
            valid_system = true
            Plist = generate_periods(star, sim_param)
            Rlist = generate_sizes(star, sim_param)
            masslist = [generate_planet_mass_from_radius(Rlist[1], sim_param)]
            eccentricity::Float64, omega_e::Float64 = generate_e_omega(sigma_ecc)::Tuple{Float64,Float64}
            ecclist, omegalist = [eccentricity], [omega_e]
        elseif num_pl>1
            Rlist::Array{Float64,1} = generate_sizes(star,sim_param, num_pl=num_pl)
            masslist::Array{Float64,1} = map(r -> generate_planet_mass_from_radius(r, sim_param), Rlist)

            ecclist::Array{Float64,1} = Array{Float64}(undef, num_pl)
            omegalist::Array{Float64,1} = Array{Float64}(undef, num_pl)
            for i in 1:num_pl
                (ecclist[i], omegalist[i]) = generate_e_omega(sigma_ecc)::Tuple{Float64,Float64}
            end

            Plist::Array{Float64,1} = Array{Float64}(undef, num_pl)

            found_good_periods = false
            attempts = 0
            max_attempts = 100 # Note: Should this be a parameter?
            while !found_good_periods && attempts < max_attempts
                attempts += 1
                Plist[1:num_pl] = generate_periods(star, sim_param; num_pl=num_pl)
                if test_stability(Plist, masslist, star.mass, sim_param; ecc=ecclist)
                    found_good_periods = true
                end
            end # while trying to draw periods

            #println("attempts for periods: ", attempts)

            if !found_good_periods
                #println("# Warning: Did not find a good set of periods, given sizes and masses for system. Re-attempting entire system...")
                attempts_system += 1
            else
                valid_system = true
            end
        end

    end # while !valid_system...

    if attempts_system > 0
        #println("attempts_system: ", attempts_system)
        if attempts_system == max_attempts_system
            println("Failed to generate a valid system after ", attempts_system, " attempts; returning just the star.")
            return PlanetarySystem(star)
        end
    end

    # To print out periods, radii, and masses (for troubleshooting):
    #=
    i_sort = sortperm(Plist)
    Plist_sorted = sort(Plist)
    if length(Plist) > 1
        ratio_list = Plist_sorted[2:end]./Plist_sorted[1:end-1]
        if minimum(ratio_list) < 1.1
            println("P: ", Plist_sorted)
            println(Rlist[i_sort])
            println(masslist[i_sort])
            println(ecclist[i_sort])
            println(omegalist[i_sort])
        end
    end
    =#

    # Now assign orbits with given periods, sizes, and masses.
    sigma_incl = deg2rad(get_real(sim_param, "sigma_incl"))
    sigma_incl_near_mmr = deg2rad(get_real(sim_param, "sigma_incl_near_mmr"))
    max_incl_sys = get_real(sim_param, "max_incl_sys")
    f_high_incl = get_real(sim_param, "f_high_incl")

    sigma_incl_use = rand() < f_high_incl ? max(sigma_incl, sigma_incl_near_mmr) : min(sigma_incl, sigma_incl_near_mmr)

    pl = Array{Planet}(undef, num_pl)
    orbit = Array{Orbit}(undef, num_pl)
    incl_sys = acos(cos(max_incl_sys*pi/180)*rand()) #acos(rand()) for isotropic distribution of system inclinations; acos(cos(X*pi/180)*rand()) gives angles from X (deg) to 90 (deg)
    idx = sortperm(Plist)       # TODO OPT: Check to see if sorting is significant time sink.  If so, could reduce redundant sortperm
    is_near_resonance = calc_if_near_resonance(Plist[idx], sim_param)
    for i in 1:num_pl
        #=
        if haskey(sim_param,"sigma_hk_one") && haskey(sim_param, "sigma_hk_multi")
            sigma_ecc = num_pl == 1 ? get_real(sim_param, "sigma_hk_one") : get_real(sim_param, "sigma_hk_multi")
        end
        =#

        #####sigma_incl_use = is_near_resonance[i] ? sigma_incl_near_mmr : sigma_incl
        incl_mut = sigma_incl_use*sqrt(randn()^2+randn()^2) # rand(Distributions.Rayleigh(sigma_incl_use))
        #
        if is_near_resonance[i]
            incl_mut = min(sigma_incl, sigma_incl_near_mmr)*sqrt(randn()^2+randn()^2)
        end
        #
        asc_node = 2pi*rand()
        mean_anom = 2pi*rand()
        incl = incl_mut!=zero(incl_mut) ? acos(cos(incl_sys)*cos(incl_mut) + sin(incl_sys)*sin(incl_mut)*cos(asc_node)) : incl_sys
        orbit[i] = Orbit(Plist[idx[i]], ecclist[idx[i]], incl, omegalist[idx[i]], asc_node, mean_anom)
        pl[i] = Planet(Rlist[idx[i]], masslist[idx[i]])
    end # for i in 1:num_pl

    return PlanetarySystem(star, pl, orbit)
end

function test_generate_planetary_system_clustered() # TODO: Update to test to not reply on generate_kepler_physical_catalog
    sim_param = setup_sim_param_model()
    cat_phys = generate_kepler_physical_catalog(sim_param)
end
