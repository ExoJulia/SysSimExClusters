if !@isdefined ExoplanetsSysSim
    using ExoplanetsSysSim
end

#=
# TODO: OPT: Try saving extracted params to a struct?
struct generate_planet_plan
generate_planet_mass_from_radius = get_function(sim_param, "generate_planet_mass_from_radius")
generate_sizes = get_function(sim_param, "generate_sizes")
generate_e_omega = get_function(sim_param, "generate_e_omega")
sigma_ecc::Float64 = haskey(sim_param, "sigma_hk") ? get_real(sim_param, "sigma_hk") : 0.0
sigma_log_radius_in_cluster = get_real(sim_param, "sigma_log_radius_in_cluster")
min_radius::Float64 = get_real(sim_param, "min_radius")
max_radius::Float64 = get_real(sim_param, "max_radius")
sigma_logperiod_per_pl_in_cluster = get_real(sim_param, "sigma_logperiod_per_pl_in_cluster")
min_period = get_real(sim_param, "min_period")
max_period = get_real(sim_param, "max_period")
log_mean_P = 0.0 # log(generate_periods_power_law(star,sim_param))
min_num_mutual_hill_radii = get_real(sim_param, "num_mutual_hill_radii")
end
=#

function generate_planet_periods_sizes_masses_eccs_in_cluster(star::StarT, sim_param::SimParam; n::Int64=1) where {StarT<:StarAbstract}  # TODO: IMPORTANT: Make this function work and test before using for science
    @assert n>=1
    generate_planet_mass_from_radius = get_function(sim_param, "generate_planet_mass_from_radius")
    generate_sizes = get_function(sim_param, "generate_sizes")
    generate_e_omega = get_function(sim_param, "generate_e_omega")
    sigma_ecc::Float64 = haskey(sim_param, "sigma_hk") ? get_real(sim_param, "sigma_hk") : 0.0

    if n==1
        R = generate_sizes(star, sim_param)
        mass = [generate_planet_mass_from_radius(R[1], sim_param)]
        P = [1.0] # generate_periods_power_law(star,sim_param)
        eccentricity::Float64, omega_e::Float64 = generate_e_omega(sigma_ecc)::Tuple{Float64,Float64}
        ecc, omega = [eccentricity], [omega_e]
        return (P, R, mass, ecc, omega)
    end

    #If reach here, then at least 2 planets in cluster

    mean_R = generate_sizes(star,sim_param)[1]
    sigma_log_radius_in_cluster = get_real(sim_param, "sigma_log_radius_in_cluster")
    #println("# mean_R = ",mean_R," sigma_log_radius_in_cluster= ",sigma_log_radius_in_cluster)

    min_radius::Float64 = get_real(sim_param, "min_radius")
    max_radius::Float64 = get_real(sim_param, "max_radius")
    Rdist = Truncated(LogNormal(log(mean_R),sigma_log_radius_in_cluster), min_radius, max_radius) #if we want clustered planet sizes
    R = rand(Rdist, n)
    #R = generate_sizes(star,sim_param, num_pl=n) #####If want non-clustered planet sizes instead

    #println("# Rp = ", R)
    mass = map(r -> generate_planet_mass_from_radius(r, sim_param), R)
    #println("# mass = ", mass)

    local ecc, omega
    ecc::Array{Float64,1} = Array{Float64}(undef, n)
    omega::Array{Float64,1} = Array{Float64}(undef, n)
    for i in 1:n
        (ecc[i], omega[i]) = generate_e_omega(sigma_ecc)::Tuple{Float64,Float64}
    end

    sigma_logperiod_per_pl_in_cluster = get_real(sim_param, "sigma_logperiod_per_pl_in_cluster")
    min_period = get_real(sim_param, "min_period")
    max_period = get_real(sim_param, "max_period")
    max_period_ratio = max_period/min_period
    log_mean_P = 0.0 # log(generate_periods_power_law(star,sim_param))
    local P
    #
    P = zeros(n)
    for i in 1:n # Draw periods one at a time
        if any(isnan.(P))
            P[i:end] .= NaN
            println("Cannot fit any more planets in cluster; returning planets that did fit.")
            break
        end
        P[i] = draw_period_lognormal_allowed_regions(P[1:i-1], mass[1:i-1], star.mass, sim_param; μ=log_mean_P, σ=n*sigma_logperiod_per_pl_in_cluster, x_min=1/sqrt(max_period_ratio), x_max=sqrt(max_period_ratio), ecc=ecc[1:i-1])
    end
    #

    #=
    # Note: Currently, drawing all periods within a cluster at once and either keeping or rejecting the whole cluster
    #       Should we instead draw periods one at a time?
    Pdist = Truncated(LogNormal(log_mean_P,sigma_logperiod_per_pl_in_cluster*n), 1/sqrt(max_period_ratio), sqrt(max_period_ratio)) #Truncated unscaled period distribution to ensure that the cluster can fit in the period range [min_period, max_period] after scaling by a period scale
    found_good_periods = false
    attempts = 0
    max_attempts = 100  # Note: Should this be a parameter?
    while !found_good_periods && attempts < max_attempts
        attempts += 1
        P = rand(Pdist, n)
        if test_stability(P, mass, star.mass, sim_param; ecc=ecc)
            found_good_periods = true
        end
    end # while trying to draw periods
    #println("attempts for cluster: ", attempts)

    if !found_good_periods
        #println("# Warning: Did not find a good set of periods, sizes and masses for one cluster.")
        return (fill(NaN,n), R, mass, ecc, omega)  # Return NaNs for periods to indicate failed
    end
    =#

    return (P, R, mass, ecc, omega)    # Note can also return earlier if only one planet in cluster or if fail to generate a good set of values
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

function generate_planetary_system_clustered(star::StarAbstract, sim_param::SimParam; verbose::Bool=false)  # TODO: Make this function work and test before using for science
    # Load functions to use for drawing parameters:
    if haskey(sim_param, "f_stars_with_planets_attempted")
        f_stars_with_planets_attempted = get_real(sim_param, "f_stars_with_planets_attempted")
        @assert 0<=f_stars_with_planets_attempted<=1
    else
        f_stars_with_planets_attempted = 1.
    end
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

    # Generate a set of periods, planet radii, and planet masses:
    attempt_system = 0
    max_attempts_system = 1 #Note: currently this should not matter; each system is always attempted just once
    local num_pl, clusteridlist, Plist, Rlist, masslist, ecclist, omegalist
    valid_system = false
    while !valid_system && attempt_system < max_attempts_system

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

        @assert num_pl_in_cluster[1] >= 1
        pl_start = 1
        pl_stop = 0
        for c in 1:num_clusters
            pl_stop += num_pl_in_cluster[c]
            clusteridlist[pl_start:pl_stop] = ones(Int64, num_pl_in_cluster[c])*c
            Plist_tmp::Array{Float64,1}, Rlist_tmp::Array{Float64,1}, masslist_tmp::Array{Float64,1}, ecclist_tmp::Array{Float64,1}, omegalist_tmp::Array{Float64,1} = generate_planet_periods_sizes_masses_eccs_in_cluster(star, sim_param, n=num_pl_in_cluster[c])
            Rlist[pl_start:pl_stop], masslist[pl_start:pl_stop], ecclist[pl_start:pl_stop], omegalist[pl_start:pl_stop] = Rlist_tmp, masslist_tmp, ecclist_tmp, omegalist_tmp
            valid_cluster = !any(isnan.(Plist_tmp)) #should this be looking for nans in Plist or Plist_tmp? Was Plist but I think it should be Plist_tmp!
            valid_period_scale = false
            attempt_period_scale = 0
            max_attempts_period_scale = 100
            while !valid_period_scale && attempt_period_scale<max_attempts_period_scale && valid_cluster
                attempt_period_scale += 1

                period_scale::Array{Float64,1} = draw_power_law(power_law_P, min_period/minimum(Plist_tmp), max_period/maximum(Plist_tmp), 1)
                #Note: this ensures that the minimum and maximum periods will be in the range [min_period, max_period]
                #Warning: not sure about the behaviour when min_period/minimum(Plist_tmp) > max_period/maximum(Plist_tmp) (i.e. when the cluster cannot fit in the given range)?
                #TODO OPT: could draw period_scale more efficiently by computing the allowed regions in [min_period, max_period] given the previous cluster draws

                Plist[pl_start:pl_stop] = Plist_tmp .* period_scale

                #if test_stability(Plist[1:pl_stop], masslist[1:pl_stop], star.mass, sim_param; ecc=ecclist[1:pl_stop])
                if test_stability(view(Plist,1:pl_stop), view(masslist,1:pl_stop), star.mass, sim_param; ecc=view(ecclist,1:pl_stop))
                    valid_period_scale = true
                end
            end  # while !valid_period_scale...

            #if attempt_period_scale > 1
                #println("attempts_period_scale: ", attempt_period_scale)
            #end

            if !valid_period_scale
                Plist[pl_start:pl_stop] .= NaN
            end
            pl_start += num_pl_in_cluster[c]
        end # for c in 1:num_clusters
        isnanPlist::Array{Bool,1} = isnan.(Plist::Array{Float64,1})
        if any(isnanPlist)  # If any loop failed to generate valid planets, it should set a NaN in the period list
            keep::Array{Bool,1} = .!(isnanPlist)  # Currently, keeping clusters that could be fit, rather than throwing out entirely and starting from scratch.  Is this a good idea?  Matthias tried the other approach in his python code.
            num_pl = sum(keep)
            clusteridlist = clusteridlist[keep]
            Plist = Plist[keep]
            Rlist = Rlist[keep]
            masslist = masslist[keep]
            ecclist = ecclist[keep]
            omegalist = omegalist[keep]
        end

        valid_system = false

        #println("P: ", Plist)
        #println("R: ", Rlist)
        #println("M: ", masslist)
        #println("ecc: ", ecclist)
        #println("omega: ", omegalist)

        #Note: This would be for drawing each cluster separately and then accepting or rejecting the whole lot. By testing for stability before adding each cluster, this last test should be unnecessary.
        if length(Plist) > 0
            if test_stability(Plist, masslist, star.mass, sim_param; ecc=ecclist)
                valid_system = true
            else
                println("Warning: re-attempting system because it fails stability test even though its clusters each pass the test.")
                #Note: this should never happen because we check for stability before adding each cluster, and unstable additions are set to NaN and then discarded
            end
        else
            valid_system = true #this else statement is to allow for systems with no planets to pass
        end

        attempt_system += 1
    end # while !valid_system...

    if attempt_system > 1
        println("attempt_system: ", attempt_system)
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
        pl[i] = Planet(Rlist[idx[i]], masslist[idx[i]], clusteridlist[idx[i]])
    end # for i in 1:num_pl

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
    attempt_system = 0
    max_attempts_system = 100 # Note: Should this be a parameter?
    local num_pl, Plist, Rlist, masslist, ecclist, omegalist
    valid_system = false
    while !valid_system && attempt_system < max_attempts_system

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
                attempt_system += 1
            else
                valid_system = true
            end
        end

    end # while !valid_system...

    if attempt_system > 0
        #println("attempt_system: ", attempt_system)
        if attempt_system == max_attempts_system
            println("Failed to generate a valid system after ", attempt_system, " attempts; returning just the star.")
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
