import DataFrames.skipmissing





##### To define functions for calculating the distances:

function calc_distance(ss1::ExoplanetsSysSim.CatalogSummaryStatistics, ss2::ExoplanetsSysSim.CatalogSummaryStatistics, return_KS_or_AD::String ; AD_mod::Bool=false, all_dist::Bool=false, save_dist::Bool=true)
    #This function calculates the total KS distance between two summary statistics (simulated observed catalogs).
    #If 'all_dist=true', the function outputs the individual distances in the distance function.
    #If 'save_dist=true', the function also saves the distances (individual and total) to a file (assuming file 'f' is open for writing).

    M_cat_obs1 = ones(Int64,0) #array to be filled with the number of transiting planets in each simulated system for ss1
    M_cat_obs2 = ones(Int64,0) #array to be filled with the number of transiting planets in each simulated system for ss2
    max_k = get_int(sim_param,"max_tranets_in_sys")
    for k in 1:max_k
        append!(M_cat_obs1, k*ones(Int64, ss1.stat["num n-tranet systems"][k]))
        append!(M_cat_obs2, k*ones(Int64, ss2.stat["num n-tranet systems"][k]))
    end
    Nmult_obs1 = ss1.stat["num n-tranet systems"]
    Nmult_obs2 = ss2.stat["num n-tranet systems"]
    planets_obs1, pairs_obs1 = sum(Nmult_obs1 .* collect(1:max_k)), sum(Nmult_obs1[2:end] .* collect(1:max_k-1)) #total number of planets, and planet pairs
    planets_obs2, pairs_obs2 = sum(Nmult_obs2 .* collect(1:max_k)), sum(Nmult_obs2[2:end] .* collect(1:max_k-1)) #total number of planets, and planet pairs

    #To handle empty arrays:
    if min(length(ss1.stat["radius_ratio_above_list"]), length(ss1.stat["radius_ratio_below_list"]), length(ss1.stat["radius_ratio_across_list"]), length(ss1.stat["duration_ratio_non_mmr_list"]), length(ss1.stat["duration_ratio_near_mmr_list"])) < 2 || min(length(ss2.stat["radius_ratio_above_list"]), length(ss2.stat["radius_ratio_below_list"]), length(ss2.stat["radius_ratio_across_list"]), length(ss2.stat["duration_ratio_non_mmr_list"]), length(ss2.stat["duration_ratio_near_mmr_list"])) < 2 #need at least 2 elements in each of these summary statistics per catalog in order to be able to compute AD distances
        println("Not enough observed multi-planet systems in one of the catalogs to compute the AD distance.")
        d = ones(Int64,17)*1e6

        println("Distances: ", d, [sum(d)])
        if save_dist
            println(f, "Mult_counts: ", Nmult_obs1, [planets_obs1, pairs_obs1], Nmult_obs2, [planets_obs2, pairs_obs2])
            println(f, "Dist_KS: ", d, [sum(d)])
            println(f, "Dist_AD: ", d, [sum(d)])
        end

        if all_dist
            if return_KS_or_AD == "KS" || return_KS_or_AD == "AD"
                return d
            elseif return_KS_or_AD == "Both"
                return [d; d[3:end]]
            end
        else
            if return_KS_or_AD == "KS" || return_KS_or_AD == "AD"
                return sum(d)
            elseif return_KS_or_AD == "Both"
                return sum([d; d[3:end]])
            end
        end
    end

    #If there are enough simulated planets to compute distances:
    max_incl_sys = get_real(sim_param,"max_incl_sys")
    cos_factor = cos(max_incl_sys*pi/180) #factor to divide the number of targets in simulation by to get the actual number of targets needed (with an isotropic distribution of system inclinations) to produce as many transiting systems for a single observer

    delta_f = abs(ss1.stat["num_tranets"]/(ss1.stat["num targets"]/cos_factor) - ss2.stat["num_tranets"]/(ss2.stat["num targets"]))
    d_mult_KS = ksstats_ints(M_cat_obs1, M_cat_obs2)[5]
    d_mult_CRPD = CRPDstats([Nmult_obs1[1:4]; sum(Nmult_obs1[5:end])], [Nmult_obs2[1:4]; sum(Nmult_obs2[5:end])])
    d_mult_CRPD_switched = CRPDstats([Nmult_obs2[1:4]; sum(Nmult_obs2[5:end])], [Nmult_obs1[1:4]; sum(Nmult_obs1[5:end])])

    #To compute the KS distances:
    d_KS = Array{Float64}(undef, 17)
    d_KS[1] = delta_f
    d_KS[2] = d_mult_KS
    d_KS[3] = d_mult_CRPD
    d_KS[4] = d_mult_CRPD_switched
    d_KS[5] = ksstats(ss1.stat["P list"], ss2.stat["P list"])[5]
    d_KS[6] = ksstats(ss1.stat["period_ratio_list"], ss2.stat["period_ratio_list"])[5]
    d_KS[7] = ksstats(ss1.stat["duration list"], ss2.stat["duration list"])[5]
    d_KS[8] = ksstats(ss1.stat["duration_ratio_list"], ss2.stat["duration_ratio_list"])[5]
    d_KS[9] = ksstats(ss1.stat["duration_ratio_non_mmr_list"], ss2.stat["duration_ratio_non_mmr_list"])[5]
    d_KS[10] = ksstats(ss1.stat["duration_ratio_near_mmr_list"], ss2.stat["duration_ratio_near_mmr_list"])[5]
    d_KS[11] = ksstats(ss1.stat["depth list"], ss2.stat["depth list"])[5]
    d_KS[12] = ksstats(ss1.stat["depth above list"], ss2.stat["depth above list"])[5]
    d_KS[13] = ksstats(ss1.stat["depth below list"], ss2.stat["depth below list"])[5]
    d_KS[14] = ksstats(ss1.stat["radius_ratio_list"], ss2.stat["radius_ratio_list"])[5]
    d_KS[15] = ksstats(ss1.stat["radius_ratio_above_list"], ss2.stat["radius_ratio_above_list"])[5]
    d_KS[16] = ksstats(ss1.stat["radius_ratio_below_list"], ss2.stat["radius_ratio_below_list"])[5]
    d_KS[17] = ksstats(ss1.stat["radius_ratio_across_list"], ss2.stat["radius_ratio_across_list"])[5]

    #To compute the AD distances:
    if AD_mod
        ADdist = ADstats_mod
    else
        ADdist = ADstats
    end

    d_AD = Array{Float64}(undef, 17)
    d_AD[1] = delta_f
    d_AD[2] = d_mult_KS
    d_AD[3] = d_mult_CRPD
    d_AD[4] = d_mult_CRPD_switched
    d_AD[5] = ADdist(ss1.stat["P list"], ss2.stat["P list"])
    d_AD[6] = ADdist(ss1.stat["period_ratio_list"], ss2.stat["period_ratio_list"])
    d_AD[7] = ADdist(ss1.stat["duration list"], ss2.stat["duration list"])
    d_AD[8] = ADdist(ss1.stat["duration_ratio_list"], ss2.stat["duration_ratio_list"])
    d_AD[9] = ADdist(ss1.stat["duration_ratio_non_mmr_list"], ss2.stat["duration_ratio_non_mmr_list"])
    d_AD[10] = ADdist(ss1.stat["duration_ratio_near_mmr_list"], ss2.stat["duration_ratio_near_mmr_list"])
    d_AD[11] = ADdist(ss1.stat["depth list"], ss2.stat["depth list"])
    d_AD[12] = ADdist(ss1.stat["depth above list"], ss2.stat["depth above list"])
    d_AD[13] = ADdist(ss1.stat["depth below list"], ss2.stat["depth below list"])
    d_AD[14] = ADdist(ss1.stat["radius_ratio_list"], ss2.stat["radius_ratio_list"])
    d_AD[15] = ADdist(ss1.stat["radius_ratio_above_list"], ss2.stat["radius_ratio_above_list"])
    d_AD[16] = ADdist(ss1.stat["radius_ratio_below_list"], ss2.stat["radius_ratio_below_list"])
    d_AD[17] = ADdist(ss1.stat["radius_ratio_across_list"], ss2.stat["radius_ratio_across_list"])

    #To print and/or write the distances to file:
    println("Multiplicity Counts: ", Nmult_obs1, [planets_obs1, pairs_obs1], Nmult_obs2, [planets_obs2, pairs_obs2])
    println("KS Distances: ", d_KS, [sum(d_KS)])
    println("AD Distances: ", d_AD, [sum(d_AD)])
    if save_dist
        println(f, "Mult_counts: ", Nmult_obs1, [planets_obs1, pairs_obs1], Nmult_obs2, [planets_obs2, pairs_obs2])
        println(f, "Dist_KS: ", d_KS, [sum(d_KS)])
        println(f, "Dist_AD: ", d_AD, [sum(d_AD)])
    end

    #To return the distances or total distance:
    if all_dist
        if return_KS_or_AD == "KS"
            return d_KS
        elseif return_KS_or_AD == "AD"
            return d_AD
        elseif return_KS_or_AD == "Both"
            return [d_KS; d_AD[3:end]]
        end
    else
        if return_KS_or_AD == "KS"
            return sum(d_KS)
        elseif return_KS_or_AD == "AD"
            return sum(d_AD)
        elseif return_KS_or_AD == "Both"
            return sum([d_KS; d_AD[3:end]])
        end
    end
end

function calc_distance_Kepler(ss1::ExoplanetsSysSim.CatalogSummaryStatistics, return_KS_or_AD::String ; AD_mod::Bool=false, all_dist::Bool=false, save_dist::Bool=true)
    #This function calculates the total KS distance between a population generated by our model and the actual Kepler population (which must be loaded in).
    #If 'all_dist=true', the function outputs the individual distances in the distance function.
    #If 'save_dist=true', the function also saves the distances (individual and total) to a file (assuming file 'f' is open for writing).

    M_cat_obs = ones(Int64,0) #array to be filled with the number of transiting planets in each simulated system
    max_k = get_int(sim_param,"max_tranets_in_sys")
    for k in 1:max_k
        append!(M_cat_obs, k*ones(Int64, ss1.stat["num n-tranet systems"][k]))
    end
    Nmult_obs = ss1.stat["num n-tranet systems"]
    planets_obs, pairs_obs = sum(Nmult_obs .* collect(1:max_k)), sum(Nmult_obs[2:end] .* collect(1:max_k-1)) #total number of planets, and planet pairs

    if min(length(ss1.stat["radius_ratio_above_list"]), length(ss1.stat["radius_ratio_below_list"]), length(ss1.stat["radius_ratio_across_list"]), length(ss1.stat["duration_ratio_non_mmr_list"]), length(ss1.stat["duration_ratio_near_mmr_list"])) < 2 #need at least 2 elements in each of these summary statistics in order to be able to compute AD distances
        println("Not enough observed multi-planet systems in the simulated catalog.")
        d = ones(Int64,17)*1e6

        println("Distances: ", d, [sum(d)])
        if save_dist
            println(f, "Mult_counts: ", Nmult_obs, [planets_obs, pairs_obs])
            println(f, "Dist_KS: ", d, [sum(d)])
            println(f, "Dist_AD: ", d, [sum(d)])
        end

        if all_dist
            if return_KS_or_AD == "KS" || return_KS_or_AD == "AD"
                return d
            elseif return_KS_or_AD == "Both"
                return [d; d[3:end]]
            end
        else
            if return_KS_or_AD == "KS" || return_KS_or_AD == "AD"
                return sum(d)
            elseif return_KS_or_AD == "Both"
                return sum([d; d[3:end]])
            end
        end
    end

    max_incl_sys = get_real(sim_param,"max_incl_sys")
    cos_factor = cos(max_incl_sys*pi/180) #factor to divide the number of targets in simulation by to get the actual number of targets needed (with an isotropic distribution of system inclinations) to produce as many transiting systems for a single observer

    delta_f = abs(ss1.stat["num_tranets"]/(ss1.stat["num targets"]/cos_factor) - length(P_confirmed)/N_Kepler_targets)
    d_mult_KS = ksstats_ints(M_cat_obs, M_confirmed)[5]
    d_mult_CRPD = CRPDstats([Nmult_obs[1:4]; sum(Nmult_obs[5:end])], [Nmult_confirmed[1:4]; sum(Nmult_confirmed[5:end])])
    d_mult_CRPD_switched = CRPDstats([Nmult_confirmed[1:4]; sum(Nmult_confirmed[5:end])], [Nmult_obs[1:4]; sum(Nmult_obs[5:end])])

    #To compute the KS distances:
    d_KS = Array{Float64}(undef, 17)
    d_KS[1] = delta_f
    d_KS[2] = d_mult_KS
    d_KS[3] = d_mult_CRPD
    d_KS[4] = d_mult_CRPD_switched
    d_KS[5] = ksstats(ss1.stat["P list"], P_confirmed)[5]
    d_KS[6] = ksstats(ss1.stat["period_ratio_list"], R_confirmed)[5]
    d_KS[7] = ksstats(ss1.stat["duration list"].*24, t_D_confirmed)[5] #transit durations in simulations are in days, while in the Kepler catalog are in hours
    d_KS[8] = ksstats(ss1.stat["duration_ratio_list"], xi_confirmed)[5]
    d_KS[9] = ksstats(ss1.stat["duration_ratio_non_mmr_list"], xi_non_mmr_confirmed)[5]
    d_KS[10] = ksstats(ss1.stat["duration_ratio_near_mmr_list"], xi_near_mmr_confirmed)[5]
    d_KS[11] = ksstats(ss1.stat["depth list"], D_confirmed)[5]
    d_KS[12] = ksstats(ss1.stat["depth above list"], D_above_confirmed)[5]
    d_KS[13] = ksstats(ss1.stat["depth below list"], D_below_confirmed)[5]
    d_KS[14] = ksstats(ss1.stat["radius_ratio_list"].^2, D_ratio_confirmed)[5] #simulations save radius ratios while we computed transit duration ratios from the Kepler catalog
    d_KS[15] = ksstats(ss1.stat["radius_ratio_above_list"].^2, D_ratio_above_confirmed)[5]
    d_KS[16] = ksstats(ss1.stat["radius_ratio_below_list"].^2, D_ratio_below_confirmed)[5]
    d_KS[17] = ksstats(ss1.stat["radius_ratio_across_list"].^2, D_ratio_across_confirmed)[5]

    #To compute the AD distances:
    if AD_mod
        ADdist = ADstats_mod
    else
        ADdist = ADstats
    end

    d_AD = Array{Float64}(undef, 17)
    d_AD[1] = delta_f
    d_AD[2] = d_mult_KS
    d_AD[3] = d_mult_CRPD
    d_AD[4] = d_mult_CRPD_switched
    d_AD[5] = ADdist(ss1.stat["P list"], P_confirmed)
    d_AD[6] = ADdist(ss1.stat["period_ratio_list"], R_confirmed)
    d_AD[7] = ADdist(ss1.stat["duration list"].*24, t_D_confirmed) #transit durations in simulations are in days, while in the Kepler catalog are in hours
    d_AD[8] = ADdist(ss1.stat["duration_ratio_list"], xi_confirmed)
    d_AD[9] = ADdist(ss1.stat["duration_ratio_non_mmr_list"], xi_non_mmr_confirmed)
    d_AD[10] = ADdist(ss1.stat["duration_ratio_near_mmr_list"], xi_near_mmr_confirmed)
    d_AD[11] = ADdist(ss1.stat["depth list"], D_confirmed)
    d_AD[12] = ADdist(ss1.stat["depth above list"], D_above_confirmed)
    d_AD[13] = ADdist(ss1.stat["depth below list"], D_below_confirmed)
    d_AD[14] = ADdist(ss1.stat["radius_ratio_list"].^2, D_ratio_confirmed) #simulations save radius ratios while we computed transit duration ratios from the Kepler catalog
    d_AD[15] = ADdist(ss1.stat["radius_ratio_above_list"].^2, D_ratio_above_confirmed)
    d_AD[16] = ADdist(ss1.stat["radius_ratio_below_list"].^2, D_ratio_below_confirmed)
    d_AD[17] = ADdist(ss1.stat["radius_ratio_across_list"].^2, D_ratio_across_confirmed)

    #To print and/or write the distances to file:
    println("Multiplicity Counts: ", Nmult_obs, [planets_obs, pairs_obs])
    println("KS Distances: ", d_KS, [sum(d_KS)])
    println("AD Distances: ", d_AD, [sum(d_AD)])
    if save_dist
        println(f, "Mult_counts: ", Nmult_obs, [planets_obs, pairs_obs])
        println(f, "Dist_KS: ", d_KS, [sum(d_KS)])
        println(f, "Dist_AD: ", d_AD, [sum(d_AD)])
    end

    #To return the distances or total distance:
    if all_dist
        if return_KS_or_AD == "KS"
            return d_KS
        elseif return_KS_or_AD == "AD"
            return d_AD
        elseif return_KS_or_AD == "Both"
            return [d_KS; d_AD[3:end]]
        end
    else
        if return_KS_or_AD == "KS"
            return sum(d_KS)
        elseif return_KS_or_AD == "AD"
            return sum(d_AD)
        elseif return_KS_or_AD == "Both"
            return sum([d_KS; d_AD[3:end]])
        end
    end
end

function target_function(active_param::Vector{Float64}, use_KS_or_AD::String, Kep_or_Sim::String ; AD_mod::Bool=false, weights::Vector{Float64}=ones(17), all_dist::Bool=false, save_dist::Bool=true)
    #This function takes in the values of the active model parameters, generates a simulated observed catalog, and computes the distance function.
    #If 'all_dist=true', the function outputs the individual distances in the distance function.
    #If 'save_dist=true', the function also saves the distances (unweighted and weighted, individual and total) to a file (assuming file 'f' is open for writing).

    println("Active parameter values: ", active_param)
    if save_dist
        println(f, "Active_params: ", active_param) #to write the params to file
    end

    global sim_param
    sim_param_here = deepcopy(sim_param)
    ExoplanetsSysSim.update_sim_param_from_vector!(active_param,sim_param_here)
    cat_phys = generate_kepler_physical_catalog(sim_param_here)
    cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param_here)
    cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param_here)
    summary_stat = calc_summary_stats_model(cat_obs,sim_param_here)

    if Kep_or_Sim == "Kep"
        dist = calc_distance_Kepler(summary_stat, use_KS_or_AD; AD_mod=AD_mod, all_dist=true, save_dist=save_dist)
    elseif Kep_or_Sim == "Sim"
        global summary_stat_ref
        dist = calc_distance(summary_stat, summary_stat_ref, use_KS_or_AD; AD_mod=AD_mod, all_dist=true, save_dist=save_dist)
    end

    weighted_dist = dist .* weights

    if weights != ones(17)
        println("Weighted distances: ", weighted_dist, [sum(weighted_dist)])
        if save_dist
            println(f, "Dist_weighted: ", weighted_dist, [sum(weighted_dist)])
        end
    end

    if all_dist
        return weighted_dist
    else
        return sum(weighted_dist)
    end
end

function compute_weights_target_fitness_std_perfect_model(num_evals::Int64, use_KS_or_AD::String ; AD_mod::Bool=false, weight::Bool=true, dists_exclude::Vector{Int64}=Int64[], save_dist::Bool=true)
    t_elapsed = @elapsed begin
        active_param_true = make_vector_of_sim_param(sim_param)
        println("# True values: ", active_param_true)
        if save_dist
            println(f, "# Format: Dist: [distances][total distance]")
        end

        dists_true = zeros(num_evals,17)
        for i in 1:num_evals
            dists_true[i,:] = target_function(active_param_true, use_KS_or_AD, "Sim"; AD_mod=AD_mod, all_dist=true, save_dist=save_dist)
        end
        mean_dists = transpose(mean(dists_true, dims=1))[:,] #array of mean distances for each individual distance
        mean_dist = mean(sum(dists_true, dims=2)) #mean total distance
        #std_dists = transpose(std(dists_true, 1))[:,] #array of std distances for each individual distance
        rms_dists = transpose(sqrt.(mean(dists_true .^2, dims=1)))[:,] #array of rms (std around 0)  distances for each individual distance
        std_dist = std(sum(dists_true, dims=2)) #std of total distance
        rms_dist = sqrt(mean(sum(dists_true, dims=2) .^2)) #rms (std around 0) of total distance; should be similar to the mean total distance

        if weight
            weights = 1 ./ rms_dists #to use the array 'rms_dists' as the weights for the individual distances
        else
            weights = ones(17)
        end
        weights[dists_exclude] .= 0. #to exclude certain distances from being used in the total distance function (but still computing and saving them during the optimization) by setting their weights to zero

        weighted_dists_true = zeros(num_evals,length(weights))
        for i in 1:num_evals
            weighted_dists_true[i,:] = dists_true[i,:]  .* weights
        end
        mean_weighted_dists = transpose(mean(weighted_dists_true, dims=1))[:,] #array of mean weighted distances for each individual distance
        mean_weighted_dist = mean(sum(weighted_dists_true, dims=2)) #mean weighted total distance
        std_weighted_dist = std(sum(weighted_dists_true, dims=2)) #std of weighted total distance
    end

    println("Mean dists: ", mean_dists)
    println("Rms dists: ", rms_dists)
    println("Weights (1/rms dists): ", weights)
    println("Mean weighted dists: ", mean_weighted_dists)
    println("Distance using true values: ", mean_dist, " +/- ", std_dist)
    println("Weighted distance using true values: ", mean_weighted_dist, " +/- ", std_weighted_dist)
    if save_dist
        println(f, "#")
        println(f, "Mean: ", mean_dists, [mean_dist])
        println(f, "Rms: ", rms_dists, [rms_dist])
        println(f, "Weights (1/rms dists): ", weights)
        println(f, "Mean weighted dists: ", mean_weighted_dists, [mean_weighted_dist])
        println(f, "# Distance using true values (default parameter values): ", mean_dist, " +/- ", std_dist)
        println(f, "# Weighted distance using true values (default parameter values): ", mean_weighted_dist, " +/- ", std_weighted_dist)
        println(f, "# elapsed time: ", t_elapsed, " seconds")
        println(f, "#")
    end

    return (active_param_true, weights, mean_weighted_dist, std_weighted_dist)
end

function compute_weights_target_fitness_std_from_file(file_name::String, use_KS_or_AD::String ; weight::Bool=true, dists_exclude::Vector{Int64}=Int64[], save_dist::Bool=true)
    t_elapsed = @elapsed begin
        active_param_true = make_vector_of_sim_param(sim_param)
        println("# True values: ", active_param_true)

        dists_true = zeros(1000,17)
        open(file_name) do f_weights
            eval_count = 1
            for (i,line) in enumerate(eachline(f_weights))
                if length(line) > 8
                    if line[1:7] == "Dist_"*use_KS_or_AD
                        dists_str, sum_str = split(line[11:end-1], "][")
                        dists_true[eval_count,:] = [parse(Float64, x) for x in split(dists_str, ", ")]
                        eval_count += 1
                    end
                end
            end
        end

        mean_dists = transpose(mean(dists_true, dims=1))[:,] #array of mean distances for each individual distance
        mean_dist = mean(sum(dists_true, dims=2)) #mean total distance
        #std_dists = transpose(std(dists_true, 1))[:,] #array of std distances for each individual distance
        rms_dists = transpose(sqrt.(mean(dists_true .^2, dims=1)))[:,] #array of rms (std around 0)  distances for each individual distance
        std_dist = std(sum(dists_true, dims=2)) #std of total distance
        rms_dist = sqrt(mean(sum(dists_true, dims=2) .^2)) #rms (std around 0) of total distance; should be similar to the mean total distance

        if weight
            weights = 1 ./ rms_dists #to use the array 'rms_dists' as the weights for the individual distances
        else
            weights = ones(17)
        end
        weights[dists_exclude] .= 0. #to exclude certain distances from being used in the total distance function (but still computing and saving them during the optimization) by setting their weights to zero

        weighted_dists_true = zeros(1000,length(weights))
        for i in 1:1000
            weighted_dists_true[i,:] = dists_true[i,:]  .* weights
        end
        mean_weighted_dists = transpose(mean(weighted_dists_true, dims=1))[:,] #array of mean weighted distances for each individual distance
        mean_weighted_dist = mean(sum(weighted_dists_true, dims=2)) #mean weighted total distance
        std_weighted_dist = std(sum(weighted_dists_true, dims=2)) #std of weighted total distance
    end

    println("# Using weights from pre-computed file:")
    println("Mean dists: ", mean_dists)
    println("Rms dists: ", rms_dists)
    println("Weights (1/rms dists): ", weights)
    println("Mean weighted dists: ", mean_weighted_dists)
    println("Distance using true values: ", mean_dist, " +/- ", std_dist)
    println("Weighted distance using true values: ", mean_weighted_dist, " +/- ", std_weighted_dist)
    if save_dist
        println(f, "#")
        println(f, "# Using weights from pre-computed file:")
        println(f, "Mean: ", mean_dists, [mean_dist])
        println(f, "Rms: ", rms_dists, [rms_dist])
        println(f, "Weights (1/rms dists): ", weights)
        println(f, "Mean weighted dists: ", mean_weighted_dists, [mean_weighted_dist])
        println(f, "# Distance using true values (default parameter values): ", mean_dist, " +/- ", std_dist)
        println(f, "# Weighted distance using true values (default parameter values): ", mean_weighted_dist, " +/- ", std_weighted_dist)
        println(f, "# elapsed time: ", t_elapsed, " seconds")
        println(f, "#")
    end

    return (active_param_true, weights, mean_weighted_dist, std_weighted_dist)
end

function map_square_to_triangle(r1::Float64, r2::Float64, A::Vector{Float64}, B::Vector{Float64}, C::Vector{Float64})
    #This function takes in a point (r1,r2) in the unit square (i.e. r1,r2 in [0,1]) and maps it to a point P=(x,y) in the triangle defined by vertices A,B,C
    #If r1,r2 are uniformly drawn in [0,1], then the point P=(x,y) is also uniformly drawn in the triangle; see http://www.cs.princeton.edu/~funk/tog02.pdf (Section 4.2) for a reference

    @assert 0. <= r1 <= 1.
    @assert 0. <= r2 <= 1.
    P = (1. - sqrt(r1)) .* A + (sqrt(r1)*(1. - r2)) .* B + (sqrt(r1)*r2) .* C
    return P
end

function target_function_transformed_params(active_param_transformed::Vector{Float64}, transformed_indices::Vector{Int64}, A::Vector{Float64}, B::Vector{Float64}, C::Vector{Float64}, use_KS_or_AD::String, Kep_or_Sim::String ; AD_mod::Bool=false, weights::Vector{Float64}=ones(17), all_dist::Bool=false, save_dist::Bool=true)

    r1, r2 = active_param_transformed[transformed_indices]
    @assert 0. <= r1 <= 1.
    @assert 0. <= r2 <= 1.
    active_param = deepcopy(active_param_transformed)
    active_param[transformed_indices] = map_square_to_triangle(r1, r2, A, B, C)

    return target_function(active_param, use_KS_or_AD, Kep_or_Sim ; AD_mod=AD_mod, weights=weights, all_dist=all_dist, save_dist=save_dist)
end
