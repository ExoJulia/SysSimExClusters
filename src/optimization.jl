using Distributed # just to compile a function
using ParallelDataTransfer # just to compile a function



"""
    draw_random_active_params(active_param_keys, active_params_box, sim_param; PT_all=ExoplanetsSysSim.ParamsTriangle[])

Draw a set of active parameter values randomly (uniformly in the given box). If there are any triangle-transformed pairs of parameters (`PT_all` is not empty), those parameters are drawn uniformly in their respective triangles. NOTE: also rewrites the active parameter values in `sim_param` with the drawn values!

# Arguments:
- `active_param_keys::Vector{String}`: list of active parameter keys.
- `active_params_box::Vector{Tuple{Float64,Float64}}`: list of tuples specifying the bounds for each active parameter.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.
- `PT_all::Vector{ExoplanetsSysSim.ParamsTriangle}=ExoplanetsSysSim.ParamsTriangle[]`: a vector of objects containing the pairs of parameters that are triangle-transformed (default is an empty list indicating no transformed parameters).
NOTE 1: the parameters in `active_param_keys` should be sorted (in alphabetical order), such that the function `make_vector_of_active_param_keys(sim_param)` returns the same list, and the tuples in `active_params_box` should correspond to these parameters.
NOTE 2: the bounds in `active_params_box` for transformed parameters should be simply (0,1).

# Returns:
- `active_param_draw::Vector{Float64}`: a list of randomly drawn active parameter values.
- `active_param_draw_r::Vector{Float64}`: a list of randomly drawn active parameter values, with transformed parameters having their transformed values (i.e. `r1`, `r2` pairs with values in (0,1)).
Note: also writes the drawn active parameter values to `sim_param`!
"""
function draw_random_active_params(active_param_keys::Vector{String}, active_params_box::Vector{Tuple{Float64,Float64}}, sim_param::SimParam; PT_all::Vector{ExoplanetsSysSim.ParamsTriangle}=ExoplanetsSysSim.ParamsTriangle[])
    @assert length(active_param_keys) == length(active_params_box)

    for (i,param_key) in enumerate(active_param_keys)
        active_param_draw = active_params_box[i][1] .+ (active_params_box[i][2] - active_params_box[i][1])*rand(1)
        add_param_active(sim_param, param_key, active_param_draw[1])
    end

    # To draw any transformed parameters:
    if length(PT_all) > 0
        params_r = zeros(length(active_param_keys))
        for (i,pt) in enumerate(PT_all)
            r1r2 = (rand(), rand())
            params_r[collect(pt.id_xy)] .= r1r2
            transformed_params = map_square_to_triangle(r1r2, pt)
            add_param_active(sim_param, active_param_keys[pt.id_xy[1]], transformed_params[1])
            add_param_active(sim_param, active_param_keys[pt.id_xy[2]], transformed_params[2])
        end

        active_param_draw = make_vector_of_sim_param(sim_param)
        active_param_draw_r = deepcopy(active_param_draw)
        active_param_draw_r[params_r .!= 0] .= params_r[params_r .!= 0]
    else # there are no transformed params
        active_param_draw = make_vector_of_sim_param(sim_param)
        active_param_draw_r = Float64[]
    end

    return (active_param_draw, active_param_draw_r)
end





"""
    simulate_catalog_and_calc_all_distances_dict(active_param, sim_param; ss_fit, AD_mod=true)

Simulate a catalog with the parameters in `active_param` and compute the distances compared to the catalog provided with `ss_fit`.

# Arguments:
- `active_param::Vector{Float64}`: a vector of values for the active model parameters.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.
- `ss_fit::CatalogSummaryStatistics`: an object containing the summary statistics of the catalog we are trying to fit to.
- `AD_mod::Bool=true`: whether to use the original (if false) or modified (if true) AD distance.

# Returns:
- `dists::Dict{String,Float64}`: a dictionary containing the individual distance terms.
- `counts::Dict{String,Any}`: a dictionary containing the multiplicities, numbers of planets, and planet pairs.
- `summary_stat::CatalogSummaryStatistics`: an object containing the summary statistics of the simulated catalog.
"""
function simulate_catalog_and_calc_all_distances_dict(active_param::Vector{Float64}, sim_param::SimParam; ss_fit::CatalogSummaryStatistics, AD_mod::Bool=true)

    # Generate a simulated catalog with the input model parameters in 'active_param':
    sim_param_here = deepcopy(sim_param)
    ExoplanetsSysSim.update_sim_param_from_vector!(active_param,sim_param_here)
    cat_phys = generate_kepler_physical_catalog(sim_param_here)
    cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param_here)
    cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param_here)
    summary_stat = calc_summary_stats_model(cat_obs,sim_param_here)

    # Compute the distances between the simulated catalog and 'ss_fit':
    dists, counts = calc_all_distances_dict(sim_param, summary_stat, ss_fit; AD_mod=AD_mod)

    return (dists, counts, summary_stat)
end



"""
    target_function(active_param, sim_param; ss_fit, dists_include, weights, AD_mod=true, all_dist=false, save_dist=true)

Evaluate the target function to return the total weighted distance (or individual weighted distance terms) between a simulated catalog with the parameters in `active_param` and the catalog provided with `ss_fit`.

# Arguments:
- `active_param::Vector{Float64}`: a vector of values for the active model parameters.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.
- `ss_fit::CatalogSummaryStatistics`: an object containing the summary statistics of the catalog we are trying to fit to.
- `dists_include::Vector{String}`: a vector of strings for the names of the distance terms to be included in the total distance.
- `weights::Dict{String,Float64}`: a dictionary containing the weights of the individual distance terms.
- `AD_mod::Bool=true`: whether to use the original (if false) or modified (if true) AD distance.
- `all_dist::Bool=false`: whether to return all the individual weighted distance terms (if true) or just the total weighted distance (if true).
- `save_dist::Bool=true`: whether to save all the weighted distances to a file.
- `f::Union{IOStream,Nothing}=nothing`: file for printing distances to, if provided.

# Returns:
The total weighted distance (a Float64) if `all_dist=false` or the individual weighted distance terms (a Vector{Float64}) if `all_dist=true` of the simulated catalog compared to the input catalog.
"""
function target_function(active_param::Vector{Float64}, sim_param::SimParam; ss_fit::CatalogSummaryStatistics, dists_include::Vector{String}, weights::Dict{String,Float64}, AD_mod::Bool=true, all_dist::Bool=false, save_dist::Bool=true, f::Union{IOStream,Nothing}=nothing)

    println("Active parameter values: ", active_param)
    if save_dist
        println(f, "Active_params: ", active_param) # to write the params to file
    end
    #=
    log_rate = active_param[2] + active_param[3]
    if log_rate > 2.5
        println("Not simulating to save time because ln(lc*lp) = $(log_rate) > 2.5 !")
        if save_dist
            println(f, "Not simulating to save time because ln(lc*lp) = $(log_rate) > 2.5 !")
        end
        return 1e6
    end
    =#

    # Generate a simulated catalog and compute the distances:
    dists, counts, summary_stat = simulate_catalog_and_calc_all_distances_dict(active_param, sim_param; ss_fit=ss_fit, AD_mod=AD_mod)

    # Compute the individual and total weighted distances:
    dists_used = Dict([(key, 0.) for key in dists_include])
    dists_used_w = Dict([(key, 0.) for key in dists_include])
    for (i,key) in enumerate(dists_include)
        dists_used[key] = dists[key]
        dists_used_w[key] = dists[key]*weights[key]
    end
    dists_used_keys = keys(dists_used_w)
    dists_used_vals = values(dists_used)
    dists_used_vals_w = values(dists_used_w)

    # Print and/or write the distances to file:
    println("Counts: ", counts["Nmult1"], [counts["n_pl1"], counts["n_pairs1"]])
    #println("d_all_keys: ", keys(dists))
    #println("d_all_vals: ", values(dists))
    println("d_used_keys: ", dists_used_keys)
    println("d_used_vals: ", dists_used_vals, [sum(dists_used_vals)])
    println("d_used_vals_w: ", dists_used_vals_w, [sum(dists_used_vals_w)])
    println("#")
    if save_dist
        println(f, "Counts: ", counts["Nmult1"], [counts["n_pl1"], counts["n_pairs1"]])
        #println(f, "d_all_keys: ", keys(dists))
        #println(f, "d_all_vals: ", values(dists))
        println(f, "d_used_keys: ", dists_used_keys)
        println(f, "d_used_vals: ", dists_used_vals, [sum(dists_used_vals)])
        println(f, "d_used_vals_w: ", dists_used_vals_w, [sum(dists_used_vals_w)])
        println(f, "#")
    end

    if all_dist
        return dists_used_vals_w
    else
        return sum(dists_used_vals_w)
    end
end



"""
    target_function_split_stars(active_param, sim_param; cssc_fit, dists_include_all, weights_all, names_samples, dists_include_samples, weights_samples, AD_mod=true, all_dist=false, save_dist=true, f=nothing)

Evaluate the target function to return the total weighted distance (or individual weighted distance terms) between two simulated catalogs splitting the stellar catalog with the parameters in `active_param` and the catalogs provided with `ss_fit_all` (split in the same way for the stellar catalog).

# Arguments:
- `active_param::Vector{Float64}`: a vector of values for the active model parameters.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.
- `cssc_fit::CatalogSummaryStatisticsCollection`: object containing a number of CatalogSummaryStatistics objects for each sample to compare to.
- `dists_include_all::Vector{String}`: vector of strings for the names of the distance terms of the combined sample to be included in the total distance.
- `weights_all::Dict{String,Float64}`: dictionary containing the weights of the individual distance terms for the combined sample.
- `names_samples::Array{String,1}`: names of the split stellar catalog files.
- `dists_include_samples::Array{Vector{String},1}`: array of vectors of strings for the names of the distance terms in the split samples to be included in the total distance.
- `weights_samples::Array{Dict{String,Float64},1}`: array of dictionaries containing the weights of the individual distance terms for the split samples.
- `AD_mod::Bool=true`: whether to use the original (if false) or modified (if true) AD distance.
- `all_dist::Bool=false`: whether to return all the individual weighted distance terms (if true) or just the total weighted distance (if true).
- `save_dist::Bool=true`: whether to save all the weighted distances to a file.
- `f::Union{IOStream,Nothing}=nothing`: file for printing distances to, if provided.

# Returns:
The total weighted distance (a Float64) if `all_dist=false` or the individual weighted distance terms (a Vector{{Vector{Float64}}}) if `all_dist=true` of the simulated catalog compared to the input catalog.
"""
function target_function_split_stars(active_param::Vector{Float64}, sim_param::SimParam; cssc_fit::CatalogSummaryStatisticsCollection, dists_include_all::Vector{String}, weights_all::Dict{String,Float64}, names_samples::Array{String,1}, dists_include_samples::Array{Vector{String},1}, weights_samples::Array{Dict{String,Float64},1}, AD_mod::Bool=true, all_dist::Bool=false, save_dist::Bool=true, f::Union{IOStream,Nothing}=nothing)
    @assert length(names_samples) == length(dists_include_samples) == length(weights_samples)

    println("Active parameter values: ", active_param)
    if save_dist
        println(f, "Active_params: ", active_param) # to write the params to file
    end

    # Generate a simulated catalog:
    ExoplanetsSysSim.update_sim_param_from_vector!(active_param,sim_param)
    cat_phys = generate_kepler_physical_catalog(sim_param)
    cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param)
    cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param)

    # Compute the summary statistics:
    cssc = calc_summary_stats_collection_model(cat_obs, names_samples, [cssc_fit.star_id_samples[name] for name in names_samples], sim_param)

    # Compute and write the distances:
    names_list = ["all"; names_samples]
    dists_include_list = [[dists_include_all]; dists_include_samples]
    weights_list = [weights_all; weights_samples]

    dists_used_vals_w_list = []
    for (n,name) in enumerate(names_list)

        # Compute the individual and total weighted distances:
        dists, counts = calc_all_distances_dict(sim_param, cssc.css_samples[name], cssc_fit.css_samples[name]; AD_mod=AD_mod)
        dists_used, dists_used_w = Dict{String,Float64}(), Dict{String,Float64}()
        for (i,key) in enumerate(dists_include_list[n])
            dists_used[key] = dists[key]
            dists_used_w[key] = dists[key]*weights_list[n][key]
        end
        dists_used_keys = keys(dists_used_w)
        dists_used_vals, dists_used_vals_w = values(dists_used), values(dists_used_w)
        push!(dists_used_vals_w_list, dists_used_vals_w)

        # Print and/or write the distances to file:
        println("[$name] Counts: ", counts["Nmult1"], [counts["n_pl1"], counts["n_pairs1"]])
        println("[$name] d_used_keys: ", dists_used_keys)
        println("[$name] d_used_vals: ", dists_used_vals, [sum(dists_used_vals)])
        println("[$name] d_used_vals_w: ", dists_used_vals_w, [sum(dists_used_vals_w)])
        if save_dist
            println(f, "[$name] Counts: ", counts["Nmult1"], [counts["n_pl1"], counts["n_pairs1"]])
            println(f, "[$name] d_used_keys: ", dists_used_keys)
            println(f, "[$name] d_used_vals: ", dists_used_vals, [sum(dists_used_vals)])
            println(f, "[$name] d_used_vals_w: ", dists_used_vals_w, [sum(dists_used_vals_w)])
        end
    end

    d_tot_w = sum([sum(x) for x in dists_used_vals_w_list])

    # Print and/or write the total distances to file:
    println("Total_dist_w: ", [d_tot_w])
    println("#")
    if save_dist
        println(f, "Total_dist_w: ", [d_tot_w])
        println(f, "#")
    end

    if all_dist
        return dists_used_vals_w_list
    else
        return d_tot_w
    end
end



"""
    target_function_transformed_params(active_param_transformed, PT_all, sim_param; ss_fit, dists_include, weights, AD_mod=true, all_dist=false, save_dist=true, f=nothing)

Same as `target_function`, but takes in a list of active parameters that have a pair of variables transformed via the triangle transformation. Evaluate `target_function` after transforming the transformed parameters back to physical model parameters.

# Arguments:
- `active_param_transformed::Vector{Float64}`: a vector of values for the transformed active model parameters.
- `PT_all::Vector{ExoplanetsSysSim.ParamsTriangle}`: a vector of objects containing the pairs of parameters that are triangle-transformed.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.
- `ss_fit::CatalogSummaryStatistics`: an object containing the summary statistics of the catalog we are trying to fit to.
- `dists_include::Vector{String}`: a vector of strings for the names of the distance terms to be included in the total distance.
- `weights::Dict{String,Float64}`: a dictionary containing the weights of the individual distance terms.
- `AD_mod::Bool=true`: whether to use the original (if false) or modified (if true) AD distance.
- `all_dist::Bool=false`: whether to return all the individual weighted distance terms (if true) or just the total weighted distance (if true).
- `save_dist::Bool=true`: whether to save all the weighted distances to a file.
- `f::Union{IOStream,Nothing}=nothing`: file for printing distances to, if provided.

# Returns:
The total weighted distance (a Float64) if `all_dist=false` or the individual weighted distance terms (a Vector{Float64}) if `all_dist=true` of the simulated catalog compared to the input catalog.
"""
function target_function_transformed_params(active_param_transformed::Vector{Float64}, PT_all::Vector{ExoplanetsSysSim.ParamsTriangle}, sim_param::SimParam; ss_fit::CatalogSummaryStatistics, dists_include::Vector{String}, weights::Dict{String,Float64}, AD_mod::Bool=true, all_dist::Bool=false, save_dist::Bool=true, f::Union{IOStream,Nothing}=nothing)

    active_param = deepcopy(active_param_transformed)
    for pt in PT_all
        id_xy = collect(pt.id_xy)
        r1r2 = Tuple(active_param_transformed[id_xy])
        active_param[id_xy] .= map_square_to_triangle(r1r2, pt)
    end

    return target_function(active_param, sim_param; ss_fit=ss_fit, dists_include=dists_include, weights=weights, AD_mod=AD_mod, all_dist=all_dist, save_dist=save_dist, f=f)
end



"""
    target_function_transformed_params_split_stars(active_param_transformed, PT_all, sim_param; cssc_fit, dists_include_all, weights_all, names_samples, dists_include_samples, weights_samples, AD_mod=true, all_dist=false, save_dist=true, f=nothing)

Same as `target_function_split_stars`, but takes in a list of active parameters that have a pair of variables transformed via the triangle transformation. Evaluate `target_function_split_stars` after transforming the transformed parameters back to physical model parameters.

# Arguments:
- `active_param_transformed::Vector{Float64}`: a vector of values for the transformed active model parameters.
- `PT_all::Vector{ExoplanetsSysSim.ParamsTriangle}`: a vector of objects containing the pairs of parameters that are triangle-transformed.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.
- `cssc_fit::CatalogSummaryStatisticsCollection`: object containing a number of CatalogSummaryStatistics objects for each sample to compare to.
- `dists_include_all::Vector{String}`: vector of strings for the names of the distance terms of the combined sample to be included in the total distance.
- `weights_all::Dict{String,Float64}`: dictionary containing the weights of the individual distance terms for the combined sample.
- `names_samples::Array{String,1}`: names of the split stellar catalog files.
- `dists_include_samples::Array{Vector{String},1}`: array of vectors of strings for the names of the distance terms in the split samples to be included in the total distance.
- `weights_samples::Array{Dict{String,Float64},1}`: array of dictionaries containing the weights of the individual distance terms for the split samples.
- `AD_mod::Bool=true`: whether to use the original (if false) or modified (if true) AD distance.
- `all_dist::Bool=false`: whether to return all the individual weighted distance terms (if true) or just the total weighted distance (if true).
- `save_dist::Bool=true`: whether to save all the weighted distances to a file.
- `f::Union{IOStream,Nothing}=nothing`: file for printing distances to, if provided.

# Returns:
The total weighted distance (a Float64) if `all_dist=false` or the individual weighted distance terms (a Vector{Float64}) if `all_dist=true` of the simulated catalog compared to the input catalog.
"""
function target_function_transformed_params_split_stars(active_param_transformed::Vector{Float64}, PT_all::Vector{ExoplanetsSysSim.ParamsTriangle}, sim_param::SimParam; cssc_fit::CatalogSummaryStatisticsCollection, dists_include_all::Vector{String}, weights_all::Dict{String,Float64}, names_samples::Array{String,1}, dists_include_samples::Array{Vector{String},1}, weights_samples::Array{Dict{String,Float64},1}, AD_mod::Bool=true, all_dist::Bool=false, save_dist::Bool=true, f::Union{IOStream,Nothing}=nothing)

    active_param = deepcopy(active_param_transformed)
    for pt in PT_all
        id_xy = collect(pt.id_xy)
        r1r2 = Tuple(active_param_transformed[id_xy])
        active_param[id_xy] .= map_square_to_triangle(r1r2, pt)
    end

    return target_function_split_stars(active_param, sim_param; cssc_fit=cssc_fit, dists_include_all=dists_include_all, weights_all=weights_all, names_samples=names_samples, dists_include_samples=dists_include_samples, weights_samples=weights_samples, AD_mod=AD_mod, all_dist=all_dist, save_dist=save_dist, f=f)
end





"""
    compute_weights_target_fitness_std_from_array(dists_vals_all, dists_keys_all; dists_include, save_dist=true, f=nothing)

Compute the weights (1/rms) for each distance term, as well as the mean and standard deviation of the total weighted distance, given an array of distance terms from multiple simulations.

# Arguments:
- `dists_vals_all::Array{Float64,2}`: an array of distance terms, where each row is one simulation.
- `dists_keys_all::Vector{String}`: a vector of strings for the names of all the distance terms (columns in `dists_vals_all`).
- `dists_include::Vector{String}`: a vector of strings for the names of the distance terms to be included in the total distance.
- `save_dist::Bool=true`: whether to save all the weighted distances to a file.
- `f::Union{IOStream,Nothing}=nothing`: file for printing distances to, if provided.
NOTE: `dists_include` specifies which distance terms are included in calculating the mean and standard deviation of the total weighted distance, but ALL the distance terms and their weights are returned regardless.

# Returns:
- `weights::Dict{String,Float64}`: a dictionary containing weights for the individual distance terms.
- `totdist_w_mean_used::Float64`: the total weighted distance including only distance terms in `dists_include`.
- `totdist_w_std_used::Float64`: the standard deviation of the total weighted distance including only distance terms in `dists_include`.
"""
function compute_weights_target_fitness_std_from_array(dists_vals_all::Array{Float64,2}, dists_keys_all::Vector{String}; dists_include::Vector{String}, save_dist::Bool=true, f::Union{IOStream,Nothing}=nothing)
    num_evals, n_dists = size(dists_vals_all)
    @assert length(dists_include) <= length(dists_keys_all) == n_dists

    dists_mean_all = transpose(mean(dists_vals_all, dims=1))[:,] # array of mean distances for each individual distance
    dists_rms_all = transpose(sqrt.(mean(dists_vals_all .^2, dims=1)))[:,] # array of rms (std around 0)  distances for each individual distance

    weights_all = 1 ./ dists_rms_all
    weights = Dict([(dists_keys_all[i], weights_all[i]) for i in 1:n_dists])

    dists_w_all = zeros(num_evals,length(weights_all))
    for i in 1:num_evals
        dists_w_all[i,:] = dists_vals_all[i,:]  .* weights_all
    end
    id_dists_used = [findfirst(dists_keys_all .== key) for key in dists_include]
    dists_w_used = view(dists_w_all, :, id_dists_used)

    dists_w_mean_all = transpose(mean(dists_w_all, dims=1))[:,] # array of mean weighted distances for each individual distance
    totdist_w_mean_all = mean(sum(dists_w_all, dims=2)) # mean weighted total distance
    totdist_w_mean_used = mean(sum(dists_w_used, dims=2))
    totdist_w_std_all = std(sum(dists_w_all, dims=2)) # std of weighted total distance
    totdist_w_std_used = std(sum(dists_w_used, dims=2))

    println("#")
    println("Dists_keys_all: ", dists_keys_all)
    println("Dists_mean_all: ", dists_mean_all)
    println("Dists_rms_all: ", dists_rms_all)
    println("Weights_all (1/dists_rms_all): ", weights_all)
    println("Dists_w_mean_all: ", dists_w_mean_all, [totdist_w_mean_all])
    println("# Total weighted distance (all terms): ", totdist_w_mean_all, " +/- ", totdist_w_std_all)
    println("# Total weighted distance (used terms): ", totdist_w_mean_used, " +/- ", totdist_w_std_used)
    if save_dist
        println(f, "#")
        println(f, "Dists_keys_all: ", dists_keys_all)
        println(f, "Dists_mean_all: ", dists_mean_all)
        println(f, "Dists_rms_all: ", dists_rms_all)
        println(f, "Weights_all (1/dists_rms_all): ", weights_all)
        println(f, "Dists_w_mean_all: ", dists_w_mean_all, [totdist_w_mean_all])
        println(f, "# Total weighted distance (all terms): ", totdist_w_mean_all, " +/- ", totdist_w_std_all)
        println(f, "# Total weighted distance (used terms): ", totdist_w_mean_used, " +/- ", totdist_w_std_used)
    end

    return (weights, totdist_w_mean_used, totdist_w_std_used)
end



"""
    compute_weights_target_fitness_std_perfect_model(num_evals, sim_param; ss_ref, AD_mod=true, save_dist=true, f=nothing)

Compute the weights (1/rms) for each distance term, as well as the mean and standard deviation of the total weighted distance, by simulating the same model many times.

# Arguments:
- `num_evals::Int64`: number of times to simulate the same model to compare to the reference model.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.
- `ss_ref::CatalogSummaryStatistics`: an object containing the summary statistics of the reference catalog.
- `AD_mod::Bool=true`: whether to use the original (if false) or modified (if true) AD distance.
- `save_dist::Bool=true`: whether to save all the weighted distances to a file.
- `f::Union{IOStream,Nothing}=nothing`: file for printing distances to, if provided.

# Returns:
- `active_param_true::Vector{Float64}`: a vector of the active model parameters assumed for the 'perfect' model.
- `weights::Dict{String,Float64}`: a dictionary containing weights for the individual distance terms.
- `totdist_w_mean::Float64`: the mean of the total weighted distance.
- `totdist_w_std::Float64`: the standard deviation of the total weighted distance.
NOTE: CURRENTLY DOES NOT SAVE OUTPUTS FOR WORKERS IF USING MULTIPLE PROCESSORS
"""
function compute_weights_target_fitness_std_perfect_model(num_evals::Int64, sim_param::SimParam; ss_ref::CatalogSummaryStatistics, AD_mod::Bool=true, save_dist::Bool=true, f::Union{IOStream,Nothing}=nothing)
    t_elapsed = @elapsed begin
        active_param_true = make_vector_of_sim_param(sim_param)
        println("# True values: ", active_param_true)
        if save_dist
            println(f, "# AD_mod: ", AD_mod)
            println(f, "# Format: Counts: [observed multiplicities][total planets, total planet pairs]")
            println(f, "# Format: d_all_keys: [names of distance terms]")
            println(f, "# Format: d_all_vals: [distance terms]")
        end

        # Simulate the same catalog once to get the total number of distance terms and their keys:
        dists, counts, summary_stat = simulate_catalog_and_calc_all_distances_dict(active_param_true, sim_param; ss_fit=ss_ref, AD_mod=AD_mod)
        dists_keys_all = collect(keys(dists))
        n_dists = length(dists_keys_all) # total number of individual distance terms

        # Simulate the same catalog for 'num_evals' times, calculating and saving the distances:
        if nprocs() == 1
            dists_vals_all = zeros(num_evals,n_dists)
            for i in 1:num_evals
                dists, counts, summary_stat = simulate_catalog_and_calc_all_distances_dict(active_param_true, sim_param; ss_fit=ss_ref, AD_mod=AD_mod)
                dists_vals_all[i,:] .= values(dists)

                # Print and/or write the distances to file:
                println("Counts: ", counts["Nmult1"], [counts["n_pl1"], counts["n_pairs1"]])
                println("d_all_keys: ", keys(dists))
                println("d_all_vals: ", values(dists))
                if save_dist
                    println(f, "Counts: ", counts["Nmult1"], [counts["n_pl1"], counts["n_pairs1"]])
                    println(f, "d_all_keys: ", keys(dists))
                    println(f, "d_all_vals: ", values(dists))
                end
            end
        else # if nprocs() > 1
            @passobj 1 workers() ss_ref # send 'ss_ref' to all workers

            dists_vals_all = SharedArray{Float64,2}(num_evals,n_dists)
            @sync @distributed for i in 1:num_evals
                dists, counts, summary_stat = simulate_catalog_and_calc_all_distances_dict(active_param_true, sim_param; ss_fit=ss_ref, AD_mod=AD_mod)
                dists_vals_all[i,:] .= values(dists)

                # Print and/or write the distances to file:
                println("Counts: ", counts["Nmult1"], [counts["n_pl1"], counts["n_pairs1"]])
                println("d_all_keys: ", keys(dists))
                println("d_all_vals: ", values(dists))
                if save_dist
                    println(f, "Counts: ", counts["Nmult1"], [counts["n_pl1"], counts["n_pairs1"]])
                    println(f, "d_all_keys: ", keys(dists))
                    println(f, "d_all_vals: ", values(dists))
                end
            end
            dists_vals_all = Array(dists_vals_all) # convert SharedArray to Array so we can pass it to the function below
        end

    end

    weights, totdist_w_mean, totdist_w_std = compute_weights_target_fitness_std_from_array(dists_vals_all, dists_keys_all; dists_include=dists_keys_all, save_dist=save_dist, f=f)

    println("#")
    println("# elapsed time: ", t_elapsed, " seconds")
    println("#")
    if save_dist
        println(f, "#")
        println(f, "# elapsed time: ", t_elapsed, " seconds")
        println(f, "#")
    end

    return (active_param_true, weights, totdist_w_mean, totdist_w_std)
end



"""
    compute_weights_target_fitness_std_perfect_model_split_stars(num_evals, sim_param; cssc_ref, names_samples, ss_refs, AD_mod=true, save_dist=true, f=nothing)

Compute the weights (1/rms) for each distance term, as well as the mean and standard deviation of the total weighted distance, by simulating the same model many times.

# Arguments:
- `num_evals::Int64`: number of times to simulate the same model to compare to the reference model.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.
- `cssc_ref::CatalogSummaryStatisticsCollection`: object containig the CatalogSummaryStatistics objects of the reference catalog for each sample.
- `names_samples::Array{String,1}`: names of the split stellar catalog files.
- `AD_mod::Bool=true`: whether to use the original (if false) or modified (if true) AD distance.
- `save_dist::Bool=true`: whether to save all the weighted distances to a file.
- `f::Union{IOStream,Nothing}=nothing`: file for printing distances to, if provided.

# Returns:
- `active_param_true::Vector{Float64}`: a vector of the active model parameters assumed for the 'perfect' model.
- `weights_dicts::Dict{String,Dict{String,Float64}}`: a dictionary of dictionaries containing weights for the individual distance terms for each sample.
- `totdist_w_mean::Float64`: the mean of the total weighted distance.
- `totdist_w_std::Float64`: the standard deviation of the total weighted distance.
"""
function compute_weights_target_fitness_std_perfect_model_split_stars(num_evals::Int64, sim_param::SimParam; cssc_ref::CatalogSummaryStatisticsCollection, names_samples::Array{String,1}, AD_mod::Bool=true, save_dist::Bool=true, f::Union{IOStream,Nothing}=nothing)
    t_elapsed = @elapsed begin
        active_param_true = make_vector_of_sim_param(sim_param)
        println("# True values: ", active_param_true)
        if save_dist
            println(f, "# AD_mod: ", AD_mod)
            println(f, "# Format: Counts: [observed multiplicities][total planets, total planet pairs]")
            println(f, "# Format: d_all_keys: [names of distance terms]")
            println(f, "# Format: d_all_vals: [distance terms]")
        end

        # Simulate a catalog once to get the total number of distance terms and their keys:
        dists, counts, summary_stat = simulate_catalog_and_calc_all_distances_dict(active_param_true, sim_param; ss_fit=cssc_ref.css_samples["all"], AD_mod=AD_mod)
        dists_keys_all = collect(keys(dists))
        n_dists = length(dists_keys_all) # total number of individual distance terms

        # Simulate the same catalog for 'num_evals' times, calculating and saving the distances:
        names_list = ["all"; names_samples]
        dists_all_vals_evals = Dict([(name, zeros(num_evals,n_dists)) for name in names_list])
        for i in 1:num_evals

            # Generate a simulated catalog:
            cat_phys = generate_kepler_physical_catalog(sim_param)
            cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param)
            cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param)

            # Compute the summary statistics:
            cssc = calc_summary_stats_collection_model(cat_obs, names_samples, [cssc_ref.star_id_samples[name] for name in names_samples], sim_param)

            # Compute and write the distances:
            for (n,name) in enumerate(names_list)

                # Compute the individual and total weighted distances:
                dists, counts = calc_all_distances_dict(sim_param, cssc.css_samples[name], cssc_ref.css_samples[name]; AD_mod=AD_mod)
                dists_all_keys, dists_all_vals = keys(dists), values(dists)

                dists_all_vals_evals[name][i,:] .= dists_all_vals

                # Print and/or write the distances to file:
                println("[$name] Counts: ", counts["Nmult1"], [counts["n_pl1"], counts["n_pairs1"]])
                println("[$name] d_all_keys: ", dists_all_keys)
                println("[$name] d_all_vals: ", dists_all_vals)
                if save_dist
                    println(f, "[$name] Counts: ", counts["Nmult1"], [counts["n_pl1"], counts["n_pairs1"]])
                    println(f, "[$name] d_all_keys: ", dists_all_keys)
                    println(f, "[$name] d_all_vals: ", dists_all_vals)
                end
            end

            println("#")
            if save_dist
                println(f, "#")
            end
        end
    end

    weights_dicts = Dict([(name, Dict{String,Float64}()) for name in names_list])
    totdist_w_mean_list = Float64[]
    totdist_w_std_list = Float64[]
    for name in names_list
        weights, totdist_w_mean, totdist_w_std = compute_weights_target_fitness_std_from_array(dists_all_vals_evals[name], dists_keys_all; dists_include=dists_keys_all, save_dist=save_dist, f=f)
        weights_dicts[name] = weights
        push!(totdist_w_mean_list, totdist_w_mean)
        push!(totdist_w_std_list, totdist_w_std)
    end

    totdist_w_mean = sum(totdist_w_mean_list)
    totdist_w_std = sqrt(sum(totdist_w_std_list.^2))

    println("#")
    println("# Total weighted distance (all terms): ", totdist_w_mean, " +/- ", totdist_w_std)
    println("# elapsed time: ", t_elapsed, " seconds")
    println("#")
    if save_dist
        println(f, "#")
        println(f, "# Total weighted distance (all terms): ", totdist_w_mean, " +/- ", totdist_w_std)
        println(f, "# elapsed time: ", t_elapsed, " seconds")
        println(f, "#")
    end

    return (active_param_true, weights_dicts, totdist_w_mean, totdist_w_std)
end



"""
    compute_weights_target_fitness_std_from_file(file_name, num_evals, sim_param; dists_include, save_dist=true, f=nothing)

Compute the weights (1/rms) for each distance term, as well as the mean and standard deviation of the total weighted distance, by loading a file with the same model simulated many times.

# Arguments:
- `file_name::String`: name of file with the precomputed distances.
- `num_evals::Int64`: number of times to simulate the same model to compare to the reference model.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.
- `dists_include::Vector{String}`: a vector of strings for the names of the distance terms to be included in the total distance.
- `save_dist::Bool=true`: whether to save all the weighted distances to a file.
- `f::Union{IOStream,Nothing}=nothing`: file for printing distances to, if provided.
NOTE: `dists_include` specifies which distance terms are included in calculating the mean and standard deviation of the total weighted distance, but ALL the distance terms and their weights are returned regardless.

# Returns:
- `active_param_true::Vector{Float64}`: a vector of the active model parameters assumed for the 'perfect' model.
- `weights::Dict{String,Float64}`: a dictionary containing weights for the individual distance terms.
- `totdist_w_mean_used::Float64`: the total weighted distance including only distance terms in `dists_include`.
- `totdist_w_std_used::Float64`: the standard deviation of the total weighted distance including only distance terms in `dists_include`.
"""
function compute_weights_target_fitness_std_from_file(file_name::String, num_evals::Int64, sim_param::SimParam; dists_include::Vector{String}, save_dist::Bool=true, f::Union{IOStream,Nothing}=nothing)
    t_elapsed = @elapsed begin
        active_param_true = make_vector_of_sim_param(sim_param)
        println("# True values: ", active_param_true)
        println("# Computing weights from pre-computed file.")
        if save_dist
            println(f, "# Computing weights from pre-computed file.")
        end

        # Read the number of distance terms and their names:
        dists_keys_all = String[]
        open(file_name) do f_weights
            for (i,line) in enumerate(eachline(f_weights))
                if length(line) > 10
                    if line[1:10] == "d_all_keys"
                        dists_keys_all_str = split(line[15:end-2], "\", \"")
                        append!(dists_keys_all, dists_keys_all_str)
                        break
                    end
                end
            end
        end

        # Read the distance terms for each simulation:
        dists_vals_all = zeros(num_evals,length(dists_keys_all))
        open(file_name) do f_weights
            eval_count = 1
            for (i,line) in enumerate(eachline(f_weights))
                if length(line) > 10
                    if line[1:10] == "d_all_vals"
                        dists_vals_all_str = split(line[15:end-1], ", ")
                        dists_vals_all[eval_count,:] = [parse(Float64, x) for x in dists_vals_all_str]
                        eval_count += 1
                    end
                end
            end
        end

    end

    weights, totdist_w_mean_used, totdist_w_std_used = compute_weights_target_fitness_std_from_array(dists_vals_all, dists_keys_all; dists_include=dists_include, save_dist=save_dist, f=f)

    println("#")
    println("# elapsed time: ", t_elapsed, " seconds")
    println("#")
    if save_dist
        println(f, "#")
        println(f, "# elapsed time: ", t_elapsed, " seconds")
        println(f, "#")
    end

    return (active_param_true, weights, totdist_w_mean_used, totdist_w_std_used)
end



"""
    compute_weights_target_fitness_std_from_file_split_samples(file_name, num_evals, sim_param; names_samples, dists_include_samples, dists_include_combined, save_dist=true, f=nothing)

Compute the weights (1/rms) for each distance term, as well as the mean and standard deviation of the total weighted distance, by loading a file with the same model simulated many times.

# Arguments:
- `file_name::String`: name of file with the precomputed distances.
- `num_evals::Int64`: number of times to simulate the same model to compare to the reference model.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.
- `names_samples::Array{String,1}`: names of the split stellar catalog files.
- `dists_include_samples::Array{Vector{String},1}`: array of vectors of strings for the names of the distance terms in the split samples to be included in the total distance.
- `dists_include_all::Vector{String}`: vector of strings for the names of the distance terms of the combined sample to be included in the total distance.
- `save_dist::Bool=true`: whether to save all the weighted distances to a file.
- `f::Union{IOStream,Nothing}=nothing`: file for printing distances to, if provided.
NOTE: `dists_include` specifies which distance terms are included in calculating the mean and standard deviation of the total weighted distance, but ALL the distance terms and their weights are returned regardless.

# Returns:
- `active_param_true::Vector{Float64}`: a vector of the active model parameters assumed for the 'perfect' model.
- `weights_dicts::Dict{String,Dict{String,Float64}}`: a dictionary of dictionaries containing weights for the individual distance terms for each sample.
- `totdist_w_mean_used::Float64`: the total weighted distance including only distance terms in `dists_include_samples` and `dists_include_combined`.
- `totdist_w_std_used::Float64`: the standard deviation of the total weighted distance including only distance terms in `dists_include_samples` and `dists_include_combined`.
"""
function compute_weights_target_fitness_std_from_file_split_samples(file_name::String, num_evals::Int64, sim_param::SimParam; names_samples::Array{String,1}, dists_include_samples::Array{Vector{String},1}, dists_include_all::Vector{String}, save_dist::Bool=true, f::Union{IOStream,Nothing}=nothing)
    t_elapsed = @elapsed begin
        active_param_true = make_vector_of_sim_param(sim_param)
        println("# True values: ", active_param_true)
        println("# Computing weights from pre-computed file.")
        if save_dist
            println(f, "# Computing weights from pre-computed file.")
        end

        names_list = ["all"; names_samples]
        dists_include_list = [[dists_include_all]; dists_include_samples]

        # Read the number of distance terms and their names:
        dists_all_keys = Dict{String,Vector{String}}()
        open(file_name) do f_weights
            for (i,line) in enumerate(eachline(f_weights))
                for name in names_list
                    chars = length(name)
                    if length(line) > 13+chars
                        if !haskey(dists_all_keys, name) && line[1:13+chars]=="[$name] d_all_keys"
                            dists_keys_all_str = split(line[18+chars:end-2], "\", \"")
                            dists_all_keys[name] = dists_keys_all_str
                        end
                    end
                end
            end
        end

        # Read the distance terms for each simulation:
        dists_all_vals_evals = Dict([(name, zeros(num_evals,length(dists_all_keys[name]))) for name in names_list])
        open(file_name) do f_weights
            eval_counts = Dict([(name, 1) for name in names_list])
            for (i,line) in enumerate(eachline(f_weights))
                for name in names_list
                    chars = length(name)
                    if length(line) > 13+chars
                        if line[1:13+chars]=="[$name] d_all_vals"
                            dists_vals_all_str = split(line[17+chars:end-1], ", ")
                            dists_all_vals_evals[name][eval_counts[name],:] = [parse(Float64, x) for x in dists_vals_all_str]
                            eval_counts[name] += 1
                        end
                    end
                end
            end
        end

    end

    weights_dicts = Dict([(name, Dict{String,Float64}()) for name in names_list])
    totdist_w_mean_used_list = Float64[]
    totdist_w_std_used_list = Float64[]
    for (i,name) in enumerate(names_list)
        weights, totdist_w_mean, totdist_w_std = compute_weights_target_fitness_std_from_array(dists_all_vals_evals[name], dists_all_keys[name]; dists_include=dists_include_list[i], save_dist=save_dist, f=f)
        weights_dicts[name] = weights
        push!(totdist_w_mean_used_list, totdist_w_mean)
        push!(totdist_w_std_used_list, totdist_w_std)
    end

    totdist_w_mean_used = sum(totdist_w_mean_used_list)
    totdist_w_std_used = sqrt(sum(totdist_w_std_used_list.^2))

    println("#")
    println("# Total weighted distance (all terms): ", totdist_w_mean_used, " +/- ", totdist_w_std_used)
    println("# elapsed time: ", t_elapsed, " seconds")
    println("#")
    if save_dist
        println(f, "#")
        println(f, "# Total weighted distance (all terms): ", totdist_w_mean_used, " +/- ", totdist_w_std_used)
        println(f, "# elapsed time: ", t_elapsed, " seconds")
        println(f, "#")
    end

    return (active_param_true, weights_dicts, totdist_w_mean_used, totdist_w_std_used)
end
