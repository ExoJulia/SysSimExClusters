using Distributed # just to compile a function



"""
    simulate_catalog_and_calc_all_distances_dict(active_param, sim_param; summary_stat_fit, AD_mod=true)

Simulate a catalog with the parameters in `active_param` and compute the distances compared to the catalog provided with `summary_stat_fit`.

# Arguments:
- `active_param::Vector{Float64}`: a vector of values for the active model parameters.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.
- `summary_stat_fit::ExoplanetsSysSim.CatalogSummaryStatistics`: an object containing the summary statistics of the catalog we are trying to fit to.
- `AD_mod::Bool=true`: whether to use the original (if false) or modified (if true) AD distance.

# Returns:
- `dists::Dict{String,Float64}`: a dictionary containing the individual distance terms.
- `counts::Dict{String,Any}`: a dictionary containing the multiplicities, numbers of planets, and planet pairs.
"""
function simulate_catalog_and_calc_all_distances_dict(active_param::Vector{Float64}, sim_param::SimParam; summary_stat_fit::ExoplanetsSysSim.CatalogSummaryStatistics, AD_mod::Bool=true)

    # Generate a simulated catalog with the input model parameters in 'active_param':
    sim_param_here = deepcopy(sim_param)
    ExoplanetsSysSim.update_sim_param_from_vector!(active_param,sim_param_here)
    cat_phys = generate_kepler_physical_catalog(sim_param_here)
    cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param_here)
    cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param_here)
    summary_stat = calc_summary_stats_model(cat_obs,sim_param_here)

    # Compute the distances between the simulated catalog and 'summary_stat_fit':
    dists, counts = calc_all_distances_dict(sim_param, summary_stat, summary_stat_fit; AD_mod=AD_mod)

    return (dists, counts)
end

"""
    target_function(active_param, sim_param; summary_stat_fit, dists_include, weights, AD_mod=true, all_dist=false, save_dist=true)

Evaluate the target function to return the total weighted distance (or individual weighted distance terms) between a simulated catalog with the parameters in `active_param` and the catalog provided with `summary_stat_fit`.

# Arguments:
- `active_param::Vector{Float64}`: a vector of values for the active model parameters.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.
- `summary_stat_fit::ExoplanetsSysSim.CatalogSummaryStatistics`: an object containing the summary statistics of the catalog we are trying to fit to.
- `dists_include::Vector{String}`: a vector of strings for the names of the distance terms to be included in the total distance.
- `weights::Dict{String,Float64}`: a dictionary containing the weights of the individual distance terms.
- `AD_mod::Bool=true`: whether to use the original (if false) or modified (if true) AD distance.
- `all_dist::Bool=false`: whether to return all the individual weighted distance terms (if true) or just the total weighted distance (if true).
- `save_dist::Bool=true`: whether to save all the weighted distances to a file.
- `f::Union{IOStream,Nothing}=nothing`: file for printing distances to, if provided.

# Returns:
The total weighted distance (a Float64) if `all_dist=false` or the individual weighted distance terms (a Vector{Float64}) if `all_dist=true` of the simulated catalog compared to the input catalog.
"""
function target_function(active_param::Vector{Float64}, sim_param::SimParam; summary_stat_fit::ExoplanetsSysSim.CatalogSummaryStatistics, dists_include::Vector{String}, weights::Dict{String,Float64}, AD_mod::Bool=true, all_dist::Bool=false, save_dist::Bool=true, f::Union{IOStream,Nothing}=nothing)

    println("Active parameter values: ", active_param)
    if save_dist && !isnothing(f)
        println(f, "Active_params: ", active_param) # to write the params to file
    end
    #=
    log_rate = active_param[2] + active_param[3]
    if log_rate > 2.5
        println("Not simulating to save time because ln(lc*lp) = $(log_rate) > 2.5 !")
        if save_dist && !isnothing(f)
            println(f, "Not simulating to save time because ln(lc*lp) = $(log_rate) > 2.5 !")
        end
        return 1e6
    end
    =#

    # Generate a simulated catalog and compute the distances:
    dists, counts = simulate_catalog_and_calc_all_distances_dict(active_param, sim_param; summary_stat_fit=summary_stat_fit, AD_mod=AD_mod)

    # Compute the individual and total weighted distances:
    dists_used = Dict([(key, 0.) for key in dists_include])
    dists_used_w = Dict([(key, 0.) for key in dists_include])
    for (i,key) in enumerate(dists_include)
        dists_used[key] = dists[key]
        dists_used_w[key] = dists[key]*weights[key]
    end
    dists_keys_list = keys(dists_used_w)
    dists_used_list = values(dists_used)
    dists_used_w_list = values(dists_used_w)

    # Print and/or write the distances to file:
    println("Counts: ", counts["Nmult1"], [counts["n_pl1"], counts["n_pairs1"]])
    #println("Dists_keys_all: ", keys(dists))
    #println("Dists_vals_all: ", values(dists))
    println("Dists_keys_used: ", dists_keys_list)
    println("Dists_vals_used: ", dists_used_list, [sum(dists_used_list)])
    println("Dists_vals_used_w: ", dists_used_w_list, [sum(dists_used_w_list)])
    if save_dist && !isnothing(f)
        println(f, "Counts: ", counts["Nmult1"], [counts["n_pl1"], counts["n_pairs1"]])
        #println(f, "Dists_keys_all: ", keys(dists))
        #println(f, "Dists_vals_all: ", values(dists))
        println(f, "Dists_keys_used: ", dists_keys_list)
        println(f, "Dists_vals_used: ", dists_used_list, [sum(dists_used_list)])
        println(f, "Dists_vals_used_w: ", dists_used_w_list, [sum(dists_used_w_list)])
    end

    if all_dist
        return dists_used_w_list
    else
        return sum(dists_used_w_list)
    end
end

function target_function_transformed_params(active_param_transformed::Vector{Float64}, transformed_indices::Vector{Int64}, A::Vector{Float64}, B::Vector{Float64}, C::Vector{Float64}, sim_param::SimParam; summary_stat_fit::ExoplanetsSysSim.CatalogSummaryStatistics, dists_include::Vector{String}, weights::Dict{String,Float64}, AD_mod::Bool=true, all_dist::Bool=false, save_dist::Bool=true, f::Union{IOStream,Nothing}=nothing)

    r1, r2 = active_param_transformed[transformed_indices]
    @assert 0. <= r1 <= 1.
    @assert 0. <= r2 <= 1.
    active_param = deepcopy(active_param_transformed)
    active_param[transformed_indices] = map_square_to_triangle(r1, r2, A, B, C)

    return target_function(active_param, sim_param; summary_stat_fit=summary_stat_fit, dists_include=dists_include, weights=weights, AD_mod=AD_mod, all_dist=all_dist, save_dist=save_dist, f=f)
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

    println("Dists_mean_all: ", dists_mean_all)
    println("Dists_rms_all: ", dists_rms_all)
    println("Weights_all (1/dists_rms_all): ", weights_all)
    println("Dists_w_mean_all: ", dists_w_mean_all, [totdist_w_mean_all])
    println("# Total weighted distance (all terms): ", totdist_w_mean_all, " +/- ", totdist_w_std_all)
    println("# Total weighted distance (used terms): ", totdist_w_mean_used, " +/- ", totdist_w_std_used)
    if save_dist && !isnothing(f)
        println(f, "#")
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
    compute_weights_target_fitness_std_perfect_model(num_evals, sim_param; dists_include, AD_mod=true, save_dist=true, f=nothing)

Compute the weights (1/rms) for each distance term, as well as the mean and standard deviation of the total weighted distance, by simulating the same model many times.

# Arguments:
- `num_evals::Int64`: number of times to simulate the same model to compare to the reference model.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.
- `dists_include::Vector{String}`: a vector of strings for the names of the distance terms to be included in the total distance.
- `AD_mod::Bool=true`: whether to use the original (if false) or modified (if true) AD distance.
- `save_dist::Bool=true`: whether to save all the weighted distances to a file.
- `f::Union{IOStream,Nothing}=nothing`: file for printing distances to, if provided.
NOTE 1: `dists_include` specifies which distance terms are included in calculating the mean and standard deviation of the total weighted distance, but ALL the distance terms and their weights are returned regardless.
NOTE 2: a `summary_stat_ref` object (summary statistics for the reference catalog) is assumed as a global variable, and must be accessible to all workers if running on multiple processors (e.g. using `@passobj 1 workers() summary_stat_ref`)!

# Returns:
- `active_param_true::Vector{Float64}`: a vector of the active model parameters assumed for the 'perfect' model.
- `weights::Dict{String,Float64}`: a dictionary containing weights for the individual distance terms.
- `totdist_w_mean_used::Float64`: the total weighted distance including only distance terms in `dists_include`.
- `totdist_w_std_used::Float64`: the standard deviation of the total weighted distance including only distance terms in `dists_include`.
NOTE: the multiprocessor version of this function is untested!
"""
function compute_weights_target_fitness_std_perfect_model(num_evals::Int64, sim_param::SimParam; dists_include::Vector{String}, AD_mod::Bool=true, save_dist::Bool=true, f::Union{IOStream,Nothing}=nothing)
    t_elapsed = @elapsed begin
        active_param_true = make_vector_of_sim_param(sim_param)
        println("# True values: ", active_param_true)
        if save_dist && !isnothing(f)
            println(f, "# AD_mod: ", AD_mod)
            println(f, "# Format: Counts: [observed multiplicities][total planets, total planet pairs]")
            println(f, "# Format: Dists_keys_all: [names of distance terms]")
            println(f, "# Format: Dists_vals_all: [distance terms]")
        end

        global summary_stat_ref

        # Simulate the same catalog once to get the total number of distance terms and their keys:
        dists, counts = simulate_catalog_and_calc_all_distances_dict(active_param_true, sim_param; summary_stat_fit=summary_stat_ref, AD_mod=AD_mod)
        dists_keys_all = collect(keys(dists))
        n_dists = length(dists_keys_all) # total number of individual distance terms

        if nprocs() == 1
            dists_vals_all = zeros(num_evals,n_dists)
            for i in 1:num_evals
                dists, counts = simulate_catalog_and_calc_all_distances_dict(active_param_true, sim_param; summary_stat_fit=summary_stat_ref, AD_mod=AD_mod)
                dists_vals_all[i,:] .= values(dists)

                # Print and/or write the distances to file:
                println("Counts: ", counts["Nmult1"], [counts["n_pl1"], counts["n_pairs1"]])
                println("Dists_keys_all: ", keys(dists))
                println("Dists_vals_all: ", values(dists))
                if save_dist && !isnothing(f)
                    println(f, "Counts: ", counts["Nmult1"], [counts["n_pl1"], counts["n_pairs1"]])
                    println(f, "Dists_keys_all: ", keys(dists))
                    println(f, "Dists_vals_all: ", values(dists))
                end
            end
        else # if nprocs() > 1
            dists_vals_all = SharedArray{Float64,2}(num_evals,n_dists)
            @sync @distributed for i in 1:num_evals
                dists, counts = simulate_catalog_and_calc_all_distances_dict(active_param_true, sim_param; summary_stat_fit=summary_stat_ref, AD_mod=AD_mod)
                dists_vals_all[i,:] .= values(dists)

                # Print and/or write the distances to file:
                println("Counts: ", counts["Nmult1"], [counts["n_pl1"], counts["n_pairs1"]])
                println("Dists_keys_all: ", keys(dists))
                println("Dists_vals_all: ", values(dists))
                if save_dist && !isnothing(f)
                    println(f, "Counts: ", counts["Nmult1"], [counts["n_pl1"], counts["n_pairs1"]])
                    println(f, "Dists_keys_all: ", keys(dists))
                    println(f, "Dists_vals_all: ", values(dists))
                end
            end
            dists_vals_all = Array(dists_vals_all) # convert SharedArray to Array so we can pass it to the function below
        end

    end

    weights, totdist_w_mean_used, totdist_w_std_used = compute_weights_target_fitness_std_from_array(dists_vals_all, dists_keys_all; dists_include=dists_include, save_dist=save_dist, f=f)

    if save_dist && !isnothing(f)
        println(f, "# elapsed time: ", t_elapsed, " seconds")
        println(f, "#")
    end

    return (active_param_true, weights, totdist_w_mean_used, totdist_w_std_used)
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
        if save_dist && !isnothing(f)
            println(f, "# Computing weights from pre-computed file.")
        end

        dists_keys_all = String[]
        open(file_name) do f_weights
            for (i,line) in enumerate(eachline(f_weights))
                if length(line) > 14
                    if line[1:14] == "Dists_keys_all"
                        dists_keys_all_str = split(line[19:end-2], "\", \"")
                        append!(dists_keys_all, dists_keys_all_str)
                        break
                    end
                end
            end
        end

        dists_vals_all = zeros(num_evals,length(dists_keys_all))
        open(file_name) do f_weights
            eval_count = 1
            for (i,line) in enumerate(eachline(f_weights))
                if length(line) > 14
                    if line[1:14] == "Dists_vals_all"
                        dists_vals_all_str = split(line[19:end-1], ", ")
                        dists_vals_all[eval_count,:] = [parse(Float64, x) for x in dists_vals_all_str]
                        eval_count += 1
                    end
                end
            end
        end

    end

    weights, totdist_w_mean_used, totdist_w_std_used = compute_weights_target_fitness_std_from_array(dists_vals_all, dists_keys_all; dists_include=dists_include, save_dist=save_dist, f=f)

    if save_dist && !isnothing(f)
        println(f, "# elapsed time: ", t_elapsed, " seconds")
        println(f, "#")
    end

    return (active_param_true, weights, totdist_w_mean_used, totdist_w_std_used)
end
