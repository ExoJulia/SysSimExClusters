dir_path = dirname(@__FILE__)

include(joinpath(dir_path, "../src/clusters.jl"))
include(joinpath(dir_path, "../src/planetary_catalog.jl"))
include(joinpath(dir_path, "../src/optimization.jl"))

##### To load model parameters found using the GP emulator and simulate catalogs if they pass a distance threshold:

save_path = "test"
file_name = "GP_train2000_meanf75.0_sigmaf2.7_lscales9.7_vol109.18_points50000_meanInf_stdInf_post-30.0.csv"
GP_points = CSV.read(joinpath(save_path, file_name), comment="#")
active_params_names = names(GP_points)[1:end-3]
active_params_best_all = GP_points[!,active_params_names]

# If transformed:
#
i_tf, j_tf = 3,4 # indices of transformed parameters
active_params_names[i_tf:j_tf] = ["log_rate_clusters", "log_rate_planets_per_cluster"]
active_params_best_all[!,i_tf], active_params_best_all[!,j_tf] = (active_params_best_all[!,i_tf] .- active_params_best_all[!,j_tf])/2., (active_params_best_all[!,i_tf] .+ active_params_best_all[!,j_tf])/2.
rename!(active_params_best_all, Symbol.(active_params_names))
#

model_name = "Clustered_P_R"
names_split = ["bluer", "redder"]
AD_mod = true
num_targs = 86760
dists_include_split = ["delta_f", "mult_CRPD_r", "periods_KS", "period_ratios_KS", "durations_KS", "duration_ratios_nonmmr_KS", "duration_ratios_mmr_KS", "depths_KS", "radius_ratios_KS"]
dists_include_all = dists_include_split

d_threshold, mean_f = 40., 75.
n_pass = 1000 # number of simulations we want to pass the distance threshold
n_save = 0 # number of simulations we want to pass the distance threshold and also save (choose a small number or else requires a lot of storage space); must not be greater than n_pass!

file_name = model_name*"_pass_GP_meanf$(mean_f)_thres$(d_threshold)_pass$(n_pass)_targs$(num_targs).txt"
f = open(joinpath(save_path, file_name), "w")





##### To split the Kepler data into redder and bluer halves:

bprp = stellar_catalog[!,:bp_rp] .- stellar_catalog[!,:e_bp_min_rp_interp]
med_bprp = median(bprp)
idx_bluer = collect(1:size(stellar_catalog,1))[bprp .< med_bprp]
idx_redder = collect(1:size(stellar_catalog,1))[bprp .>= med_bprp]
star_id_split = [idx_bluer, idx_redder]

cssck = calc_summary_stats_collection_Kepler(stellar_catalog, planet_catalog, names_split, star_id_split, sim_param)





##### To load a file with the weights:

active_param_true, weights, target_fitness, target_fitness_std = compute_weights_target_fitness_std_from_file_split_samples(joinpath(dir_path, "../src/Clustered_P_R_split_stars_weights_ADmod_$(AD_mod)_targs88912_evals100_all_pairs.txt"), 4950, sim_param; names_samples=names_split, dists_include_samples=[dists_include_split, dists_include_split], dists_include_all=dists_include_all, f=f)
weights_all = weights["all"]
weights_split = [weights["bluer"], weights["redder"]]





##### To simulate a catalog with each set of params from the GP emulator that passed the threshold, saving all the distances, and saving the catalog if the true distance also passes the threshold:

sim_param = setup_sim_param_model()
add_param_fixed(sim_param,"num_targets_sim_pass_one", num_targs)

println(f, "# Active parameters: ", String.(active_params_names))
println(f, "# AD_mod: ", AD_mod)
println(f, "# Distances used (split): ", dists_include_split)
println(f, "# Distances used (all): ", dists_include_all)
println(f, "#")
println(f, "# Format: Active_params: [active parameter values]")
println(f, "# Format: [sample] Counts: [observed multiplicities][total planets, total planet pairs]")
println(f, "# Format: [sample] d_used_keys: [names of distance terms]")
println(f, "# Format: [sample] d_used_vals: [distance terms][sum of distance terms]")
println(f, "# Format: [sample] d_used_vals_w: [weighted distance terms][sum of weighted distance terms]")
println(f, "#")

sim_count = 0
pass_count = 0
save_count = 0
summary_array = Array{Float64,2}(undef, 0, size(GP_points,2)+1)

t_elapsed = @elapsed begin
    while pass_count < n_pass && sim_count < size(active_params_best_all,1)
        global sim_count, pass_count, save_count, summary_array
        sim_count += 1
        println("# Generating simulated catalog ", sim_count)

        # To set up the model parameters:
        for (i,param_name) in enumerate(active_params_names)
            add_param_active(sim_param, string(param_name), active_params_best_all[sim_count, param_name])
        end
        println(f, "Active_params: ", collect(active_params_best_all[sim_count,:])) # to write the params to file

        # Generate a simulated catalog:
        cat_phys = generate_kepler_physical_catalog(sim_param)
        cat_phys_copy = deepcopy(cat_phys) # need to deepcopy to save later, since cat_phys_cut overwrites cat_phys too
        cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param)
        cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param)

        # Compute the summary statistics:
        cssc = calc_summary_stats_collection_model(cat_obs, names_split, [cssck.star_id_samples[name] for name in names_split], sim_param)
        summary_stat = cssc.css_samples["all"]

        # Compute and write the distances:
        names_list = ["all"; names_split]
        dists_include_list = [dists_include_all, dists_include_split, dists_include_split]
        weights_list = [weights_all; weights_split]

        dists_used_vals_w_list = []
        for (n,name) in enumerate(names_list)

            # Compute the individual and total weighted distances:
            dists, counts = calc_all_distances_dict(sim_param, cssc.css_samples[name], cssck.css_samples[name]; AD_mod=AD_mod)
            dists_used, dists_used_w = Dict{String,Float64}(), Dict{String,Float64}()
            for (i,key) in enumerate(dists_include_list[n])
                dists_used[key] = dists[key]
                dists_used_w[key] = dists[key]*weights_list[n][key]
            end
            dists_used_keys = keys(dists_used_w)
            dists_used_vals, dists_used_vals_w = values(dists_used), values(dists_used_w)
            push!(dists_used_vals_w_list, dists_used_vals_w)

            # Write the distances to file:
            println(f, "[$name] Counts: ", counts["Nmult1"], [counts["n_pl1"], counts["n_pairs1"]])
            println(f, "[$name] d_used_keys: ", dists_used_keys)
            println(f, "[$name] d_used_vals: ", dists_used_vals, [sum(dists_used_vals)])
            println(f, "[$name] d_used_vals_w: ", dists_used_vals_w, [sum(dists_used_vals_w)])
        end

        d_tot_w = sum([sum(x) for x in dists_used_vals_w_list])

        # Write the total distances to file:
        println(f, "Total_dist_w: ", [d_tot_w])
        println(f, "#")



        # To save the catalog if the weighted distance passes the distance threshold:
        if d_tot_w <= d_threshold
            pass_count += 1
            if save_count < n_save
                save_count += 1
                println("$sim_count: d_tot_w = $d_tot_w; catalog saved (pass_count = $pass_count, save_count = $save_count)")

                save_physical_catalog_given_cat_phys(cat_phys_copy, sim_param; save_path=save_path, run_number=save_count)
                save_observed_catalog_given_cat_phys_obs(cat_phys, cat_obs, summary_stat, sim_param; save_path=save_path, run_number=save_count)
            else
                println("$sim_count: d_tot_w = $d_tot_w; (pass_count = $pass_count, max save_count reached)")
            end
        else
            println("$sim_count: d_tot_w = $d_tot_w")
        end

        summary_array = vcat(summary_array, reshape([[GP_points[sim_count,j] for j in 1:size(GP_points,2)]; d_tot_w - mean_f], (1,size(GP_points,2)+1)))
    end
end

println(f, "# elapsed time: ", t_elapsed, " seconds")
close(f)

summary_table = DataFrame(summary_array, [names(GP_points); :dist_tot_weighted])
CSV.write(joinpath(save_path, "Simulate_GP_points_summary.txt"), summary_table)
