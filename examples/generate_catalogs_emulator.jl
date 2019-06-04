dir_path = dirname(@__FILE__)

include(joinpath(dir_path, "../src/clusters.jl"))
include(joinpath(dir_path, "../src/planetary_catalog.jl"))
include(joinpath(dir_path, "../src/optimization.jl"))

##### To load model parameters found using the GP emulator and simulate catalogs if they pass a distance threshold:

save_path = "test"
file_name = "GP_train2000_meanf75.0_sigmaf1.0_lscales36.37_vol273.72_points10000_meanInf_stdInf_post-45.0.csv"
GP_points = CSV.read(joinpath(save_path, file_name), comment="#", allowmissing=:none)
active_params_names = names(GP_points)[1:end-3]
active_params_best_all = GP_points[active_params_names]

# If transformed:
#
active_params_names[2:3] = [Symbol("log_rate_clusters"), Symbol("log_rate_planets_per_cluster")]
active_params_best_all[:,2], active_params_best_all[:,3] = (active_params_best_all[:,2] .- active_params_best_all[:,3])/2., (active_params_best_all[:,2] .+ active_params_best_all[:,3])/2.
names!(active_params_best_all, active_params_names)
#

use_KS_or_AD = "KS"
AD_mod = true
dists_exclude = [2,3,8,12,13,15,16,17]

active_param_true, weights, target_fitness, target_fitness_std = compute_weights_target_fitness_std_from_file(joinpath(dir_path, "../src/Clustered_P_R_broken_R_weights_ADmod_$(AD_mod)_targs399675_evals1000.txt"), 1000, use_KS_or_AD ; weight=true, dists_exclude=dists_exclude, save_dist=false)





sim_param = setup_sim_param_model()
add_param_fixed(sim_param,"num_targets_sim_pass_one", 79935)

d_threshold, mean_f = 30., 75.
n_keep = 100

sim_count = 0
save_count = 0
summary_array = Array{Float64,2}(undef, 0, size(GP_points,2)+1)
while save_count < n_keep && sim_count < size(active_params_best_all,1)
    global sim_count, save_count, summary_array
    sim_count += 1
    println("Generating simulated catalog ", sim_count)

    # To set up the model parameters:
    for (i,param_name) in enumerate(active_params_names)
        add_param_active(sim_param, string(param_name), active_params_best_all[sim_count, param_name])
    end

    # To simulate a catalog:
    @time cat_phys = generate_kepler_physical_catalog(sim_param)
    cat_phys_copy = deepcopy(cat_phys)
    @time begin
        cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param)
        cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param)
        summary_stat = calc_summary_stats_model(cat_obs,sim_param)
    end

    # To compute the weighted distance:
    dists = calc_distance_Kepler(summary_stat, use_KS_or_AD; AD_mod=AD_mod, all_dist=true, save_dist=false)
    dists_w = dists .* weights
    dist_w = sum(dists_w)

    # To save the catalog if the weighted distance passes the distance threshold:
    if dist_w <= d_threshold
        println("### Saving catalog $sim_count: dist_w = $dist_w")
        save_count += 1
        save_physical_catalog_given_cat_phys(cat_phys_copy, sim_param; save_path=save_path, run_number=save_count)
        save_observed_catalog_given_cat_phys_obs(cat_phys, cat_obs, summary_stat, sim_param; save_path=save_path, run_number=save_count)
    else
        println("### Not saving catalog $sim_count: dist_w = $dist_w")
    end

    summary_array = vcat(summary_array, reshape([[GP_points[sim_count,j] for j in 1:size(GP_points,2)]; dist_w - mean_f], (1,size(GP_points,2)+1)))
end

summary_table = DataFrame(summary_array, [names(GP_points); :dist_tot_weighted])
CSV.write(joinpath(save_path, "Simulate_GP_points_summary.txt"), summary_table)
