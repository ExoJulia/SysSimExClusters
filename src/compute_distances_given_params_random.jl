import DataFrames.skipmissing

include("clusters.jl")
include("planetary_catalog.jl")
include("optimization.jl")

sim_param = setup_sim_param_model()





##### To start saving the model iterations in the optimization into a file:

model_name = "Clustered_P_R_broken_R"
run_number = "_random"*ARGS[1]
use_KS_or_AD = "KS" #'KS' or 'AD' or 'Both' (need to be careful counting indices for 'dists_exclude'!!!)
AD_mod = true
Kep_or_Sim = "Kep" #'Kep' or 'Sim'
num_targs = 400030
evals = 2000
dists_exclude = [3,4,8,12,13,15,16,17] #Int64[] if want to include all distances

file_name = model_name*run_number*"_targs"*string(num_targs)*"_evals"*string(evals)*".txt"
f = open(file_name, "w")
println(f, "# All initial parameters:")
write_model_params(f, sim_param)





##### To run the same model multiple times to see how it compares to a simulated catalog with the same parameters:

using Random
Random.seed!(1234) #to have the same reference catalog and simulated catalogs for calculating the weights

#To generate a simulated catalog to fit to:
cat_phys = generate_kepler_physical_catalog(sim_param)
cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param)
cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param)
summary_stat_ref = calc_summary_stats_model(cat_obs,sim_param)

#To simulate more observed planets for the subsequent model generations:
add_param_fixed(sim_param,"num_targets_sim_pass_one", num_targs)
add_param_fixed(sim_param,"max_incl_sys", 0.0) #degrees; 0 (deg) for isotropic system inclinations; set closer to 90 (deg) for more transiting systems

active_param_true, weights, target_fitness, target_fitness_std = compute_weights_target_fitness_std_from_file("Weights1000_targs200015_maxincl60.txt", use_KS_or_AD ; weight=true, dists_exclude=dists_exclude, save_dist=true)





##### To draw the active parameters randomly (uniformly) within a search range:

active_param_keys = ["f_high_incl", "log_rate_clusters", "log_rate_planets_per_cluster", "power_law_P", "power_law_r1", "power_law_r2", "sigma_hk", "sigma_incl", "sigma_incl_near_mmr", "sigma_log_radius_in_cluster", "sigma_logperiod_per_pl_in_cluster"]
    #["break_radius", "f_high_incl", "log_rate_clusters", "log_rate_planets_per_cluster", "mr_power_index", "num_mutual_hill_radii", "power_law_P", "power_law_r1", "power_law_r2", "sigma_hk", "sigma_incl", "sigma_incl_near_mmr", "sigma_log_radius_in_cluster", "sigma_logperiod_per_pl_in_cluster"]
active_params_box = [(0., 1.), (log(0.5), log(5.)), (log(0.5), log(5.)), (-2., 1.), (-6., 0.), (-6., 0.), (0., 0.1), (10., 90.), (0., 10.), (0., 0.5), (0., 0.3)] #search ranges for all of the active parameters

println(f, "# Active parameters: ", make_vector_of_active_param_keys(sim_param))
println(f, "# Active parameters bounds: ", active_params_box)
println(f, "# Format: Active_params: [active parameter values]")
println(f, "# Format: Dist: [distances][total distance]")
println(f, "# Format: Dist_weighted: [weighted distances][total weighted distance]")
println(f, "# Distances used: ", use_KS_or_AD)
println(f, "# AD_mod: ", AD_mod)
println(f, "#")

Random.seed!()

t_elapsed = @elapsed begin
    for i in 1:evals
        for (i,param_key) in enumerate(active_param_keys)
            active_param_draw = active_params_box[i][1] .+ (active_params_box[i][2] - active_params_box[i][1])*rand(1)
            add_param_active(sim_param,param_key,active_param_draw[1])
        end
        active_params_draw = make_vector_of_sim_param(sim_param)

        target_function(active_params_draw, use_KS_or_AD, Kep_or_Sim ; AD_mod=AD_mod, weights=weights, all_dist=false, save_dist=true)
    end
end

println(f, "# elapsed time: ", t_elapsed, " seconds")
close(f)
