import DataFrames.skipmissing

include("clusters.jl")
include("planetary_catalog.jl")
include("optimization.jl")

sim_param = setup_sim_param_model()





##### To start saving the model iterations in the optimization into a file:

add_param_fixed(sim_param,"num_targets_sim_pass_one",78005)

model_name = "Clustered_P_R_broken_R_optimization_blocks"
optimization_number = "_random"*ARGS[1] #if want to run on the cluster with random initial active parameters: "_random"*ARGS[1]
use_KS_or_AD = "KS" #'KS' or 'AD'
AD_mod = true
Kep_or_Sim = "Sim" #'Kep' or 'Sim'
max_evals = 300
num_evals_weights = 20
dists_exclude = [7,8,13] #Int64[] if want to include all distances
cycles = 3
Pop_per_param = 4

file_name = model_name*optimization_number*"_targs"*string(get_int(sim_param,"num_targets_sim_pass_one"))*"_evals"*string(max_evals)*"_cycles"*string(cycles)*".txt"
f = open(file_name, "w")
println(f, "# All initial parameters:")
write_model_params(f, sim_param)





##### To run the same model multiple times to see how it compares to a simulated catalog with the same parameters:

srand(1234) #to have the same reference catalog and simulated catalogs for calculating the weights

#To generate a simulated catalog to fit to:
cat_phys = generate_kepler_physical_catalog(sim_param)
cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param)
cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param)
summary_stat_ref = calc_summary_stats_model(cat_obs,sim_param)

#To simulate more observed planets for the subsequent model generations:
add_param_fixed(sim_param,"max_incl_sys",80.0) #degrees; 0 (deg) for isotropic system inclinations; set closer to 90 (deg) for more transiting systems

active_param_true, weights, target_fitness, target_fitness_std = compute_weights_target_fitness_std_perfect_model(num_evals_weights, use_KS_or_AD ; AD_mod=AD_mod, weight=true, dists_exclude=dists_exclude, save_dist=true)





##### To set up a block-optimization routine:

#To specify the blocks of random parameters, their search ranges, and to draw initial values:
active_param_keys = [["log_rate_clusters", "log_rate_planets_per_cluster", "num_mutual_hill_radii", "power_law_P", "sigma_logperiod_per_pl_in_cluster"], ["break_radius", "mr_power_index", "power_law_r1", "power_law_r2", "sigma_log_radius_in_cluster"], ["sigma_hk", "sigma_incl", "sigma_incl_near_mmr"]]
    #[["log_rate_clusters", "log_rate_planets_per_cluster", "num_mutual_hill_radii", "power_law_P", "sigma_logperiod_per_pl_in_cluster"], ["break_radius", "mr_power_index", "power_law_r1", "power_law_r2", "sigma_log_radius_in_cluster"], ["sigma_hk", "sigma_incl", "sigma_incl_near_mmr"]]
active_params_box = [[(log(1.), log(5.)), (log(1.), log(5.)), (3., 20.), (-0.5, 1.5), (0., 0.3)], [(0.5*ExoplanetsSysSim.earth_radius, 10.*ExoplanetsSysSim.earth_radius),(1., 4.), (-6., 0.), (-6., 0.), (0.1, 1.0)], [(0., 0.1), (0., 5.), (0., 5.)]] #search ranges for all of the active parameters
    #[[(log(1.), log(5.)), (log(1.), log(5.)), (3., 20.), (-0.5, 1.5), (0., 0.3)], [(0.5*ExoplanetsSysSim.earth_radius, 10.*ExoplanetsSysSim.earth_radius),(1., 4.), (-6., 0.), (-6., 0.), (0.1, 1.0)], [(0., 0.1), (0., 5.), (0., 5.)]] #search ranges for all of the active parameters

active_param_draws = [Float64[], Float64[], Float64[]] #[Float64[], Float64[], Float64[]]
for (i,active_param_block) in enumerate(active_param_keys)
    for (j,param_key) in enumerate(active_param_block)
        active_param_draw = active_params_box[i][j][1] + (active_params_box[i][j][2] - active_params_box[i][j][1])*rand(1)
        append!(active_param_draws[i], active_param_draw)
        add_param_fixed(sim_param, param_key, active_param_draw[1])
    end
end

println("# All active parameters: ", active_param_keys)
println(f, "# All active parameters: ", active_param_keys)
println(f, "# Starting active parameter values: ", active_param_draws)
println(f, "# Optimization active parameters search bounds: ", active_params_box)
println(f, "# Method: adaptive_de_rand_1_bin_radiuslimited")
println(f, "# PopulationSize: ", Pop_per_param, " (per param)")
println(f, "# Format: Active_params: [active parameter values]")
println(f, "# Format: Dist: [distances][total distance]")
println(f, "# Format: Dist_weighted: [weighted distances][total weighted distance]")
println(f, "# Distances used: ", use_KS_or_AD)
println(f, "# AD_mod: ", AD_mod)

#To start the block-optimization:

# Pkg.add("BlackBoxOptim")       # only need to do these once
# Pkg.checkout("BlackBoxOptim")  # needed to get the lastest version
using BlackBoxOptim              # see https://github.com/robertfeldt/BlackBoxOptim.jl for documentation

tic()
for c in 1:cycles
    println("# Cycle: ", c)
    println(f, "#")
    println(f, "# Cycle: ", c)
    println(f, "#")

    tic()
    for (i,active_param_block) in enumerate(active_param_keys)
        println("# Active parameters: ", active_param_block)
        println(f, "# Active parameters: ", active_param_block)

        #To make the current block of parameters active:
        for (j,param_key) in enumerate(active_param_block)
            add_param_active(sim_param, param_key, active_param_draws[i][j])
        end

        tic()
        opt_result = bboptimize(active_params -> target_function(active_params, use_KS_or_AD, Kep_or_Sim ; AD_mod=AD_mod, weights=weights, all_dist=false, save_dist=true); SearchRange = active_params_box[i], NumDimensions = length(active_param_block), Method = :adaptive_de_rand_1_bin_radiuslimited, PopulationSize = length(active_param_block)*Pop_per_param, MaxFuncEvals = max_evals, TargetFitness = target_fitness, FitnessTolerance = target_fitness_std, TraceMode = :verbose)
        t_block = toc()
        active_params_best = best_candidate(opt_result)
        best_distance = best_fitness(opt_result)

        println(f, "# Best parameters block: ", active_params_best)
        println(f, "# Best fitness block: ", best_distance)
        println(f, "# Block elapsed time: ", t_block, " seconds")
        println(f, "#")

        #To make the block of parameters fixed at the best values found:
        for (j,param_key) in enumerate(active_param_block)
            add_param_fixed(sim_param, param_key, active_params_best[j])
        end
    end
    t_cycle = toc()

    println(f, "# Cycle elapsed time: ", t_cycle, " seconds")
end
t_optim = toc()

println(f, "#")
println(f, "# Total elapsed time: ", t_optim, " seconds")
close(f)
