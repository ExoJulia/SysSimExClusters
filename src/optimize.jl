include("clusters.jl")
include("planetary_catalog.jl")
include("optimization.jl")
using Random

sim_param = setup_sim_param_model()





##### To start saving the model iterations in the optimization into a file:

model_name = "Clustered_P_R_optimization"
optimization_number = "_random"*ARGS[1] # if want to run on the cluster with random initial active parameters: "_random"*ARGS[1]
AD_mod = true
num_targs = 79935
max_incl_sys = 0.
max_evals = 5000
dists_include = ["delta_f", "mult_CRPD_r", "periods_KS", "period_ratios_KS", "durations_KS", "duration_ratios_nonmmr_KS", "duration_ratios_mmr_KS", "depths_KS", "radius_ratios_KS"]
Pop_per_param = 4

file_name = model_name*optimization_number*"_targs$(num_targs)_evals$(max_evals).txt"
f = open(file_name, "w")
println(f, "# All initial parameters:")
write_model_params(f, sim_param)





##### To load a file with the weights, and simulate a reference catalog if we want to fit to a model:

# To generate a reference catalog, if we want to fit to a simulated catalog instead of data:
#=
Random.seed!(1234) # to have the same reference catalog
cat_phys = generate_kepler_physical_catalog(sim_param)
cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param)
cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param)
summary_stat_ref = calc_summary_stats_model(cat_obs,sim_param)
=#

# To simulate more observed planets for the subsequent model generations:
add_param_fixed(sim_param,"num_targets_sim_pass_one", num_targs)
add_param_fixed(sim_param,"max_incl_sys", max_incl_sys)

# To load and compute the weights, target distance, and target distance std from a precomputed file:
active_param_true, weights, target_fitness, target_fitness_std = compute_weights_target_fitness_std_from_file("Clustered_P_R_weights_ADmod_$(AD_mod)_targs399675_evals1000.txt", 1000, sim_param; dists_include=dists_include, f=f)





##### To draw the initial values of the active parameters randomly within a search range:

ParamsTriangles_all = [ExoplanetsSysSim.ParamsTriangle((9,10), (0.,0.), (30.,30.), (30.,0.))] #[]
active_param_keys = ["f_high_incl", "f_stars_with_planets_attempted", "log_rate_clusters", "log_rate_planets_per_cluster", "power_law_P", "power_law_r1", "power_law_r2", "sigma_hk", "sigma_incl", "sigma_incl_near_mmr", "sigma_log_radius_in_cluster", "sigma_logperiod_per_pl_in_cluster"]
active_params_box = [(0., 1.), (0., 1.), (log(0.5), log(5.)), (log(0.5), log(5.)), (-2., 2.), (-4., 2.), (-6., 0.), (0., 0.1), (0., 1.), (0., 1.), (0., 0.5), (0., 0.3)] #search ranges for all of the active parameters

Random.seed!() # to have a random set of initial parameters and optimization run

active_param_start, active_param_transformed_start = draw_random_active_params(active_param_keys, active_params_box, sim_param; PT_all=ParamsTriangles_all)



PopSize = length(active_param_true)*Pop_per_param

println("# Active parameters: ", make_vector_of_active_param_keys(sim_param))
println(f, "# Active parameters: ", make_vector_of_active_param_keys(sim_param))
println(f, "# Starting active parameter values: ", active_param_start)
println(f, "# Optimization active parameters search bounds: ", active_params_box)
if length(ParamsTriangles_all) > 0
    for pt in ParamsTriangles_all
        println(f, "# Transformed active parameters: ", make_vector_of_active_param_keys(sim_param)[collect(pt.id_xy)], "; Triangle bounds: $(pt.A), $(pt.B), $(pt.C)")
    end
end
println(f, "# Method: adaptive_de_rand_1_bin_radiuslimited")
println(f, "# PopulationSize: ", PopSize)
println(f, "# AD_mod: ", AD_mod)
println(f, "# Distances used: ", dists_include)
println(f, "#")
println(f, "# Format: Active_params: [active parameter values]")
println(f, "# Format: Counts: [observed multiplicities][total planets, total planet pairs]")
println(f, "# Format: d_used_keys: [names of distance terms]")
println(f, "# Format: d_used_vals: [distance terms][sum of distance terms]")
println(f, "# Format: d_used_vals_w: [weighted distance terms][sum of weighted distance terms]")
println(f, "#")

target_function(active_param_start, sim_param; ss_fit=ssk, dists_include=dists_include, weights=weights, AD_mod=AD_mod, f=f)
target_function_transformed_params(active_param_transformed_start, ParamsTriangles_all, sim_param; ss_fit=ssk, dists_include=dists_include, weights=weights, AD_mod=AD_mod, f=f)




##### To use automated optimization routines to optimize the active model parameters:

# Pkg.add("BlackBoxOptim")       # only need to do these once
# Pkg.checkout("BlackBoxOptim")  # needed to get the lastest version
using BlackBoxOptim              # see https://github.com/robertfeldt/BlackBoxOptim.jl for documentation

t_elapsed = @elapsed begin
    #opt_result = bboptimize(active_params -> target_function(active_params, sim_param; ss_fit=ssk, dists_include=dists_include, weights=weights, AD_mod=AD_mod, f=f); SearchRange = active_params_box, NumDimensions = length(active_param_true), Method = :adaptive_de_rand_1_bin_radiuslimited, PopulationSize = PopSize, MaxFuncEvals = max_evals, TargetFitness = target_fitness, FitnessTolerance = target_fitness_std, TraceMode = :verbose)

    opt_result = bboptimize(active_params -> target_function_transformed_params(active_params, ParamsTriangles_all, sim_param; ss_fit=ssk, dists_include=dists_include, weights=weights, AD_mod=AD_mod, f=f); SearchRange = active_params_box, NumDimensions = length(active_param_true), Method = :adaptive_de_rand_1_bin_radiuslimited, PopulationSize = PopSize, MaxFuncEvals = max_evals, TargetFitness = target_fitness, FitnessTolerance = target_fitness_std, TraceMode = :verbose)
end

println(f, "# best_candidate: ", best_candidate(opt_result))
println(f, "# best_fitness: ", best_fitness(opt_result))
println(f, "# elapsed time: ", t_elapsed, " seconds")
close(f)

