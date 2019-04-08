#Pkg.add("ParallelDataTransfer")
using Distributed

addprocs(3) #number of additional processors

@everywhere using ParallelDataTransfer
@everywhere using ExoplanetsSysSim

import DataFrames.skipmissing

@everywhere include("clusters.jl")
@everywhere include("planetary_catalog.jl")
@everywhere include("optimization.jl")

@everywhere sim_param = setup_sim_param_model()





##### To start saving the model iterations in the optimization into a file:

@everywhere add_param_fixed(sim_param,"num_targets_sim_pass_one",80006)

model_name = "Clustered_P_R_broken_R_optimization"
optimization_number = "_random"*ARGS[1] #if want to run on the cluster with random initial active parameters: "_random"*ARGS[1]
use_KS_or_AD = "KS" #'KS' or 'AD'
AD_mod = true
Kep_or_Sim = "Kep" #'Kep' or 'Sim'
max_evals = 10000
num_evals_weights = 20
dists_exclude = [7,8,13] #Int64[] if want to include all distances
Pop_per_param = 4

file_name = model_name*optimization_number*"_targs"*string(get_int(sim_param,"num_targets_sim_pass_one"))*"_evals"*string(max_evals)

sendto(workers(), max_evals=max_evals, file_name=file_name)

@everywhere f = open(file_name*"_worker"*string(myid())*".txt", "w")
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

@passobj 1 workers() summary_stat_ref #to send the 'summary_stat_ref' object to all workers

#To simulate more observed planets for the subsequent model generations:
@everywhere add_param_fixed(sim_param,"max_incl_sys",80.0) #degrees; 0 (deg) for isotropic system inclinations; set closer to 90 (deg) for more transiting systems

active_param_true, weights, target_fitness, target_fitness_std = compute_weights_target_fitness_std_perfect_model(num_evals_weights, use_KS_or_AD ; AD_mod=AD_mod, weight=true, dists_exclude=dists_exclude, save_dist=true)





##### To draw the initial values of the active parameters randomly within a search range:

transformed_indices = []
active_param_keys = ["log_rate_clusters", "log_rate_planets_per_cluster", "power_law_P", "power_law_r1", "power_law_r2", "sigma_hk", "sigma_incl", "sigma_incl_near_mmr", "sigma_log_radius_in_cluster", "sigma_logperiod_per_pl_in_cluster"]
    #["break_radius", "log_rate_clusters", "log_rate_planets_per_cluster", "mr_power_index", "num_mutual_hill_radii", "power_law_P", "power_law_r1", "power_law_r2", "sigma_hk", "sigma_incl", "sigma_incl_near_mmr", "sigma_log_radius_in_cluster", "sigma_logperiod_per_pl_in_cluster"]
active_params_box = [(log(1.), log(5.)), (log(1.), log(5.)), (-0.5, 1.5), (-6., 0.), (-6., 0.), (0., 0.1), (0., 5.), (0., 5.), (0., 0.5), (0., 0.3)] #search ranges for all of the active parameters
    #[(0.5*ExoplanetsSysSim.earth_radius, 10.*ExoplanetsSysSim.earth_radius), (log(1.), log(5.)), (log(1.), log(5.)), (1., 4.), (3., 20.), (-0.5, 1.5), (-6., 0.), (-6., 0.), (0., 0.1), (0., 5.), (0., 5.), (0.1, 1.0), (0., 0.3)] #search ranges for all of the active parameters
#active_params_transformed_box = [(log(1.), log(5.)), (log(1.), log(5.)), (-0.5, 1.5), (-6., 0.), (-6., 0.), (0., 0.1), (0., 1.), (0., 1.), (0.1, 1.0), (0., 0.3)] #search ranges for all of the active parameters
#transformed_triangle = [[0., 0.], [5., 5.], [5., 0.]] #vertices (x,y) of the triangle for the transformed params

#To randomly draw (uniformly) a value for each active model parameter within its search range:

Random.seed!() #to have a random set of initial parameters and optimization run

for (i,param_key) in enumerate(active_param_keys)
    if ~in(i,transformed_indices)
        active_param_draw = active_params_box[i][1] .+ (active_params_box[i][2] - active_params_box[i][1])*rand(1)
        add_param_active(sim_param,param_key,active_param_draw[1])
    end
end
#r1_r2 = rand(2)
#transformed_params = map_square_to_triangle(r1_r2[1], r1_r2[2], transformed_triangle[1], transformed_triangle[2], transformed_triangle[3])
#add_param_active(sim_param, active_param_keys[transformed_indices[1]], transformed_params[1])
#add_param_active(sim_param, active_param_keys[transformed_indices[2]], transformed_params[2])
active_param_start = make_vector_of_sim_param(sim_param)
#active_param_transformed_start = deepcopy(active_param_start)
#active_param_transformed_start[transformed_indices] = r1_r2

PopSize = length(active_param_true)*Pop_per_param

println("# Active parameters: ", make_vector_of_active_param_keys(sim_param))
println(f, "# Active parameters: ", make_vector_of_active_param_keys(sim_param))
println(f, "# Starting active parameter values: ", active_param_start)
println(f, "# Optimization active parameters search bounds: ", active_params_box)
println(f, "# Transformed active parameters: ", make_vector_of_active_param_keys(sim_param)[transformed_indices])
#println(f, "# Transformed active parameters search triangle vertices: ", transformed_triangle)
println(f, "# Method: dxnes")
println(f, "# PopulationSize: ", PopSize)
println(f, "# Format: Active_params: [active parameter values]")
println(f, "# Format: Dist: [distances][total distance]")
println(f, "# Format: Dist_weighted: [weighted distances][total weighted distance]")
println(f, "# Distances used: ", use_KS_or_AD)
println(f, "# AD_mod: ", AD_mod)

target_function(active_param_start, use_KS_or_AD, Kep_or_Sim ; AD_mod=AD_mod, weights=weights, all_dist=false, save_dist=true) #to simulate the model once with the drawn parameters before starting the optimization
#target_function_transformed_params(active_param_transformed_start, transformed_indices, transformed_triangle[1], transformed_triangle[2], transformed_triangle[3], use_KS_or_AD, Kep_or_Sim ; AD_mod=AD_mod, weights=weights, all_dist=false, save_dist=true) #to simulate the model once with the drawn parameters before starting the optimization





##### To use automated optimization routines to optimize the active model parameters:

# Pkg.add("BlackBoxOptim")       # only need to do these once
# Pkg.checkout("BlackBoxOptim")  # needed to get the lastest version
using BlackBoxOptim              # see https://github.com/robertfeldt/BlackBoxOptim.jl for documentation

t_elapsed = @elapsed begin
    opt_result = bboptimize(active_params -> target_function(active_params, use_KS_or_AD, Kep_or_Sim ; AD_mod=AD_mod, weights=weights, all_dist=false, save_dist=true); SearchRange = active_params_box, NumDimensions = length(active_param_true), Method = :dxnes, PopulationSize = PopSize, MaxFuncEvals = max_evals, TargetFitness = target_fitness, FitnessTolerance = target_fitness_std, TraceMode = :verbose, Workers = workers())
end

println(f, "# best_candidate: ", best_candidate(opt_result))
println(f, "# best_fitness: ", best_fitness(opt_result))
println(f, "# elapsed time: ", t_elapsed, " seconds")
@everywhere close(f)
