#Pkg.add("ParallelDataTransfer")
using Distributed

addprocs(2) # number of additional processors

@everywhere using SharedArrays
@everywhere using ParallelDataTransfer
@everywhere using ExoplanetsSysSim

@everywhere include("clusters.jl")
@everywhere include("planetary_catalog.jl")
@everywhere include("optimization.jl")

@everywhere sim_param = setup_sim_param_model()





##### To start saving the model iterations in the optimization into a file:

use_KS_or_AD = "KS" # doesn't actually matter for this script
AD_mod = true
num_targs = 1000
num_evals_weights = 20
dists_exclude = [2,4,8,12,13,15,16,17] # Int64[] if want to include all distances in weighted sum; all distances are saved regardless

file_name = "Clustered_P_R_broken_R_weights_ADmod_$(AD_mod)_targs$(num_targs)_evals$(num_evals_weights)"

sendto(workers(), num_targs=num_targs, file_name=file_name)

@everywhere f = open(file_name*"_worker"*string(myid())*".txt", "w")
println(f, "# All initial parameters:")
write_model_params(f, sim_param)





##### To run the same model multiple times to see how it compares to a simulated catalog with the same parameters:

using Random
Random.seed!(1234) # to have the same reference catalog and simulated catalogs for calculating the weights

#To generate a simulated catalog to fit to:
cat_phys = generate_kepler_physical_catalog(sim_param)
cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param)
cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param)
summary_stat_ref = calc_summary_stats_model(cat_obs,sim_param)

@passobj 1 workers() summary_stat_ref # to send the 'summary_stat_ref' object to all workers

# To simulate more observed planets for the subsequent model generations:
@everywhere add_param_fixed(sim_param,"num_targets_sim_pass_one", num_targs)
@everywhere add_param_fixed(sim_param,"max_incl_sys", 0.) # degrees; 0 (deg) for isotropic system inclinations; set closer to 90 (deg) for more transiting systems

active_param_true, weights, target_fitness, target_fitness_std = compute_weights_target_fitness_std_perfect_model(num_evals_weights, use_KS_or_AD ; AD_mod=AD_mod, weight=true, dists_exclude=dists_exclude, save_dist=true)

@everywhere close(f)
