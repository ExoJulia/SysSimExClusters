#Pkg.add("ParallelDataTransfer")
using Distributed

#addprocs(10) # number of additional processors

@everywhere using SharedArrays
@everywhere using ParallelDataTransfer
@everywhere using ExoplanetsSysSim

@everywhere include("clusters.jl")
@everywhere include("planetary_catalog.jl")
@everywhere include("optimization.jl")

@everywhere sim_param = setup_sim_param_model()





##### To start saving the model iterations in the optimization into a file:

AD_mod = true
num_targs = 79935*5
max_incl_sys = 0.
num_evals_weights = 1000

file_name_base = "Clustered_P_R_weights_ADmod_$(AD_mod)_targs$(num_targs)_evals$(num_evals_weights)"

sendto(workers(), num_targs=num_targs, max_incl_sys=max_incl_sys, file_name_base=file_name_base)

@everywhere f = open("$(file_name_base)_worker$(myid()).txt", "w")
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

# To simulate more observed planets for the subsequent model generations:
@everywhere add_param_fixed(sim_param,"num_targets_sim_pass_one", num_targs)
@everywhere add_param_fixed(sim_param,"max_incl_sys", max_incl_sys)

active_param_true, weights, target_fitness, target_fitness_std = compute_weights_target_fitness_std_perfect_model(num_evals_weights, sim_param; ss_ref=summary_stat_ref, AD_mod=AD_mod, f=f)

@everywhere close(f)





# Test the loading of the weights file:

#dists_include = ["delta_f", "mult_CRPD_r", "periods_KS", "period_ratios_KS", "durations_KS", "duration_ratios_KS", "duration_ratios_nonmmr_KS", "duration_ratios_mmr_KS", "depths_KS", "radius_ratios_KS"]
#active_param_true2, weights2, target_fitness2, target_fitness_std2 = compute_weights_target_fitness_std_from_file("$(file_name_base)_worker$(myid()).txt", num_evals_weights, sim_param; dists_include=dists_include, save_dist=false)
