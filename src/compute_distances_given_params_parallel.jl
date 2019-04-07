#Pkg.add("ParallelDataTransfer")
using Distributed

addprocs(5) #number of additional processors

@everywhere using ParallelDataTransfer
@everywhere using ExoplanetsSysSim

import DataFrames.skipmissing

@everywhere include("clusters.jl")
@everywhere include("planetary_catalog.jl")
@everywhere include("optimization.jl")

@everywhere sim_param = setup_sim_param_model()





##### To start saving the model iterations in the optimization into a file:

model_name = "Clustered_P_R_broken_R"
use_KS_or_AD = "KS" #'KS' or 'AD' or 'Both' (need to be careful counting indices for 'dists_exclude'!!!)
AD_mod = true
Kep_or_Sim = "Kep" #'Kep' or 'Sim'
num_targs = 400030
dists_exclude = [3,4,8,12,13,15,16,17] #Int64[] if want to include all distances

data_table = CSV.read("Active_params_distances_table_best100000_every10.txt", delim=" ", allowmissing=:none)
n_params = length(make_vector_of_active_param_keys(sim_param))
params_keys = names(data_table)[1:n_params]
@assert all(make_vector_of_active_param_keys(sim_param) .== String.(params_keys))

params_array = convert(Matrix, data_table[1:end, params_keys])

file_name = model_name*"_recompute_optim_best100000_every10_targs"*string(num_targs)*".txt"

sendto(workers(), num_targs=num_targs, file_name=file_name)

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
@everywhere add_param_fixed(sim_param,"num_targets_sim_pass_one", num_targs)
@everywhere add_param_fixed(sim_param,"max_incl_sys", 0.0) #degrees; 0 (deg) for isotropic system inclinations; set closer to 90 (deg) for more transiting systems

active_param_true, weights, target_fitness, target_fitness_std = compute_weights_target_fitness_std_from_file("Weights1000_targs200015_maxincl60.txt", use_KS_or_AD ; weight=true, dists_exclude=dists_exclude, save_dist=true)





##### To recompute the model with the parameters in the table:

println(f, "# Active parameters: ", String.(params_keys))
println(f, "# Format: Active_params: [active parameter values]")
println(f, "# Format: Dist: [distances][total distance]")
println(f, "# Format: Dist_weighted: [weighted distances][total weighted distance]")
println(f, "# Distances used: ", use_KS_or_AD)
println(f, "# AD_mod: ", AD_mod)
println(f, "#")

Random.seed!()

t_elapsed = @elapsed begin
    @sync @distributed for i in 1:size(params_array,1)
        target_function(params_array[i,:], use_KS_or_AD, Kep_or_Sim ; AD_mod=AD_mod, weights=weights, all_dist=false, save_dist=true)
    end
end

println(f, "# elapsed time: ", t_elapsed, " seconds")
@everywhere close(f)
