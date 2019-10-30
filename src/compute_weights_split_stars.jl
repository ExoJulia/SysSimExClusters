include("clusters.jl")
include("planetary_catalog.jl")
include("optimization.jl")

sim_param = setup_sim_param_model()





##### To start saving the model iterations in the optimization into a file:

names_split = ["bluer", "redder"]
AD_mod = true
num_targs = 79935*5
max_incl_sys = 0.
num_evals_weights = 1000

file_name = "Clustered_P_R_split_stars_weights_ADmod_$(AD_mod)_targs$(num_targs)_evals$(num_evals_weights).txt"

f = open(file_name, "w")
println(f, "# All initial parameters:")
write_model_params(f, sim_param)





##### To split the Kepler data into redder and bluer halves:

bprp = stellar_catalog[:bp_rp]
med_bprp = median(bprp)
idx_bluer = collect(1:size(stellar_catalog,1))[bprp .< med_bprp]
idx_redder = collect(1:size(stellar_catalog,1))[bprp .>= med_bprp]
star_id_split = [idx_bluer, idx_redder]

cssck = calc_summary_stats_collection_Kepler(stellar_catalog, planet_catalog, names_split, star_id_split, sim_param)





##### To run the same model multiple times to see how it compares to a simulated catalog with the same parameters:

using Random
Random.seed!(1234) # to have the same reference catalog and simulated catalogs for calculating the weights

# To generate a simulated catalog to fit to:
add_param_fixed(sim_param,"num_targets_sim_pass_one", 79935)
cat_phys = generate_kepler_physical_catalog(sim_param)
cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param)
cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param)

cssc_ref = calc_summary_stats_collection_model(cat_obs, names_split, star_id_split, sim_param)

# To simulate more observed planets for the subsequent model generations:
add_param_fixed(sim_param,"num_targets_sim_pass_one", num_targs)
add_param_fixed(sim_param,"max_incl_sys", max_incl_sys)

active_param_true, weights, target_fitness, target_fitness_std = compute_weights_target_fitness_std_perfect_model_split_stars(num_evals_weights, sim_param; cssc_ref=cssc_ref, names_samples=names_split, AD_mod=AD_mod, f=f)

close(f)





# Test the loading of the weights file:

#dists_include_split = ["delta_f", "mult_CRPD_r", "period_ratios_KS", "durations_KS", "depths_KS", "radius_ratios_KS"]
#dists_include_all = ["delta_f", "mult_CRPD_r", "periods_KS", "period_ratios_KS", "durations_KS", "duration_ratios_KS", "duration_ratios_nonmmr_KS", "duration_ratios_mmr_KS", "depths_KS", "radius_ratios_KS"]
#active_param_true2, weights2, target_fitness2, target_fitness_std2 = compute_weights_target_fitness_std_from_file_split_samples(file_name, num_evals_weights, sim_param; names_samples=names_split, dists_include_samples=[dists_include_split, dists_include_split], dists_include_all=dists_include_all, save_dist=false)
