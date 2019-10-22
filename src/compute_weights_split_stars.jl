include("clusters.jl")
include("planetary_catalog.jl")
include("optimization.jl")

sim_param = setup_sim_param_model()





##### To start saving the model iterations in the optimization into a file:

sc_bluer = "q1q17_dr25_gaia_fgk_bluer.jld2"
sc_redder = "q1q17_dr25_gaia_fgk_redder.jld2"
AD_mod = true
num_targs = 79935*5
max_incl_sys = 0.
num_evals_weights = 1000

file_name = "Clustered_P_R_split_stars_weights_ADmod_$(AD_mod)_targs$(num_targs)_evals$(num_evals_weights).txt"

f = open(file_name, "w")
println(f, "# All initial parameters:")
write_model_params(f, sim_param)





##### To run the same model multiple times to see how it compares to a simulated catalog with the same parameters:

using Random
Random.seed!(1234) # to have the same reference catalog and simulated catalogs for calculating the weights

# To generate a simulated catalog to fit to:
add_param_fixed(sim_param,"num_targets_sim_pass_one", div(79935,2))

# Set up the summary statistics for the planets around the bluer half:
add_param_fixed(sim_param,"stellar_catalog", sc_bluer)
stellar_catalog = ExoplanetsSysSim.StellarTable.setup_star_table(sim_param; force_reread=true)
cat_phys = generate_kepler_physical_catalog(sim_param)
cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param)
cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param)

ss_ref_bluer = calc_summary_stats_model(cat_obs,sim_param)

# Set up the summary statistics for the planets around the redder half:
add_param_fixed(sim_param,"stellar_catalog", sc_redder)
stellar_catalog = ExoplanetsSysSim.StellarTable.setup_star_table(sim_param; force_reread=true)
cat_phys = generate_kepler_physical_catalog(sim_param)
cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param)
cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param)

ss_ref_redder = calc_summary_stats_model(cat_obs,sim_param)

# To simulate more observed planets for the subsequent model generations:
add_param_fixed(sim_param,"num_targets_sim_pass_one", div(num_targs,2))
add_param_fixed(sim_param,"max_incl_sys", max_incl_sys)

active_param_true, weights1, weights2, weightsc, target_fitness, target_fitness_std = compute_weights_target_fitness_std_perfect_model_split_stars(num_evals_weights, sim_param; stellar_catalog_all=[sc_bluer, sc_redder], ss_refs=[ss_ref_bluer, ss_ref_redder], AD_mod=AD_mod, f=f)

close(f)





# Test the loading of the weights file:

#dists_include_split = ["delta_f", "mult_CRPD_r", "pratios_KS", "durations_KS", "depths_KS", "rratios_KS"]
#dists_include_combined = ["delta_f", "mult_CRPD_r", "periods_KS", "pratios_KS", "durations_KS", "xis_KS", "xis_nonmmr_KS", "xis_mmr_KS", "depths_KS", "rratios_KS"]
#active_param_true, weights1, weights2, weightsc, target_fitness, target_fitness_std = compute_weights_target_fitness_std_from_file_split_samples(file_name, num_evals_weights, sim_param; stellar_catalog_all=[sc_bluer, sc_redder], dists_include_all=[dists_include_split, dists_include_split], dists_include_combined=dists_include_combined)
