dir_path = dirname(@__FILE__)

include(joinpath(dir_path, "../src/clusters.jl"))
include(joinpath(dir_path, "../src/planetary_catalog.jl"))

##### To generate one physical and observed catalog:

sim_param = setup_sim_param_model()
add_param_fixed(sim_param,"num_targets_sim_pass_one", 88912*5)

using Random
rng_seed = rand(1:10000)
println("Random seed: ", rng_seed)
Random.seed!(rng_seed)

cat_phys, cat_phys_cut, cat_obs, summary_stat = generate_and_save_physical_and_observed_catalogs(sim_param)
println("Finished generating physical and observed catalogs.")
