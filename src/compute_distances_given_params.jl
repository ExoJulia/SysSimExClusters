include("clusters.jl")
include("planetary_catalog.jl")
include("optimization.jl")
using Random

sim_param = setup_sim_param_model()





##### To start saving the model iterations in the optimization into a file:

model_name = "Clustered_P_R"
names_split = ["bluer", "redder"]
AD_mod = true
num_targs = 86760
max_incl_sys = 0.
dists_include_split = ["delta_f", "mult_CRPD_r", "periods_KS", "period_ratios_KS", "durations_KS", "duration_ratios_nonmmr_KS", "duration_ratios_mmr_KS", "depths_KS", "radius_ratios_KS"]
dists_include_all = ["delta_f", "mult_CRPD_r", "periods_KS", "period_ratios_KS", "durations_KS", "duration_ratios_nonmmr_KS", "duration_ratios_mmr_KS", "depths_KS", "radius_ratios_KS"]

data_table = CSV.read("../emulator/GP_files/durations_KS/Active_params_distances_table_best100000_every10.txt", delim=" ", allowmissing=:none)
n_params = length(make_vector_of_active_param_keys(sim_param))
params_keys = names(data_table)[1:n_params]
@assert all(make_vector_of_active_param_keys(sim_param) .== String.(params_keys))

params_array = convert(Matrix, data_table[1:end, params_keys])

run_number, runs = parse(Int64, ARGS[1]), parse(Int64, ARGS[2])
evals = Int(size(params_array,1)/runs)
start, stop = 1+(run_number-1)*evals, run_number*evals

file_name = model_name*"_recompute_optim_best100000_every10_evals$(start)to$(stop)_targs$(num_targs).txt"
f = open("durations_KS/"*file_name, "w")
println(f, "# All initial parameters:")
write_model_params(f, sim_param)





##### To split the Kepler data into redder and bluer halves:

bprp = stellar_catalog[!,:bp_rp] .- stellar_catalog[!,:e_bp_min_rp_interp]
med_bprp = median(bprp)
idx_bluer = collect(1:size(stellar_catalog,1))[bprp .< med_bprp]
idx_redder = collect(1:size(stellar_catalog,1))[bprp .>= med_bprp]
star_id_split = [idx_bluer, idx_redder]

cssck = calc_summary_stats_collection_Kepler(stellar_catalog, planet_catalog, names_split, star_id_split, sim_param)





##### To load a file with the weights:

# To simulate more observed planets for the subsequent model generations:
add_param_fixed(sim_param,"num_targets_sim_pass_one", num_targs)
add_param_fixed(sim_param,"max_incl_sys", max_incl_sys)

# To load and compute the weights, target distance, and target distance std from a precomputed file:
active_param_true, weights, target_fitness, target_fitness_std = compute_weights_target_fitness_std_from_file_split_samples("Clustered_P_R_split_stars_weights_ADmod_$(AD_mod)_targs88912_evals100_all_pairs.txt", 4950, sim_param; names_samples=names_split, dists_include_samples=[dists_include_split, dists_include_split], dists_include_all=dists_include_all, f=f)
weights_split = [weights["bluer"], weights["redder"]]





##### To recompute the model with the parameters in the table:

println(f, "# Active parameters: ", String.(params_keys))
println(f, "# AD_mod: ", AD_mod)
println(f, "# Distances used (split): ", dists_include_split)
println(f, "# Distances used (all): ", dists_include_all)
println(f, "#")
println(f, "# Format: Active_params: [active parameter values]")
println(f, "# Format: [sample] Counts: [observed multiplicities][total planets, total planet pairs]")
println(f, "# Format: [sample] d_used_keys: [names of distance terms]")
println(f, "# Format: [sample] d_used_vals: [distance terms][sum of distance terms]")
println(f, "# Format: [sample] d_used_vals_w: [weighted distance terms][sum of weighted distance terms]")
println(f, "#")

Random.seed!()

t_elapsed = @elapsed begin
    for i in start:stop #1:size(params_array,1)
        target_function_split_stars(params_array[i,:], sim_param; cssc_fit=cssck, dists_include_all=dists_include_all, weights_all=weights["all"], names_samples=names_split, dists_include_samples=[dists_include_split, dists_include_split], weights_samples=weights_split, AD_mod=AD_mod, f=f)
    end
end

println(f, "# elapsed time: ", t_elapsed, " seconds")
close(f)
