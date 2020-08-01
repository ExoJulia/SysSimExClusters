dir_path = dirname(@__FILE__)

include(joinpath(dir_path, "../src/clusters.jl"))
include(joinpath(dir_path, "../src/planetary_catalog.jl"))

sim_param = setup_sim_param_model()
add_param_fixed(sim_param,"num_targets_sim_pass_one", 86760*5)
add_param_fixed(sim_param,"max_incl_sys", 0.)





##### To load a file with all the best active parameters from a set of optimization runs:

files_directory = "Sims/"
file_name = "Active_params_table.txt"

active_params_best_all = CSV.read(files_directory*file_name, delim=" ")
active_params_names = names(active_params_best_all)[2:end]





##### To simulate a catalog for each set of best active parameters:

for num_cat in 1:10 #1:size(active_params_best_all)[1]
    println("Generating simulated catalog ", num_cat)

    # To set up the model parameters:
    for (i,param_name) in enumerate(active_params_names)
        add_param_active(sim_param, string(param_name), active_params_best_all[num_cat, param_name])
    end

    run_number = active_params_best_all[num_cat, :run_number]

    @time cat_phys, cat_phys_cut, cat_obs, summary_stat = generate_and_save_physical_and_observed_catalogs(sim_param; save_path=files_directory, run_number=run_number)
end
