dir_path = dirname(@__FILE__)

include(joinpath(dir_path, "../src/clusters.jl"))

sim_param = setup_sim_param_model()
add_param_fixed(sim_param,"num_targets_sim_pass_one", 139232*5)
add_param_fixed(sim_param,"max_incl_sys", 0.) #degrees; 0 (deg) for isotropic system inclinations; set closer to 90 (deg) for more transiting systems





##### To load a file with all the best active parameters from a set of optimization runs:

files_directory = "Sims/"
file_name = "Active_params_table.txt"

active_params_best_all = CSV.read(files_directory*file_name, delim=" ")
active_params_names = names(active_params_best_all)[2:end]





##### To simulate a catalog for each set of best active parameters:

for num_cat in 1:10 #1:size(active_params_best_all)[1]
    println("Generating simulated catalog ", num_cat)
    tic()

    #To set up the model parameters:
    for (i,param_name) in enumerate(active_params_names)
        add_param_active(sim_param, string(param_name), active_params_best_all[num_cat, param_name])
    end

    save_name_end = string(active_params_best_all[num_cat, :run_number])*".out"


    #To generate the underlying systems:
    cat_phys = generate_kepler_physical_catalog(sim_param)

    #To save the underlying/true planets/systems:
    f = open(files_directory*"periods_all"*save_name_end, "w")
    write_model_params(f, sim_param)
    for (i,targ) in enumerate(cat_phys.target)
        if length(targ.sys[1].orbit) > 0
            periods_sys = Array{Float64}(undef, length(targ.sys[1].orbit))
            for (j,planet) in enumerate(targ.sys[1].orbit)
                periods_sys[j] = planet.P #days
            end
            println(f, periods_sys)
        end
    end
    close(f)

    f = open(files_directory*"eccentricities_all"*save_name_end, "w")
    write_model_params(f, sim_param)
    for (i,targ) in enumerate(cat_phys.target)
        if length(targ.sys[1].orbit) > 0
            ecc_sys = Array{Float64}(undef, length(targ.sys[1].orbit))
            for (j,planet) in enumerate(targ.sys[1].orbit)
                ecc_sys[j] = planet.ecc
            end
            println(f, ecc_sys)
        end
    end
    close(f)

    f = open(files_directory*"radii_all"*save_name_end, "w")
    write_model_params(f, sim_param)
    for (i,targ) in enumerate(cat_phys.target)
        if length(targ.sys[1].planet) > 0
            radii_sys = Array{Float64}(undef, length(targ.sys[1].planet))
            for (j,planet) in enumerate(targ.sys[1].planet)
                radii_sys[j] = planet.radius #solar radii
            end
            println(f, radii_sys)
        end
    end
    close(f)

    f = open(files_directory*"masses_all"*save_name_end, "w")
    write_model_params(f, sim_param)
    for (i,targ) in enumerate(cat_phys.target)
        if length(targ.sys[1].planet) > 0
            masses_sys = Array{Float64}(undef, length(targ.sys[1].planet))
            for (j,planet) in enumerate(targ.sys[1].planet)
                masses_sys[j] = planet.mass #solar masses
            end
            println(f, masses_sys)
        end
    end
    close(f)

    #To save the stellar properties of all the systems:
    f = open(files_directory*"stellar_masses_with_planets"*save_name_end, "w")
    write_model_params(f, sim_param)
    for (i,targ) in enumerate(cat_phys.target)
        if length(targ.sys[1].planet) > 0
            star_mass = targ.sys[1].star.mass #solar masses
            println(f, star_mass)
        end
    end
    close(f)

    f = open(files_directory*"stellar_radii_with_planets"*save_name_end, "w")
    write_model_params(f, sim_param)
    for (i,targ) in enumerate(cat_phys.target)
        if length(targ.sys[1].planet) > 0
            star_radius = targ.sys[1].star.radius #solar radii
            println(f, star_radius)
        end
    end
    close(f)



    #To generate the observed systems:
    cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param)
    cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param)
    summary_stat = calc_summary_stats_model(cat_obs,sim_param)

    #To save the observed planets/systems:
    f = open(files_directory*"periods"*save_name_end, "w")
    write_model_params(f, sim_param)
    for num_pl_in_sys in 1:length(summary_stat.cache["idx_n_tranets"])
        num_targets = length(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
        period_array = Array{Float64}(undef, num_pl_in_sys, num_targets)
        for (i,j) in enumerate(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
            period_array[:,i] = map(ExoplanetsSysSim.period, cat_obs.target[j].obs)[1:num_pl_in_sys]
        end
        println(f,"# Periods of systems with ", num_pl_in_sys, " detected planets.")
        if length(period_array) > 0
            println(f,period_array')
        end
    end
    close(f)

    f = open(files_directory*"depths"*save_name_end, "w")
    write_model_params(f, sim_param)
    for num_pl_in_sys in 1:length(summary_stat.cache["idx_n_tranets"])
        num_targets = length(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
        depths_array = Array{Float64}(undef, num_pl_in_sys, num_targets)
        for (i,j) in enumerate(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
            depths_array[:,i] = map(ExoplanetsSysSim.depth, cat_obs.target[j].obs)[1:num_pl_in_sys]
        end
        println(f,"# Transit depths of systems with ", num_pl_in_sys, " detected planets.")
        if length(depths_array) > 0
            println(f,depths_array')
        end
    end
    close(f)

    f = open(files_directory*"durations"*save_name_end, "w")
    write_model_params(f, sim_param)
    for num_pl_in_sys in 1:length(summary_stat.cache["idx_n_tranets"])
        num_targets = length(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
        duration_array = Array{Float64}(undef, num_pl_in_sys, num_targets)
        for (i,j) in enumerate(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
            duration_array[:,i] = map(ExoplanetsSysSim.duration, cat_obs.target[j].obs)[1:num_pl_in_sys]
        end
        println(f,"# Transit durations of systems with ", num_pl_in_sys, " detected planets.")
        if length(duration_array) > 0
            println(f,duration_array')
        end
    end
    close(f)

    #To save the stellar properties of the systems with observed planets:
    f = open(files_directory*"stellar_masses_obs"*save_name_end, "w")
    write_model_params(f, sim_param)
    for num_pl_in_sys in 1:length(summary_stat.cache["idx_n_tranets"])
        num_targets = length(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
        stellar_mass_array = Array{Float64}(undef, num_targets)
        for (i,j) in enumerate(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
            stellar_mass_array[i] = cat_obs.target[j].star.mass
        end
        println(f,"# Stellar masses of systems with ", num_pl_in_sys, " detected planets.")
        if length(stellar_mass_array) > 0
            println(f, stellar_mass_array)
        end
    end
    close(f)

    f = open(files_directory*"stellar_radii_obs"*save_name_end, "w")
    write_model_params(f, sim_param)
    for num_pl_in_sys in 1:length(summary_stat.cache["idx_n_tranets"])
        num_targets = length(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
        stellar_radius_array = Array{Float64}(undef, num_targets)
        for (i,j) in enumerate(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
            stellar_radius_array[i] = cat_obs.target[j].star.radius
        end
        println(f,"# Stellar radii of systems with ", num_pl_in_sys, " detected planets.")
        if length(stellar_radius_array) > 0
            println(f, stellar_radius_array)
        end
    end
    close(f)

    toc()
end
