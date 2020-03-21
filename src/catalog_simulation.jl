##### To define functions for simulating physical catalogs and saving them:

function save_physical_catalog(cat_phys::KeplerPhysicalCatalog, sim_param::SimParam; save_path::String="", run_number::Union{String,Int64}="")

    f = open(joinpath(save_path, "physical_catalog$run_number.csv"), "w")
    write_model_params(f, sim_param)
    println(f, "target_id,star_id,planet_mass,planet_radius,clusterid,period,ecc,incl_mut,star_mass,star_radius")
    for (i,targ) in enumerate(cat_phys.target)
        if length(targ.sys) > 1 #this should never happen
            println("There is more than one system for a given target? Check index: ", i)
        end
        sys = targ.sys[1]
        if length(sys.planet) > 0
            incl_ref = sys.system_plane.incl
            Ω_ref = sys.system_plane.asc_node
            for (j,planet) in enumerate(sys.planet)
                incl_pl, Ω_pl = sys.orbit[j].incl, sys.orbit[j].asc_node
                inclmut_pl = calc_incl_spherical_cosine_law(incl_ref, incl_pl, Ω_pl-Ω_ref)
                println(f, join([i, sys.star.id, planet.mass, planet.radius, planet.id, sys.orbit[j].P, sys.orbit[j].ecc, inclmut_pl, sys.star.mass, sys.star.radius], ","))
            end
        end
    end
    close(f)
end

function save_physical_catalog_stars_only(cat_phys::KeplerPhysicalCatalog, sim_param::SimParam; save_path::String="", run_number::Union{String,Int64}="")

    f = open(joinpath(save_path, "physical_catalog_stars$run_number.csv"), "w")
    write_model_params(f, sim_param)
    println(f, "target_id,star_id,star_mass,star_radius,num_planets")
    for (i,targ) in enumerate(cat_phys.target)
        if length(targ.sys) > 1 #this should never happen
            println("There is more than one system for a given target? Check index: ", i)
        end
        sys = targ.sys[1]
        println(f, join([i, sys.star.id, sys.star.mass, sys.star.radius, length(sys.planet)], ","))
    end
    close(f)
end

function save_clusterids_all(cat_phys::KeplerPhysicalCatalog, sim_param::SimParam; save_path::String="", run_number::Union{String,Int64}="")

    f = open(joinpath(save_path, "clusterids_all$run_number.out"), "w")
    write_model_params(f, sim_param)
    for (i,targ) in enumerate(cat_phys.target)
        if length(targ.sys[1].planet) > 0
            clusterids_sys = Array{Int64}(undef, length(targ.sys[1].planet))
            for (j,planet) in enumerate(targ.sys[1].planet)
                clusterids_sys[j] = planet.id
            end
            println(f, clusterids_sys)
        end
    end
    close(f)
end

function save_periods_all(cat_phys::KeplerPhysicalCatalog, sim_param::SimParam; save_path::String="", run_number::Union{String,Int64}="")

    f = open(joinpath(save_path, "periods_all$run_number.out"), "w")
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
end

function save_eccentricities_all(cat_phys::KeplerPhysicalCatalog, sim_param::SimParam; save_path::String="", run_number::Union{String,Int64}="")

    f = open(joinpath(save_path, "eccentricities_all$run_number.out"), "w")
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
end

function save_mutualinclinations_all(cat_phys::KeplerPhysicalCatalog, sim_param::SimParam; save_path::String="", run_number::Union{String,Int64}="")

    f = open(joinpath(save_path, "mutualinclinations_all$run_number.out"), "w")
    write_model_params(f, sim_param)
    for (i,targ) in enumerate(cat_phys.target)
        if length(targ.sys[1].orbit) > 0
            incl_ref = targ.sys[1].system_plane.incl
            Ω_ref = targ.sys[1].system_plane.asc_node
            inclmut_sys = Array{Float64}(undef, length(targ.sys[1].orbit))
            for (j,planet) in enumerate(targ.sys[1].orbit)
                incl_pl, Ω_pl = planet.incl, planet.asc_node
                inclmut_sys[j] = calc_incl_spherical_cosine_law(incl_ref, incl_pl, Ω_pl-Ω_ref)
            end
            println(f, inclmut_sys)
        end
    end
    close(f)
end

function save_radii_all(cat_phys::KeplerPhysicalCatalog, sim_param::SimParam; save_path::String="", run_number::Union{String,Int64}="")

    f = open(joinpath(save_path, "radii_all$run_number.out"), "w")
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
end

function save_masses_all(cat_phys::KeplerPhysicalCatalog, sim_param::SimParam; save_path::String="", run_number::Union{String,Int64}="")

    f = open(joinpath(save_path, "masses_all$run_number.out"), "w")
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
end

function save_stellar_masses_with_planets_all(cat_phys::KeplerPhysicalCatalog, sim_param::SimParam; save_path::String="", run_number::Union{String,Int64}="")

    f = open(joinpath(save_path, "stellar_masses_with_planets$run_number.out"), "w")
    write_model_params(f, sim_param)
    for (i,targ) in enumerate(cat_phys.target)
        if length(targ.sys[1].planet) > 0
            star_mass = targ.sys[1].star.mass #solar masses
            println(f, star_mass)
        end
    end
    close(f)
end

function save_stellar_radii_with_planets_all(cat_phys::KeplerPhysicalCatalog, sim_param::SimParam; save_path::String="", run_number::Union{String,Int64}="")

    f = open(joinpath(save_path, "stellar_radii_with_planets$run_number.out"), "w")
    write_model_params(f, sim_param)
    for (i,targ) in enumerate(cat_phys.target)
        if length(targ.sys[1].planet) > 0
            star_radius = targ.sys[1].star.radius #solar radii
            println(f, star_radius)
        end
    end
    close(f)
end

function save_physical_catalog_given_cat_phys(cat_phys::KeplerPhysicalCatalog, sim_param::SimParam; save_path::String="", run_number::Union{String,Int64}="")

    save_physical_catalog(cat_phys, sim_param; save_path=save_path, run_number=run_number)
    save_physical_catalog_stars_only(cat_phys, sim_param; save_path=save_path, run_number=run_number)

    save_clusterids_all(cat_phys, sim_param; save_path=save_path, run_number=run_number)
    save_periods_all(cat_phys, sim_param; save_path=save_path, run_number=run_number)
    save_eccentricities_all(cat_phys, sim_param; save_path=save_path, run_number=run_number)
    save_mutualinclinations_all(cat_phys, sim_param; save_path=save_path, run_number=run_number)
    save_radii_all(cat_phys, sim_param; save_path=save_path, run_number=run_number)
    save_masses_all(cat_phys, sim_param; save_path=save_path, run_number=run_number)

    save_stellar_radii_with_planets_all(cat_phys, sim_param; save_path=save_path, run_number=run_number)
    save_stellar_masses_with_planets_all(cat_phys, sim_param; save_path=save_path, run_number=run_number)
end

function generate_and_save_physical_catalog(sim_param::SimParam; save_path::String="", run_number::Union{String,Int64}="")

    @time cat_phys = generate_kepler_physical_catalog(sim_param)

    save_physical_catalog_given_cat_phys(cat_phys, sim_param; save_path=save_path, run_number=run_number)

    return cat_phys
end





##### To define functions for simulating observed catalogs and saving them:

function save_observed_catalog(cat_phys::KeplerPhysicalCatalog, cat_obs::KeplerObsCatalog, summary_stat::CatalogSummaryStatistics, sim_param::SimParam; save_path::String="", run_number::Union{String,Int64}="")

    f = open(joinpath(save_path, "observed_catalog$run_number.csv"), "w")
    write_model_params(f, sim_param)
    println(f, "target_id,star_id,period,depth,duration,star_mass,star_radius")
    for num_pl_in_sys in 1:length(summary_stat.cache["idx_n_tranets"])
        num_targets = length(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
        for (i,j) in enumerate(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
            for k in 1:num_pl_in_sys
                println(f, join([j, cat_phys.target[j].sys[1].star.id, map(ExoplanetsSysSim.period, cat_obs.target[j].obs)[k], map(ExoplanetsSysSim.depth, cat_obs.target[j].obs)[k], map(ExoplanetsSysSim.duration, cat_obs.target[j].obs)[k], cat_obs.target[j].star.mass, cat_obs.target[j].star.radius], ","))
            end
        end
    end
    close(f)
end

function save_observed_catalog_stars_only(cat_phys::KeplerPhysicalCatalog, cat_obs::KeplerObsCatalog, summary_stat::CatalogSummaryStatistics, sim_param::SimParam; save_path::String="", run_number::Union{String,Int64}="")

    f = open(joinpath(save_path, "observed_catalog_stars$run_number.csv"), "w")
    write_model_params(f, sim_param)
    println(f, "target_id,star_id,star_mass,star_radius,num_obs_planets")
    for num_pl_in_sys in 1:length(summary_stat.cache["idx_n_tranets"])
        num_targets = length(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
        for (i,j) in enumerate(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
            println(f, join([j, cat_phys.target[j].sys[1].star.id, cat_obs.target[j].star.mass, cat_obs.target[j].star.radius, length(cat_obs.target[j].obs)], ","))
        end
    end
    close(f)
end

function save_periods(cat_obs::KeplerObsCatalog, summary_stat::CatalogSummaryStatistics, sim_param::SimParam; save_path::String="", run_number::Union{String,Int64}="")

    f = open(joinpath(save_path, "periods$run_number.out"), "w")
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
end

function save_depths(cat_obs::KeplerObsCatalog, summary_stat::CatalogSummaryStatistics, sim_param::SimParam; save_path::String="", run_number::Union{String,Int64}="")

    f = open(joinpath(save_path, "depths$run_number.out"), "w")
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
end

function save_durations(cat_obs::KeplerObsCatalog, summary_stat::CatalogSummaryStatistics, sim_param::SimParam; save_path::String="", run_number::Union{String,Int64}="")

    f = open(joinpath(save_path, "durations$run_number.out"), "w")
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
end

function save_eccentricities(cat_obs::KeplerObsCatalog, summary_stat::CatalogSummaryStatistics, sim_param::SimParam; save_path::String="", run_number::Union{String,Int64}="")

    f = open(joinpath(save_path, "eccentricities$run_number.out"), "w")
    for num_pl_in_sys in 1:length(summary_stat.cache["idx_n_tranets"])
        num_targets = length(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
        ecc_array = Array{Float64}(num_pl_in_sys,num_targets)
        for (i,j) in enumerate(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
            ecc_array[:,i] = map(ExoplanetsSysSim.eccentricity, cat_obs.target[j].obs)[1:num_pl_in_sys]
        end
        println(f,"# Eccentricities of systems with ", num_pl_in_sys, " detected planets.")
        if length(ecc_array) > 0
            println(f,ecc_array')
        end
    end
    close(f)
end

function save_stellar_masses_with_planets_obs(cat_obs::KeplerObsCatalog, summary_stat::CatalogSummaryStatistics, sim_param::SimParam; save_path::String="", run_number::Union{String,Int64}="")

    f = open(joinpath(save_path, "stellar_masses_obs$run_number.out"), "w")
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
end

function save_stellar_radii_with_planets_obs(cat_obs::KeplerObsCatalog, summary_stat::CatalogSummaryStatistics, sim_param::SimParam; save_path::String="", run_number::Union{String,Int64}="")

    f = open(joinpath(save_path, "stellar_radii_obs$run_number.out"), "w")
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
end

function save_observed_catalog_given_cat_phys_obs(cat_phys::KeplerPhysicalCatalog, cat_obs::KeplerObsCatalog, summary_stat::CatalogSummaryStatistics, sim_param::SimParam; save_path::String="", run_number::Union{String,Int64}="")

    save_observed_catalog(cat_phys, cat_obs, summary_stat, sim_param; save_path=save_path, run_number=run_number)
    save_observed_catalog_stars_only(cat_phys, cat_obs, summary_stat, sim_param; save_path=save_path, run_number=run_number)

    save_periods(cat_obs, summary_stat, sim_param; save_path=save_path, run_number=run_number)
    save_depths(cat_obs, summary_stat, sim_param; save_path=save_path, run_number=run_number)
    save_durations(cat_obs, summary_stat, sim_param; save_path=save_path, run_number=run_number)
    #save_eccentricities(cat_obs, summary_stat, sim_param; save_path=save_path, run_number=run_number)

    save_stellar_radii_with_planets_obs(cat_obs, summary_stat, sim_param; save_path=save_path, run_number=run_number)
    save_stellar_masses_with_planets_obs(cat_obs, summary_stat, sim_param; save_path=save_path, run_number=run_number)
end

function generate_and_save_observed_catalog_from_physical(cat_phys::KeplerPhysicalCatalog, sim_param::SimParam; save_path::String="", run_number::Union{String,Int64}="")

    @time begin
        cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param)
        cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param)
        summary_stat = calc_summary_stats_model(cat_obs,sim_param)
    end

    save_observed_catalog_given_cat_phys_obs(cat_phys, cat_obs, summary_stat, sim_param; save_path=save_path, run_number=run_number)

    return cat_phys_cut, cat_obs, summary_stat
end





function generate_and_save_physical_and_observed_catalogs(sim_param::SimParam; save_path::String="", run_number::Union{String,Int64}="")

    cat_phys = generate_and_save_physical_catalog(sim_param; save_path=save_path, run_number=run_number)
    cat_phys_cut, cat_obs, summary_stat = generate_and_save_observed_catalog_from_physical(cat_phys, sim_param; save_path=save_path, run_number=run_number)

    return cat_phys, cat_phys_cut, cat_obs, summary_stat
end
