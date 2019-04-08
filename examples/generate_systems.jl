dir_path = dirname(@__FILE__)

include(joinpath(dir_path, "../src/clusters.jl"))
sim_param = setup_sim_param_model()





##### To generate the underlying systems:

add_param_fixed(sim_param,"num_targets_sim_pass_one", 80006)
@time cat_phys = generate_kepler_physical_catalog(sim_param)

##### For saving the underlying/true planets/systems:

f = open("physical_catalog_planets.txt", "w")
write_model_params(f, sim_param)
println(f, "target_id star_id planet_mass planet_radius period ecc")
for (i,targ) in enumerate(cat_phys.target)
    if length(targ.sys) > 1 #this should never happen
        println("There is more than one system for a given target? Check index: ", i)
    end
    if length(targ.sys[1].planet) > 0
        for (j,planet) in enumerate(targ.sys[1].planet)
            println(f, join([i, targ.sys[1].star.id, planet.mass, planet.radius, targ.sys[1].orbit[j].P, targ.sys[1].orbit[j].ecc], " "))
        end
    end
end
close(f)


f = open("physical_catalog_stars.txt", "w")
write_model_params(f, sim_param)
println(f, "target_id star_id star_mass star_radius num_planets")
for (i,targ) in enumerate(cat_phys.target)
    if length(targ.sys) > 1 #this should never happen
        println("There is more than one system for a given target? Check index: ", i)
    end
    println(f, join([i, targ.sys[1].star.id, targ.sys[1].star.mass, targ.sys[1].star.radius, length(targ.sys[1].planet)], " "))
end
close(f)



f = open("periods_all.out", "w")
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

f = open("eccentricities_all.out", "w")
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

f = open("radii_all.out", "w")
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

f = open("masses_all.out", "w")
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

##### For saving the stellar properties of all the systems:

f = open("stellar_masses_with_planets.out", "w")
write_model_params(f, sim_param)
for (i,targ) in enumerate(cat_phys.target)
    if length(targ.sys[1].planet) > 0
        star_mass = targ.sys[1].star.mass #solar masses
        println(f, star_mass)
    end
end
close(f)

f = open("stellar_radii_with_planets.out", "w")
write_model_params(f, sim_param)
for (i,targ) in enumerate(cat_phys.target)
    if length(targ.sys[1].planet) > 0
        star_radius = targ.sys[1].star.radius #solar radii
        println(f, star_radius)
    end
end
close(f)



##### To generate the observed systems:

@time begin
    cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param)
    cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param)
    summary_stat = calc_summary_stats_model(cat_obs,sim_param)
end

##### For saving the observed planets/systems:

f = open("observed_catalog_planets.txt", "w")
write_model_params(f, sim_param)
println(f, "target_id star_id period depth duration")
for num_pl_in_sys in 1:length(summary_stat.cache["idx_n_tranets"])
    num_targets = length(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
    for (i,j) in enumerate(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
        for k in 1:num_pl_in_sys
            println(f, join([j, cat_phys.target[j].sys[1].star.id, map(ExoplanetsSysSim.period, cat_obs.target[j].obs)[k], map(ExoplanetsSysSim.depth, cat_obs.target[j].obs)[k], map(ExoplanetsSysSim.duration, cat_obs.target[j].obs)[k]], " "))
        end
    end
end
close(f)


f = open("observed_catalog_stars.txt", "w")
write_model_params(f, sim_param)
println(f, "target_id star_id star_mass star_radius num_obs_planets")
for num_pl_in_sys in 1:length(summary_stat.cache["idx_n_tranets"])
    num_targets = length(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
    for (i,j) in enumerate(summary_stat.cache["idx_n_tranets"][num_pl_in_sys])
        println(f, join([j, cat_phys.target[j].sys[1].star.id, cat_obs.target[j].star.mass, cat_obs.target[j].star.radius, length(cat_obs.target[j].obs)], " "))
    end
end
close(f)



f = open("periods.out", "w")
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

f = open("depths.out", "w")
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

f = open("durations.out", "w")
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

#=
f = open("eccentricities.out", "w")
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
=#

##### For saving the stellar properties of the systems with observed planets:

f = open("stellar_masses_obs.out", "w")
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

f = open("stellar_radii_obs.out", "w")
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
