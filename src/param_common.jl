if !@isdefined ExoplanetsSysSim
    using ExoplanetsSysSim
end

import Compat #: UTF8String, ASCIIString

##### Simulation_parameters:
function setup_sim_param_model(args::Vector{String} = Array{String}(undef, 0)) # Allow this to take a list of parameter (e.g., from command line)
    sim_param = SimParam()

    # How many targets to generate:
    add_param_fixed(sim_param,"num_targets_sim_pass_one", 79935) # Note this is used for the number of stars in the simulations, not necessarily related to number of Kepler targets
    add_param_fixed(sim_param,"num_kepler_targets", 79935) # Note this is used for the number of Kepler targets for the observational catalog

    # For generating target star properties:
    add_param_fixed(sim_param,"star_table_setup", ExoplanetsSysSim.StellarTable.setup_star_table)
    add_param_fixed(sim_param,"stellar_catalog", "q1q17_dr25_gaia_fgk.jld2")
    add_param_fixed(sim_param,"generate_kepler_target", ExoplanetsSysSim.generate_kepler_target_from_table)
    add_param_fixed(sim_param,"window_function", "DR25topwinfuncs.jld2")
    add_param_fixed(sim_param,"osd_file","dr25fgk_small_osds.jld2")
    #add_param_fixed(sim_param,"osd_file","dr25fgk_relaxcut_osds.jld2") # WARNING: need 8gb of memory to read this file

    # For generating planetary system properties:
    add_param_fixed(sim_param,"generate_planetary_system", generate_planetary_system_clustered) # For Non-clustered model: "generate_planetary_system_non_clustered"

    add_param_active(sim_param,"f_stars_with_planets_attempted_color_slope", 0.6)
    add_param_active(sim_param,"f_stars_with_planets_attempted_at_med_color", 0.6)
    add_param_fixed(sim_param,"med_color", 0.95)
    #add_param_active(sim_param,"f_stars_with_planets_attempted", 0.6)
    add_param_fixed(sim_param,"generate_num_clusters", generate_num_clusters_ZTP)
    add_param_fixed(sim_param,"generate_num_planets_in_cluster", generate_num_planets_in_cluster_ZTP)
    add_param_active(sim_param,"log_rate_clusters", log(2.0))
    add_param_fixed(sim_param,"max_clusters_in_sys", 10)
    add_param_active(sim_param,"log_rate_planets_per_cluster", log(2.3))
    add_param_fixed(sim_param,"max_planets_in_cluster", 10)

    # Generate_num_planets_in_cluster currently calls:
    add_param_fixed(sim_param,"generate_periods", ExoplanetsSysSim.generate_periods_power_law)
    add_param_fixed(sim_param,"generate_sizes", ExoplanetsSysSim.generate_sizes_broken_power_law) # To choose the way we draw planetary radii; if "generate_sizes_power_law", then takes "power_law_r"; if "generate_sizes_broken_power_law", then takes "power_law_r1", "power_law_r2", and "break_radius"
    add_param_active(sim_param,"power_law_P", 0.)
    #add_param_fixed(sim_param,"power_law_P1", 2.0)
    #add_param_fixed(sim_param,"power_law_P2", 0.5)
    #add_param_fixed(sim_param,"power_law_r", -2.5)
    add_param_active(sim_param,"power_law_r1", -1.4)
    add_param_active(sim_param,"power_law_r2", -4.8)
    add_param_fixed(sim_param,"min_period", 3.0)
    add_param_fixed(sim_param,"max_period", 300.0)
    #add_param_fixed(sim_param,"break_period", 10.0)
    add_param_fixed(sim_param,"min_radius", 0.5*ExoplanetsSysSim.earth_radius)
    add_param_fixed(sim_param,"max_radius", 10.0*ExoplanetsSysSim.earth_radius)
    add_param_fixed(sim_param,"break_radius", 3.0*ExoplanetsSysSim.earth_radius)

    # Generate_num_planets_in_cluster currently use these for the inclination distribution:
    add_param_fixed(sim_param,"resonance_width", 0.05)
    add_param_fixed(sim_param,"period_ratios_mmr", [2.0, 1.5, 4/3, 5/4])
    #add_param_active(sim_param,"f_high_incl", 0.4) # fraction of systems with higher mutual inclinations
    #add_param_active(sim_param,"sigma_incl", 50.) # degrees; 0 = coplanar w/ generate_kepler_target_simple; ignored by generate_planetary_system_uncorrelated_incl
    #add_param_active(sim_param,"sigma_incl_near_mmr", 1.3)

    #add_param_fixed(sim_param,"max_incl_sys", 0.0) # degrees; gives system inclinations from "max_incl_sys" (deg) to 90 (deg), so set to 0 for isotropic distribution of system inclinations; NOTE: make sure the difference between this and 90 (deg) is at least greater than "sigma_incl" and "sigma_incl_near_mmr"!

    # Generate_num_planets_in_cluster currently use these for the eccentricity distribution:
    add_param_fixed(sim_param,"generate_e_omega", ExoplanetsSysSim.generate_e_omega_rayleigh)
    #add_param_active(sim_param,"sigma_hk_color_slope", 0.04)
    #add_param_active(sim_param,"sigma_hk_at_med_color", 0.02)
    add_param_active(sim_param,"sigma_hk", 0.018)

    # Generate_num_planets_in_cluster currently use these for the stability tests:
    add_param_active(sim_param,"num_mutual_hill_radii", 8.0)
    add_param_fixed(sim_param,"f_amd_crit", 1.0) # fraction of critical AMD to distribute

    add_param_fixed(sim_param,"generate_planet_mass_from_radius", generate_planet_mass_from_radius_Ning2018_table) # "ExoplanetsSysSim.generate_planet_mass_from_radius_powerlaw" or "generate_planet_mass_from_radius_Ning2018" or "generate_planet_mass_from_radius_Ning2018_table"
    add_param_active(sim_param,"sigma_log_radius_in_cluster", 0.3)
    add_param_active(sim_param,"sigma_logperiod_per_pl_in_cluster", 0.2)

    # Functions to calculate observables from physical system properties:
    add_param_fixed(sim_param,"calc_target_obs_single_obs", ExoplanetsSysSim.calc_target_obs_single_obs)
    add_param_fixed(sim_param,"max_tranets_in_sys", 8) # SysSim ignores some planets in any systems with more than this many transiting planets to avoid wasting time on unphysical parameter values
    add_param_fixed(sim_param,"transit_noise_model", ExoplanetsSysSim.transit_noise_model_diagonal)
    #add_param_fixed(sim_param,"rng_seed",1234) # If you want to be able to reproduce simulations
    add_param_fixed(sim_param,"read_target_obs", ExoplanetsSysSim.simulated_read_kepler_observations) # Read Kepler observations to compare to from disk

    # Read any customizations, so its easy to keep them separate from the essentials:
    #=
    if isfile("param_custom.jl")
        include("param_custom.jl")
    end
    =#
    if @isdefined add_param_custom
        sim_param = add_param_custom(sim_param)
    end
    return sim_param
end

function test_setup_sim_param()
    setup_sim_param_model()
end



function write_model_params(f, sim_param::SimParam)
    # This function writes all the model parameters to a file f as a header
    println(f, "# num_targets_sim_pass_one: ", get_int(sim_param,"num_targets_sim_pass_one"))
    println(f, "# max_incl_sys: ", get_real(sim_param,"max_incl_sys"))
    #println(f, "# f_stars_with_planets_attempted: ", get_real(sim_param,"f_stars_with_planets_attempted"))
    println(f, "# f_stars_with_planets_attempted_at_med_color: ", get_real(sim_param,"f_stars_with_planets_attempted_at_med_color"))
    println(f, "# f_stars_with_planets_attempted_color_slope: ", get_real(sim_param,"f_stars_with_planets_attempted_color_slope"))
    println(f, "# generate_num_clusters: ", string(get_function(sim_param,"generate_num_clusters")))
    println(f, "# log_rate_clusters: ", get_real(sim_param,"log_rate_clusters"))
    println(f, "# max_clusters_in_sys: ", get_int(sim_param,"max_clusters_in_sys"))
    if string(get_function(sim_param,"generate_planetary_system")) == "generate_planetary_system_clustered"
        println(f, "# generate_num_planets_in_cluster: ", string(get_function(sim_param,"generate_num_planets_in_cluster")))
        println(f, "# log_rate_planets_per_cluster: ", get_real(sim_param,"log_rate_planets_per_cluster"))
        println(f, "# max_planets_in_clusters: ", get_int(sim_param,"max_planets_in_cluster"))
    end

    if string(get_function(sim_param,"generate_periods")) == "ExoplanetsSysSim.generate_periods_power_law"
        println(f, "# power_law_P: ", get_real(sim_param,"power_law_P"))
    elseif string(get_function(sim_param,"generate_periods")) == "ExoplanetsSysSim.generate_periods_broken_power_law"
        println(f, "# power_law_P1: ", get_real(sim_param,"power_law_P1"))
        println(f, "# power_law_P2: ", get_real(sim_param,"power_law_P2"))
        println(f, "# break_period: ", get_real(sim_param,"break_period"))
    end
    println(f, "# min_period: ", get_real(sim_param,"min_period"))
    println(f, "# max_period: ", get_real(sim_param,"max_period"))

    if string(get_function(sim_param,"generate_sizes")) == "ExoplanetsSysSim.generate_sizes_power_law"
        println(f, "# power_law_r: ", get_real(sim_param,"power_law_r"))
    elseif string(get_function(sim_param,"generate_sizes")) == "ExoplanetsSysSim.generate_sizes_broken_power_law"
        println(f, "# power_law_r1: ", get_real(sim_param,"power_law_r1"))
        println(f, "# power_law_r2: ", get_real(sim_param,"power_law_r2"))
        println(f, "# break_radius (R_earth): ", get_real(sim_param,"break_radius")/ExoplanetsSysSim.earth_radius)
    end
    println(f, "# min_radius (R_earth): ", get_real(sim_param,"min_radius")/ExoplanetsSysSim.earth_radius)
    println(f, "# max_radius (R_earth): ", get_real(sim_param,"max_radius")/ExoplanetsSysSim.earth_radius)

    #println(f, "# f_high_incl: ", get_real(sim_param,"f_high_incl"))
    #println(f, "# sigma_incl: ", get_real(sim_param,"sigma_incl"))
    #println(f, "# sigma_incl_near_mmr: ", get_real(sim_param,"sigma_incl_near_mmr"))
    println(f, "# sigma_hk: ", get_real(sim_param,"sigma_hk"))
    #println(f, "# sigma_hk_at_med_color: ", get_real(sim_param,"sigma_hk_at_med_color"))
    #println(f, "# sigma_hk_color_slope: ", get_real(sim_param,"sigma_hk_color_slope"))
    println(f, "# num_mutual_hill_radii: ", get_real(sim_param,"num_mutual_hill_radii"))
    println(f, "# f_amd_crit: ", get_real(sim_param,"f_amd_crit"))

    if string(get_function(sim_param,"generate_planet_mass_from_radius")) == "ExoplanetsSysSim.generate_planet_mass_from_radius_powerlaw"
        println(f, "# mr_power_index: ", get_real(sim_param,"mr_power_index"))
        println(f, "# mr_max_mass (M_earth): ", get_real(sim_param,"mr_max_mass")/ExoplanetsSysSim.earth_mass)
    elseif string(get_function(sim_param,"generate_planet_mass_from_radius")) == "generate_planet_mass_from_radius_Ning2018"
        println(f, "# mr_model: Ning2018")
    elseif string(get_function(sim_param,"generate_planet_mass_from_radius")) == "generate_planet_mass_from_radius_Ning2018_table"
        println(f, "# mr_model: Ning2018_table")
    end

    if string(get_function(sim_param,"generate_planetary_system")) == "generate_planetary_system_clustered"
        println(f, "# sigma_log_radius_in_cluster: ", get_real(sim_param,"sigma_log_radius_in_cluster"))
        println(f, "# sigma_logperiod_per_pl_in_cluster: ", get_real(sim_param,"sigma_logperiod_per_pl_in_cluster"))
    end
    println(f, "#")
end
