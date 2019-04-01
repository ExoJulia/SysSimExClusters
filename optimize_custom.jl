import DataFrames.DataFrame
import DataArrays.skipmissing

include("clusters.jl")

sim_param = setup_sim_param_model()
add_param_fixed(sim_param,"num_targets_sim_pass_one",150061) #9)   # For "observed" data, use a realistic number of targets (after any cuts you want to perform)
add_param_fixed(sim_param,"max_incl_sys",80.0) #degrees; 0 (deg) for isotropic system inclinations; set closer to 90 (deg) for more transiting systems

const max_incl_sys = get_real(sim_param,"max_incl_sys")
cos_factor = cos(max_incl_sys*pi/180) #factor to divide the number of targets in simulation by to get the actual number of targets needed (with an isotropic distribution of system inclinations) to produce as many transiting systems for a single observer





##### To load the DR25 catalog and compute arrays of multiplicities, periods, period ratios, transit durations, transit depths, period-normalized transit duration ratios (xi), and transit depth ratios:

Q1Q17_DR25 = CSV.read("q1_q17_dr25_koi.tab_selectcols_new.csv", header=19, nullable=true)
Q1Q17_DR25_stellar = CSV.read("q1_q17_dr25_stellar_koi.tab_selectcols.csv", header=30, nullable=true)
Q1Q17_DR25_stellar_all = CSV.read("q1_q17_dr25_stellar_koi.tab_all.csv", header=1, nullable=true)

N_Kepler_targets = sum((Q1Q17_DR25_stellar_all[:teff] .> 4000.) .& (Q1Q17_DR25_stellar_all[:teff] .< 7000.) .& (Q1Q17_DR25_stellar_all[:logg] .> 4.))
println("Total number of Kepler targets satisfying our cuts: ", N_Kepler_targets)

table_confirmed = Q1Q17_DR25[(Q1Q17_DR25[:koi_disposition] .== "CONFIRMED") .| (Q1Q17_DR25[:koi_disposition] .== "CANDIDATE"), :] #Table containing only the confirmed and candidate objects
table_stellar = Q1Q17_DR25_stellar

table_confirmed = table_confirmed[(table_confirmed[:koi_period] .> get_real(sim_param,"min_period")) .& (table_confirmed[:koi_period] .< get_real(sim_param,"max_period")), :] #to make additional cuts in period P to be comparable to our simulated sample
table_confirmed = table_confirmed[(table_confirmed[:koi_prad] .> get_real(sim_param,"min_radius")/ExoplanetsSysSim.earth_radius) .& (table_confirmed[:koi_prad] .< get_real(sim_param,"max_radius")/ExoplanetsSysSim.earth_radius) .& (.~ismissing.(table_confirmed[:koi_prad])), :] #to make additional cuts in planetary radii to be comparable to our simulated sample

#To make cuts based on stellar properties of T_eff and logg:
teff_confirmed = zeros(Int64, size(table_confirmed,1)) #list to be filled with the T_eff (K) for the objects
logg_confirmed = zeros(Float64, size(table_confirmed,1)) #list to be filled with the logg(cgs) for the objects
cdpp5_confirmed = zeros(Float64, size(table_confirmed,1)) #list to be filled with the RMS CDPP 5h values for the objects
cdpp6_confirmed = zeros(Float64, size(table_confirmed,1)) #list to be filled with the RMS CDPP 6h values for the objects
for (i,KepID) in enumerate(table_confirmed[:kepid])
    teff_confirmed[i] = table_stellar[:teff][table_stellar[:kepid] .== KepID, :][1]
    logg_confirmed[i] = table_stellar[:logg][table_stellar[:kepid] .== KepID, :][1]
    cdpp5_confirmed[i] = table_stellar[:rrmscdpp05p0][table_stellar[:kepid] .== KepID, :][1]
    cdpp6_confirmed[i] = table_stellar[:rrmscdpp06p0][table_stellar[:kepid] .== KepID, :][1]
end

cdpp_cut = 250.
println("Fraction of CONFIRMED and CANDIDATE planets after CDPP cut: ", size(table_confirmed[(teff_confirmed .> 4000.) .& (teff_confirmed .< 7000.) .& (logg_confirmed .> 4.) .& (cdpp5_confirmed .< cdpp_cut), :], 1), "/", size(table_confirmed[(teff_confirmed .> 4000.) .& (teff_confirmed .< 7000.) .& (logg_confirmed .> 4.), :], 1))

table_confirmed = table_confirmed[(teff_confirmed .> 4000.) .& (teff_confirmed .< 7000.) .& (logg_confirmed .> 4.) .& (cdpp5_confirmed .< cdpp_cut), :]


KOI_systems = [x[1:6] for x in table_confirmed[:kepoi_name]]
checked_bools = zeros(size(table_confirmed,1)) #0's denote KOI that were not checked yet; 1's denote already checked KOI

M_confirmed = Int64[] #list to be filled with the planet multiplicities of the systems
R_confirmed = Float64[] #list to be filled with period ratios of adjacent planet pairs
xi_confirmed = Float64[] #list to be filled with the period-normalized transit duration ratios of adjacent planet pairs
D_ratio_confirmed = Float64[] #list to be filled with the transit depth ratios of adjacent planet pairs
P_confirmed = table_confirmed[:koi_period] #array of the periods (days)
P_confirmed = collect(skipmissing(P_confirmed))
t_D_confirmed = table_confirmed[:koi_duration] #array of the transit durations (hrs)
t_D_confirmed = collect(skipmissing(t_D_confirmed))
D_confirmed = table_confirmed[:koi_depth]/(1e6) #array of the transit depths (fraction)
D_confirmed = D_confirmed[ismissing.(D_confirmed) .== false] #to get rid of NA values

for i in 1:length(KOI_systems)
    if checked_bools[i] == 0 #if the KOI has not been checked (included while looking at another planet in the same system)
        system_i = (1:length(KOI_systems))[KOI_systems .== KOI_systems[i]]
        checked_bools[system_i] = 1

        #To get the periods and transit durations in this system:
        system_P = table_confirmed[:koi_period][system_i] #periods of all the planets in this system
        system_t_D = table_confirmed[:koi_duration][system_i] #transit durations of all the planets in this system
        system_D = table_confirmed[:koi_depth][system_i] #transit depths (in ppm) of all the planets in this system
        system_sort_i = sortperm(system_P) #indices that would sort the periods of the planets in this system
        system_P = system_P[system_sort_i] #periods of all the planets in this system, sorted
        system_t_D = system_t_D[system_sort_i] #transit durations of all the planets in this system, sorted by period
        system_D = system_D[system_sort_i] #transit depths of all the planets in this system, sorted by period

        #To count the total number of planets in this system:
        push!(M_confirmed, length(system_P))

        #To compute the period ratios, period-normalized transit duration ratios, and transit depth ratios in this system:
        system_R = system_P[2:end]./system_P[1:end-1] #period ratios of all the adjacent planet pairs in this system
        system_D_ratio = system_D[2:end]./system_D[1:end-1] #transit depth ratios of all the adjacent planet pairs in this system
        system_xi = (system_t_D[1:end-1]./system_t_D[2:end]).*(system_P[2:end]./system_P[1:end-1]).^(1//3) #period-normalized transit duration ratios of all the adjacent planet pairs in this system

        append!(R_confirmed, system_R)
        append!(D_ratio_confirmed, system_D_ratio)
        append!(xi_confirmed, system_xi)
    end
end





##### To run a set of simulations on a grid of two parameters, saving KS statistics:

#=
#To run our model on a 2D grid of 'power_law_P' and 'power_law_r':

power_law_P_range = linspace(-0.5, 0.5, 11)
power_law_r_range = linspace(-3., -1.2, 11)

model_fit_stats = zeros(6,length(power_law_P_range),length(power_law_r_range)) #3D array to be filled with all the statistics of the model fits: number of observed planets, and KS distances for multiplicities, periods, period ratios, transit depths, and period-normalized transit duration ratios
model_fit_stats_fields = ["Number of observed planets", "KS multiplicities", "KS periods", "KS period ratios", "KS transit depths", "KS period-normalized transit duration ratios"] #names of the fields associated with 'model_fit_stats' array

sim_param_default = setup_sim_param_model() #To set up all the default parameters

for (i,power_law_P) in enumerate(power_law_P_range)
    for (j,power_law_r) in enumerate(power_law_r_range)
        tic()
        sim_param_here = deepcopy(sim_param_default) #To set up all the default parameters
        add_param_fixed(sim_param_here,"power_law_P",power_law_P) #To replace the 'power_law_P' parameter
        add_param_fixed(sim_param_here,"power_law_r",power_law_r) #To replace the 'power_law_r' parameter

        cat_phys = generate_kepler_physical_catalog(sim_param_here)
        cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param_here)
        cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param_here)
        summary_stat = calc_summary_stats_model(cat_obs,sim_param_here)

        M_cat_obs = ones(Int64,0) #array to be filled with the number of transiting planets in each simulated system
        for k in 1:get_int(sim_param_here,"max_tranets_in_sys")
            append!(M_cat_obs, k*ones(Int64, summary_stat.stat["num n-tranet systems"][k]))
        end

        model_fit_stats[1,i,j] = summary_stat.stat["num_tranets"]
        model_fit_stats[2,i,j] = ksstats_ints(M_cat_obs, M_confirmed)[5]
        model_fit_stats[3,i,j] = ksstats(summary_stat.stat["P list"], P_confirmed)[5]
        model_fit_stats[4,i,j] = ksstats(summary_stat.stat["period_ratio_list"], R_confirmed)[5]
        model_fit_stats[5,i,j] = ksstats(summary_stat.stat["depth list"], D_confirmed)[5]
        model_fit_stats[6,i,j] = ksstats(summary_stat.stat["duration_ratio_list"], xi_confirmed)[5]

        println("# power_law_P, power_law_r: ", get_real(sim_param_here,"power_law_P"), ", ",  get_real(sim_param_here,"power_law_r"))
        toc()
    end
end
=#

#=
f = open("Model_stats_P_r_grid.out", "w")
println(f,"# Model statistics on a grid of (power_law_P, power_law_r) = (", power_law_P_range, ", ",  power_law_r_range, ")")
println(f,"# num_targets_sim_pass_one: ", get_int(sim_param_default,"num_targets_sim_pass_one"))
println(f,"# log_rate_clusters: ", get_real(sim_param_default,"log_rate_clusters"))
println(f,"# log_rate_planets_per_cluster: ", get_real(sim_param_default,"log_rate_planets_per_cluster"))
println(f,"# sigma_incl: ", get_real(sim_param_default,"sigma_incl"))
println(f,"# sigma_incl_near_mmr: ", get_real(sim_param_default,"sigma_incl_near_mmr"))
println(f,"# sigma_hk: ", get_real(sim_param_default,"sigma_hk"))
println(f,"# num_mutual_hill_radii: ", get_real(sim_param_default,"num_mutual_hill_radii"))
println(f,"# mr_power_index: ", get_real(sim_param_default,"mr_power_index"))
println(f,"# mr_max_mass: ", get_real(sim_param_default,"mr_max_mass"))
println(f,"# sigma_logperiod_per_pl_in_cluster: ", get_real(sim_param_default,"sigma_logperiod_per_pl_in_cluster"))
for i in 1:size(model_fit_stats)[1]
    println(f,"# ", model_fit_stats_fields[i], ":")
    println(f,model_fit_stats[i,:,:]')
end
close(f)
=#





#=
#To run our model on a 2D grid of 'num_mutual_hill_radii' and 'sigma_logperiod_per_pl_in_cluster':

num_mutual_hill_radii_range = linspace(5., 20., 11)
sigma_logperiod_per_pl_in_cluster_range = linspace(0.05, 0.5, 11)

model_fit_stats = zeros(8,length(num_mutual_hill_radii_range),length(sigma_logperiod_per_pl_in_cluster_range)) #3D array to be filled with all the statistics of the model fits: number of observed planets, and KS distances for multiplicities, periods, period ratios, transit durations, period-normalized transit duration ratios, transit depths, and transit depth ratios
model_fit_stats_fields = ["Number of observed planets", "KS multiplicities", "KS periods", "KS period ratios", "KS transit durations", "KS period-normalized transit duration ratios", "KS transit depths", "KS transit depth ratios"] #names of the fields associated with 'model_fit_stats' array

for (i,num_mutual_hill_radii) in enumerate(num_mutual_hill_radii_range)
    for (j,sigma_logperiod_per_pl_in_cluster) in enumerate(sigma_logperiod_per_pl_in_cluster_range)
        tic()
        sim_param_here = deepcopy(sim_param) #To set up all the default parameters
        add_param_fixed(sim_param_here,"num_mutual_hill_radii",num_mutual_hill_radii) #To replace the 'num_mutual_hill_radii' parameter
        add_param_fixed(sim_param_here,"sigma_logperiod_per_pl_in_cluster",sigma_logperiod_per_pl_in_cluster) #To replace the 'sigma_logperiod_per_pl_in_cluster' parameter

        cat_phys = generate_kepler_physical_catalog(sim_param_here)
        cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param_here)
        cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param_here)
        summary_stat = calc_summary_stats_model(cat_obs,sim_param_here)

        M_cat_obs = ones(Int64,0) #array to be filled with the number of transiting planets in each simulated system
        for k in 1:get_int(sim_param_here,"max_tranets_in_sys")
            append!(M_cat_obs, k*ones(Int64, summary_stat.stat["num n-tranet systems"][k]))
        end

        model_fit_stats[1,i,j] = summary_stat.stat["num_tranets"]
        model_fit_stats[2,i,j] = ksstats_ints(M_cat_obs, M_confirmed)[5]
        model_fit_stats[3,i,j] = ksstats(summary_stat.stat["P list"], P_confirmed)[5]
        model_fit_stats[4,i,j] = ksstats(summary_stat.stat["period_ratio_list"], R_confirmed)[5]
        model_fit_stats[5,i,j] = ksstats(summary_stat.stat["duration list"].*24, t_D_confirmed)[5] #transit durations in simulations are in days, while in the Kepler catalog are in hours
        model_fit_stats[6,i,j] = ksstats(summary_stat.stat["duration_ratio_list"], xi_confirmed)[5]
        model_fit_stats[7,i,j] = ksstats(summary_stat.stat["depth list"], D_confirmed)[5]
        model_fit_stats[8,i,j] = ksstats(summary_stat.stat["radius_ratio_list"].^2, D_ratio_confirmed)[5] #simulations save radius ratios while we computed transit duration ratios from the Kepler catalog

        println("# num_mutual_hill_radii, sigma_logperiod_per_pl_in_cluster: ", get_real(sim_param_here,"num_mutual_hill_radii"), ", ",  get_real(sim_param_here,"sigma_logperiod_per_pl_in_cluster"))
        toc()
    end
end
=#

#=
f = open("Model_stats_Nhill_sigmaP_grid.out", "w")
println(f,"# Model statistics on a grid of (num_mutual_hill_radii, sigma_logperiod_per_pl_in_cluster) = (", num_mutual_hill_radii_range, ", ",  sigma_logperiod_per_pl_in_cluster_range, ")")
println(f, "# num_targets_sim_pass_one: ", get_int(sim_param,"num_targets_sim_pass_one"))
println(f, "# max_incl_sys: ", get_real(sim_param,"max_incl_sys"))
println(f, "# log_rate_clusters: ", get_real(sim_param,"log_rate_clusters"))
println(f, "# max_clusters_in_sys: ", get_int(sim_param,"max_clusters_in_sys"))
println(f, "# log_rate_planets_per_cluster: ", get_real(sim_param,"log_rate_planets_per_cluster"))
println(f, "# max_planets_in_clusters: ", get_int(sim_param,"max_planets_in_cluster"))
println(f, "# power_law_P: ", get_real(sim_param,"power_law_P"))
println(f, "# min_period: ", get_real(sim_param,"min_period"))
println(f, "# max_period: ", get_real(sim_param,"max_period"))
println(f, "# power_law_r1: ", get_real(sim_param,"power_law_r1"))
println(f, "# power_law_r2: ", get_real(sim_param,"power_law_r2"))
println(f, "# min_radius (R_earth): ", get_real(sim_param,"min_radius")/ExoplanetsSysSim.earth_radius)
println(f, "# max_radius (R_earth): ", get_real(sim_param,"max_radius")/ExoplanetsSysSim.earth_radius)
println(f, "# break_radius (R_earth): ", get_real(sim_param,"break_radius")/ExoplanetsSysSim.earth_radius)
println(f, "# sigma_incl: ", get_real(sim_param,"sigma_incl"))
println(f, "# sigma_incl_near_mmr: ", get_real(sim_param,"sigma_incl_near_mmr"))
println(f, "# sigma_hk: ", get_real(sim_param,"sigma_hk"))
#####println(f, "# num_mutual_hill_radii: ", get_real(sim_param,"num_mutual_hill_radii"))
println(f, "# mr_power_index: ", get_real(sim_param,"mr_power_index"))
println(f, "# mr_max_mass (M_earth): ", get_real(sim_param,"mr_max_mass")/ExoplanetsSysSim.earth_mass)
println(f, "# sigma_log_radius_in_cluster: ", get_real(sim_param,"sigma_log_radius_in_cluster"))
#####println(f, "# sigma_logperiod_per_pl_in_cluster: ", get_real(sim_param,"sigma_logperiod_per_pl_in_cluster"))
for i in 1:size(model_fit_stats)[1]
    println(f,"# ", model_fit_stats_fields[i], ":")
    println(f,model_fit_stats[i,:,:]')
end
close(f)
=#





#=
#To run our model on a 2D grid of 'sigma_incl' and 'sigma_hk':

sigma_incl_range = linspace(0., 5., 11)
sigma_hk_range = linspace(0., 0.2, 11)

model_fit_stats = zeros(8,length(sigma_incl_range),length(sigma_hk_range)) #3D array to be filled with all the statistics of the model fits: number of observed planets, and KS distances for multiplicities, periods, period ratios, transit durations, period-normalized transit duration ratios, transit depths, and transit depth ratios
model_fit_stats_fields = ["Number of observed planets", "KS multiplicities", "KS periods", "KS period ratios", "KS transit durations", "KS period-normalized transit duration ratios", "KS transit depths", "KS transit depth ratios"] #names of the fields associated with 'model_fit_stats' array

for (i,sigma_incl) in enumerate(sigma_incl_range)
    for (j,sigma_hk) in enumerate(sigma_hk_range)
        tic()
        sim_param_here = deepcopy(sim_param) #To set up all the default parameters
        add_param_fixed(sim_param_here,"sigma_incl",sigma_incl) #To replace the 'sigma_incl' parameter
        add_param_fixed(sim_param_here,"sigma_hk",sigma_hk) #To replace the 'sigma_hk' parameter

        cat_phys = generate_kepler_physical_catalog(sim_param_here)
        cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param_here)
        cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param_here)
        summary_stat = calc_summary_stats_model(cat_obs,sim_param_here)

        M_cat_obs = ones(Int64,0) #array to be filled with the number of transiting planets in each simulated system
        for k in 1:get_int(sim_param_here,"max_tranets_in_sys")
            append!(M_cat_obs, k*ones(Int64, summary_stat.stat["num n-tranet systems"][k]))
        end

        model_fit_stats[1,i,j] = summary_stat.stat["num_tranets"]
        model_fit_stats[2,i,j] = ksstats_ints(M_cat_obs, M_confirmed)[5]
        model_fit_stats[3,i,j] = ksstats(summary_stat.stat["P list"], P_confirmed)[5]
        model_fit_stats[4,i,j] = ksstats(summary_stat.stat["period_ratio_list"], R_confirmed)[5]
        model_fit_stats[5,i,j] = ksstats(summary_stat.stat["duration list"].*24, t_D_confirmed)[5] #transit durations in simulations are in days, while in the Kepler catalog are in hours
        model_fit_stats[6,i,j] = ksstats(summary_stat.stat["duration_ratio_list"], xi_confirmed)[5]
        model_fit_stats[7,i,j] = ksstats(summary_stat.stat["depth list"], D_confirmed)[5]
        model_fit_stats[8,i,j] = ksstats(summary_stat.stat["radius_ratio_list"].^2, D_ratio_confirmed)[5] #simulations save radius ratios while we computed transit duration ratios from the Kepler catalog

        println("# sigma_incl, sigma_hk: ", get_real(sim_param_here,"sigma_incl"), ", ",  get_real(sim_param_here,"sigma_hk"))
        toc()
    end
end
=#

#=
f = open("Model_stats_sigma_incl_sigma_hk_grid.out", "w")
println(f,"# Model statistics on a grid of (sigma_incl, sigma_hk) = (", sigma_incl_range, ", ",  sigma_hk_range, ")")
println(f, "# num_targets_sim_pass_one: ", get_int(sim_param,"num_targets_sim_pass_one"))
println(f, "# max_incl_sys: ", get_real(sim_param,"max_incl_sys"))
println(f, "# log_rate_clusters: ", get_real(sim_param,"log_rate_clusters"))
println(f, "# max_clusters_in_sys: ", get_int(sim_param,"max_clusters_in_sys"))
println(f, "# log_rate_planets_per_cluster: ", get_real(sim_param,"log_rate_planets_per_cluster"))
println(f, "# max_planets_in_clusters: ", get_int(sim_param,"max_planets_in_cluster"))
println(f, "# power_law_P: ", get_real(sim_param,"power_law_P"))
println(f, "# min_period: ", get_real(sim_param,"min_period"))
println(f, "# max_period: ", get_real(sim_param,"max_period"))
println(f, "# power_law_r1: ", get_real(sim_param,"power_law_r1"))
println(f, "# power_law_r2: ", get_real(sim_param,"power_law_r2"))
println(f, "# min_radius (R_earth): ", get_real(sim_param,"min_radius")/ExoplanetsSysSim.earth_radius)
println(f, "# max_radius (R_earth): ", get_real(sim_param,"max_radius")/ExoplanetsSysSim.earth_radius)
println(f, "# break_radius (R_earth): ", get_real(sim_param,"break_radius")/ExoplanetsSysSim.earth_radius)
#####println(f, "# sigma_incl: ", get_real(sim_param,"sigma_incl"))
println(f, "# sigma_incl_near_mmr: ", get_real(sim_param,"sigma_incl_near_mmr"))
#####println(f, "# sigma_hk: ", get_real(sim_param,"sigma_hk"))
println(f, "# num_mutual_hill_radii: ", get_real(sim_param,"num_mutual_hill_radii"))
println(f, "# mr_power_index: ", get_real(sim_param,"mr_power_index"))
println(f, "# mr_max_mass (M_earth): ", get_real(sim_param,"mr_max_mass")/ExoplanetsSysSim.earth_mass)
println(f, "# sigma_log_radius_in_cluster: ", get_real(sim_param,"sigma_log_radius_in_cluster"))
println(f, "# sigma_logperiod_per_pl_in_cluster: ", get_real(sim_param,"sigma_logperiod_per_pl_in_cluster"))
for i in 1:size(model_fit_stats)[1]
    println(f,"# ", model_fit_stats_fields[i], ":")
    println(f,model_fit_stats[i,:,:]')
end
close(f)
=#





#=
#To run our model on a 2D grid of 'mr_power_index' and 'sigma_log_radius_in_cluster':

mr_power_index_range = linspace(1.5, 3.5, 11)
sigma_log_radius_in_cluster_range = linspace(0.1, 1.0, 11)

model_fit_stats = zeros(8,length(mr_power_index_range),length(sigma_log_radius_in_cluster_range)) #3D array to be filled with all the statistics of the model fits: number of observed planets, and KS distances for multiplicities, periods, period ratios, transit durations, period-normalized transit duration ratios, transit depths, and transit depth ratios
model_fit_stats_fields = ["Number of observed planets", "KS multiplicities", "KS periods", "KS period ratios", "KS transit durations", "KS period-normalized transit duration ratios", "KS transit depths", "KS transit depth ratios"] #names of the fields associated with 'model_fit_stats' array

for (i,max_radius) in enumerate(mr_power_index_range)
    for (j,mr_power_index) in enumerate(sigma_log_radius_in_cluster_range)
        tic()
        sim_param_here = deepcopy(sim_param) #To set up all the default parameters
        add_param_fixed(sim_param_here,"mr_power_index",mr_power_index) #To replace the 'mr_power_index' parameter
        add_param_fixed(sim_param_here,"sigma_log_radius_in_cluster",sigma_log_radius_in_cluster) #To replace the 'sigma_log_radius_in_cluster' parameter

        cat_phys = generate_kepler_physical_catalog(sim_param_here)
        cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param_here)
        cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param_here)
        summary_stat = calc_summary_stats_model(cat_obs,sim_param_here)

        M_cat_obs = ones(Int64,0) #array to be filled with the number of transiting planets in each simulated system
        for k in 1:get_int(sim_param_here,"max_tranets_in_sys")
            append!(M_cat_obs, k*ones(Int64, summary_stat.stat["num n-tranet systems"][k]))
        end

        model_fit_stats[1,i,j] = summary_stat.stat["num_tranets"]
        model_fit_stats[2,i,j] = ksstats_ints(M_cat_obs, M_confirmed)[5]
        model_fit_stats[3,i,j] = ksstats(summary_stat.stat["P list"], P_confirmed)[5]
        model_fit_stats[4,i,j] = ksstats(summary_stat.stat["period_ratio_list"], R_confirmed)[5]
        model_fit_stats[5,i,j] = ksstats(summary_stat.stat["duration list"].*24, t_D_confirmed)[5] #transit durations in simulations are in days, while in the Kepler catalog are in hours
        model_fit_stats[6,i,j] = ksstats(summary_stat.stat["duration_ratio_list"], xi_confirmed)[5]
        model_fit_stats[7,i,j] = ksstats(summary_stat.stat["depth list"], D_confirmed)[5]
        model_fit_stats[8,i,j] = ksstats(summary_stat.stat["radius_ratio_list"].^2, D_ratio_confirmed)[5] #simulations save radius ratios while we computed transit duration ratios from the Kepler catalog

        println("# mr_power_index, sigma_log_radius_in_cluster: ", get_real(sim_param_here,"mr_power_index"), ", ",  get_real(sim_param_here,"sigma_log_radius_in_cluster"))
        toc()
    end
end
=#

#=
f = open("Model_stats_mr_sigmaR_grid.out", "w")
println(f,"# Model statistics on a grid of (mr_power_index,sigma_log_radius_in_cluster) = (", mr_power_index_range, ", ",  sigma_log_radius_in_cluster_range, ")")
println(f, "# num_targets_sim_pass_one: ", get_int(sim_param,"num_targets_sim_pass_one"))
println(f, "# max_incl_sys: ", get_real(sim_param,"max_incl_sys"))
println(f, "# log_rate_clusters: ", get_real(sim_param,"log_rate_clusters"))
println(f, "# max_clusters_in_sys: ", get_int(sim_param,"max_clusters_in_sys"))
println(f, "# log_rate_planets_per_cluster: ", get_real(sim_param,"log_rate_planets_per_cluster"))
println(f, "# max_planets_in_clusters: ", get_int(sim_param,"max_planets_in_cluster"))
println(f, "# power_law_P: ", get_real(sim_param,"power_law_P"))
println(f, "# min_period: ", get_real(sim_param,"min_period"))
println(f, "# max_period: ", get_real(sim_param,"max_period"))
println(f, "# power_law_r1: ", get_real(sim_param,"power_law_r1"))
println(f, "# power_law_r2: ", get_real(sim_param,"power_law_r2"))
println(f, "# min_radius (R_earth): ", get_real(sim_param,"min_radius")/ExoplanetsSysSim.earth_radius)
println(f, "# max_radius (R_earth): ", get_real(sim_param,"max_radius")/ExoplanetsSysSim.earth_radius)
println(f, "# break_radius (R_earth): ", get_real(sim_param,"break_radius")/ExoplanetsSysSim.earth_radius)
println(f, "# sigma_incl: ", get_real(sim_param,"sigma_incl"))
println(f, "# sigma_incl_near_mmr: ", get_real(sim_param,"sigma_incl_near_mmr"))
println(f, "# sigma_hk: ", get_real(sim_param,"sigma_hk"))
println(f, "# num_mutual_hill_radii: ", get_real(sim_param,"num_mutual_hill_radii"))
#####println(f, "# mr_power_index: ", get_real(sim_param,"mr_power_index"))
println(f, "# mr_max_mass (M_earth): ", get_real(sim_param,"mr_max_mass")/ExoplanetsSysSim.earth_mass)
#####println(f, "# sigma_log_radius_in_cluster: ", get_real(sim_param,"sigma_log_radius_in_cluster"))
println(f, "# sigma_logperiod_per_pl_in_cluster: ", get_real(sim_param,"sigma_logperiod_per_pl_in_cluster"))
for i in 1:size(model_fit_stats)[1]
    println(f,"# ", model_fit_stats_fields[i], ":")
    println(f,model_fit_stats[i,:,:]')
end
close(f)
=#





#=
#To run our model on a 2D grid of 'log_rate_clusters' and 'log_rate_planets_per_cluster':

log_rate_clusters_range = log.(linspace(1., 5., 11))
log_rate_planets_per_cluster_range = log.(linspace(1., 5., 11))

model_fit_stats = zeros(9,length(log_rate_clusters_range),length(log_rate_planets_per_cluster_range)) #3D array to be filled with all the statistics of the model fits: number of observed planets, and KS distances for multiplicities, periods, period ratios, transit durations, period-normalized transit duration ratios, transit depths, and transit depth ratios
model_fit_stats_fields = ["Number of observed planets", "Number of systems with observed singles", "KS multiplicities", "KS periods", "KS period ratios", "KS transit durations", "KS period-normalized transit duration ratios", "KS transit depths", "KS transit depth ratios"] #names of the fields associated with 'model_fit_stats' array

for (i,log_rate_clusters) in enumerate(log_rate_clusters_range)
    for (j,log_rate_planets_per_cluster) in enumerate(log_rate_planets_per_cluster_range)
        tic()
        sim_param_here = deepcopy(sim_param) #To set up all the default parameters
        add_param_fixed(sim_param_here,"log_rate_clusters",log_rate_clusters) #To replace the 'log_rate_clusters' parameter
        add_param_fixed(sim_param_here,"log_rate_planets_per_cluster",log_rate_planets_per_cluster) #To replace the 'log_rate_planets_per_cluster' parameter

        cat_phys = generate_kepler_physical_catalog(sim_param_here)
        cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param_here)
        cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param_here)
        summary_stat = calc_summary_stats_model(cat_obs,sim_param_here)

        M_cat_obs = ones(Int64,0) #array to be filled with the number of transiting planets in each simulated system
        for k in 1:get_int(sim_param_here,"max_tranets_in_sys")
            append!(M_cat_obs, k*ones(Int64, summary_stat.stat["num n-tranet systems"][k]))
        end

        model_fit_stats[1,i,j] = summary_stat.stat["num_tranets"]
        model_fit_stats[2,i,j] = summary_stat.stat["num n-tranet systems"][1]
        model_fit_stats[3,i,j] = ksstats_ints(M_cat_obs, M_confirmed)[5]
        model_fit_stats[4,i,j] = ksstats(summary_stat.stat["P list"], P_confirmed)[5]
        model_fit_stats[5,i,j] = ksstats(summary_stat.stat["period_ratio_list"], R_confirmed)[5]
        model_fit_stats[6,i,j] = ksstats(summary_stat.stat["duration list"].*24, t_D_confirmed)[5] #transit durations in simulations are in days, while in the Kepler catalog are in hours
        model_fit_stats[7,i,j] = ksstats(summary_stat.stat["duration_ratio_list"], xi_confirmed)[5]
        model_fit_stats[8,i,j] = ksstats(summary_stat.stat["depth list"], D_confirmed)[5]
        model_fit_stats[9,i,j] = ksstats(summary_stat.stat["radius_ratio_list"].^2, D_ratio_confirmed)[5] #simulations save radius ratios while we computed transit duration ratios from the Kepler catalog

        println("# exp(log_rate_clusters), exp(log_rate_planets_per_cluster): ", exp(get_real(sim_param_here,"log_rate_clusters")), ", ",  exp(get_real(sim_param_here,"log_rate_planets_per_cluster")))
        toc()
    end
end
=#

#=
f = open("Model_stats_Nc_Np_grid.out", "w")
println(f,"# Model statistics on a grid of (log_rate_clusters, log_rate_planets_per_cluster) = (", log_rate_clusters_range, ", ",  log_rate_planets_per_cluster_range, ")")
println(f, "# num_targets_sim_pass_one: ", get_int(sim_param,"num_targets_sim_pass_one"))
println(f, "# max_incl_sys: ", get_real(sim_param,"max_incl_sys"))
#####println(f, "# log_rate_clusters: ", get_real(sim_param,"log_rate_clusters"))
println(f, "# max_clusters_in_sys: ", get_int(sim_param,"max_clusters_in_sys"))
#####println(f, "# log_rate_planets_per_cluster: ", get_real(sim_param,"log_rate_planets_per_cluster"))
println(f, "# max_planets_in_clusters: ", get_int(sim_param,"max_planets_in_cluster"))
println(f, "# power_law_P: ", get_real(sim_param,"power_law_P"))
println(f, "# min_period: ", get_real(sim_param,"min_period"))
println(f, "# max_period: ", get_real(sim_param,"max_period"))
println(f, "# power_law_r1: ", get_real(sim_param,"power_law_r1"))
println(f, "# power_law_r2: ", get_real(sim_param,"power_law_r2"))
println(f, "# min_radius (R_earth): ", get_real(sim_param,"min_radius")/ExoplanetsSysSim.earth_radius)
println(f, "# max_radius (R_earth): ", get_real(sim_param,"max_radius")/ExoplanetsSysSim.earth_radius)
println(f, "# break_radius (R_earth): ", get_real(sim_param,"break_radius")/ExoplanetsSysSim.earth_radius)
println(f, "# sigma_incl: ", get_real(sim_param,"sigma_incl"))
println(f, "# sigma_incl_near_mmr: ", get_real(sim_param,"sigma_incl_near_mmr"))
println(f, "# sigma_hk: ", get_real(sim_param,"sigma_hk"))
println(f, "# num_mutual_hill_radii: ", get_real(sim_param,"num_mutual_hill_radii"))
println(f, "# mr_power_index: ", get_real(sim_param,"mr_power_index"))
println(f, "# mr_max_mass (M_earth): ", get_real(sim_param,"mr_max_mass")/ExoplanetsSysSim.earth_mass)
println(f, "# sigma_log_radius_in_cluster: ", get_real(sim_param,"sigma_log_radius_in_cluster"))
println(f, "# sigma_logperiod_per_pl_in_cluster: ", get_real(sim_param,"sigma_logperiod_per_pl_in_cluster"))
for i in 1:size(model_fit_stats)[1]
    println(f,"# ", model_fit_stats_fields[i], ":")
    println(f,model_fit_stats[i,:,:]')
end
close(f)
=#





##### The following optimization code is for exploring our updated model with a broken power law for the planet radii distribution instead of a single power law:

#=
#To run our model on a 2D grid of 'power_law_r1' and 'power_law_r2':

power_law_r1_range = linspace(-6., 3., 11)
power_law_r2_range = linspace(-6., 3., 11)

model_fit_stats = zeros(8,length(power_law_r1_range),length(power_law_r2_range)) #3D array to be filled with all the statistics of the model fits: number of observed planets, and KS distances for multiplicities, periods, period ratios, transit durations, period-normalized transit duration ratios, transit depths, and transit depth ratios
model_fit_stats_fields = ["Number of observed planets", "KS multiplicities", "KS periods", "KS period ratios", "KS transit durations", "KS period-normalized transit duration ratios", "KS transit depths", "KS transit depth ratios"] #names of the fields associated with 'model_fit_stats' array

for (i,power_law_r1) in enumerate(power_law_r1_range)
    for (j,power_law_r2) in enumerate(power_law_r2_range)
        tic()
        sim_param_here = deepcopy(sim_param) #To set up all the default parameters
        add_param_fixed(sim_param_here,"power_law_r1",power_law_r1) #To replace the 'power_law_r1' parameter
        add_param_fixed(sim_param_here,"power_law_r2",power_law_r2) #To replace the 'power_law_r2' parameter

        cat_phys = generate_kepler_physical_catalog(sim_param_here)
        cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param_here)
        cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param_here)
        summary_stat = calc_summary_stats_model(cat_obs,sim_param_here)

        M_cat_obs = ones(Int64,0) #array to be filled with the number of transiting planets in each simulated system
        for k in 1:get_int(sim_param_here,"max_tranets_in_sys")
            append!(M_cat_obs, k*ones(Int64, summary_stat.stat["num n-tranet systems"][k]))
        end

        model_fit_stats[1,i,j] = summary_stat.stat["num_tranets"]
        model_fit_stats[2,i,j] = ksstats_ints(M_cat_obs, M_confirmed)[5]
        model_fit_stats[3,i,j] = ksstats(summary_stat.stat["P list"], P_confirmed)[5]
        model_fit_stats[4,i,j] = ksstats(summary_stat.stat["period_ratio_list"], R_confirmed)[5]
        model_fit_stats[5,i,j] = ksstats(summary_stat.stat["duration list"].*24, t_D_confirmed)[5] #transit durations in simulations are in days, while in the Kepler catalog are in hours
        model_fit_stats[6,i,j] = ksstats(summary_stat.stat["duration_ratio_list"], xi_confirmed)[5]
        model_fit_stats[7,i,j] = ksstats(summary_stat.stat["depth list"], D_confirmed)[5]
        model_fit_stats[8,i,j] = ksstats(summary_stat.stat["radius_ratio_list"].^2, D_ratio_confirmed)[5] #simulations save radius ratios while we computed transit duration ratios from the Kepler catalog

        println("# power_law_r1, power_law_r2: ", get_real(sim_param_here,"power_law_r1"), ", ",  get_real(sim_param_here,"power_law_r2"))
        toc()
    end
end
=#

#=
f = open("Model_stats_r1_r2_rbreak2_grid.out", "w")
println(f,"# Model statistics on a grid of (power_law_r1, power_law_r2) = (", power_law_r1_range, ", ",  power_law_r2_range, ")")
println(f, "# num_targets_sim_pass_one: ", get_int(sim_param,"num_targets_sim_pass_one"))
println(f, "# max_incl_sys: ", get_real(sim_param,"max_incl_sys"))
println(f, "# log_rate_clusters: ", get_real(sim_param,"log_rate_clusters"))
println(f, "# max_clusters_in_sys: ", get_int(sim_param,"max_clusters_in_sys"))
println(f, "# log_rate_planets_per_cluster: ", get_real(sim_param,"log_rate_planets_per_cluster"))
println(f, "# max_planets_in_clusters: ", get_int(sim_param,"max_planets_in_cluster"))
println(f, "# power_law_P: ", get_real(sim_param,"power_law_P"))
println(f, "# min_period: ", get_real(sim_param,"min_period"))
println(f, "# max_period: ", get_real(sim_param,"max_period"))
#####println(f, "# power_law_r1: ", get_real(sim_param,"power_law_r1"))
#####println(f, "# power_law_r2: ", get_real(sim_param,"power_law_r2"))
println(f, "# min_radius (R_earth): ", get_real(sim_param,"min_radius")/ExoplanetsSysSim.earth_radius)
println(f, "# max_radius (R_earth): ", get_real(sim_param,"max_radius")/ExoplanetsSysSim.earth_radius)
println(f, "# break_radius (R_earth): ", get_real(sim_param,"break_radius")/ExoplanetsSysSim.earth_radius)
println(f, "# sigma_incl: ", get_real(sim_param,"sigma_incl"))
println(f, "# sigma_incl_near_mmr: ", get_real(sim_param,"sigma_incl_near_mmr"))
println(f, "# sigma_hk: ", get_real(sim_param,"sigma_hk"))
println(f, "# num_mutual_hill_radii: ", get_real(sim_param,"num_mutual_hill_radii"))
println(f, "# mr_power_index: ", get_real(sim_param,"mr_power_index"))
println(f, "# mr_max_mass (M_earth): ", get_real(sim_param,"mr_max_mass")/ExoplanetsSysSim.earth_mass)
println(f, "# sigma_log_radius_in_cluster: ", get_real(sim_param,"sigma_log_radius_in_cluster"))
println(f, "# sigma_logperiod_per_pl_in_cluster: ", get_real(sim_param,"sigma_logperiod_per_pl_in_cluster"))
for i in 1:size(model_fit_stats)[1]
    println(f,"# ", model_fit_stats_fields[i], ":")
    println(f,model_fit_stats[i,:,:]')
end
close(f)
=#





##### The following optimization code is for exploring our updated model with clustered (lognormal) planet radii, but still a single power law:

#=
#To run our model on a 2D grid of 'power_law_r' and 'sigma_log_radius_in_cluster':

power_law_r_range = linspace(-3., -1.2, 11)
sigma_log_radius_in_cluster_range = linspace(0.05, 0.5, 11)

model_fit_stats = zeros(7,length(power_law_r_range),length(sigma_log_radius_in_cluster_range)) #3D array to be filled with all the statistics of the model fits: number of observed planets, number of systems with observed singles, and KS distances for multiplicities, periods, period ratios, transit depths, and period-normalized transit duration ratios
model_fit_stats_fields = ["Number of observed planets", "Number of systems with observed singles", "KS multiplicities", "KS periods", "KS period ratios", "KS transit depths", "KS period-normalized transit duration ratios"] #names of the fields associated with 'model_fit_stats' array

sim_param_default = setup_sim_param_model() #To set up all the default parameters

for (i,power_law_r) in enumerate(power_law_r_range)
    for (j,sigma_log_radius_in_cluster) in enumerate(sigma_log_radius_in_cluster_range)
        tic()
        sim_param_here = deepcopy(sim_param_default) #To set up all the default parameters
        add_param_fixed(sim_param_here,"power_law_r",power_law_r) #To replace the 'power_law_r' parameter
        add_param_fixed(sim_param_here,"sigma_log_radius_in_cluster",sigma_log_radius_in_cluster) #To replace the 'sigma_log_radius_in_cluster' parameter

        cat_phys = generate_kepler_physical_catalog(sim_param_here)
        cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys,sim_param_here)
        cat_obs = observe_kepler_targets_single_obs(cat_phys_cut,sim_param_here)
        summary_stat = calc_summary_stats_model(cat_obs,sim_param_here)

        M_cat_obs = ones(Int64,0) #array to be filled with the number of transiting planets in each simulated system
        for k in 1:get_int(sim_param_here,"max_tranets_in_sys")
            append!(M_cat_obs, k*ones(Int64, summary_stat.stat["num n-tranet systems"][k]))
        end

        model_fit_stats[1,i,j] = summary_stat.stat["num_tranets"]
        model_fit_stats[2,i,j] = summary_stat.stat["num n-tranet systems"][1]
        model_fit_stats[3,i,j] = ksstats_ints(M_cat_obs, M_confirmed)[5]
        model_fit_stats[4,i,j] = ksstats(summary_stat.stat["P list"], P_confirmed)[5]
        model_fit_stats[5,i,j] = ksstats(summary_stat.stat["period_ratio_list"], R_confirmed)[5]
        model_fit_stats[6,i,j] = ksstats(summary_stat.stat["depth list"], D_confirmed)[5]
        model_fit_stats[7,i,j] = ksstats(summary_stat.stat["duration_ratio_list"], xi_confirmed)[5]

        println("# power_law_r, sigma_log_radius_in_cluster: ", get_real(sim_param_here,"power_law_r"), ", ",  get_real(sim_param_here,"sigma_log_radius_in_cluster"))
        toc()
    end
end
=#

#=
f = open("Model_stats_r_sigma_r_grid.out", "w")
println(f,"# Model statistics on a grid of (power_law_r, sigma_log_radius_in_cluster) = (", power_law_r_range, ", ",  sigma_log_radius_in_cluster_range, ")")
println(f,"# num_targets_sim_pass_one: ", get_int(sim_param_default,"num_targets_sim_pass_one"))
println(f,"# log_rate_clusters: ", get_real(sim_param_default,"log_rate_clusters"))
println(f,"# log_rate_planets_per_cluster: ", get_real(sim_param_default,"log_rate_planets_per_cluster"))
println(f,"# power_law_P: ", get_real(sim_param_default,"power_law_P"))
println(f,"# max_radius: ", get_real(sim_param_default,"max_radius"))
println(f,"# sigma_incl: ", get_real(sim_param_default,"sigma_incl"))
println(f,"# sigma_incl_near_mmr: ", get_real(sim_param_default,"sigma_incl_near_mmr"))
println(f,"# sigma_hk: ", get_real(sim_param_default,"sigma_hk"))
println(f,"# num_mutual_hill_radii: ", get_real(sim_param_default,"num_mutual_hill_radii"))
println(f,"# mr_power_index: ", get_real(sim_param_default,"mr_power_index"))
println(f,"# mr_max_mass: ", get_real(sim_param_default,"mr_max_mass"))
println(f,"# sigma_logperiod_per_pl_in_cluster: ", get_real(sim_param_default,"sigma_logperiod_per_pl_in_cluster"))
for i in 1:size(model_fit_stats)[1]
    println(f,"# ", model_fit_stats_fields[i], ":")
    println(f,model_fit_stats[i,:,:]')
end
close(f)
=#


