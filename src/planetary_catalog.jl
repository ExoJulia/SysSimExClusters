using JLD
using JLD2
using FileIO
using Plots
using Dates

sim_param = setup_sim_param_model()

include("misc_functions.jl")





##### To load the DR25 KOI catalog, update the planet radii using Gaia stellar radii, apply the necessary cuts:

#planet_catalog = CSV.read(joinpath(dirname(pathof(ExoplanetsSysSim)), "../data/q1_q17_dr25_koi.csv"), header=157, allowmissing=:all)
planet_catalog = load(joinpath(dirname(pathof(ExoplanetsSysSim)), "../data/q1_q17_dr25_koi.jld2"))["koi"]

stellar_catalog = ExoplanetsSysSim.StellarTable.setup_star_table(sim_param)

function keep_planet_candidates_given_sim_param(planet_catalog::DataFrame; sim_param::SimParam, stellar_catalog::DataFrame, recompute_radii::Bool=true)

    planets_keep = planet_catalog[(planet_catalog[:koi_disposition] .== "CONFIRMED") .| (planet_catalog[:koi_disposition] .== "CANDIDATE"), :] # table containing only the confirmed and candidate objects
    println("Candidate and confirmed planets: ", size(planets_keep, 1))

    in_stellar_catalog = [] # will be filled with booleans indicating whether each koi in 'planets_keep' is found in the 'stellar_catalog' or not
    for i in 1:length(planets_keep[:kepid])
        if any(x->x==planets_keep[i,:kepid], stellar_catalog[:kepid])
            push!(in_stellar_catalog, true)
        else
            push!(in_stellar_catalog, false)
        end
    end
    in_stellar_catalog_indices = findall(in_stellar_catalog)

    planets_keep = planets_keep[in_stellar_catalog_indices, :]
    println("After removing planets not around stars in stellar catalog: ", size(planets_keep, 1))

    planets_keep = planets_keep[(ismissing.(planets_keep[:koi_duration]) .== false) .& (planets_keep[:koi_duration] .>= 0), :]
    planets_keep = planets_keep[(ismissing.(planets_keep[:koi_depth]) .== false) .& (planets_keep[:koi_depth] .> 0), :]
    println("After removing planets with missing or negative transit durations or depths: ", size(planets_keep, 1))

    if recompute_radii
        for (i,kepid) in enumerate(planets_keep[:kepid]) # to replace the stellar and planetary radii in 'planets_keep' with the more reliable values as derived from the stellar properties in 'stellar_catalog'
            stellar_radii_new = stellar_catalog[stellar_catalog[:kepid] .== kepid, :radius][1]
            stellar_mass = stellar_catalog[stellar_catalog[:kepid] .== kepid, :mass][1]
            planets_keep[i, :koi_srad] = stellar_radii_new
            planets_keep[i, :koi_smass] = stellar_mass # to save stellar masses because for some reason, currently the stellar mass for planets after the first planet in multi-planet systems are set to nan
            planets_keep[i, :koi_prad] = (1 ./ExoplanetsSysSim.earth_radius)*stellar_radii_new*sqrt(planets_keep[i, :koi_depth]/(1e6))
        end
    end

    planets_keep = planets_keep[(planets_keep[:koi_period] .> get_real(sim_param,"min_period")) .& (planets_keep[:koi_period] .< get_real(sim_param,"max_period")), :] # to make additional cuts in period P to be comparable to our simulated sample
    planets_keep = planets_keep[(planets_keep[:koi_prad] .> get_real(sim_param,"min_radius")/ExoplanetsSysSim.earth_radius) .& (planets_keep[:koi_prad] .< get_real(sim_param,"max_radius")/ExoplanetsSysSim.earth_radius) .& (.~ismissing.(planets_keep[:koi_prad])), :] # to make additional cuts in planetary radii to be comparable to our simulated sample
    println("After applying our period and radius cuts (final count): ", size(planets_keep, 1))
    return planets_keep
end

@time planets_cleaned = keep_planet_candidates_given_sim_param(planet_catalog; sim_param=sim_param, stellar_catalog=stellar_catalog, recompute_radii=true)

# If we want to write the cleaned planetary catalog and the stellar catalog to a csv file, keeping only the columns we need:
#CSV.write("q1_q17_dr25_gaia_fgk_koi_cleaned.csv", planets_cleaned[[:kepid, :kepoi_name, :koi_disposition, :koi_pdisposition, :koi_score, :koi_period, :koi_duration, :koi_depth, :koi_prad, :koi_steff, :koi_slogg, :koi_srad, :koi_smass]])
#CSV.write("q1_q17_dr25_gaia_fgk_cleaned.csv", stellar_catalog[[:kepid, :mass, :radius, :teff]])





##### To compute arrays of the observables (multiplicities, periods, period ratios, transit durations, transit depths, period-normalized transit duration ratios (xi), and transit depth ratios) from the remaining sample of planets:

KOI_systems = [x[1:6] for x in planets_cleaned[:kepoi_name]]
checked_bools = zeros(size(planets_cleaned,1)) #0's denote KOI that were not checked yet; 1's denote already checked KOI

M_confirmed = Int64[] #list to be filled with the planet multiplicities of the systems
R_confirmed = Float64[] #list to be filled with period ratios of adjacent planet pairs
xi_confirmed = Float64[] #list to be filled with the period-normalized transit duration ratios of adjacent planet pairs
xi_non_mmr_confirmed = Float64[] #list to be filled with the period-normalized transit duration ratios of adjacent planet pairs not near any resonances
xi_near_mmr_confirmed = Float64[] #list to be filled with the period-normalized transit duration ratios of adjacent planet pairs near a resonance
D_ratio_confirmed = Float64[] #list to be filled with the transit depth ratios of adjacent planet pairs
P_confirmed = planets_cleaned[:koi_period] #array of the periods (days)
P_confirmed = collect(skipmissing(P_confirmed))
t_D_confirmed = planets_cleaned[:koi_duration] #array of the transit durations (hrs)
t_D_confirmed = collect(skipmissing(t_D_confirmed))
D_confirmed = planets_cleaned[:koi_depth]/(1e6) #array of the transit depths (fraction)
D_confirmed = D_confirmed[ismissing.(D_confirmed) .== false] #to get rid of NA values

D_above_confirmed = Float64[] #list to be filled with the transit depths of planets above the photoevaporation boundary in Carrera et al 2018
D_below_confirmed = Float64[] #list to be filled with the transit depths of planets below the boundary
D_ratio_above_confirmed = Float64[] #list to be filled with the transit depth ratios of adjacent planet pairs, both above the boundary
D_ratio_below_confirmed = Float64[] #list to be filled with the transit depth ratios of adjacent planet pairs, both below the boundary
D_ratio_across_confirmed = Float64[] #list to be filled with the transit depth ratios of adjacent planet pairs, across the boundary

for i in 1:length(KOI_systems)
    if checked_bools[i] == 0 #if the KOI has not been checked (included while looking at another planet in the same system)
        system_i = (1:length(KOI_systems))[KOI_systems .== KOI_systems[i]]
        checked_bools[system_i] .= 1

        #To get the periods and transit durations in this system:
        system_P = planets_cleaned[:koi_period][system_i] #periods of all the planets in this system
        system_t_D = planets_cleaned[:koi_duration][system_i] #transit durations of all the planets in this system
        system_D = planets_cleaned[:koi_depth][system_i] #transit depths (in ppm) of all the planets in this system
        system_radii = planets_cleaned[:koi_prad][system_i] #radii of all the planets in this system
        system_sort_i = sortperm(system_P) #indices that would sort the periods of the planets in this system
        system_P = system_P[system_sort_i] #periods of all the planets in this system, sorted
        system_t_D = system_t_D[system_sort_i] #transit durations of all the planets in this system, sorted by period
        system_D = system_D[system_sort_i]/(1e6) #transit depths of all the planets in this system, sorted by period
        system_radii = system_radii[system_sort_i] #radii of all the planets in this system, sorted by period

        #To count the total number of planets in this system:
        push!(M_confirmed, length(system_P))

        #To compute the period ratios, period-normalized transit duration ratios, and transit depth ratios in this system:
        system_R = system_P[2:end] ./ system_P[1:end-1] #period ratios of all the adjacent planet pairs in this system
        system_D_ratio = system_D[2:end] ./ system_D[1:end-1] #transit depth ratios of all the adjacent planet pairs in this system
        system_xi = (system_t_D[1:end-1] ./ system_t_D[2:end]) .* (system_R .^(1//3)) #period-normalized transit duration ratios of all the adjacent planet pairs in this system

        append!(R_confirmed, system_R)
        append!(D_ratio_confirmed, system_D_ratio)
        append!(xi_confirmed, system_xi)
        for (j,period_ratio) in enumerate(system_R)
            if is_period_ratio_near_resonance(period_ratio, sim_param)
                append!(xi_near_mmr_confirmed, system_xi[j])
            else
                append!(xi_non_mmr_confirmed, system_xi[j])
            end
        end

        #To separate the planets in the system as above and below the boundary:
        system_above_bools = [photoevap_boundary_Carrera2018(system_radii[x], system_P[x]) for x in 1:length(system_P)]
        #if length(system_above_bools) > 1 println(system_above_bools) end

        #To record the transit depths of the planets above and below the boundary:
        for (j,D) in enumerate(system_D)
            if system_above_bools[j] == 1
                append!(D_above_confirmed, D)
            elseif system_above_bools[j] == 0
                append!(D_below_confirmed, D)
            end
        end

        #To record the transit depth ratios of the planets above, below, and across the boundary:
        for (j,D_ratio) in enumerate(system_D_ratio)
            if system_above_bools[j] + system_above_bools[j+1] == 2 #both planets are above the boundary
                append!(D_ratio_above_confirmed, D_ratio)
            elseif system_above_bools[j] + system_above_bools[j+1] == 1 #one planet is above, the other is below the boundary
                append!(D_ratio_across_confirmed, D_ratio)
            elseif system_above_bools[j] + system_above_bools[j+1] == 0 #both planets are below the boundary
                append!(D_ratio_below_confirmed, D_ratio)
            end
        end
    end
end

Nmult_confirmed = [sum(M_confirmed .== k) for k in 1:maximum(M_confirmed)]
