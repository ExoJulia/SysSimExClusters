include("misc_functions.jl")
include("complexity_stats.jl")



"""
    calc_summary_stats_idx_n_tranets!(css, cat_obs, sim_param)

Compile a list of indices for systems with each observed multiplicity in an observed catalog `cat_obs` and add it to the summary statistics (`css.cache`).

# Arguments:
- `css::CatalogSummaryStatistics`: object containing the summary statistics.
- `cat_obs::KeplerObsCatalog`: an observed catalog.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.

# Returns:
- `idx_n_tranets::Vector{Vector{Int64}}`: list of indices for the systems with 1,2,...,8 observed planets in `cat_obs`. Also writes to `css.cache`.
"""
function calc_summary_stats_idx_n_tranets!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, sim_param::SimParam)
    if haskey(css.cache,"idx_n_tranets")
        return css.cache["idx_n_tranets"]
    end

    max_tranets_in_sys = get_int(sim_param, "max_tranets_in_sys")
    idx_n_tranets = Vector{Int64}[Int64[] for m in 1:max_tranets_in_sys]
    for n in 1:max_tranets_in_sys-1
        idx_n_tranets[n] = findall(x::KeplerTargetObs -> length(x.obs)==n, cat_obs.target)
    end
    idx_n_tranets[max_tranets_in_sys] = findall(x::KeplerTargetObs -> length(x.obs)>=max_tranets_in_sys, cat_obs.target)
    css.cache["idx_n_tranets"] = idx_n_tranets
    return idx_n_tranets
end



# Count total number of tranets using lists of indices for N-tranet systems
"""
    calc_summary_stats_num_tranets!(css, cat_obs, sim_param)

Compute the total number of observed transiting planets in `cat_obs` and add it to the summary statistics (`css.stat`).

# Arguments:
- `css::CatalogSummaryStatistics`: object containing the summary statistics.
- `cat_obs::KeplerObsCatalog`: an observed catalog.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.

# Returns:
- `num_tranets::Int64`: total number of observed transiting planets. Also writes to `css.stat`.
"""
function calc_summary_stats_num_tranets!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, sim_param::SimParam)
    if haskey(css.stat,"num_tranets")
        return css.stat["num_tranets"]
    elseif haskey(css.cache,"num_tranets")
        css.stat["num_tranets"] = css.cache["num_tranets"]
        return css.stat["num_tranets"]
    end

    idx_n_tranets = calc_summary_stats_idx_n_tranets!(css, cat_obs, sim_param)
    num_tranets = 0
    for n in 1:length(idx_n_tranets)
        num_tranets += n*length(idx_n_tranets[n])
    end
    css.stat["num_tranets"] = num_tranets
    return num_tranets
end



"""
    calc_summary_stats_num_targets!(css, cat_obs, sim_param)

Read the total number of observed/simulated targets from `sim_param` and add it to the summary statistics (`css.stat`).

# Arguments:
- `css::CatalogSummaryStatistics`: object containing the summary statistics.
- `cat_obs::KeplerObsCatalog`: an observed catalog.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.

# Returns:
The total number of targets in the observed catalog.
"""
function calc_summary_stats_num_targets!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, sim_param::SimParam)
    css.stat["num_targets"] = length(cat_obs.target)
end



"""
    calc_summary_stats_num_n_tranet_systems!(css, cat_obs, sim_param)

Compute the observed number of systems with each planet multiplicity and add it to the summary statistics (`css.stat`).

# Arguments:
- `css::CatalogSummaryStatistics`: object containing the summary statistics.
- `cat_obs::KeplerObsCatalog`: an observed catalog.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.

# Returns:
- `num_n_tranet_systems::Vector{Int64}`: list of the number of systems with each observed planet multiplicity. Also writes to `css.stat`.
"""
function calc_summary_stats_num_n_tranet_systems!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, sim_param::SimParam)
    idx_n_tranets = calc_summary_stats_idx_n_tranets!(css, cat_obs, sim_param)
    #max_tranets_in_sys = get_int(sim_param, "max_tranets_in_sys")
    num_n_tranet_systems = map(n -> length(idx_n_tranets[n]), 1:length(idx_n_tranets))
    css.stat["num_n-tranet_systems"] = num_n_tranet_systems
end



"""
    calc_summary_stats_duration_ratios_neighbors!(css, cat_obs, sim_param)

Compute the (period-normalized) transit duration ratios for adjacent planet pairs in the observed catalog and add it to the summary statistics (`css.stat`).

# Arguments:
- `css::CatalogSummaryStatistics`: object containing the summary statistics.
- `cat_obs::KeplerObsCatalog`: an observed catalog.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.

# Returns:
- `duration_ratio_list::Vector{Float64}`: list of transit duration ratios.
- `duration_ratio_non_mmr_list::Vector{Float64}`: list of transit duration ratios for planet pairs not near a (first order) mean-motion resonance (MMR).
- `duration_ratio_near_mmr_list::Vector{Float64}`: list of transit duration ratios for planet pairs near a (first order) MMR.
Also writes these to `css.stat`.

Note: see `sim_param` for what is considered "near" and "not near" a first order MMR, and which MMRs are included in the model.
"""
function calc_summary_stats_duration_ratios_neighbors!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, sim_param::SimParam)
    if haskey(css.stat,"duration_ratios")
        return css.stat["duration_ratios"]
    elseif haskey(css.cache,"duration_ratios")
        return css.cache["duration_ratios"]
    end
    idx_n_tranets = calc_summary_stats_idx_n_tranets!(css, cat_obs, sim_param)
    @assert length(idx_n_tranets) >= 1

    # Calculate how many duration ratios there will be & allocate storage:
    num_ratios = 0
    for i in 2:length(idx_n_tranets)
        num_ratios += length(idx_n_tranets[i])*(i-1)
    end
    duration_ratio_list = Array{Float64}(undef, num_ratios)
    duration_ratio_non_mmr_list = Float64[]
    duration_ratio_near_mmr_list = Float64[]

    k = 0
    for n in 2:length(idx_n_tranets) # Loop over number of tranets in system
        for i in idx_n_tranets[n] # Loop over systems with n tranets
            period_in_sys = Array{Float64}(undef, n)
            duration_in_sys = Array{Float64}(undef, n)
            for j in 1:n # Loop over periods within a system
                period_in_sys[j] = cat_obs.target[i].obs[j].period
                duration_in_sys[j] = cat_obs.target[i].obs[j].duration
            end
            perm = sortperm(period_in_sys)
            for j in 1:(n-1) # Loop over period ratios within a system
                period_ratio = period_in_sys[perm[j+1]]/period_in_sys[perm[j]]
                if 1<period_ratio<Inf
                    k += 1
                    xi = duration_in_sys[perm[j]]/duration_in_sys[perm[j+1]] * period_ratio^(1//3)
                    duration_ratio_list[k] = xi

                    if is_period_ratio_near_resonance(period_ratio, sim_param)
                        append!(duration_ratio_near_mmr_list, xi)
                    else
                        append!(duration_ratio_non_mmr_list, xi)
                    end
                end
            end
        end
    end
    resize!(duration_ratio_list, k)
    css.stat["duration_ratios"] = duration_ratio_list
    css.stat["duration_ratios_nonmmr"] = duration_ratio_non_mmr_list
    css.stat["duration_ratios_mmr"] = duration_ratio_near_mmr_list

    return (duration_ratio_list, duration_ratio_non_mmr_list, duration_ratio_near_mmr_list)
end



"""
    calc_summary_stats_period_radius_ratios_neighbors_internal!(css, cat_obs, sim_param)

Compute the period ratios and radius ratios for adjacent planet pairs in the observed catalog, and add them to the cached summary statistics (`css.cache`).

# Arguments:
- `css::CatalogSummaryStatistics`: object containing the summary statistics.
- `cat_obs::KeplerObsCatalog`: an observed catalog.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.

# Returns:
- `period_ratio_list::Vector{Float64}`: list of period ratios.
- `radius_ratio_list::Vector{Float64}`: list of radius ratios.
- `radius_ratio_above_list::Vector{Float64}`: list of radius ratios where both planets are observed to be above the photo-evaporation boundary defined by Carrera et al. (2018).
- `radius_ratio_below_list::Vector{Float64}`: list of radius ratios where both planets are observed to be below the photo-evaporation boundary.
- `radius_ratio_across_list::Vector{Float64}`: list of radius ratios where the planets straddle the photo-evaporation boundary.
Also writes these to `css.cache`.
"""
function calc_summary_stats_period_radius_ratios_neighbors_internal!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, sim_param::SimParam)
    #=
    if haskey(css.stat,"period_ratios") && haskey(css.stat,"radius_ratios")
        return (css.stat["period_ratios"], css.stat["radius_ratios"])
    elseif haskey(css.cache,"period_ratios") && haskey(css.cache,"radius_ratios")
        return (css.cache["period_ratios"], css.cache["radius_ratios"])
    end
    =#
    idx_n_tranets = calc_summary_stats_idx_n_tranets!(css, cat_obs, sim_param)
    @assert length(idx_n_tranets) >= 1

    # Calculate how many period ratios there will be & allocate storage
    num_ratios = 0
    for i in 2:length(idx_n_tranets)
        num_ratios += length(idx_n_tranets[i])*(i-1)
    end
    period_ratio_list = Array{Float64}(undef, num_ratios)
    radius_ratio_list = Array{Float64}(undef, num_ratios)

    radius_ratio_above_list = Float64[] # list to be filled with the radius ratios of adjacent planet pairs, both above the photoevaporation boundary
    radius_ratio_below_list = Float64[] # list to be filled with the radius ratios of adjacent planet pairs, both below the boundary
    radius_ratio_across_list = Float64[] # list to be filled with the radius ratios of adjacent planet pairs, across the boundary

    k = 0
    for n in 2:length(idx_n_tranets) # Loop over number of tranets in system
        period_in_sys = Array{Float64}(undef, n)
        radius_in_sys = Array{Float64}(undef, n)
        depth_in_sys = Array{Float64}(undef, n)
        for i in idx_n_tranets[n] # Loop over systems with n tranets
            for j in 1:n # Loop over periods within a system
                period_in_sys[j] = cat_obs.target[i].obs[j].period
                radius_in_sys[j] = sqrt(cat_obs.target[i].obs[j].depth)*cat_obs.target[i].star.radius
                depth_in_sys[j] = cat_obs.target[i].obs[j].depth
            end
            perm = sortperm(period_in_sys)
            for j in 1:(n-1) # Loop over period ratios within a system
                k = k+1
                radius_in_earths, radius_out_earths = radius_in_sys[perm[j]]/ExoplanetsSysSim.earth_radius, radius_in_sys[perm[j+1]]/ExoplanetsSysSim.earth_radius
                period_in, period_out = period_in_sys[perm[j]], period_in_sys[perm[j+1]]

                period_ratio_list[k] = period_out/period_in
                #radius_ratio_list[k] = radius_out_earths/radius_in_earths
                radius_ratio_list[k] = sqrt(depth_in_sys[perm[j+1]]/depth_in_sys[perm[j]])

                if photoevap_boundary_Carrera2018(radius_in_earths, period_in) + photoevap_boundary_Carrera2018(radius_out_earths, period_out) == 2
                    append!(radius_ratio_above_list, sqrt(depth_in_sys[perm[j+1]]/depth_in_sys[perm[j]]))
                elseif photoevap_boundary_Carrera2018(radius_in_earths, period_in) + photoevap_boundary_Carrera2018(radius_out_earths, period_out) == 1
                    append!(radius_ratio_across_list, sqrt(depth_in_sys[perm[j+1]]/depth_in_sys[perm[j]]))
                elseif photoevap_boundary_Carrera2018(radius_in_earths, period_in) + photoevap_boundary_Carrera2018(radius_out_earths, period_out) == 0
                    append!(radius_ratio_below_list, sqrt(depth_in_sys[perm[j+1]]/depth_in_sys[perm[j]]))
                end
            end
        end
    end
    css.cache["period_ratios"] = period_ratio_list
    css.cache["radius_ratios"] = radius_ratio_list

    css.cache["radius_ratios_above"] = radius_ratio_above_list
    css.cache["radius_ratios_below"] = radius_ratio_below_list
    css.cache["radius_ratios_across"] = radius_ratio_across_list

    return (period_ratio_list, radius_ratio_list, radius_ratio_above_list, radius_ratio_below_list, radius_ratio_across_list)
end



"""
    calc_summary_stats_period_radius_ratios_neighbors!(css, cat_obs, sim_param)

Compute the period ratios and radius ratios for adjacent planet pairs in the observed catalog by calling `calc_summary_stats_period_radius_ratios_neighbors_internal!`, and add them to the summary statistics (`css.stat`).

# Arguments:
- `css::CatalogSummaryStatistics`: object containing the summary statistics.
- `cat_obs::KeplerObsCatalog`: an observed catalog.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.

# Returns:
- `period_ratio_list::Vector{Float64}`: list of period ratios.
- `radius_ratio_list::Vector{Float64}`: list of radius ratios.
Also writes these to `css.stat`.
"""
function calc_summary_stats_period_radius_ratios_neighbors!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, sim_param::SimParam)
    period_ratio_list, radius_ratio_list = calc_summary_stats_period_radius_ratios_neighbors_internal!(css, cat_obs, sim_param)[1:2]
    css.stat["period_ratios"] = period_ratio_list
    css.stat["radius_ratios"] = radius_ratio_list
    return (period_ratio_list, radius_ratio_list)
end



"""
    calc_summary_stats_period_ratios_neighbors!(css, cat_obs, sim_param)

Compute the period ratios for adjacent planet pairs in the observed catalog by calling `calc_summary_stats_period_radius_ratios_neighbors_internal!`, and add it to the summary statistics (`css.stat`).

# Arguments:
- `css::CatalogSummaryStatistics`: object containing the summary statistics.
- `cat_obs::KeplerObsCatalog`: an observed catalog.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.

# Returns:
- `period_ratio_list::Vector{Float64}`: list of period ratios. Also writes to `css.stat`.
"""
function calc_summary_stats_period_ratios_neighbors!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, sim_param::SimParam)
    period_ratio_list = calc_summary_stats_period_radius_ratios_neighbors_internal!(css, cat_obs, sim_param)[1]
    css.stat["period_ratios"] = period_ratio_list
    return period_ratio_list
end



"""
    calc_summary_stats_radius_ratios_neighbors!(css, cat_obs, sim_param)

Compute the radius ratios for adjacent planet pairs in the observed catalog by calling `calc_summary_stats_period_radius_ratios_neighbors_internal!`, and add it to the summary statistics (`css.stat`).

# Arguments:
- `css::CatalogSummaryStatistics`: object containing the summary statistics.
- `cat_obs::KeplerObsCatalog`: an observed catalog.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.

# Returns:
- `radius_ratio_list::Vector{Float64}`: list of radius ratios. Also writes to `css.stat`.
"""
function calc_summary_stats_radius_ratios_neighbors!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, sim_param::SimParam)
    radius_ratio_list = calc_summary_stats_period_radius_ratios_neighbors_internal!(css, cat_obs, sim_param)[2]
    css.stat["radius_ratios"] = radius_ratio_list
    return radius_ratio_list
end



"""
    calc_summary_stats_radius_ratios_neighbors_photoevap_boundary_Carrera2018!(css, cat_obs, sim_param)

Compute the radius ratios for adjacent planet pairs in the observed catalog, split into three catagories based on whether the planets are both above, both below, or straddling across the photo-evaporation boundary defined by Carrera et al. (2018), by calling `calc_summary_stats_period_radius_ratios_neighbors_internal!`, and add them to the cached summary statistics (`css.cache`).

# Arguments:
- `css::CatalogSummaryStatistics`: object containing the summary statistics.
- `cat_obs::KeplerObsCatalog`: an observed catalog.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.

# Returns:
- `radius_ratio_above_list::Vector{Float64}`: list of radius ratios where both planets are observed to be above the photo-evaporation boundary defined by Carrera et al. (2018).
- `radius_ratio_below_list::Vector{Float64}`: list of radius ratios where both planets are observed to be below the photo-evaporation boundary.
- `radius_ratio_across_list::Vector{Float64}`: list of radius ratios where the planets straddle the photo-evaporation boundary.
Also writes these to `css.stat`.
"""
function calc_summary_stats_radius_ratios_neighbors_photoevap_boundary_Carrera2018!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, sim_param::SimParam)
    radius_ratio_above_list, radius_ratio_below_list, radius_ratio_across_list = calc_summary_stats_period_radius_ratios_neighbors_internal!(css, cat_obs, sim_param)[3:5]
    css.stat["radius_ratios_above"] = radius_ratio_above_list
    css.stat["radius_ratios_below"] = radius_ratio_below_list
    css.stat["radius_ratios_across"] = radius_ratio_across_list
    return (radius_ratio_above_list, radius_ratio_below_list, radius_ratio_across_list)
end



# Untested and unused function:
function calc_summary_stats_mean_std_log_period_depth!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, sim_param::SimParam)
    # Allocate arrays to store values for each tranet
    num_tranets  = calc_summary_stats_num_tranets!(css, cat_obs, sim_param)
    period_list = zeros(num_tranets)
    depth_list = zeros(num_tranets)
    #weight_list = ones(num_tranets)

    idx_n_tranets = calc_summary_stats_idx_n_tranets!(css, cat_obs, sim_param)
    max_tranets_in_sys = get_int(sim_param, "max_tranets_in_sys")
    @assert max_tranets_in_sys >= 1
    i = 1 # tranet id
    for targ in cat_obs.target # For each target
        for j in 1:min(length(targ.obs),max_tranets_in_sys) # For each tranet around that target (but truncated if too many tranets in one system)
            #println("# i= ",i," j= ",j)
            period_list[i] = targ.obs[j].period
            depth_list[i] = targ.obs[j].depth
            #weight_list[i] = 1.0
            i = i+1
        end
    end

    css.cache["periods"] = period_list
    css.cache["depths"] = depth_list
    #css.cache["weight list"] = weight_list

    idx_good = Bool[ period_list[i]>0.0 && depth_list[i]>0.0 for i in 1:length(period_list) ]
    log_period_list = log10(period_list[idx_good])
    log_depth_list = log10(depth_list[idx_good])
    css.stat["mean log10 P"] = mean_log_P = mean(log_period_list)
    css.stat["mean log10 depth"] = mean_log_depth = mean(log_depth_list)
    css.stat["std log10 P"] = std_log_P = stdm(log_period_list,mean_log_P)
    css.stat["std log10 depth"] = std_log_depth = stdm(log_depth_list,mean_log_depth)

    return (mean_log_P, std_log_P, mean_log_depth, std_log_depth)
end



"""
    calc_summary_stats_cuml_period_depth_duration!(css, cat_obs, sim_param)

Compute the lists of periods, transit depths, transit durations, and normalized transit durations in the observed catalog and add them to the summary statistics (`css.stat`).

# Arguments:
- `css::CatalogSummaryStatistics`: object containing the summary statistics.
- `cat_obs::KeplerObsCatalog`: an observed catalog.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.

# Returns:
- `period_list::Vector{Float64}`: list of observed periods.
- `depth_list::Vector{Float64}`: list of transit depths.
- `duration_list::Vector{Float64}`: list of transit durations.
- `duration_norm_circ_list::Vector{Float64}`: list of transit durations normalized by the circular, central durations.
- `depth_above_list::Vector{Float64}`: list of transit depths above the photo-evaporation boundary defined by Carrera et al. (2018).
- `depth_below_list::Vector{Float64}`: list of transit depths below the photo-evaporation boundary.
Also writes the period, depth, duration, and normalized duration lists to `css.stat`, and writes the depth lists above and below the photo-evaporation boundary to `css.cache`.
"""
function calc_summary_stats_cuml_period_depth_duration!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, sim_param::SimParam)
    # Allocate arrays to store values for each tranet:
    num_tranets = calc_summary_stats_num_tranets!(css, cat_obs, sim_param)
    period_list = zeros(num_tranets)
    depth_list = zeros(num_tranets)
    duration_list = zeros(num_tranets)
    duration_norm_circ_list = zeros(num_tranets)
    duration_norm_circ_singles_list = Float64[]
    duration_norm_circ_multis_list = Float64[]
    #weight_list = ones(num_tranets)

    depth_above_list = Float64[] # list to be filled with the transit depths of planets above the photoevaporation boundary in Carrera et al 2018
    depth_below_list = Float64[] # list to be filled with the transit depths of planets below the boundary

    idx_n_tranets = calc_summary_stats_idx_n_tranets!(css, cat_obs, sim_param)
    max_tranets_in_sys = get_int(sim_param, "max_tranets_in_sys")
    @assert max_tranets_in_sys >= 1
    i = 0 # tranet id
    for targ in cat_obs.target # For each target
        for j in 1:min(length(targ.obs), max_tranets_in_sys) # For each tranet around that target (but truncated if too many tranets in one system)
            i = i+1
            #println("# i= ",i," j= ",j)
            period_list[i] = targ.obs[j].period
            depth_list[i] = targ.obs[j].depth
            duration_list[i] = targ.obs[j].duration
            duration_norm_circ = targ.obs[j].duration/calc_transit_duration_central_circ_obs(targ, j)
            duration_norm_circ_list[i] = duration_norm_circ
            if length(targ.obs) == 1
                append!(duration_norm_circ_singles_list, duration_norm_circ)
            elseif length(targ.obs) > 1
                append!(duration_norm_circ_multis_list, duration_norm_circ)
            end
            #weight_list[i] = 1.0

            radius_earths, period = (sqrt(targ.obs[j].depth)*targ.star.radius) / ExoplanetsSysSim.earth_radius, targ.obs[j].period
            if photoevap_boundary_Carrera2018(radius_earths, period) == 1
                append!(depth_above_list, targ.obs[j].depth)
            elseif photoevap_boundary_Carrera2018(radius_earths, period) == 0
                append!(depth_below_list, targ.obs[j].depth)
            end
        end
    end
    resize!(period_list, i)
    resize!(depth_list, i)
    resize!(duration_list, i)
    resize!(duration_norm_circ_list, i)
    css.stat["periods"] = period_list
    css.stat["depths"] = depth_list
    css.stat["durations"] = duration_list
    css.stat["durations_norm_circ"] = duration_norm_circ_list
    css.stat["durations_norm_circ_singles"] = duration_norm_circ_singles_list
    css.stat["durations_norm_circ_multis"] = duration_norm_circ_multis_list
    #css.cache["weight list"] = weight_list

    css.cache["depths_above"] = depth_above_list
    css.cache["depths_below"] = depth_below_list

    #println("# P = ",period_list)
    return (period_list, depth_list, duration_list, duration_norm_circ_list, depth_above_list, depth_below_list)
end



"""
    calc_summary_stats_depths_photoevap_boundary_Carrera2018!(css, cat_obs, sim_param)

Compute the transit depths lists in the observed catalog, split by planets above and below the photo-evaporation boundary defined by Carrera et al. (2018), and add them to the summary statistics (`css.stat`).

# Arguments:
- `css::CatalogSummaryStatistics`: object containing the summary statistics.
- `cat_obs::KeplerObsCatalog`: an observed catalog.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.

# Returns:
- `depth_above_list::Vector{Float64}`: list of transit depths above the photo-evaporation boundary defined by Carrera et al. (2018).
- `depth_below_list::Vector{Float64}`: list of transit depths below the photo-evaporation boundary.
Also writes these to `css.stat`.
"""
function calc_summary_stats_depths_photoevap_boundary_Carrera2018!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, sim_param::SimParam)
    depth_above_list, depth_below_list = calc_summary_stats_cuml_period_depth_duration!(css, cat_obs, sim_param)[5:6]
    css.stat["depths_above"] = depth_above_list
    css.stat["depths_below"] = depth_below_list
    return (depth_above_list, depth_below_list)
end



# Untested and unused function:
function calc_summary_stats_obs_binned_rates!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, sim_param::SimParam)
    num_tranets  = calc_summary_stats_num_tranets!(css, cat_obs, sim_param)
    idx_n_tranets = calc_summary_stats_idx_n_tranets!(css, cat_obs, sim_param)

    idx_tranets = findall(x::KeplerTargetObs-> length(x.obs) > 0, cat_obs.target)::Array{Int64,1} # Find indices of systems with at least 1 tranet = potentially detectable transiting planet
    css.cache["idx_tranets"] = idx_tranets # We can save lists of indices to summary stats for pass 2, even though we will not use these for computing a distance or probability

    #=
    max_tranets_in_sys = get_int(sim_param, "max_tranets_in_sys") # Demo that simulation parameters can specify how to evaluate models, too
    @assert max_tranets_in_sys >= 1

    # Count total number of tranets and compile indices for N-tranet systems
    num_tranets = 0
    idx_n_tranets = Vector{Int64}[ Int64[] for m = 1:max_tranets_in_sys]
    for n in 1:max_tranets_in_sys-1
    idx_n_tranets[n] = findall(x::KeplerTargetObs-> length(x.obs) == n, cat_obs.target )
    num_tranets += n*length(idx_n_tranets[n])
    end
    idx_n_tranets[max_tranets_in_sys] = findall(x::KeplerTargetObs-> length(x.obs) >= max_tranets_in_sys, cat_obs.target )
    css.cache["idx_n_tranets"] = idx_n_tranets

    num_tranets += max_tranets_in_sys*length(idx_n_tranets[max_tranets_in_sys])  # WARNING: this means we need to ignore planets w/ indices > max_tranets_in_sys
    if ( length( findall(x::KeplerTargetObs-> length(x.obs) > max_tranets_in_sys, cat_obs.target ) ) > 0) # Make sure max_tranets_in_sys is at least big enough for observed systems
    warn("Observational data has more transiting planets in one systems than max_tranets_in_sys allows.")
    end
    num_tranets  = convert(Int64,num_tranets) # TODO OPT: Figure out why is not this already an Int. I may be doing something that prevents some optimizations
    css.cache["num_tranets"] = num_tranets

    num_sys_tranets = zeros(max_tranets_in_sys) # Since observed data, don't need to calculate probabilities.
    for n in 1:max_tranets_in_sys # Make histogram of N-tranet systems
    num_sys_tranets[n] = length(idx_n_tranets[n])
    end
    css.stat["num_sys_tranets"] = num_sys_tranets
    css.stat["planets detected"] = num_tranets
    =#

    period_list = zeros(num_tranets)
    radius_list = zeros(num_tranets)
    weight_list = ones(num_tranets)

    n = 1 # tranet id
    for i in 1:length(cat_obs.target)
    for j in 1:num_planets(cat_obs.target[i])
        period_list[n] = cat_obs.target[i].obs[j].period
        radius_list[n] = sqrt(cat_obs.target[i].obs[j].depth)*cat_obs.target[i].star.radius
        n = n+1
    end
    end

    limitP::Array{Float64,1} = get_any(sim_param, "p_lim_arr", Array{Float64,1})
    limitRp::Array{Float64,1} = get_any(sim_param, "r_lim_arr", Array{Float64,1})
    @assert length(limitP)>=2 && length(limitRp)>=2

    np_bin = zeros((length(limitP)-1) * (length(limitRp)-1))
    np_bin_idx = 1
    for i in 1:(length(limitP)-1)
        P_match = findall(x -> ((x > limitP[i]) && (x < limitP[i+1])), period_list)
        for j in 1:(length(limitRp)-1)
            R_match = findall(x -> ((x > limitRp[j]) && (x < limitRp[j+1])), radius_list)

            bin_match = intersect(P_match, R_match)

            np_bin[np_bin_idx] = sum(weight_list[bin_match])
            np_bin_idx += 1
        end
    end

    #css.stat["planets detected"] = sum(np_bin)
    css.stat["planets table"] = np_bin

    return css
end



"""
    calc_summary_stats_complexity_stats_GF2020!(css, cat_obs, sim_param)

Compute the lists of system level metrics taken from or inspired by Gilbert & Fabrycky (2020) in the observed catalog and add them to the summary statistics (`css.stat`).

# Arguments:
- `css::CatalogSummaryStatistics`: object containing the summary statistics.
- `cat_obs::KeplerObsCatalog`: an observed catalog.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.

# Returns:
- `radii_partitioning_list::Vector{Float64}`: list of "radius partitioning" (disequilibirum normalized to (0,1)), for systems with 2+ observed planets.
- `radii_monotonicity_list::Vector{Float64}`: list of "radius monotonicity", for systems with 2+ observed planets.
- `gap_complexity_list::Vector{Float64}`: list of "gap complexity", for systems with 3+ observed planets.
Also writes these lists to `css.stat`.
"""
function calc_summary_stats_complexity_stats_GF2020!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, sim_param::SimParam)
    # Allocate arrays to store values for each tranet:
    num_n_tranet_systems = calc_summary_stats_num_n_tranet_systems!(css, cat_obs, sim_param)
    radii_partitioning_list = zeros(sum(num_n_tranet_systems[2:end]))
    radii_monotonicity_list = zeros(sum(num_n_tranet_systems[2:end]))
    gap_complexity_list = zeros(sum(num_n_tranet_systems[3:end]))

    max_tranets_in_sys = get_int(sim_param, "max_tranets_in_sys")
    @assert max_tranets_in_sys >= 1
    i2 = 0
    i3 = 0
    for targ in cat_obs.target # For each target
        radii_earths = map(j -> (sqrt(targ.obs[j].depth)*targ.star.radius) / ExoplanetsSysSim.earth_radius, 1:length(targ.obs))
        periods = map(j -> targ.obs[j].period, 1:length(targ.obs))

        if length(targ.obs) >= 2
            i2 += 1
            radii_partitioning_list[i2] = partitioning(radii_earths)
            radii_monotonicity_list[i2] = monotonicity_GF2020(radii_earths)
            if length(targ.obs) >= 3
                i3 += 1
                gap_complexity_list[i3] = gap_complexity_GF2020(periods)
            end
        end
    end
    css.stat["radii_partitioning"] = radii_partitioning_list
    css.stat["radii_monotonicity"] = radii_monotonicity_list
    css.stat["gap_complexity"] = gap_complexity_list

    return (radii_partitioning_list, radii_monotonicity_list, gap_complexity_list)
end



"""
    calc_summary_stats_model(cat_obs, sim_param)

Compute all the summary statistics and compile them into a `CatalogSummaryStatistics` object.

# Arguments:
- `cat_obs::KeplerObsCatalog`: an observed catalog.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.

# Returns:
- `css::CatalogSummaryStatistics`: object containing the summary statistics.
"""
function calc_summary_stats_model(cat_obs::KeplerObsCatalog, sim_param::SimParam)
    css = CatalogSummaryStatistics()

    calc_summary_stats_num_targets!(css, cat_obs, sim_param)
    calc_summary_stats_num_tranets!(css, cat_obs, sim_param)
    calc_summary_stats_num_n_tranet_systems!(css, cat_obs, sim_param)
    calc_summary_stats_cuml_period_depth_duration!(css, cat_obs, sim_param)
    calc_summary_stats_depths_photoevap_boundary_Carrera2018!(css, cat_obs, sim_param) # to also compute arrays of depths above and below photoevaporation boundary
    calc_summary_stats_period_radius_ratios_neighbors!(css, cat_obs, sim_param)
    calc_summary_stats_radius_ratios_neighbors_photoevap_boundary_Carrera2018!(css, cat_obs, sim_param) # to also compute arrays of radius ratios above, below, and across photoevaporation boundary
    calc_summary_stats_duration_ratios_neighbors!(css, cat_obs, sim_param)
    calc_summary_stats_complexity_stats_GF2020!(css, cat_obs, sim_param)

    return css
end



"""
    test_summary_stats()

Test the summary statistics functions by simulating a physical and observed catalog and computing the summary statistics.

# Returns:
- `summary_stat::CatalogSummaryStatistics`: object containing the summary statistics.
"""
function test_summary_stats()
    sim_param = setup_sim_param_model()
    cat_phys = generate_kepler_physical_catalog(sim_param)
    cat_obs = observe_kepler_targets_single_obs(cat_phys, sim_param)
    summary_stat = calc_summary_stats_model(cat_obs, sim_param)
end



"""
    combine_summary_stats(ss1, ss2)

Combine two summary statistics, assuming they have the same keys.

# Arguments:
- `ss1::CatalogSummaryStatistics`: object containing the summary statistics of a catalog.
- `ss2::CatalogSummaryStatistics`: object containing the summary statistics of another catalog.
NOTE: these two summary statistics must have the same keys!

# Returns:
- `ssc::CatalogSummaryStatistics`: object containing the summary statistics of the combined catalog.
"""
function combine_summary_stats(ss1::CatalogSummaryStatistics, ss2::CatalogSummaryStatistics)
    @assert all(keys(ss1.stat) .== keys(ss2.stat))
    keys_all = keys(ss1.stat)
    stat_combined = Dict{String,Any}()
    for key in keys_all
        if any(key .== ["num_targets", "num_tranets", "num_n-tranet_systems"])
            # Sum these summary statistics
            stat_combined[key] = ss1.stat[key] .+ ss2.stat[key]
        else
            # Concatenate these summary statistics
            stat_combined[key] = [ss1.stat[key]; ss2.stat[key]]
        end
    end
    return ssc = CatalogSummaryStatistics(stat_combined, Dict{String,Any}())
end





mutable struct CatalogSummaryStatisticsCollection
    star_id_samples::Dict{String,Vector{Int64}} # for storing a collection of stellar samples (lists of star ids)
    css_samples::Dict{String,CatalogSummaryStatistics} # for storing a collection of CatalogSummaryStatistics objects
end

# Constructor function for a `CatalogSummaryStatisticsCollection` object
function CatalogSummaryStatisticsCollection()
    CatalogSummaryStatisticsCollection( Dict{String,Vector{Int64}}(), Dict{String,CatalogSummaryStatistics}() )
end



"""
    calc_summary_stats_collection_model(cat_obs, names_samples, star_id_samples, sim_param)

Compute all the summary statistics for an observed catalog, and for the same catalog split into a number of sub-samples given by lists of star ids in `star_id_samples`. Compiles these summary statistics into a `CatalogSummaryStatisticsCollection` object.

# Arguments:
- `cat_obs::KeplerObsCatalog`: an observed catalog.
- `names_samples::Vector{String}`: list of names for each sample (keys for each `CatalogSummaryStatistic` object in `css_samples`).
- `star_id_samples::Vector{Vector{Int64}}`: list of vectors of the star ids (i.e. row numbers of the input stellar catalog) in each sample, corresponding to `names_samples`.
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.

# Returns:
- `cssc::CatalogSummaryStatisticsCollection`: object containing a collection of `CatalogSummaryStatistics` objects, including the summary statistics for the full `cat_obs`, and also for sub-samples of `cat_obs` corresponding to the stars in `star_id_samples`.
"""
function calc_summary_stats_collection_model(cat_obs::KeplerObsCatalog, names_samples::Vector{String}, star_id_samples::Vector{Vector{Int64}}, sim_param::SimParam)
    @assert length(names_samples) == length(star_id_samples)

    cssc = CatalogSummaryStatisticsCollection()
    cssc.star_id_samples["all"] = collect(1:length(cat_obs.target))
    cssc.css_samples["all"] = calc_summary_stats_model(cat_obs, sim_param)

    for (i,name) in enumerate(names_samples)
        idx_sample = [(cat_obs.target[j].star.id in star_id_samples[i]) for j in 1:length(cat_obs.target)]
        cat_obs_sample = KeplerObsCatalog(cat_obs.target[idx_sample])
        cssc.star_id_samples[name] = star_id_samples[i]
        cssc.css_samples[name] = calc_summary_stats_model(cat_obs_sample, sim_param)
    end

    return cssc
end
