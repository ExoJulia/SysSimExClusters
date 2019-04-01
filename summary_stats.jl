#using ExoplanetsSysSim
#using StatsFuns
#using JLD
#using DataFrames
#import ExoplanetsSysSim.StellarTable.df
#import ExoplanetsSysSim.StellarTable.usable

include("misc_functions.jl")

## summary_statistics borrowed from multiple_planets example, eventually should be merged into main branch

# Compile indices for N-tranet systems
function calc_summary_stats_idx_n_tranets!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  if haskey(css.cache,"idx_n_tranets")
     return css.cache["idx_n_tranets"]
  end
     max_tranets_in_sys = get_int(param,"max_tranets_in_sys")    
     idx_n_tranets = Vector{Int64}[ Int64[] for m = 1:max_tranets_in_sys]
     for n in 1:max_tranets_in_sys-1
       idx_n_tranets[n] = findall(x::KeplerTargetObs-> length(x.obs) == n, cat_obs.target )
     end
     idx_n_tranets[max_tranets_in_sys] = findall(x::KeplerTargetObs-> length(x.obs) >= max_tranets_in_sys, cat_obs.target )
     css.cache["idx_n_tranets"] = idx_n_tranets
  return idx_n_tranets 
end

# Count total number of tranets using lists of indices for N-tranet systems
function calc_summary_stats_num_tranets!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  if haskey(css.stat,"num_tranets")
     return css.stat["num_tranets"]
  elseif haskey(css.cache,"num_tranets")
     css.stat["num_tranets"] = css.cache["num_tranets"]
     return css.stat["num_tranets"]
  end
     idx_n_tranets = calc_summary_stats_idx_n_tranets!(css,cat_obs,param)
     num_tranets = 0
     for n in 1:length(idx_n_tranets)
         num_tranets += n*length(idx_n_tranets[n])
     end
     css.stat["num_tranets"] = num_tranets
  return num_tranets 
end

function calc_summary_stats_num_targets!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam ; trueobs_cat::Bool = false)
  if !trueobs_cat
    css.stat["num targets"] = get_int(param,"num_targets_sim_pass_one")
  else
    css.stat["num targets"] = get_int(param,"num_kepler_targets")
  end
end

function calc_summary_stats_num_n_tranet_systems!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  idx_n_tranets = calc_summary_stats_idx_n_tranets!(css,cat_obs,param)
  #max_tranets_in_sys = get_int(param,"max_tranets_in_sys")    
  num_n_tranet_systems = map(n->length(idx_n_tranets[n]), 1:length(idx_n_tranets) )
  #for n in 1:length(idx_n_tranets)
  #  num_n_tranet_systems[n] = length(idx_n_tranets[n])
  #end
  css.stat["num n-tranet systems"] = num_n_tranet_systems
end

function calc_summary_stats_duration_ratios_neighbors!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  if haskey(css.stat,"duration_ratio_list")
     return css.stat["duration_ratio_list"]
  elseif haskey(css.cache,"duration_ratio_list")
     return css.cache["duration_ratio_list"]
  end
  idx_n_tranets = calc_summary_stats_idx_n_tranets!(css,cat_obs,param)
  @assert length(idx_n_tranets) >= 1 

  # Calculate how many duration ratios there will be & allocate storage
  num_ratios = 0
  for i in 2:length(idx_n_tranets)
     num_ratios += length(idx_n_tranets[i])*(i-1)
  end
  duration_ratio_list = Array{Float64}(undef, num_ratios)
  duration_ratio_non_mmr_list = Float64[]
  duration_ratio_near_mmr_list = Float64[]

  k = 0
  for n in 2:length(idx_n_tranets)         # Loop over number of tranets in system
    for i in idx_n_tranets[n]              # Loop over systems with n tranets
       period_in_sys = Array{Float64}(undef, n)
       duration_in_sys = Array{Float64}(undef, n)
       for j in 1:n                        # Loop over periods within a system
         period_in_sys[j] = cat_obs.target[i].obs[j].period
         duration_in_sys[j] = cat_obs.target[i].obs[j].duration
       end
       perm = sortperm(period_in_sys)
       for j in 1:(n-1)                       # Loop over period ratios within a system
          period_ratio = period_in_sys[perm[j+1]]/period_in_sys[perm[j]]
          if 1<period_ratio<Inf
             k += 1
             xi = duration_in_sys[perm[j]]/duration_in_sys[perm[j+1]] * period_ratio^(1//3)
             duration_ratio_list[k] = xi

             if is_period_ratio_near_resonance(period_ratio, param)
                append!(duration_ratio_near_mmr_list, xi)
             else
                append!(duration_ratio_non_mmr_list, xi)
             end
          end
       end
    end
  end
  resize!(duration_ratio_list,k)
  css.stat["duration_ratio_list"] = duration_ratio_list
  css.stat["duration_ratio_non_mmr_list"] = duration_ratio_non_mmr_list
  css.stat["duration_ratio_near_mmr_list"] = duration_ratio_near_mmr_list

  return duration_ratio_list, duration_ratio_non_mmr_list, duration_ratio_near_mmr_list
end

function calc_summary_stats_period_radius_ratios_neighbors_internal!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  #=
  if haskey(css.stat,"period_ratio_list") && haskey(css.stat,"radius_ratio_list")
     return (css.stat["period_ratio_list"], css.stat["radius_ratio_list"])
  elseif haskey(css.cache,"period_ratio_list") && haskey(css.cache,"radius_ratio_list")
     return (css.cache["period_ratio_list"], css.cache["radius_ratio_list"])
  end
  =#
  idx_n_tranets = calc_summary_stats_idx_n_tranets!(css,cat_obs,param)
  @assert length(idx_n_tranets) >= 1 

  # Calculate how many period ratios there will be & allocate storage
  num_ratios = 0
  for i in 2:length(idx_n_tranets)
     num_ratios += length(idx_n_tranets[i])*(i-1)
  end
  period_ratio_list = Array{Float64}(undef, num_ratios)
  radius_ratio_list = Array{Float64}(undef, num_ratios)

  radius_ratio_above_list = Float64[] #list to be filled with the radius ratios of adjacent planet pairs, both above the photoevaporation boundary
  radius_ratio_below_list = Float64[] #list to be filled with the radius ratios of adjacent planet pairs, both below the boundary
  radius_ratio_across_list = Float64[] #list to be filled with the radius ratios of adjacent planet pairs, across the boundary

  k = 0
  for n in 2:length(idx_n_tranets)         # Loop over number of tranets in system
    period_in_sys = Array{Float64}(undef, n)
    radius_in_sys = Array{Float64}(undef, n)
    depth_in_sys = Array{Float64}(undef, n)
    for i in idx_n_tranets[n]              # Loop over systems with n tranets
       for j in 1:n                        # Loop over periods within a system
         period_in_sys[j] = cat_obs.target[i].obs[j].period
         radius_in_sys[j] = sqrt(cat_obs.target[i].obs[j].depth)*cat_obs.target[i].star.radius
         depth_in_sys[j] = cat_obs.target[i].obs[j].depth
       end
       perm = sortperm(period_in_sys)
       for j in 1:(n-1)                       # Loop over period ratios within a system
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
  css.cache["period_ratio_list"] = period_ratio_list
  css.cache["radius_ratio_list"] = radius_ratio_list

  css.cache["radius_ratio_above_list"] = radius_ratio_above_list
  css.cache["radius_ratio_below_list"] = radius_ratio_below_list
  css.cache["radius_ratio_across_list"] = radius_ratio_across_list

  return (period_ratio_list, radius_ratio_list, radius_ratio_above_list, radius_ratio_below_list, radius_ratio_across_list)
end


function calc_summary_stats_period_radius_ratios_neighbors!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  (period_ratio_list,radius_ratio_list) = calc_summary_stats_period_radius_ratios_neighbors_internal!(css,cat_obs,param)[1:2]
  css.stat["period_ratio_list"] = period_ratio_list
  css.stat["radius_ratio_list"] = radius_ratio_list
  return (period_ratio_list, radius_ratio_list)
end

function calc_summary_stats_period_ratios_neighbors!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  period_ratio_list = calc_summary_stats_period_radius_ratios_neighbors_internal!(css,cat_obs,param)[1]
  css.stat["period_ratio_list"] = period_ratio_list
  return period_ratio_list
end

function calc_summary_stats_radius_ratios_neighbors!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  radius_ratio_list = calc_summary_stats_period_radius_ratios_neighbors_internal!(css,cat_obs,param)[2]
  css.stat["radius_ratio_list"] = radius_ratio_list
  return radius_ratio_list
end

function calc_summary_stats_radius_ratios_neighbors_photoevap_boundary_Carrera2018!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
    (radius_ratio_above_list, radius_ratio_below_list, radius_ratio_across_list) = calc_summary_stats_period_radius_ratios_neighbors_internal!(css,cat_obs,param)[3:5]
    css.stat["radius_ratio_above_list"] = radius_ratio_above_list
    css.stat["radius_ratio_below_list"] = radius_ratio_below_list
    css.stat["radius_ratio_across_list"] = radius_ratio_across_list
    return (radius_ratio_above_list, radius_ratio_below_list, radius_ratio_across_list)
end

function calc_summary_stats_mean_std_log_period_depth!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  # Allocate arrays to store values for each tranet
  num_tranets  = calc_summary_stats_num_tranets!(css, cat_obs, param)
  period_list = zeros(num_tranets)
  depth_list = zeros(num_tranets)
  #weight_list = ones(num_tranets)

  idx_n_tranets = calc_summary_stats_idx_n_tranets!(css, cat_obs, param)
  max_tranets_in_sys = get_int(param,"max_tranets_in_sys") 
  @assert max_tranets_in_sys >= 1
   i = 1   # tranet id
   for targ in cat_obs.target                        # For each target 
     for j in 1:min(length(targ.obs),max_tranets_in_sys)          # For each tranet around that target (but truncated if too many tranets in one system)
         #println("# i= ",i," j= ",j)
         period_list[i] = targ.obs[j].period
         depth_list[i] = targ.obs[j].depth
         #weight_list[i] = 1.0
         i = i+1
      end
   end

  css.cache["P list"] = period_list                                     # We can store whole lists, e.g., if we want to compute K-S distances
  css.cache["depth list"] = depth_list
  #css.cache["weight list"] = weight_list

  idx_good = Bool[ period_list[i]>0.0 && depth_list[i]>0.0 for i in 1:length(period_list) ]
  log_period_list = log10(period_list[idx_good])
  log_depth_list = log10(depth_list[idx_good])
  css.stat["mean log10 P"]  =  mean_log_P = mean(log_period_list)
  css.stat["mean log10 depth"]  =  mean_log_depth = mean(log_depth_list)
  css.stat["std log10 P"]  = std_log_P = stdm(log_period_list,mean_log_P)
  css.stat["std log10 depth"]  = std_log_depth = stdm(log_depth_list,mean_log_depth)

  return (mean_log_P, std_log_P, mean_log_depth, std_log_depth) 
end

function calc_summary_stats_cuml_period_depth_duration!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
  # Allocate arrays to store values for each tranet
  num_tranets  = calc_summary_stats_num_tranets!(css, cat_obs, param)
  period_list = zeros(num_tranets)
  depth_list = zeros(num_tranets)
  duration_list = zeros(num_tranets)
  #weight_list = ones(num_tranets)

  depth_above_list = Float64[] #list to be filled with the transit depths of planets above the photoevaporation boundary in Carrera et al 2018
  depth_below_list = Float64[] #list to be filled with the transit depths of planets below the boundary

  idx_n_tranets = calc_summary_stats_idx_n_tranets!(css, cat_obs, param)
  max_tranets_in_sys = get_int(param,"max_tranets_in_sys") 
  @assert max_tranets_in_sys >= 1
  i = 0   # tranet id
  for targ in cat_obs.target                        # For each target
      for j in 1:min(length(targ.obs),max_tranets_in_sys)          # For each tranet around that target (but truncated if too many tranets in one system)
          i = i+1
          #println("# i= ",i," j= ",j)
          period_list[i] = targ.obs[j].period
          depth_list[i] = targ.obs[j].depth
          duration_list[i] = targ.obs[j].duration
          #weight_list[i] = 1.0

          radius_earths, period = (sqrt(targ.obs[j].depth)*targ.star.radius)/ExoplanetsSysSim.earth_radius, targ.obs[j].period
          if photoevap_boundary_Carrera2018(radius_earths, period) == 1
             append!(depth_above_list, targ.obs[j].depth)
          elseif photoevap_boundary_Carrera2018(radius_earths, period) == 0
             append!(depth_below_list, targ.obs[j].depth)
          end
      end
  end
  resize!(period_list,i)
  resize!(depth_list,i)
  resize!(duration_list,i)
  css.stat["P list"] = period_list                                     # We can store whole lists, e.g., if we want to compute K-S distances
  css.stat["depth list"] = depth_list
  css.stat["duration list"] = duration_list
  #css.cache["weight list"] = weight_list

  css.cache["depth above list"] = depth_above_list
  css.cache["depth below list"] = depth_below_list

  #println("# P list = ",period_list)
  return (period_list, depth_list, duration_list, depth_above_list, depth_below_list)
end

function calc_summary_stats_depths_photoevap_boundary_Carrera2018!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam)
    (depth_above_list, depth_below_list) = calc_summary_stats_cuml_period_depth_duration!(css, cat_obs, param)[4:5]
    css.stat["depth above list"] = depth_above_list
    css.stat["depth below list"] = depth_below_list
    return (depth_above_list, depth_below_list)
end

function calc_summary_stats_obs_binned_rates!(css::CatalogSummaryStatistics, cat_obs::KeplerObsCatalog, param::SimParam; trueobs_cat::Bool = false)
  num_tranets  = calc_summary_stats_num_tranets!(css, cat_obs, param)
  idx_n_tranets = calc_summary_stats_idx_n_tranets!(css, cat_obs, param)

  idx_tranets = findall(x::KeplerTargetObs-> length(x.obs) > 0, cat_obs.target)::Array{Int64,1}             # Find indices of systems with at least 1 tranet = potentially detectable transiting planet
  css.cache["idx_tranets"] = idx_tranets                                   # We can save lists of indices to summary stats for pass 2, even though we will not use these for computing a distance or probability

  #=
  max_tranets_in_sys = get_int(param,"max_tranets_in_sys")    # Demo that simulation parameters can specify how to evalute models, too
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
  if ( length( findall(x::KeplerTargetObs-> length(x.obs) > max_tranets_in_sys, cat_obs.target ) ) > 0)   # Make sure max_tranets_in_sys is at least big enough for observed systems
    warn("Observational data has more transiting planets in one systems than max_tranets_in_sys allows.")
  end
  num_tranets  = convert(Int64,num_tranets)            # TODO OPT: Figure out why is not this already an Int.  I may be doing something that prevents some optimizations
  css.cache["num_tranets"] = num_tranets                                   

  num_sys_tranets = zeros(max_tranets_in_sys)                           # Since observed data, don't need to calculate probabilities.
  for n in 1:max_tranets_in_sys                                         # Make histogram of N-tranet systems
    num_sys_tranets[n] = length(idx_n_tranets[n])
  end
  css.stat["num_sys_tranets"] = num_sys_tranets
  css.stat["planets detected"] = num_tranets 
  =#

  period_list = zeros(num_tranets)
  radius_list = zeros(num_tranets)
  weight_list = ones(num_tranets)

  n = 1    # tranet id
  for i in 1:length(cat_obs.target)
    for j in 1:num_planets(cat_obs.target[i])
      period_list[n] = cat_obs.target[i].obs[j].period
      radius_list[n] = sqrt(cat_obs.target[i].obs[j].depth)*cat_obs.target[i].star.radius
      n = n+1
    end
  end

  limitP::Array{Float64,1} = get_any(param, "p_lim_arr", Array{Float64,1})
  limitRp::Array{Float64,1} = get_any(param, "r_lim_arr", Array{Float64,1})
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


function calc_summary_stats_model(cat_obs::KeplerObsCatalog, param::SimParam; trueobs_cat::Bool = false)
  css = CatalogSummaryStatistics()
  calc_summary_stats_num_targets!(css,cat_obs,param,trueobs_cat=trueobs_cat)
  calc_summary_stats_num_tranets!(css,cat_obs,param)
  calc_summary_stats_num_n_tranet_systems!(css,cat_obs,param)
  calc_summary_stats_cuml_period_depth_duration!(css,cat_obs,param)
  calc_summary_stats_depths_photoevap_boundary_Carrera2018!(css,cat_obs,param) #to also compute arrays of depths above and below photoevaporation boundary
  #calc_summary_stats_obs_binned_rates!(css,cat_obs,param)
  #calc_summary_stats_mean_std_log_period_depth!(css,cat_obs,param)
  calc_summary_stats_period_radius_ratios_neighbors!(css,cat_obs,param)
  calc_summary_stats_radius_ratios_neighbors_photoevap_boundary_Carrera2018!(css,cat_obs,param) #to also compute arrays of radius ratios above, below, and across photoevaporation boundary
  calc_summary_stats_duration_ratios_neighbors!(css,cat_obs,param)
  return css
end

## test functions in this file, requires model and function from ExoplanetsSysSim
function test_summary_stats()
  sim_param = setup_sim_param_model()
  cat_phys = generate_kepler_physical_catalog(sim_param)
  cat_obs = observe_kepler_targets_single_obs(cat_phys,sim_param)
  summary_stat = calc_summary_stats_model(cat_obs,sim_param)
end
