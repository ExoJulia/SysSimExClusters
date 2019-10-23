#using ExoplanetsSysSim
#using StatsFuns
#using JLD
#using DataFrames
#import ExoplanetsSysSim.StellarTable.df
#import ExoplanetsSysSim.StellarTable.usable
#import Compat: UTF8String, ASCIIString

function calc_distance_model_ks(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  d1 = calc_distance_num_planets(summary1,summary2,sim_param,verbose=verbose)
  d2 = calc_distance_num_n_tranet_systems(summary1,summary2,sim_param,verbose=verbose)
  #d3 = calc_distance_num_planets_binned(summary1,summary2,sim_param,verbose=verbose)
  #d4 = calc_distance_mean_std_log_period_depth(summary1,summary2,sim_param,verbose=verbose)
  d5 = calc_distance_ks_period(summary1,summary2,sim_param,verbose=verbose)
  d6 = calc_distance_ks_depth(summary1,summary2,sim_param,verbose=verbose)
  d7 = calc_distance_ks_period_ratios(summary1,summary2,sim_param,verbose=verbose)
  d8 = calc_distance_ks_duration_ratios(summary1,summary2,sim_param,verbose=verbose)
  return vcat(d1, d2, d5, d6, d7, d8)
end

function calc_distance_model_kl(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  d1 = calc_distance_kl_num_planets(summary1,summary2,sim_param,verbose=verbose)
  d2 = calc_distance_kl_num_n_tranet_systems(summary1,summary2,sim_param,verbose=verbose)
  #d2 = calc_distance_hellinger_num_n_tranet_systems(summary1,summary2,sim_param,verbose=verbose)
  #d3 = calc_distance_num_planets_binned(summary1,summary2,sim_param,verbose=verbose)
  #d4 = calc_distance_mean_std_log_period_depth(summary1,summary2,sim_param,verbose=verbose)
  d5 = calc_distance_kl_period(summary1,summary2,sim_param,verbose=verbose)
  d6 = calc_distance_kl_depth(summary1,summary2,sim_param,verbose=verbose)
  d7 = calc_distance_kl_period_ratios(summary1,summary2,sim_param,verbose=verbose)
  d8 = calc_distance_kl_duration_ratios(summary1,summary2,sim_param,verbose=verbose)
  return vcat(d1, d2, d5, d6, d7, d8)
end
calc_distance_model(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false) = calc_distance_model_kl(summary1,summary2,sim_param,verbose=verbose)

function calc_distance_num_planets(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
    np1 = haskey(summary1.stat,"num_tranets") ? summary1.stat["num_tranets"] : summary1.stat["expected planets detected"]
    np2 = haskey(summary2.stat,"num_tranets") ? summary2.stat["num_tranets"] : summary2.stat["expected planets detected"]
    #println("np1 = ",np1,", np2 = ",np2)

    dist_np = dist_L1_abs(np1/summary1.stat["num targets"], np2/summary2.stat["num targets"])
    #println("np1 (normalized) = ",np1/summary1.stat["num targets"],", np2 (normalized) = ",np2/summary2.stat["num targets"],", d = ",dist_np)
  return dist_np
end

function calc_distance_num_planets_binned(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
    np1 = haskey(summary1.stat,"planets table") ? summary1.stat["planets table"] : summary1.stat["expected planets table"]
    np2 = haskey(summary2.stat,"planets table") ? summary2.stat["planets table"] : summary2.stat["expected planets table"]
    #println("np1 = ",np1,", np2 = ",np2)

    dist_np_bin = zeros(length(np1))
    for n in 1:length(np1)
      dist_np_bin[n] = dist_L1_abs(np1[n]/summary1.stat["num targets"], np2[n]/summary2.stat["num targets"])
      #println("True # [Bin ", n,"] = ",np1[n],", Expected # [Bin ", n,"] = ",np2[n])
    end
    #println("np1 (normalized) = ",np1/summary1.stat["num targets"],", np2 (normalized) = ",np2/summary2.stat["num targets"],", dist = ",dist_np_bin)
  return dist_np_bin
end

function calc_distance_num_n_tranet_systems(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  max_tranets_in_sys = get_int(sim_param,"max_tranets_in_sys")
  d = zeros(max_tranets_in_sys)
  for n in 1:max_tranets_in_sys
    d[n] = dist_L1_abs(summary1.stat["num n-tranet systems"][n]/summary1.stat["num targets"], summary2.stat["num n-tranet systems"][n]/summary2.stat["num targets"])
  end
  return d
end

function calc_distance_mean_std_log_period_depth(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
    d1 = dist_L1_abs(summary1.stat["mean log10 P"], summary2.stat["mean log10 P"])
    d2 = dist_L1_abs(summary1.stat["mean log10 depth"], summary2.stat["mean log10 depth"])
    d3 = dist_L1_abs(summary1.stat["std log10 P"], summary2.stat["std log10 P"])
    d4 = dist_L1_abs(summary1.stat["std log10 depth"], summary2.stat["std log10 depth"])
    return [d1, d2, d3, d4]
end

function calc_distance_ks_period(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  samp1 = summary1.stat["P list"]
  samp2 = summary2.stat["P list"]
  return ExoplanetsSysSim.ksstats(samp1,samp2)[5]
end

function calc_distance_ks_depth(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  samp1 = summary1.stat["depth list"]
  samp2 = summary2.stat["depth list"]
  return ExoplanetsSysSim.ksstats(samp1,samp2)[5]
end

function calc_distance_ks_period_ratios(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  samp1 = summary1.stat["period_ratio_list"]
  samp2 = summary2.stat["period_ratio_list"]
  return ExoplanetsSysSim.ksstats(samp1,samp2)[5]
end

function calc_distance_ks_duration_ratios(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  samp1 = summary1.stat["duration_ratio_list"]
  samp2 = summary2.stat["duration_ratio_list"]
  return ExoplanetsSysSim.ksstats(samp1,samp2)[5]
end

# Function for Relative Entropy / K-L divergence
include("kde.jl")

function calc_distance_kl_num_planets(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
    np1 = (haskey(summary1.stat,"num_tranets") ? summary1.stat["num_tranets"] : summary1.stat["expected planets detected"])
    np2 = (haskey(summary2.stat,"num_tranets") ? summary2.stat["num_tranets"] : summary2.stat["expected planets detected"])
    ntarg1 = summary1.stat["num targets"]
    ntarg2 = summary2.stat["num targets"]
    mintarg = min(ntarg1,ntarg2)
    alpha1 = 1+np1
    beta1 = 1+ntarg1
    alpha2 = 1+np2
    beta2 = 1+ntarg2
    #= Approximate Estimated Rate Posterior as a Normal Distribution
    mu1 = alpha1/beta1
    sigma1 = mu1 /sqrt(alpha1)
    mu2 = alpha2/beta2
    sigma2 = mu2 /sqrt(alpha1)
    dist = 0.5*( (sigma2/sigma1)^2 - 1 + 2*log(sigma1/sigma2) + (mu2-mu1)^2/sigma1^2 )  # Normal approximation, substitute denom for difference in rates
    =#
    #= Approximate Estimated Rate Posterior as a Gamma Distribution
    dist = abs( (alpha2-alpha1)*digamma(alpha1) - lgamma(alpha1) + lgamma(alpha2)  ) # Gamma approximation, same beta, why need abs?
    dist += abs( (alpha1-alpha2)*digamma(alpha2) - lgamma(alpha2) + lgamma(alpha1)  )  # Gamma approximation, same beta, why need abs?
    print("# a1=",alpha1, " a2=",alpha2, " digamma=",digamma(alpha1), " ", digamma(alpha2)," gamma term=",lgamma(alpha1)-lgamma(alpha2))
    if ntarg1!=ntarg2
       print("# ntarg=",ntarg1," ",ntarg2," dist(before betas)=", dist)
       dist += alpha2*(log(beta1)-log(beta2)) + alpha1*(beta2-beta1)/beta1  # Gammma approximation, beta terms
       dist += alpha1*(log(beta2)-log(beta1)) + alpha2*(beta1-beta2)/beta2  # Gammma approximation, beta terms
       println(" dist (after betas)=",dist)
    end
    =#
    ## #= Approximate Estimated Rate Posterior as a Poisson Distribution
    #dist = abs( alpha2-alpha1+alpha1*log(alpha1/alpha2) )
    #dist += abs( alpha2-alpha1+alpha1*log(alpha1/alpha2) )
    rate1 = alpha1/beta1*mintarg
    rate2 = alpha2/beta2*mintarg
    dist = (rate1-rate2)*log(rate1/rate2)
    ## =#
    dist /= 2
    return dist
end

function calc_distance_kl_num_n_tranet_systems(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  max_tranets_in_sys = get_int(sim_param,"max_tranets_in_sys")
  #= categorical distribution
  f1sum = sum(summary1.stat["num n-tranet systems"]) # summary1.stat["num targets"]
  f2sum = sum(summary2.stat["num n-tranet systems"]) # summary2.stat["num targets"]
  if !(f1sum>0 && f2sum>0) return 0.0  end
  d = zeros(max_tranets_in_sys)
  for n in 1:max_tranets_in_sys
    f1 = summary1.stat["num n-tranet systems"][n]/f1sum
    f2 = summary2.stat["num n-tranet systems"][n]/f2sum
    m = (f1+f2)/2
    if m>zero(m)
    if f1>zero(f1)
       d[n] += 0.5*f1*log(f1/m)
    end
    if f2>zero(f2)
       d[n] += 0.5*f2*log(f2/m)
    end
#   else
#       d += Inf
    end
  end
  =#
  # Poisson distributions for each
    ntarg1 = summary1.stat["num targets"]
    ntarg2 = summary2.stat["num targets"]
    #mintarg = min(ntarg1,ntarg2)
  d = zeros(max_tranets_in_sys)
  for n in 1:max_tranets_in_sys
    f1 = (1+summary1.stat["num n-tranet systems"][n])/(1+ntarg1) # *(1+mintarg)
    f2 = (1+summary2.stat["num n-tranet systems"][n])/(1+ntarg2) # *(1+mintarg)
    d[n] = (f1-f2)*log(f1/f2)
  end
  return d
end

function calc_distance_hellinger_num_n_tranet_systems(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  max_tranets_in_sys = get_int(sim_param,"max_tranets_in_sys")
  f1sum = sum(summary1.stat["num n-tranet systems"]) # summary1.stat["num targets"]
  f2sum = sum(summary2.stat["num n-tranet systems"]) # summary2.stat["num targets"]
  d = 1
  for n in 1:max_tranets_in_sys
    f1 = summary1.stat["num n-tranet systems"][n]/f1sum
    f2 = summary2.stat["num n-tranet systems"][n]/f2sum
    d -= sqrt(f1*f2)
  end
  #d = sqrt(d)
  return d
end

function calc_distance_kl_period(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  #println("# P list 1: n=",length(summary1.stat["P list"])," min=",minimum(summary1.stat["P list"]), " max=",maximum(summary1.stat["P list"]))
  #println("# P list 2: n=",length(summary2.stat["P list"])," min=",minimum(summary2.stat["P list"]), " max=",maximum(summary2.stat["P list"]))
  samp1 = log.(summary1.stat["P list"] )
  samp2 = log.(summary2.stat["P list"] )
  calc_kl_distance_ab(samp1,samp2,log(0.5),log(320.) )
end

function calc_distance_kl_depth(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  samp1 = log.(summary1.stat["depth list"] )
  samp2 = log.(summary2.stat["depth list"] )
  calc_kl_distance_ab(samp1,samp2,log(0.000025),log(0.025) )
end

function calc_distance_kl_period_ratios(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  min_ratios_to_compute_distances = 3
  distance_when_not_enough_ratios = 10.0
  samp1 = summary1.stat["period_ratio_list"]
  samp2 = summary2.stat["period_ratio_list"]
  if length(samp1)<min_ratios_to_compute_distances || length(samp2)<min_ratios_to_compute_distances
     return distance_when_not_enough_ratios
  end
  calc_kl_distance_ab(samp1,samp2,0.0,1.0)
end

function calc_distance_kl_duration_ratios(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, sim_param::SimParam ; verbose::Bool = false)
  min_ratios_to_compute_distances = 3
  distance_when_not_enough_ratios = 10.0
  samp1 = log.( summary1.stat["duration_ratio_list"] )
  samp2 = log.( summary2.stat["duration_ratio_list"] )
  if length(samp1)<min_ratios_to_compute_distances || length(samp2)<min_ratios_to_compute_distances
     return distance_when_not_enough_ratios
  end
  calc_kl_distance_ab(samp1,samp2,-log(10.0),log(10.0) )
end


## test functions in this file, requires model, summary stats and functions from ExoplanetsSysSim
function test_distance()
  sim_param = setup_sim_param_model()
  cat_phys1 = generate_kepler_physical_catalog(sim_param)
  cat_ref = observe_kepler_targets_single_obs(cat_phys1,sim_param)
  summary_stat_ref = calc_summary_stats_model(cat_ref,sim_param)
  cat_phys2 = generate_kepler_physical_catalog(sim_param)
  cat_obs = observe_kepler_targets_single_obs(cat_phys2,sim_param)
  summary_stat_obs = calc_summary_stats_model(cat_obs,sim_param)
  dist = calc_distance_model(summary_stat_ref,summary_stat_obs,sim_param)
end





"""
    calc_all_distances_dict(sim_param, ss1, ss2; AD_mod=true)

Compute and return the distances between two planet catalogs as a dictionary.

# Arguments:
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.
- `ss1::ExoplanetsSysSim.CatalogSummaryStatistics`: an object containing the summary statistics of a catalog.
- `ss2::ExoplanetsSysSim.CatalogSummaryStatistics`: an object containing the summary statistics of another catalog.
- `AD_mod::Bool=true`: whether to use the original (if false) or modified (if true) AD distance.

# Returns:
- `dists::Dict{String,Float64}`: a dictionary containing the individual distance terms.
- `counts::Dict{String,Any}`: a dictionary containing the multiplicities, numbers of planets, and planet pairs.

NOTE: individual distance terms are set to `Inf` if they cannot be computed (i.e. if not enough planets)!
"""
function calc_all_distances_dict(sim_param::SimParam, ss1::ExoplanetsSysSim.CatalogSummaryStatistics, ss2::ExoplanetsSysSim.CatalogSummaryStatistics; AD_mod::Bool=true)

    max_incl_sys = get_real(sim_param, "max_incl_sys")
    cos_factor = cos(max_incl_sys*pi/180)

    Nmult1 = ss1.stat["num n-tranet systems"]
    Nmult2 = ss2.stat["num n-tranet systems"]
    M_obs1 = Int64[] # array to be filled with the number of transiting planets in each simulated system for ss1
    M_obs2 = Int64[] # array to be filled with the number of transiting planets in each simulated system for ss2
    max_k1, max_k2 = length(Nmult1), length(Nmult2)
    for k in 1:max_k1
        append!(M_obs1, k*ones(Int64, Nmult1[k]))
    end
    for k in 1:max_k2
        append!(M_obs2, k*ones(Int64, Nmult2[k]))
    end
    n_pl1, n_pairs1 = sum(Nmult1 .* collect(1:max_k1)), sum(Nmult1[2:end] .* collect(1:max_k1-1)) # total numbers of planets, and planet pairs
    n_pl2, n_pairs2 = sum(Nmult2 .* collect(1:max_k2)), sum(Nmult2[2:end] .* collect(1:max_k2-1)) # total numbers of planets, and planet pairs

    # Make a dictionary for the planet counts (multiplicity, total planets, total planet pairs):
    counts = Dict("Nmult1"=>Nmult1, "Nmult2"=>Nmult2, "n_pl1"=>n_pl1, "n_pl2"=>n_pl2, "n_pairs1"=>n_pairs1, "n_pairs2"=>n_pairs2)

    # Initialize a dictionary for all the distances:
    dists = Dict{String,Float64}()

    # Compute distances for rates of planets:
    dists["delta_f"] = abs(ss1.stat["num_tranets"]/(ss1.stat["num targets"]/cos_factor) - ss2.stat["num_tranets"]/(ss2.stat["num targets"]))
    dists["mult_KS"] = ExoplanetsSysSim.ksstats_ints(M_obs1, M_obs2)[5]
    dists["mult_CRPD"] = ExoplanetsSysSim.CRPDstats([Nmult1[1:4]; sum(Nmult1[5:end])], [Nmult2[1:4]; sum(Nmult2[5:end])]) # NOTE: binning 5+ planet systems together
    dists["mult_CRPD_r"] = ExoplanetsSysSim.CRPDstats([Nmult2[1:4]; sum(Nmult2[5:end])], [Nmult1[1:4]; sum(Nmult1[5:end])]) # NOTE: binning 5+ planet systems together

    # Compute KS and AD distances for the marginals:
    ADdist = AD_mod ? ExoplanetsSysSim.ADstats_mod : ExoplanetsSysSim.ADstats
    obs = ["periods", "pratios", "durations", "xis", "xis_nonmmr", "xis_mmr", "depths", "depths_above", "depths_below", "rratios", "rratios_above", "rratios_below", "rratios_across"]
    for (i,key) in enumerate(obs)
        # NOTE: if the KS or AD distance cannot be computed (i.e. not enough observed planets), assign Inf
        dists[key*"_KS"] = min(length(ss1.stat[key]), length(ss2.stat[key]))>0 ? ExoplanetsSysSim.ksstats(ss1.stat[key], ss2.stat[key])[5] : Inf
        dists[key*"_AD"] = min(length(ss1.stat[key]), length(ss2.stat[key]))>1 ? ADdist(ss1.stat[key], ss2.stat[key]) : Inf
    end

    return (dists, counts)
end
