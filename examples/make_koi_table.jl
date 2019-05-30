##### For writing simulated KOI table
using Printf
using DataFrames

function calc_snr_if_transit(t::ExoplanetsSysSim.KeplerTarget,s::Integer,p::Integer,sim_param::ExoplanetsSysSim.SimParam)
    period = t.sys[s].orbit[p].P
    size_ratio = t.sys[s].planet[p].radius/t.sys[s].star.radius
    depth = ExoplanetsSysSim.calc_transit_depth(t,s,p)
    duration = ExoplanetsSysSim.calc_transit_duration(t,s,p)
    b = ExoplanetsSysSim.calc_impact_parameter(t.sys[s],p)
    snr_correction = ExoplanetsSysSim.calc_depth_correction_for_grazing_transit(b,size_ratio)
    depth *= snr_correction
    kepid = StellarTable.star_table(t.sys[s].star.id,:kepid)
    osd_duration = ExoplanetsSysSim.get_durations_searched_Kepler(period,duration)
    osd = WindowFunction.interp_OSD_from_table(kepid,period,osd_duration)
    if osd_duration > duration
        osd *= osd_duration/duration
    end
    #ntr = calc_expected_num_transits(t,s,p,sim_param)
    snr = depth/osd*1.0e6
    return snr
end

function calc_snr_if_transit_nontransiting(t::ExoplanetsSysSim.KeplerTarget,s::Integer,p::Integer,sim_param::ExoplanetsSysSim.SimParam)
    period = t.sys[s].orbit[p].P
    size_ratio = t.sys[s].planet[p].radius/t.sys[s].star.radius
    depth = ExoplanetsSysSim.calc_transit_depth(t,s,p)
    duration_central = ExoplanetsSysSim.calc_transit_duration_central(t,s,p)
    b = rand()
    duration = duration_central * sqrt(1-b*b)
    snr_correction = ExoplanetsSysSim.calc_depth_correction_for_grazing_transit(b,size_ratio)
    depth *= snr_correction
    kepid = StellarTable.star_table(t.sys[s].star.id,:kepid)
    osd_duration = ExoplanetsSysSim.get_durations_searched_Kepler(period,duration)
    osd = WindowFunction.interp_OSD_from_table(kepid,period,osd_duration)
    if osd_duration > duration
        osd *= osd_duration/duration
    end
    #ntr = calc_expected_num_transits(t,s,p,sim_param)
    snr = depth/osd*1.0e6
    return snr
end

function make_koi_list(cat_phys,cat_obs; output_nondetections::Bool = true)
    df = DataFrame()
    # Entries following Jason Rowe's tables for Arch 3 paper
    df.KIC = Int64[]
    df.KOI = String[]
    df.Period = Float64[]
    df.Period_e = Float64[]
    df.Epoch = Float64[]
    df.RpRs = Float64[]
    df.RpRs_ep = Float64[]
    df.RpRs_em = Float64[]
    df.b = Float64[]
    df.b_ep = Float64[]
    df.b_em = Float64[]
    df.tdepth = Float64[]
    df.tdepth_e = Float64[]
    df.tdur = Float64[]
    df.tdur_e = Float64[]
    df.radius = Float64[]
    df.radius_ep = Float64[]
    df.radius_em = Float64[]
    df.rhostarm = Float64[]
    df.rhostarm_ep = Float64[]
    df.rhostarm_em = Float64[]
    df.teff = Float64[]
    df.teff_e = Float64[]
    df.rstar = Float64[]
    df.rstar_ep = Float64[]
    df.rstar_em = Float64[]
    df.TTVflag = Int64[]
    df.SNR = Float64[]
    df.ModelSource = Int64[]
    df.DR24_Status = String[]
    df.DR25_Status = String[]
    df.DR25_supStatus = String[]
    # Extras added for Dan Fabrycky
    df.mass_true = Float64[]
    df.Period_true = Float64[]
    df.ecc_true = Float64[]
    df.incl_true = Float64[]
    df.omega_true = Float64[]
    df.asc_node_true = Float64[]
    df.mean_anom_true = Float64[]
    df.b_true = Float64[]
    df.syssim_targ_id = Int64[]
    df.syssim_obs_pl_id = Int64[]
    df.syssim_phys_pl_id = Int64[]
    koi_star = 0
    #println(f, "target_id star_id planet_mass planet_radius period ecc")
    for (i,targ) in enumerate(cat_obs.target)
        if length(targ.obs) >= 1
            star = StellarTable.star_table(targ.star.id)
            rn = randn()
            rstar_obs = star[:radius] + rn * ( rn>0 ? star[:radius_err1] : abs(star[:radius_err2]))
            mstar_obs = star[:mass] + rn * ( rn>0 ? star[:mass_err1] : abs(star[:mass_err2]))
            rhostar_obs = star[:dens] + rn * ( rn>0 ? star[:dens_err1] : abs(star[:dens_err2]))
            sys = cat_phys.target[i].sys[1]
            phys_to_obs_id = zeros(length(sys.orbit))
            koi_star += 1
            koi_pl = 0
            for (j,obs) in enumerate(targ.obs)
                koi_pl += 1
                sigma = targ.sigma[j]
                entry = Dict()

                # Find corresponding physical planet
                k = findmin(abs.(obs.period .- map(k->sys.orbit[k].P,1:length(sys.orbit))))[2]
                phys_to_obs_id[k] = j
                entry[:mass_true] = sys.planet[k].mass
                entry[:Period_true] = sys.orbit[k].P
                entry[:ecc_true] = sys.orbit[k].ecc
                entry[:incl_true] = sys.orbit[k].incl
                entry[:omega_true] = sys.orbit[k].omega
                entry[:asc_node_true] = sys.orbit[k].asc_node
                entry[:mean_anom_true] = sys.orbit[k].mean_anom
                entry[:b_true] = ExoplanetsSysSim.calc_impact_parameter(cat_phys.target[i].sys[1],k)

                # Create entry in KOI table
                entry[:KIC] = star[:kepid]
                entry[:KOI] = @sprintf("%04d.%02d",koi_star,koi_pl)
                entry[:Period] = obs.period
                entry[:Period_e] = sigma.period
                entry[:Epoch] = obs.t0
                entry[:RpRs] = sqrt(obs.depth) # TODO: account for limb darkening
                entry[:RpRs_ep] = sqrt(0.5) * sigma.depth/obs.depth * entry[:RpRs]
                entry[:RpRs_em] = -entry[:RpRs_ep]
                D_central_circ = ExoplanetsSysSim.calc_transit_duration_central_circ(cat_phys.target[i],1,k)
                b = (obs.duration<D_central_circ) ? sqrt(1-(obs.duration/D_central_circ)^2) : 0.0
                entry[:b] = b # ExoplanetsSysSim.calc_impact_parameter(cat_phys.target[i].sys[1],k)
                entry[:b_ep] = min(1+entry[:RpRs]-entry[:b],obs.duration * sigma.duration / (max(entry[:b_true],entry[:RpRs])*D_central_circ^2)) # TODO: deal with b \approx 0 case
                entry[:b_em] = -min(entry[:b_ep],entry[:b])
                entry[:tdepth] = 1e6*obs.depth
                entry[:tdepth_e] = 1e6*sigma.depth
                entry[:tdur] = 24*obs.duration
                entry[:tdur_e] = 24*sigma.duration
                entry[:radius] = entry[:RpRs] * rstar_obs / ExoplanetsSysSim.earth_radius
                entry[:radius_ep] = (0.5 * (sigma.depth/obs.depth)^2 + (star[:radius_err1]/rstar_obs)^2) * entry[:radius]
                entry[:radius_em] = -(0.5 * (sigma.depth/obs.depth)^2 + (star[:radius_err2]/rstar_obs)^2) * entry[:radius]
                entry[:rhostarm] = rhostar_obs
                entry[:rhostarm_ep] = star[:dens_err1]
                entry[:rhostarm_em] = star[:dens_err2]
                entry[:teff] = star[:teff]
                entry[:teff_e] = 0.0
                entry[:rstar] = rstar_obs
                entry[:rstar_ep] = star[:radius_err1]
                entry[:rstar_em] = star[:radius_err2]
                entry[:TTVflag] = 0
                entry[:SNR] = calc_snr_if_transit(cat_phys.target[i],1,k,sim_param)
                entry[:ModelSource] = 0
                entry[:DR24_Status] = "NA"
                entry[:DR25_Status] = "Candidate"
                entry[:DR25_supStatus] = "NA"
                entry[:syssim_targ_id] = i
                entry[:syssim_obs_pl_id] = j
                entry[:syssim_phys_pl_id] = k

                push!(df,entry)
            end
            # Add non-transiting planets
            if !output_nondetections || (length(sys.orbit) == length(targ.obs)) # No undetected planets
                continue
            end
            for j in 1:length(sys.orbit)
                if phys_to_obs_id[j] != 0  # matches an observed planet
                    continue
                end
                koi_pl += 1
                entry = Dict()
                entry[:mass_true] = cat_phys.target[i].sys[1].planet[j].mass
                entry[:Period_true] = cat_phys.target[i].sys[1].orbit[j].P
                entry[:ecc_true] = cat_phys.target[i].sys[1].orbit[j].ecc
                entry[:incl_true] = cat_phys.target[i].sys[1].orbit[j].incl
                entry[:omega_true] = cat_phys.target[i].sys[1].orbit[j].omega
                entry[:asc_node_true] = cat_phys.target[i].sys[1].orbit[j].asc_node
                entry[:mean_anom_true] = cat_phys.target[i].sys[1].orbit[j].mean_anom
                entry[:b_true] = ExoplanetsSysSim.calc_impact_parameter(cat_phys.target[i].sys[1],j)

                entry[:KIC] = star[:kepid]
                entry[:KOI] = @sprintf("%04d.%02d",koi_star,koi_pl)
                entry[:Period] = 0.0
                entry[:Period_e] = 0.0
                entry[:Epoch] = 0.0
                entry[:RpRs] = 0.0
                entry[:RpRs_ep] = 0.0
                entry[:RpRs_em] = 0.0
                entry[:b] = 0.0
                entry[:b_ep] = 0.0
                entry[:b_em] = 0.0
                entry[:tdepth] = 0.0
                entry[:tdepth_e] = 0.0
                entry[:tdur] = 0.0
                entry[:tdur_e] = 0.0
                entry[:radius] = 0.0
                entry[:radius_ep] = 0.0
                entry[:radius_em] = 0.0
                entry[:rhostarm] = rhostar_obs
                entry[:rhostarm_ep] = star[:dens_err1]
                entry[:rhostarm_em] = star[:dens_err2]
                entry[:teff] = star[:teff]
                entry[:teff_e] = 0.0
                entry[:rstar] = rstar_obs
                entry[:rstar_ep] = star[:radius_err1]
                entry[:rstar_em] = star[:radius_err2]
                entry[:TTVflag] = 0
                entry[:SNR] = calc_snr_if_transit_nontransiting(cat_phys.target[i],1,j,sim_param)
                entry[:ModelSource] = 0
                entry[:DR24_Status] = "NA"
                entry[:DR25_Status] = "Not Detected"
                entry[:DR25_supStatus] = "NA"
                entry[:syssim_targ_id] = i
                entry[:syssim_obs_pl_id] = 0
                entry[:syssim_phys_pl_id] = j
                push!(df,entry)

            end
        end
    end
    return df
end
