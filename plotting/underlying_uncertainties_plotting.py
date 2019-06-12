from ExoplanetsSysSim_functions import *





savefigures = False
loadfiles_directory = 'ACI/Simulated_Data/Julia_v0.7/Kepler_catalog_optimization/q1q17_dr25_gaia_fgk_stars79935/Clustered_P_R/f_high_incl_low_incl_mmr/Fit_rate_mult_P_Pratios_D_Dratios_dur_durratios_mmr/Some11_params_CRPDr_KS/Fixed_Rbreak3_Ncrit8/lc_0p2_5_lp_0p5_10_alphaP_-2_2_alphaR1_-4_2_alphaR2_-6_0_ecc_0_0p1_incl_inclmmr_0_90_sigmaR_0_0p5_sigmaP_0_0p3/targs79935_maxincl0_maxiters5000/sigma_i_greater_sigma_i_mmr/GP_med/'
savefigures_directory = 'ExoplanetsSysSim_Clusters/Figures/'
run_number = ''
model_name = 'Clustered_P_R_Model' + run_number





##### To load the underlying populations:

# To first read the number of simulated targets and bounds for the periods and radii:
N_sim, cos_factor, P_min, P_max, radii_min, radii_max = read_targets_period_radius_bounds(loadfiles_directory + 'periods%s.out' % run_number)

# To read the simulation parameters from the file:
param_vals_all = read_sim_params(loadfiles_directory + 'periods%s.out' % run_number)

# To load and process the simulated physical catalog of stars and planets:

#star_phys = load_star_phys(loadfiles_directory + 'physical_catalog_stars%s.txt' % run_number)
#cat_phys = load_cat_phys(loadfiles_directory + 'physical_catalog_planets%s.txt' % run_number)
#sssp_per_sys, sssp = compute_summary_stats_from_cat_phys(cat_phys=cat_phys, star_phys=star_phys)

sssp_per_sys, sssp = compute_summary_stats_from_cat_phys(file_name_path=loadfiles_directory, run_number=run_number)





##### To plot the simulated catalog as marginal distributions:

subdirectory = 'Paper_Figures/Models/Underlying/Clustered_P_R/' #'Paper_Figures/'; 'Talk_Figures/'

fig_size = (8,3) #size of each panel (figure)
fig_lbrt = [0.15, 0.3, 0.95, 0.925]

n_bins = 100
lw = 3 #linewidth

afs = 20 #axes labels font size
tfs = 20 #text labels font size
lfs = 16 #legend labels font size





##### To load and compute the same statistics for a large number of models, computing the confidence intervals for each bin:

loadfiles_directory = 'ACI/Simulated_Data/Julia_v0.7/Kepler_catalog_optimization/q1q17_dr25_gaia_fgk_stars79935/Clustered_P_R/f_high_incl_low_incl_mmr/Fit_rate_mult_P_Pratios_D_Dratios_dur_durratios_mmr/Some11_params_CRPDr_KS/Fixed_Rbreak3_Ncrit8/lc_0p2_5_lp_0p5_10_alphaP_-2_2_alphaR1_-4_2_alphaR2_-6_0_ecc_0_0p1_incl_inclmmr_0_90_sigmaR_0_0p5_sigmaP_0_0p3/targs79935_maxincl0_maxiters5000/sigma_i_greater_sigma_i_mmr/GP_best_models/'
runs = 100

Mtot_bins = np.arange(22)-0.5
Mtot_bins_mid = (Mtot_bins[:-1] + Mtot_bins[1:])/2.
Mtot_counts_all = []
Mtot_earth_counts_all = []

clustertot_bins = np.arange(12)-0.5 #includes 0 bin but will not count them
clustertot_bins_mid = (clustertot_bins[:-1] + clustertot_bins[1:])/2.
clustertot_counts_all = []

pl_per_cluster_bins = np.arange(17)-0.5 #includes 0 bin but will not count them
pl_per_cluster_bins_mid = (pl_per_cluster_bins[:-1] + pl_per_cluster_bins[1:])/2.
pl_per_cluster_counts_all = []

P_bins = np.logspace(np.log10(P_min), np.log10(P_max), n_bins+1)
P_bins_mid = (P_bins[:-1] + P_bins[1:])/2.
P_counts_all = []

Rm_bins = np.logspace(np.log10(1.), np.log10(10.), n_bins+1)
Rm_bins_mid = (Rm_bins[:-1] + Rm_bins[1:])/2.
Rm_counts_all = []

e_bins = np.linspace(0., 0.1, n_bins+1)
e_bins_mid = (e_bins[:-1] + e_bins[1:])/2.
e_counts_all = []

mass_bins = np.logspace(np.log10(0.07), np.log10(1e3), n_bins+1)
mass_bins_mid = (mass_bins[:-1] + mass_bins[1:])/2.
mass_counts_all = []

radii_bins = np.logspace(np.log10(radii_min), np.log10(radii_max), n_bins+1)
radii_bins_mid = (radii_bins[:-1] + radii_bins[1:])/2.
radii_counts_all = []

radii_ratio_bins = np.logspace(-1., 1., n_bins+1)
radii_ratio_bins_mid = (radii_ratio_bins[:-1] + radii_ratio_bins[1:])/2.
radii_ratio_counts_all = []

N_mH_bins = np.logspace(np.log10(8.), np.log10(200.), n_bins+1)
N_mH_bins_mid = (N_mH_bins[:-1] + N_mH_bins[1:])/2.
N_mH_counts_all = []

Rstar_bins = np.linspace(0.5, 2.5, n_bins+1)
Rstar_bins_mid = (Rstar_bins[:-1] + Rstar_bins[1:])/2.
Rstar_counts_all = []

for i in range(1,runs+1): #range(1,runs+1)
    run_number = i
    N_sim_i = read_targets_period_radius_bounds(loadfiles_directory + 'periods%s.out' % run_number)[0]
    sssp_per_sys_i, sssp_i = compute_summary_stats_from_cat_phys(file_name_path=loadfiles_directory, run_number=run_number)
    
    # Multiplicities:
    counts, bins = np.histogram(sssp_i['Mtot_all'], bins=Mtot_bins)
    counts[0] = N_sim_i - len(sssp_i['Mtot_all'])
    Mtot_counts_all.append(counts/float(np.sum(counts)))
    
    # Multiplicities for Earth-sized planets:
    Earth_bools_per_sys = (sssp_per_sys_i['radii_all'] > 0.75) & (sssp_per_sys_i['radii_all'] < 1.25)
    Earth_counts_per_sys = np.sum(Earth_bools_per_sys, axis=1)
    counts, bins = np.histogram(Earth_counts_per_sys, bins=Mtot_bins)
    counts[0] = N_sim_i - len(Earth_counts_per_sys)
    Mtot_earth_counts_all.append(counts/float(np.sum(counts)))
    
    # Numbers of clusters:
    counts, bins = np.histogram(sssp_i['clustertot_all'], bins=clustertot_bins)
    clustertot_counts_all.append(counts/float(np.sum(counts)))
    
    # Numbers of planets per cluster:
    counts, bins = np.histogram(sssp_i['pl_per_cluster_all'], bins=pl_per_cluster_bins)
    pl_per_cluster_counts_all.append(counts/float(np.sum(counts)))
    
    # Periods:
    counts, bins = np.histogram(sssp_i['P_all'], bins=P_bins)
    P_counts_all.append(counts/float(np.sum(counts)))
    
    # Period ratios:
    counts, bins = np.histogram(sssp_i['Rm_all'], bins=Rm_bins)
    Rm_counts_all.append(counts/float(np.sum(counts)))
    
    # Eccentricities:
    counts, bins = np.histogram(sssp_i['e_all'], bins=e_bins)
    e_counts_all.append(counts/float(np.sum(counts)))
    
    # Planet masses:
    counts, bins = np.histogram(sssp_i['mass_all'], bins=mass_bins)
    mass_counts_all.append(counts/float(np.sum(counts)))
    
    # Planet radii:
    counts, bins = np.histogram(sssp_i['radii_all'], bins=radii_bins)
    radii_counts_all.append(counts/float(np.sum(counts)))
    
    # Planet radii ratios:
    counts, bins = np.histogram(sssp_i['radii_ratio_all'], bins=radii_ratio_bins)
    radii_ratio_counts_all.append(counts/float(np.sum(counts)))
    
    # Separations:
    counts, bins = np.histogram(sssp_i['N_mH_all'], bins=N_mH_bins)
    N_mH_counts_all.append(counts/float(np.sum(counts)))
    
    # Stellar radii:
    counts, bins = np.histogram(sssp_i['Rstar_all'], bins=Rstar_bins)
    Rstar_counts_all.append(counts/float(np.sum(counts)))

Mtot_counts_all = np.array(Mtot_counts_all)
Mtot_earth_counts_all = np.array(Mtot_earth_counts_all)
clustertot_counts_all = np.array(clustertot_counts_all)
pl_per_cluster_counts_all = np.array(pl_per_cluster_counts_all)
P_counts_all = np.array(P_counts_all)
Rm_counts_all = np.array(Rm_counts_all)
e_counts_all = np.array(e_counts_all)
mass_counts_all = np.array(mass_counts_all)
radii_counts_all = np.array(radii_counts_all)
radii_ratio_counts_all = np.array(radii_ratio_counts_all)
N_mH_counts_all = np.array(N_mH_counts_all)
Rstar_counts_all = np.array(Rstar_counts_all)



Mtot_counts_16, Mtot_counts_84 = np.zeros(len(Mtot_bins_mid)), np.zeros(len(Mtot_bins_mid))
clustertot_counts_16, clustertot_counts_84 = np.zeros(len(clustertot_bins_mid)), np.zeros(len(clustertot_bins_mid))
pl_per_cluster_counts_16, pl_per_cluster_counts_84 = np.zeros(len(pl_per_cluster_bins_mid)), np.zeros(len(pl_per_cluster_bins_mid))
for b in range(len(Mtot_bins_mid)):
    counts_bin_sorted = np.sort(Mtot_counts_all[:,b])
    Mtot_counts_16[b], Mtot_counts_84[b] = counts_bin_sorted[16], counts_bin_sorted[84]
for b in range(len(clustertot_bins_mid)):
    counts_bin_sorted = np.sort(clustertot_counts_all[:,b])
    clustertot_counts_16[b], clustertot_counts_84[b] = counts_bin_sorted[16], counts_bin_sorted[84]
for b in range(len(pl_per_cluster_bins_mid)):
    counts_bin_sorted = np.sort(pl_per_cluster_counts_all[:,b])
    pl_per_cluster_counts_16[b], pl_per_cluster_counts_84[b] = counts_bin_sorted[16], counts_bin_sorted[84]

P_counts_16, P_counts_84 = np.zeros(n_bins), np.zeros(n_bins)
Rm_counts_16, Rm_counts_84 = np.zeros(n_bins), np.zeros(n_bins)
e_counts_16, e_counts_84 = np.zeros(n_bins), np.zeros(n_bins)
mass_counts_16, mass_counts_84 = np.zeros(n_bins), np.zeros(n_bins)
radii_counts_16, radii_counts_84 = np.zeros(n_bins), np.zeros(n_bins)
radii_ratio_counts_16, radii_ratio_counts_84 = np.zeros(n_bins), np.zeros(n_bins)
N_mH_counts_16, N_mH_counts_84 = np.zeros(n_bins), np.zeros(n_bins)
Rstar_counts_16, Rstar_counts_84 = np.zeros(n_bins), np.zeros(n_bins)
for b in range(n_bins):
    # Periods:
    counts_bin_sorted = np.sort(P_counts_all[:,b])
    P_counts_16[b], P_counts_84[b] = counts_bin_sorted[16], counts_bin_sorted[84]
    
    # Period ratios:
    counts_bin_sorted = np.sort(Rm_counts_all[:,b])
    Rm_counts_16[b], Rm_counts_84[b] = counts_bin_sorted[16], counts_bin_sorted[84]
    
    # Eccentricities:
    counts_bin_sorted = np.sort(e_counts_all[:,b])
    e_counts_16[b], e_counts_84[b] = counts_bin_sorted[16], counts_bin_sorted[84]
    
    # Planet masses:
    counts_bin_sorted = np.sort(mass_counts_all[:,b])
    mass_counts_16[b], mass_counts_84[b] = counts_bin_sorted[16], counts_bin_sorted[84]
    
    # Planet radii:
    counts_bin_sorted = np.sort(radii_counts_all[:,b])
    radii_counts_16[b], radii_counts_84[b] = counts_bin_sorted[16], counts_bin_sorted[84]
    
    # Planet radii ratios:
    counts_bin_sorted = np.sort(radii_ratio_counts_all[:,b])
    radii_ratio_counts_16[b], radii_ratio_counts_84[b] = counts_bin_sorted[16], counts_bin_sorted[84]
    
    # Separations:
    counts_bin_sorted = np.sort(N_mH_counts_all[:,b])
    N_mH_counts_16[b], N_mH_counts_84[b] = counts_bin_sorted[16], counts_bin_sorted[84]
    
    # Stellar radii:
    counts_bin_sorted = np.sort(Rstar_counts_all[:,b])
    Rstar_counts_16[b], Rstar_counts_84[b] = counts_bin_sorted[16], counts_bin_sorted[84]

#####





#'''
# Multiplicities:
ax = setup_fig_single(fig_size, fig_lbrt[0], fig_lbrt[1], fig_lbrt[2], fig_lbrt[3])
x = sssp['Mtot_all']
counts, bins = np.histogram(x, bins=np.max(x)+1, range=(-0.5, np.max(x)+0.5))
counts[0] = N_sim - len(x) #to compute the number of systems with no planets
bins_mid = (bins[:-1] + bins[1:])/2.
#plt.plot(bins_mid, counts/float(np.sum(counts)), 'o-', color='k', lw=lw, label='Simulated catalog')
plt.plot(bins_mid, counts/float(np.sum(counts)), drawstyle='steps-mid', color='k', lw=lw, label='Simulated catalog')
#plt.plot(Mtot_bins_mid, Mtot_counts_16, color='r', lw=1, ls='--', label=r'16% and 84%')
#plt.plot(Mtot_bins_mid, Mtot_counts_84, color='r', lw=1, ls='--')
plt.plot(Mtot_bins_mid, Mtot_counts_16, drawstyle='steps-mid', color='r', lw=1, ls='--', label=r'16% and 84%')
plt.plot(Mtot_bins_mid, Mtot_counts_84, drawstyle='steps-mid', color='r', lw=1, ls='--')
ax.tick_params(axis='both', labelsize=afs)
plt.xlim([0, 12]) #[0, np.max(x)]
plt.xlabel('Planets per system', fontsize=tfs)
plt.ylabel('Fraction', fontsize=tfs)
plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=lfs) #show the legend
if savefigures:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_underlying_multiplicities.pdf')
    plt.close()

# Number of clusters:
plot_fig_pdf_simple(fig_size, [sssp['clustertot_all']], [], x_min=0.5, x_max=len(clustertot_bins_mid)+0.5, n_bins=len(clustertot_bins_mid), log_y=True, lw=lw, xlabel_text=r'Clusters per system $N_c$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt)
plt.plot(clustertot_bins_mid, clustertot_counts_16, drawstyle='steps-mid', color='r', lw=1, ls='--', label='16')
plt.plot(clustertot_bins_mid, clustertot_counts_84, drawstyle='steps-mid', color='r', lw=1, ls='--', label='84')
plt.xlim([0.5, 7])
if savefigures:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_underlying_clusters.pdf')
    plt.close()

# Number of planets per cluster:
plot_fig_pdf_simple(fig_size, [sssp['pl_per_cluster_all']], [], x_min=0.5, x_max=len(pl_per_cluster_bins_mid)+0.5, y_max=0.3, n_bins=len(pl_per_cluster_bins_mid), lw=lw, xlabel_text=r'Planets per cluster $N_p$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt)
plt.plot(pl_per_cluster_bins_mid, pl_per_cluster_counts_16, drawstyle='steps-mid', color='r', lw=1, ls='--', label='16')
plt.plot(pl_per_cluster_bins_mid, pl_per_cluster_counts_84, drawstyle='steps-mid', color='r', lw=1, ls='--', label='84')
plt.xlim([0.5, 11])
if savefigures:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_underlying_planets_per_cluster.pdf')
    plt.close()

# Periods:
plot_fig_pdf_simple(fig_size, [sssp['P_all']], [], x_min=P_min, x_max=P_max, n_bins=n_bins, log_x=True, lw=lw, xticks_custom=[3,10,30,100,300], xlabel_text=r'$P$ (days)', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt)
plt.plot(P_bins_mid, P_counts_16, drawstyle='steps-mid', color='r', lw=1, ls='--', label='16')
plt.plot(P_bins_mid, P_counts_84, drawstyle='steps-mid', color='r', lw=1, ls='--', label='84')
if savefigures:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_underlying_periods.pdf')
    plt.close()

# Period ratios (all):
plot_fig_pdf_simple(fig_size, [sssp['Rm_all']], [], x_min=1., x_max=10., y_max=0.05, n_bins=n_bins, log_x=True, lw=lw, xticks_custom=[1,2,3,4,5,10], xlabel_text=r'$P_{i+1}/P_i$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt)
plt.plot(Rm_bins_mid, Rm_counts_16, drawstyle='steps-mid', color='r', lw=1, ls='--', label='16')
plt.plot(Rm_bins_mid, Rm_counts_84, drawstyle='steps-mid', color='r', lw=1, ls='--', label='84')
plt.minorticks_off()
if savefigures:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_underlying_periodratios.pdf')
    plt.close()

# Period ratios (< 5):
#plot_fig_pdf_simple(fig_size, [sssp['Rm_all'][sssp['Rm_all'] < 5]], [], x_min=1., x_max=5., n_bins=n_bins, log_x=True, lw=lw, xticks_custom=[1,2,3,4,5], xlabel_text=r'$P_{i+1}/P_i$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt)
#plt.minorticks_off()

# Eccentricities:
plot_fig_pdf_simple(fig_size, [sssp['e_all']], [], x_min=0., x_max=0.1, y_max=0.05, n_bins=n_bins, lw=lw, xlabel_text=r'$e$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt)
plt.plot(e_bins_mid, e_counts_16, drawstyle='steps-mid', color='r', lw=1, ls='--', label='16')
plt.plot(e_bins_mid, e_counts_84, drawstyle='steps-mid', color='r', lw=1, ls='--', label='84')
if savefigures:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_underlying_eccentricities.pdf')
    plt.close()

# Planet masses:
plot_fig_pdf_simple(fig_size, [sssp['mass_all']], [], x_min=0.07, x_max=1e3, n_bins=n_bins, log_x=True, lw=lw, xlabel_text=r'$M_p$ ($M_\oplus$)', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt)
plt.plot(mass_bins_mid, mass_counts_16, drawstyle='steps-mid', color='r', lw=1, ls='--', label='16')
plt.plot(mass_bins_mid, mass_counts_84, drawstyle='steps-mid', color='r', lw=1, ls='--', label='84')
if savefigures:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_underlying_masses.pdf')
    plt.close()

# Planet radii:
plot_fig_pdf_simple(fig_size, [sssp['radii_all']], [], x_min=radii_min, x_max=radii_max, y_max=0.025, n_bins=n_bins, log_x=True, lw=lw, xticks_custom=[0.5,1,2,4,10], xlabel_text=r'$R_p$ ($R_\oplus$)', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt)
plt.plot(radii_bins_mid, radii_counts_16, drawstyle='steps-mid', color='r', lw=1, ls='--', label='16')
plt.plot(radii_bins_mid, radii_counts_84, drawstyle='steps-mid', color='r', lw=1, ls='--', label='84')
if savefigures:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_underlying_radii.pdf')
    plt.close()

# Planet radii ratios:
plot_fig_pdf_simple(fig_size, [sssp['radii_ratio_all']], [], x_min=1e-1, x_max=10., y_max=0.06, n_bins=n_bins, log_x=True, lw=lw, xlabel_text=r'$R_{p,i+1}/R_{p,i}$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt)
plt.plot(radii_ratio_bins_mid, radii_ratio_counts_16, drawstyle='steps-mid', color='r', lw=1, ls='--', label='16')
plt.plot(radii_ratio_bins_mid, radii_ratio_counts_84, drawstyle='steps-mid', color='r', lw=1, ls='--', label='84')
if savefigures:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_underlying_radii_ratios.pdf')
    plt.close()

# Separations in mutual Hill radii:
plot_fig_pdf_simple(fig_size, [sssp['N_mH_all']], [], x_min=8., x_max=200., y_max=0.03, n_bins=n_bins, log_x=True, lw=lw, xlabel_text=r'$\Delta$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt)
plt.plot(N_mH_bins_mid, N_mH_counts_16, drawstyle='steps-mid', color='r', lw=1, ls='--', label='16')
plt.plot(N_mH_bins_mid, N_mH_counts_84, drawstyle='steps-mid', color='r', lw=1, ls='--', label='84')
if savefigures:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_underlying_deltas.pdf')
    plt.close()

# Stellar radii:
plot_fig_pdf_simple(fig_size, [sssp['Rstar_all']], [], x_min=0.5, x_max=2.5, n_bins=n_bins, lw=lw, xlabel_text=r'$R_\star (R_\odot)$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt)
plt.plot(Rstar_bins_mid, Rstar_counts_16, drawstyle='steps-mid', color='r', lw=1, ls='--', label='16')
plt.plot(Rstar_bins_mid, Rstar_counts_84, drawstyle='steps-mid', color='r', lw=1, ls='--', label='84')
if savefigures:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_underlying_stellar_radii.pdf')
    plt.close()

plt.show()
plt.close()
#'''





##### To compute some occurrence rates statistics:

# Fraction of stars with no planets/at least 1 planet:
fswnp, fswnp_16, fswnp_84 = np.median(Mtot_counts_all[:,0]), Mtot_counts_16[0], Mtot_counts_84[0]
fswp, fswp_16, fswp_84 = 1.-fswnp, 1.-fswnp_84, 1.-fswnp_16

# Mean number of planets per star:
pl_tot_all = np.sum((Mtot_counts_all*N_sim_i)*np.arange(len(Mtot_bins_mid)), axis=1)
pl_per_star_all = np.sort(pl_tot_all/np.float(N_sim_i))
mps, mps_16, mps_84 = np.median(pl_per_star_all), pl_per_star_all[16], pl_per_star_all[84]

# Mean number of planets per planetary system (i.e. stars with at least 1 planet):
n_1plus_all = N_sim_i - (Mtot_counts_all*N_sim_i)[:,0]
pl_per_1plus_all = np.sort(pl_tot_all/n_1plus_all)
mpps, mpps_16, mpps_84 = np.median(pl_per_1plus_all), pl_per_1plus_all[16], pl_per_1plus_all[84]

# Fraction of stars with Earth-sized planets:
fswne, fswne_16, fswne_84 = np.median(Mtot_earth_counts_all[:,0]), np.sort(Mtot_earth_counts_all[:,0])[16], np.sort(Mtot_earth_counts_all[:,0])[84]
fswe, fswe_16, fswe_84 = 1.-fswne, 1.-fswne_84, 1.-fswne_16

# Mean number of Earth-sized planets per star:
ep_tot_all = np.sum((Mtot_earth_counts_all*N_sim_i)*np.arange(len(Mtot_bins_mid)), axis=1)
ep_per_star_all = np.sort(ep_tot_all/np.float(N_sim_i))
meps, meps_16, meps_84 = np.median(ep_per_star_all), ep_per_star_all[16], ep_per_star_all[84]

# Mean number of Earth-sized planets per planetary system (i.e. stars with at least 1 planet):
ep_per_1plus_all = np.sort(ep_tot_all/n_1plus_all)
mepps, mepps_16, mepps_84 = np.median(ep_per_1plus_all), ep_per_1plus_all[16], ep_per_1plus_all[84]

# Fraction of planetary systems with at least one Earth-sized planet (i.e. probability that a planet-hosting star contains an Earth-sized planet):
fswpwe_all = np.sort((1.-Mtot_earth_counts_all[:,0])/(1.-Mtot_counts_all[:,0]))
fswpwe, fswpwe_16, fswpwe_84 = np.median(fswpwe_all), fswpwe_all[16], fswpwe_all[84]

