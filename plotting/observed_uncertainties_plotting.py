from ExoplanetsSysSim_functions import *





savefigures = False
#loadfiles_directory = '/Users/hematthi/Documents/GradSchool/Eric_Ford_Research/ExoplanetsSysSim_Clusters/SysSimExClusters/examples/'
loadfiles_directory = 'ACI/Simulated_Data/Julia_v0.7/Kepler_catalog_optimization/q1q17_dr25_gaia_fgk_stars79935/Clustered_P_R/f_high_incl_low_incl_mmr/Fit_rate_mult_P_Pratios_D_Dratios_dur_durratios_mmr/Some11_params_CRPDr_KS/Fixed_Rbreak3_Ncrit8/lc_0p2_5_lp_0p5_10_alphaP_-2_2_alphaR1_-4_2_alphaR2_-6_0_ecc_0_0p1_incl_inclmmr_0_90_sigmaR_0_0p5_sigmaP_0_0p3/targs79935_maxincl0_maxiters5000/sigma_i_greater_sigma_i_mmr/GP_med/'
savefigures_directory = 'ExoplanetsSysSim_Clusters/Figures/'
run_number = ''
model_name = 'Clustered_P_R_Model' + run_number #'Non_Clustered_Model', 'Clustered_P_Model', 'Clustered_P_R_Model'

AD_mod = 'true' # 'true' or 'false'
dists_exclude = np.array([2,3,8,12,13,15,16,17]) - 1





##### To load the files with the systems with observed planets:

# To first read the number of simulated targets and bounds for the periods and radii:
N_sim, cos_factor, P_min, P_max, radii_min, radii_max = read_targets_period_radius_bounds(loadfiles_directory + 'periods%s.out' % run_number)

# To read the simulation parameters from the file:
param_vals_all = read_sim_params(loadfiles_directory + 'periods%s.out' % run_number)

# To load and process the simulated observed catalog of stars and planets:
#star_obs = load_star_obs(loadfiles_directory + 'observed_catalog_stars%s.txt' % run_number)
#cat_obs = load_cat_obs(loadfiles_directory + 'observed_catalog_planets%s.txt' % run_number)
#sss_per_sys, sss = compute_summary_stats_from_cat_obs(cat_obs=cat_obs, star_obs=star_obs)

sss_per_sys, sss = compute_summary_stats_from_cat_obs(file_name_path=loadfiles_directory, run_number=run_number)

# To load and process the observed Kepler catalog and compare with our simulated catalog:
ssk_per_sys, ssk = compute_summary_stats_from_Kepler_catalog(P_min, P_max, radii_min, radii_max)

#dists_misc, dists_misc_w, dists_KS, dists_KS_w, dists_AD, dists_AD_w = compute_distances_sim_Kepler(loadfiles_directory + 'observed_catalog_planets%s.txt' % run_number, cat_obs, star_obs, P_min, P_max, radii_min, radii_max)
dists_misc, dists_misc_w, dists_KS, dists_KS_w, dists_AD, dists_AD_w = compute_distances_sim_Kepler(loadfiles_directory, None, None, P_min, P_max, radii_min, radii_max, AD_mod=AD_mod, dists_exclude=dists_exclude, run_number=run_number)





#'''
##### To plot the simulated and Kepler catalogs as marginal distributions:

subdirectory = 'Paper_Figures/Models/Observed/Clustered_P_R/' #'Paper_Figures/'; 'Talk_Figures/'

fig_size = (8,3) #size of each panel (figure)
fig_lbrt = [0.15, 0.3, 0.95, 0.925]

n_bins = 100
lw = 3 #linewidth
alpha = 0.2 #transparency of histograms

afs = 20 #axes labels font size
tfs = 20 #text labels font size
lfs = 16 #legend labels font size





##### To load and compute the same statistics for a large number of models, computing the confidence intervals for each bin:

loadfiles_directory = 'ACI/Simulated_Data/Julia_v0.7/Kepler_catalog_optimization/q1q17_dr25_gaia_fgk_stars79935/Clustered_P_R/f_high_incl_low_incl_mmr/Fit_rate_mult_P_Pratios_D_Dratios_dur_durratios_mmr/Some11_params_CRPDr_KS/Fixed_Rbreak3_Ncrit8/lc_0p2_5_lp_0p5_10_alphaP_-2_2_alphaR1_-4_2_alphaR2_-6_0_ecc_0_0p1_incl_inclmmr_0_90_sigmaR_0_0p5_sigmaP_0_0p3/targs79935_maxincl0_maxiters5000/sigma_i_greater_sigma_i_mmr/GP_best_models/'
runs = 100

Mtot_bins = np.arange(10)-0.5
Mtot_bins_mid = (Mtot_bins[:-1] + Mtot_bins[1:])/2.
Mtot_counts_all = []

P_bins = np.logspace(np.log10(P_min), np.log10(P_max), n_bins+1)
P_bins_mid = (P_bins[:-1] + P_bins[1:])/2.
P_counts_all = []

Rm_bins = np.logspace(np.log10(1.), np.log10(30.), n_bins+1)
Rm_bins_mid = (Rm_bins[:-1] + Rm_bins[1:])/2.
Rm_counts_all = []

tdur_bins = np.linspace(0., 1000., n_bins+1)
tdur_bins_mid = (tdur_bins[:-1] + tdur_bins[1:])/2.
tdur_counts_all = []

D_bins = np.logspace(-5., -1.5, n_bins+1)
D_bins_mid = (D_bins[:-1] + D_bins[1:])/2.
D_counts_all = []

radii_bins = np.linspace(radii_min, radii_max, n_bins+1)
radii_bins_mid = (radii_bins[:-1] + radii_bins[1:])/2.
radii_counts_all = []

Rstar_bins = np.linspace(0.5, 2.5, n_bins+1)
Rstar_bins_mid = (Rstar_bins[:-1] + Rstar_bins[1:])/2.
Rstar_counts_all = []

D_ratio_bins = np.logspace(-1.5, 1.5, n_bins+1)
D_ratio_bins_mid = (D_ratio_bins[:-1] + D_ratio_bins[1:])/2.
D_ratio_counts_all = []

xi_bins = np.linspace(-0.5, 0.5, n_bins+1)
xi_bins_mid = (xi_bins[:-1] + xi_bins[1:])/2.
xi_counts_all = []

xi_res_bins = np.linspace(-0.5, 0.5, n_bins+1)
xi_res_bins_mid = (xi_res_bins[:-1] + xi_res_bins[1:])/2.
xi_res_counts_all = []

xi_nonres_bins = np.linspace(-0.5, 0.5, n_bins+1)
xi_nonres_bins_mid = (xi_nonres_bins[:-1] + xi_nonres_bins[1:])/2.
xi_nonres_counts_all = []

for i in range(1,runs+1): #range(1,runs+1)
    run_number = i
    sss_per_sys_i, sss_i = compute_summary_stats_from_cat_obs(file_name_path=loadfiles_directory, run_number=run_number)
    dists_misc_i, dists_misc_w_i, dists_KS_i, dists_KS_w_i, dists_AD_i, dists_AD_w_i = compute_distances_sim_Kepler(loadfiles_directory, None, None, P_min, P_max, radii_min, radii_max, AD_mod=AD_mod, dists_exclude=dists_exclude, run_number=run_number)

    # Multiplicities:
    counts, bins = np.histogram(sss_i['Mtot_obs'], bins=Mtot_bins)
    Mtot_counts_all.append(counts/float(np.sum(counts)))
    
    # Periods:
    counts, bins = np.histogram(sss_i['P_obs'], bins=P_bins)
    P_counts_all.append(counts/float(np.sum(counts)))

    # Period ratios:
    counts, bins = np.histogram(sss_i['Rm_obs'], bins=Rm_bins)
    Rm_counts_all.append(counts/float(np.sum(counts)))

    # Durations:
    counts, bins = np.histogram(sss_i['tdur_obs'], bins=tdur_bins)
    tdur_counts_all.append(counts/float(np.sum(counts)))

    # Depths:
    counts, bins = np.histogram(sss_i['D_obs'], bins=D_bins)
    D_counts_all.append(counts/float(np.sum(counts)))

    # Planet radii:
    counts, bins = np.histogram(sss_i['radii_obs'], bins=radii_bins)
    radii_counts_all.append(counts/float(np.sum(counts)))

    # Stellar radii:
    counts, bins = np.histogram(sss_i['Rstar_obs'], bins=Rstar_bins)
    Rstar_counts_all.append(counts/float(np.sum(counts)))

    # Depth ratios:
    counts, bins = np.histogram(sss_i['D_ratio_obs'], bins=D_ratio_bins)
    D_ratio_counts_all.append(counts/float(np.sum(counts)))

    # Log(xi):
    counts, bins = np.histogram(np.log10(sss_i['xi_obs']), bins=xi_bins)
    xi_counts_all.append(counts/float(np.sum(counts)))

    # Log(xi) (res):
    counts, bins = np.histogram(np.log10(sss_i['xi_res_obs']), bins=xi_res_bins)
    xi_res_counts_all.append(counts/float(np.sum(counts)))

    # Log(xi) (non-res):
    counts, bins = np.histogram(np.log10(sss_i['xi_nonres_obs']), bins=xi_nonres_bins)
    xi_nonres_counts_all.append(counts/float(np.sum(counts)))

Mtot_counts_all = np.array(Mtot_counts_all)
P_counts_all = np.array(P_counts_all)
Rm_counts_all = np.array(Rm_counts_all)
tdur_counts_all = np.array(tdur_counts_all)
D_counts_all = np.array(D_counts_all)
radii_counts_all = np.array(radii_counts_all)
Rstar_counts_all = np.array(Rstar_counts_all)
D_ratio_counts_all = np.array(D_ratio_counts_all)
xi_counts_all = np.array(xi_counts_all)
xi_res_counts_all = np.array(xi_res_counts_all)
xi_nonres_counts_all = np.array(xi_nonres_counts_all)



Mtot_counts_16, Mtot_counts_84 = np.zeros(len(Mtot_bins_mid)), np.zeros(len(Mtot_bins_mid))
for b in range(len(Mtot_bins_mid)):
    counts_bin_sorted = np.sort(Mtot_counts_all[:,b])
    Mtot_counts_16[b], Mtot_counts_84[b] = counts_bin_sorted[16], counts_bin_sorted[84]

P_counts_16, P_counts_84 = np.zeros(n_bins), np.zeros(n_bins)
Rm_counts_16, Rm_counts_84 = np.zeros(n_bins), np.zeros(n_bins)
tdur_counts_16, tdur_counts_84 = np.zeros(n_bins), np.zeros(n_bins)
D_counts_16, D_counts_84 = np.zeros(n_bins), np.zeros(n_bins)
radii_counts_16, radii_counts_84 = np.zeros(n_bins), np.zeros(n_bins)
Rstar_counts_16, Rstar_counts_84 = np.zeros(n_bins), np.zeros(n_bins)
D_ratio_counts_16, D_ratio_counts_84 = np.zeros(n_bins), np.zeros(n_bins)
xi_counts_16, xi_counts_84 = np.zeros(n_bins), np.zeros(n_bins)
xi_res_counts_16, xi_res_counts_84 = np.zeros(n_bins), np.zeros(n_bins)
xi_nonres_counts_16, xi_nonres_counts_84 = np.zeros(n_bins), np.zeros(n_bins)
for b in range(n_bins):
    # Periods:
    counts_bin_sorted = np.sort(P_counts_all[:,b])
    P_counts_16[b], P_counts_84[b] = counts_bin_sorted[16], counts_bin_sorted[84]

    # Period ratios:
    counts_bin_sorted = np.sort(Rm_counts_all[:,b])
    Rm_counts_16[b], Rm_counts_84[b] = counts_bin_sorted[16], counts_bin_sorted[84]

    # Durations:
    counts_bin_sorted = np.sort(tdur_counts_all[:,b])
    tdur_counts_16[b], tdur_counts_84[b] = counts_bin_sorted[16], counts_bin_sorted[84]

    # Depths:
    counts_bin_sorted = np.sort(D_counts_all[:,b])
    D_counts_16[b], D_counts_84[b] = counts_bin_sorted[16], counts_bin_sorted[84]
    
    # Planet radii:
    counts_bin_sorted = np.sort(radii_counts_all[:,b])
    radii_counts_16[b], radii_counts_84[b] = counts_bin_sorted[16], counts_bin_sorted[84]
    
    # Stellar radii:
    counts_bin_sorted = np.sort(Rstar_counts_all[:,b])
    Rstar_counts_16[b], Rstar_counts_84[b] = counts_bin_sorted[16], counts_bin_sorted[84]
    
    # Depth ratios:
    counts_bin_sorted = np.sort(D_ratio_counts_all[:,b])
    D_ratio_counts_16[b], D_ratio_counts_84[b] = counts_bin_sorted[16], counts_bin_sorted[84]
    
    # Log(xi):
    counts_bin_sorted = np.sort(xi_counts_all[:,b])
    xi_counts_16[b], xi_counts_84[b] = counts_bin_sorted[16], counts_bin_sorted[84]
    
    # Log(xi) (res):
    counts_bin_sorted = np.sort(xi_res_counts_all[:,b])
    xi_res_counts_16[b], xi_res_counts_84[b] = counts_bin_sorted[16], counts_bin_sorted[84]
    
    # Log(xi) (non-res):
    counts_bin_sorted = np.sort(xi_nonres_counts_all[:,b])
    xi_nonres_counts_16[b], xi_nonres_counts_84[b] = counts_bin_sorted[16], counts_bin_sorted[84]

#####





# To make a 'plot' listing the model parameters:
fig = plt.figure(figsize=fig_size)
plot = GridSpec(1,1,left=fig_lbrt[0],bottom=fig_lbrt[1],right=fig_lbrt[2],top=fig_lbrt[3],wspace=0.1,hspace=0.1)
nrows = 8
for i in range(len(param_keys_all)): #range(len(param_keys_all))
    plt.figtext(x=0.05+0.3*int(i/float(nrows)), y=0.875-0.1*(i%nrows), s=r'%s = %s' % (param_keys_all[i][1], param_vals_all[i]), fontsize=lfs)
if savefigures == True:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_sim_params.pdf')
    plt.close()

# Multiplicities:
ax = setup_fig_single(fig_size, fig_lbrt[0], fig_lbrt[1], fig_lbrt[2], fig_lbrt[3])
x = sss['Mtot_obs'][sss['Mtot_obs'] > 0]
max_M = np.max((np.max(sss['Mtot_obs']), np.max(ssk['Mtot_obs'])))
counts, bins = np.histogram(x, bins=max_M+1, range=(-0.5, max_M+0.5))
bins_mid = (bins[:-1] + bins[1:])/2.
plt.plot(bins_mid, counts/float(np.sum(counts)), drawstyle='steps-mid', color='k', lw=lw, label='Simulated') #label='%s simulated observed systems' % len(x)
plt.plot(Mtot_bins_mid, Mtot_counts_16, drawstyle='steps-mid', color='r', lw=1, ls='--', label=r'16% and 84%')
plt.plot(Mtot_bins_mid, Mtot_counts_84, drawstyle='steps-mid', color='r', lw=1, ls='--')
counts, bins = np.histogram(ssk['Mtot_obs'], bins=bins)
plt.scatter(bins_mid, counts/float(np.sum(counts)), marker='x', s=50, color='k', label='Kepler') #label='%s Kepler systems' % len(ssk['Mtot_obs'])
plt.gca().set_yscale("log")
ax.tick_params(axis='both', labelsize=afs)
plt.xlim([0.5, max_M+0.5])
plt.ylim([1e-4, 1.])
plt.xlabel('Number of planets', fontsize=tfs)
plt.ylabel('Fraction', fontsize=tfs)
plt.legend(loc='lower left', bbox_to_anchor=(0.01,0.01), ncol=1, frameon=False, fontsize=lfs) #show the legend
if savefigures:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_multiplicities_compare.pdf')
    plt.close()

# Periods:
plot_fig_pdf_simple(fig_size, [sss['P_obs']], [ssk['P_obs']], x_min=P_min, x_max=P_max, y_min=1e-3, y_max=0.1, n_bins=n_bins, log_x=True, log_y=True, lw=lw, xticks_custom=[3,10,30,100,300], xlabel_text=r'$P$ (days)', afs=afs, tfs=tfs, lfs=lfs, legend=True, fig_lbrt=fig_lbrt)
plt.plot(P_bins_mid, P_counts_16, drawstyle='steps-mid', color='r', lw=1, ls='--', label='16')
plt.plot(P_bins_mid, P_counts_84, drawstyle='steps-mid', color='r', lw=1, ls='--', label='84')
if savefigures:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_periods_compare.pdf')
    plt.close()

# Period ratios (all, with some upper cut-off):
R_max_cut = 30. #upper cut-off for plotting period ratios; np.max(sss['Rm_obs'])
plot_fig_pdf_simple(fig_size, [sss['Rm_obs'][sss['Rm_obs'] < R_max_cut]], [ssk['Rm_obs'][ssk['Rm_obs'] < R_max_cut]], x_min=1., x_max=R_max_cut, n_bins=n_bins, log_x=True, lw=lw, xticks_custom=[1,2,3,4,5,10,20], xlabel_text=r'$P_{i+1}/P_i$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt)
plt.plot(Rm_bins_mid, Rm_counts_16, drawstyle='steps-mid', color='r', lw=1, ls='--', label='16')
plt.plot(Rm_bins_mid, Rm_counts_84, drawstyle='steps-mid', color='r', lw=1, ls='--', label='84')
if savefigures:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_periodratios_compare.pdf')
    plt.close()

# Transit durations:
plot_fig_pdf_simple(fig_size, [sss['tdur_obs']], [ssk['tdur_obs']*60.], x_min=0., x_max=1000., n_bins=n_bins, lw=lw, xlabel_text=r'$t_{\rm dur}$ (mins)', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt)
plt.plot(tdur_bins_mid, tdur_counts_16, drawstyle='steps-mid', color='r', lw=1, ls='--', label='16')
plt.plot(tdur_bins_mid, tdur_counts_84, drawstyle='steps-mid', color='r', lw=1, ls='--', label='84')
if savefigures:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_durations_compare.pdf')
    plt.close()

# Transit depths:
plot_fig_pdf_simple(fig_size, [sss['D_obs']], [ssk['D_obs']], x_min=np.min(D_bins), x_max=np.max(D_bins), log_x=True, lw=lw, xlabel_text=r'$\delta$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt)
plt.plot(D_bins_mid, D_counts_16, drawstyle='steps-mid', color='r', lw=1, ls='--', label='16')
plt.plot(D_bins_mid, D_counts_84, drawstyle='steps-mid', color='r', lw=1, ls='--', label='84')
if savefigures:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_depths_compare.pdf')
    plt.close()

# Planet radii:
plot_fig_pdf_simple(fig_size, [sss['radii_obs']], [ssk['radii_obs']], x_min=radii_min, x_max=radii_max, n_bins=n_bins, lw=lw, xlabel_text=r'$R_p (R_\oplus)$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt)
plt.plot(radii_bins_mid, radii_counts_16, drawstyle='steps-mid', color='r', lw=1, ls='--', label='16')
plt.plot(radii_bins_mid, radii_counts_84, drawstyle='steps-mid', color='r', lw=1, ls='--', label='84')
if savefigures:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_radii_compare.pdf')
    plt.close()

# Stellar radii:
plot_fig_pdf_simple(fig_size, [sss['Rstar_obs']], [ssk['Rstar_obs']], x_min=0.5, x_max=2.5, n_bins=n_bins, lw=lw, xlabel_text=r'$R_\star (R_\odot)$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt)
plt.plot(Rstar_bins_mid, Rstar_counts_16, drawstyle='steps-mid', color='r', lw=1, ls='--', label='16')
plt.plot(Rstar_bins_mid, Rstar_counts_84, drawstyle='steps-mid', color='r', lw=1, ls='--', label='84')
if savefigures:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_stellar_radii_compare.pdf')
    plt.close()

# Transit depth ratios:
plot_fig_pdf_simple(fig_size, [sss['D_ratio_obs']], [ssk['D_ratio_obs']], x_min=np.min(D_ratio_bins), x_max=np.max(D_ratio_bins), n_bins=n_bins, log_x=True, lw=lw, xlabel_text=r'$\delta_{i+1}/\delta_i$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt)
plt.plot(D_ratio_bins_mid, D_ratio_counts_16, drawstyle='steps-mid', color='r', lw=1, ls='--', label='16')
plt.plot(D_ratio_bins_mid, D_ratio_counts_84, drawstyle='steps-mid', color='r', lw=1, ls='--', label='84')
if savefigures:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_depthratios_compare.pdf')
    plt.close()

# Log(xi):
plot_fig_pdf_simple(fig_size, [np.log10(sss['xi_obs'])], [np.log10(ssk['xi_obs'])], x_min=-0.5, x_max=0.5, n_bins=n_bins, lw=lw, xlabel_text=r'$\log{\xi}$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt)
plt.plot(xi_bins_mid, xi_counts_16, drawstyle='steps-mid', color='r', lw=1, ls='--', label='16')
plt.plot(xi_bins_mid, xi_counts_84, drawstyle='steps-mid', color='r', lw=1, ls='--', label='84')
if savefigures:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_logxi_all_compare.pdf')
    plt.close()

# Log(xi) by res/non-res:
plot_fig_pdf_simple(fig_size, [np.log10(sss['xi_res_obs']), np.log10(sss['xi_nonres_obs'])], [np.log10(ssk['xi_res_obs']), np.log10(ssk['xi_nonres_obs'])], x_min=-0.5, x_max=0.5, n_bins=n_bins, c_sim=['m','g'], c_Kep=['m','g'], ls_sim=['-','-'], ls_Kep=['-','-'], lw=lw, labels_sim=['Near MMR', 'Not near MMR'], labels_Kep=[None, None], xlabel_text=r'$\log{\xi}$', afs=afs, tfs=tfs, lfs=lfs, legend=True, fig_lbrt=fig_lbrt)
plt.plot(xi_res_bins_mid, xi_res_counts_16, drawstyle='steps-mid', color='m', lw=1, ls='--', label='16')
plt.plot(xi_res_bins_mid, xi_res_counts_84, drawstyle='steps-mid', color='m', lw=1, ls='--', label='84')
plt.plot(xi_nonres_bins_mid, xi_nonres_counts_16, drawstyle='steps-mid', color='g', lw=1, ls='--', label='16')
plt.plot(xi_nonres_bins_mid, xi_nonres_counts_84, drawstyle='steps-mid', color='g', lw=1, ls='--', label='84')
if savefigures:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_logxi_compare.pdf')
    plt.close()

plt.show()
plt.close()
