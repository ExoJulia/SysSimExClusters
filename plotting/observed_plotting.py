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



'''
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
counts, bins = np.histogram(ssk['Mtot_obs'], bins=bins)
plt.scatter(bins_mid, counts/float(np.sum(counts)), marker='x', s=50, color='k', label='Kepler') #label='%s Kepler systems' % len(ssk['Mtot_obs'])
plt.gca().set_yscale("log")
ax.tick_params(axis='both', labelsize=afs)
plt.xlim([0.5, max_M+0.5])
plt.xlabel('Observed planets per system', fontsize=tfs)
plt.ylabel('Fraction', fontsize=tfs)
plt.legend(loc='upper right', bbox_to_anchor=(0.99,0.99), ncol=1, frameon=False, fontsize=lfs) #show the legend
if savefigures:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_multiplicities_compare.pdf')
    plt.close()

# Periods:
plot_fig_pdf_simple(fig_size, [sss['P_obs']], [ssk['P_obs']], x_min=3., x_max=300., y_min=1e-3, y_max=0.1, n_bins=n_bins, log_x=True, log_y=True, lw=lw, xticks_custom=[3,10,30,100,300], xlabel_text=r'$P$ (days)', afs=afs, tfs=tfs, lfs=lfs, legend=True, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + model_name + '_periods_compare.pdf', save_fig=savefigures)

# Period ratios (all, with some upper cut-off):
R_max_cut = 30. #upper cut-off for plotting period ratios; np.max(sss['Rm_obs'])
plot_fig_pdf_simple(fig_size, [sss['Rm_obs'][sss['Rm_obs'] < R_max_cut]], [ssk['Rm_obs'][ssk['Rm_obs'] < R_max_cut]], x_min=1., x_max=R_max_cut, n_bins=n_bins, log_x=True, lw=lw, xticks_custom=[1,2,3,4,5,10,20], xlabel_text=r'$P_{i+1}/P_i$', afs=afs, tfs=tfs, lfs=lfs, legend=True, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + model_name + '_periodratios_compare.pdf', save_fig=savefigures)

# Period ratios (< 5):
plot_fig_pdf_simple(fig_size, [sss['Rm_obs'][sss['Rm_obs'] < 5.]], [ssk['Rm_obs'][ssk['Rm_obs'] < 5.]], x_min=1., x_max=5., n_bins=n_bins, lw=lw, xlabel_text=r'$P_{i+1}/P_i$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + model_name + '_periodratios_less5_compare.pdf', save_fig=savefigures)

# Transit durations:
plot_fig_pdf_simple(fig_size, [sss['tdur_obs']], [ssk['tdur_obs']*60.], x_max=1000., n_bins=n_bins, lw=lw, xlabel_text=r'$t_{\rm dur}$ (mins)', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + model_name + '_durations_compare.pdf', save_fig=savefigures)

# Transit depths:
plot_fig_pdf_simple(fig_size, [sss['D_obs']], [ssk['D_obs']], log_x=True, lw=lw, xlabel_text=r'$\delta$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + model_name + '_depths_compare.pdf', save_fig=savefigures)

# Transit depths (above and below the photoevaporation boundary):
plot_fig_pdf_simple(fig_size, [sss['D_above_obs'], sss['D_below_obs']], [ssk['D_above_obs'], ssk['D_below_obs']], n_bins=n_bins, log_x=True, c_sim=['b','r'], c_Kep=['b','r'], ls_sim=['-','-'], ls_Kep=['-','-'], lw=lw, labels_sim=['Above', 'Below'], labels_Kep=[None, None], xlabel_text=r'$\delta$', afs=afs, tfs=tfs, lfs=lfs, legend=True, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + model_name + '_depths_photoevap_compare.pdf', save_fig=savefigures)

# Planet radii:
plot_fig_pdf_simple(fig_size, [sss['radii_obs']], [ssk['radii_obs']], n_bins=n_bins, lw=lw, xlabel_text=r'$R_p (R_\oplus)$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + model_name + '_radii_compare.pdf', save_fig=savefigures)

# Stellar radii:
plot_fig_pdf_simple(fig_size, [sss['Rstar_obs']], [ssk['Rstar_obs']], n_bins=n_bins, lw=lw, xlabel_text=r'$R_\star (R_\odot)$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + model_name + '_stellar_radii_compare.pdf', save_fig=savefigures)

# Transit depth ratios:
plot_fig_pdf_simple(fig_size, [sss['D_ratio_obs']], [ssk['D_ratio_obs']], n_bins=n_bins, log_x=True, lw=lw, xlabel_text=r'$\delta_{i+1}/\delta_i$', afs=afs, tfs=tfs, lfs=lfs, legend=True, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + model_name + '_depthratios_compare.pdf', save_fig=savefigures)

# Transit depth ratios (above, below, and across the photoevaporation boundary):
plot_fig_pdf_simple(fig_size, [sss['D_ratio_above_obs'], sss['D_ratio_below_obs'], sss['D_ratio_across_obs']], [ssk['D_ratio_above_obs'], ssk['D_ratio_below_obs'], ssk['D_ratio_across_obs']], n_bins=n_bins, log_x=True, c_sim=['b','r','k'], c_Kep=['b','r','k'], ls_sim=['-','-','-'], ls_Kep=['-','-','-'], lw=lw, labels_sim=['Above', 'Below', 'Across'], labels_Kep=[None, None, None], xlabel_text=r'$\delta_{i+1}/\delta_i$', afs=afs, tfs=tfs, lfs=lfs, legend=True, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + model_name + '_depthratios_photoevap_compare.pdf', save_fig=savefigures)

# Log(xi):
plot_fig_pdf_simple(fig_size, [np.log10(sss['xi_obs'])], [np.log10(ssk['xi_obs'])], x_min=-0.5, x_max=0.5, n_bins=n_bins, lw=lw, xlabel_text=r'$\log{\/xi}$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + model_name + '_logxi_all_compare.pdf', save_fig=savefigures)

# Log(xi) by res/non-res:
plot_fig_pdf_simple(fig_size, [np.log10(sss['xi_res_obs']), np.log10(sss['xi_nonres_obs'])], [np.log10(ssk['xi_res_obs']), np.log10(ssk['xi_nonres_obs'])], x_min=-0.5, x_max=0.5, n_bins=n_bins, c_sim=['m','g'], c_Kep=['m','g'], ls_sim=['-','-'], ls_Kep=['-','-'], lw=lw, labels_sim=['Near MMR', 'Not near MMR'], labels_Kep=[None, None], xlabel_text=r'$\log{\/xi}$', afs=afs, tfs=tfs, lfs=lfs, legend=True, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + model_name + '_logxi_compare.pdf', save_fig=savefigures)

# Log(xi) within res:
plot_fig_pdf_simple(fig_size, [np.log10(sss['xi_res32_obs']), np.log10(sss['xi_res21_obs'])], [np.log10(ssk['xi_res32_obs']), np.log10(ssk['xi_res21_obs'])], x_min=-0.5, x_max=0.5, n_bins=n_bins, c_sim=['r','b'], c_Kep=['r','b'], ls_sim=['-','-'], ls_Kep=['-','-'], lw=lw, labels_sim=['Near 3:2 MMR', 'Near 2:1 MMR'], labels_Kep=[None, None], xlabel_text=r'$\log{\/xi}$', afs=afs, tfs=tfs, lfs=lfs, legend=True, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + model_name + '_logxi_res_compare.pdf', save_fig=savefigures)

plt.show()
plt.close()
'''


# Multiplicity CDFs:
ax = setup_fig_single(fig_size, fig_lbrt[0], fig_lbrt[1], fig_lbrt[2], fig_lbrt[3])
x = sss['Mtot_obs'][sss['Mtot_obs'] > 0]
counts_cumu = np.array([sum(x <= xi) for xi in range(1, np.max(x)+1)])
plt.plot(range(1, np.max(x)+1), counts_cumu/np.float(len(x)), drawstyle='steps-post', color='k', lw=lw, label='Simulated')
counts_cumu = np.array([sum(ssk['Mtot_obs'] <= xi) for xi in range(1, np.max(ssk['Mtot_obs'])+1)])
plt.plot(range(1, np.max(ssk['Mtot_obs'])+1), counts_cumu/np.float(len(ssk['Mtot_obs'])), drawstyle='steps-post', color='k', ls='--', label='Kepler')
plt.annotate(s='', xy=(dists_KS['M_pos'] + 0.5, sum(x <= dists_KS['M_pos'])/np.float(len(x))), xytext=(dists_KS['M_pos'] + 0.5, sum(ssk['Mtot_obs'] <= dists_KS['M_pos'])/np.float(len(ssk['Mtot_obs']))), arrowprops=dict(arrowstyle='<->'))
plt.figtext(0.925, 0.45, r'$D_f = %s$' % np.round(dists_misc['delta_f'], 4), ha='right', fontsize=lfs)
#plt.figtext(0.925, 0.35, r'$\mathcal{D}_{\rm KS} = %s$' % np.round(dists_KS['M'], 3), ha='right', fontsize=lfs)
plt.figtext(0.925, 0.35, r'$\rho_{\rm CRPD} = %s$' % np.round(dists_misc['d_mult_CRPD_switched'], 3), ha='right', fontsize=lfs)
ax.tick_params(axis='both', labelsize=afs)
plt.xlim([1, np.max((np.max(x), np.max(ssk['Mtot_obs'])))+1])
plt.ylim([0.6,1])
plt.xlabel(r'Number of planets', fontsize=tfs)
plt.ylabel('CDF', fontsize=tfs)
plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, frameon=False, fontsize=lfs) #show the legend
if savefigures:
    plt.savefig(savefigures_directory + subdirectory + model_name + '_multiplicities_CDFs.pdf')
    plt.close()

# Periods CDFs:
plot_fig_cdf_simple(fig_size, [sss['P_obs']], [ssk['P_obs']], x_min=P_min, x_max=P_max, log_x=True, lw=lw, xticks_custom=[3,10,30,100,300], xlabel_text=r'$P$ (days)', afs=afs, tfs=tfs, lfs=lfs, label_dist=True, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + model_name + '_periods_CDFs.pdf', save_fig=savefigures)

# Period ratios CDFs:
plot_fig_cdf_simple(fig_size, [sss['Rm_obs']], [ssk['Rm_obs']], log_x=True, lw=lw, xticks_custom=[1,2,3,4,5,10,20,40], xlabel_text=r'$P_{i+1}/P_i$', afs=afs, tfs=tfs, lfs=lfs, label_dist=True, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + model_name + '_periodratios_CDFs.pdf', save_fig=savefigures)

# Transit durations CDFs:
plot_fig_cdf_simple(fig_size, [sss['tdur_obs']], [ssk['tdur_obs']*60.], x_min=0., x_max=1000., lw=lw, xlabel_text=r'$t_{\rm dur}$ (mins)', afs=afs, tfs=tfs, lfs=lfs, label_dist=True, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + model_name + '_durations_CDFs.pdf', save_fig=savefigures)

# Transit depths CDFs:
plot_fig_cdf_simple(fig_size, [sss['D_obs']], [ssk['D_obs']], log_x=True, lw=lw, xlabel_text=r'$\delta$', afs=afs, tfs=tfs, lfs=lfs, label_dist=True, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + model_name + '_depths_CDFs.pdf', save_fig=savefigures)

# Transit depths CDFs (above and below the photoevaporation boundary):
plot_fig_cdf_simple(fig_size, [sss['D_above_obs'], sss['D_below_obs']], [ssk['D_above_obs'], ssk['D_below_obs']], log_x=True, c_sim=['b','r'], c_Kep=['b','r'], ls_sim=['-','-'], ls_Kep=['--','--'], lw=lw, labels_sim=['Above', 'Below'], labels_Kep=[None, None], xlabel_text=r'$\delta$', afs=afs, tfs=tfs, lfs=lfs, legend=True, label_dist=True, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + model_name + '_depths_photoevap_CDFs.pdf', save_fig=savefigures)

# Planet radii CDFs:
plot_fig_cdf_simple(fig_size, [sss['radii_obs']], [ssk['radii_obs']], x_min=radii_min, x_max=radii_max, lw=lw, xlabel_text=r'$R_p (R_\oplus)$', afs=afs, tfs=tfs, lfs=lfs, label_dist=True, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + model_name + '_radii_CDFs.pdf', save_fig=savefigures)

# Transit depth ratios CDFs:
plot_fig_cdf_simple(fig_size, [sss['D_ratio_obs']], [ssk['D_ratio_obs']], x_min=10.**-1.5, x_max=10.**1.5, log_x=True, lw=lw, xlabel_text=r'$\delta_{i+1}/\delta_i$', afs=afs, tfs=tfs, lfs=lfs, label_dist=True, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + model_name + '_depthratios_CDFs.pdf', save_fig=savefigures)

# Transit depth ratios CDFs (above, below, and across the photoevaporation boundary):
plot_fig_cdf_simple(fig_size, [sss['D_ratio_above_obs'], sss['D_ratio_below_obs'], sss['D_ratio_across_obs']], [ssk['D_ratio_above_obs'], ssk['D_ratio_below_obs'], ssk['D_ratio_across_obs']], log_x=True, c_sim=['b','r','k'], c_Kep=['b','r','k'], ls_sim=['-','-','-'], ls_Kep=['--','--','--'], lw=lw, labels_sim=['Above', 'Below', 'Across'], labels_Kep=[None, None, None], xlabel_text=r'$\delta_{i+1}/\delta_i$', afs=afs, tfs=tfs, lfs=lfs, legend=True, label_dist=True, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + model_name + '_depthratios_photoevap_CDFs.pdf', save_fig=savefigures)

# Log(xi) CDFs:
plot_fig_cdf_simple(fig_size, [np.log10(sss['xi_obs'])], [np.log10(ssk['xi_obs'])], x_min=-0.5, x_max=0.5, lw=lw, xlabel_text=r'$\log{\xi}$', afs=afs, tfs=tfs, lfs=lfs, label_dist=True, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + model_name + '_logxi_all_CDFs.pdf', save_fig=savefigures)

# Log(xi) CDFs by res/non-res:
plot_fig_cdf_simple(fig_size, [np.log10(sss['xi_res_obs']), np.log10(sss['xi_nonres_obs'])], [np.log10(ssk['xi_res_obs']), np.log10(ssk['xi_nonres_obs'])], x_min=-0.5, x_max=0.5, c_sim=['m','g'], c_Kep=['m','g'], ls_sim=['-','-'], ls_Kep=['--','--'], lw=lw, labels_sim=['Near MMR', 'Not near MMR'], labels_Kep=[None, None], xlabel_text=r'$\log{\xi}$', afs=afs, tfs=tfs, lfs=lfs, legend=True, label_dist=True, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + model_name + '_logxi_CDFs.pdf', save_fig=savefigures)

# Log(xi) CDFs within res:
plot_fig_cdf_simple(fig_size, [np.log10(sss['xi_res32_obs']), np.log10(sss['xi_res21_obs'])], [np.log10(ssk['xi_res32_obs']), np.log10(ssk['xi_res21_obs'])], x_min=-0.5, x_max=0.5, c_sim=['r','b'], c_Kep=['r','b'], ls_sim=['-','-'], ls_Kep=['--','--'], lw=lw, labels_sim=['Near 3:2 MMR', 'Near 2:1 MMR'], labels_Kep=[None, None], xlabel_text=r'$\log{\xi}$', afs=afs, tfs=tfs, lfs=lfs, legend=True, label_dist=True, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + model_name + '_logxi_res_CDFs.pdf', save_fig=savefigures)

plt.show()
plt.close()
#'''





##### To plot the observed multi-systems by period to visualize the systems (similar to Fig 1 in Fabrycky et al. 2014):

#load_cat_obs_and_plot_figs_multis_gallery(loadfiles_directory, run_number=run_number, save_name_base=savefigures_directory + subdirectory + model_name + '_multis', save_fig=savefigures)
plt.show()
plt.close()
