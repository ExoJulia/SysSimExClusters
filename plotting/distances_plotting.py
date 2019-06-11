from ExoplanetsSysSim_functions import *





savefigures = False

savefigures_directory = 'ExoplanetsSysSim_Clusters/Figures/'
save_name = 'Distances_Compare'

AD_mod = 'true' # 'true' or 'false'
dists_exclude = np.array([2,3,8,12,13,15,16,17]) - 1





##### To load the files with the systems with observed planets:

runs = 100

# Model 1:
loadfiles_directory1 = 'ACI/Simulated_Data/Julia_v0.7/Kepler_catalog_optimization/q1q17_dr25_gaia_fgk_stars79935/Clustered_P_R/f_high_incl_low_incl_mmr/Fit_rate_mult_P_Pratios_D_Dratios_dur_durratios_mmr/Some11_params_CRPDr_KS/Fixed_Rbreak3_Ncrit8/lc_0p2_5_lp_0p5_10_alphaP_-2_2_alphaR1_-4_2_alphaR2_-6_0_ecc_0_0p1_incl_inclmmr_0_90_sigmaR_0_0p5_sigmaP_0_0p3/targs79935_maxincl0_maxiters5000/sigma_i_greater_sigma_i_mmr/GP_best_models/'

# To first read the number of simulated targets and bounds for the periods and radii:
N_sim, cos_factor, P_min, P_max, radii_min, radii_max = read_targets_period_radius_bounds(loadfiles_directory1 + 'periods1.out')

# Model 2:
loadfiles_directory2 = 'ACI/Simulated_Data/Julia_v0.7/Kepler_catalog_optimization/q1q17_dr25_gaia_fgk_stars79935/Clustered_P/f_high_incl_low_incl_mmr/Fit_rate_mult_P_Pratios_D_Dratios_dur_durratios_mmr/Some10_params_CRPDr_KS/Fixed_Rbreak3_Ncrit8/lc_0p2_5_lp_0p5_10_alphaP_-2_2_alphaR1_-4_2_alphaR2_-6_0_ecc_0_0p1_incl_inclmmr_0_90_sigmaP_0_0p3/targs79935_maxincl0_maxiters5000/sigma_i_greater_sigma_i_mmr/GP_best_models/'

# Model 3:
loadfiles_directory3 = 'ACI/Simulated_Data/Julia_v0.7/Kepler_catalog_optimization/q1q17_dr25_gaia_fgk_stars79935/Non_Clustered/f_high_incl_low_incl_mmr/Fit_rate_mult_P_Pratios_D_Dratios_dur_durratios_mmr/Some8_params_CRPDr_KS/Fixed_Rbreak3_Ncrit8/lc_1_8_alphaP_-2_2_alphaR1_-4_2_alphaR2_-6_0_ecc_0_0p1_incl_inclmmr_0_90/targs79935_maxincl0_maxiters5000/sigma_i_greater_sigma_i_mmr/GP_best_models/'



model_dirs = [loadfiles_directory1, loadfiles_directory2, loadfiles_directory3]
model_names = ['Clustered P+R', 'Clustered P', 'Non-clustered'] # Make sure this matches the models loaded!
model_linestyles = ['-', '-', '-']
model_colors = ['b', 'g', 'r']
n_models = len(model_names)

# To load and process the observed Kepler catalog and compare with our simulated catalog:
ssk_per_sys, ssk = compute_summary_stats_from_Kepler_catalog(P_min, P_max, radii_min, radii_max)





##### To load and compute the distances for a large number of models:

dist_tot_w_KS_all = np.zeros((runs,n_models))
dist_tot_w_AD_all = np.zeros((runs,n_models))
dist_delta_f_all = np.zeros((runs,n_models))
dist_CRPDr_all = np.zeros((runs,n_models))

dist_P_KS_all = np.zeros((runs,n_models))
dist_Rm_KS_all = np.zeros((runs,n_models))
dist_tdur_KS_all = np.zeros((runs,n_models))
dist_logxi_res_KS_all = np.zeros((runs,n_models))
dist_logxi_nonres_KS_all = np.zeros((runs,n_models))
dist_D_KS_all = np.zeros((runs,n_models))
dist_D_ratio_KS_all = np.zeros((runs,n_models))

dist_P_AD_all = np.zeros((runs,n_models))
dist_Rm_AD_all = np.zeros((runs,n_models))
dist_tdur_AD_all = np.zeros((runs,n_models))
dist_logxi_res_AD_all = np.zeros((runs,n_models))
dist_logxi_nonres_AD_all = np.zeros((runs,n_models))
dist_D_AD_all = np.zeros((runs,n_models))
dist_D_ratio_AD_all = np.zeros((runs,n_models))

for m,load_dir in enumerate(model_dirs):
    for i in range(runs): #range(1,runs+1)
        run_number = i+1
        sss_per_sys_i, sss_i = compute_summary_stats_from_cat_obs(file_name_path=load_dir, run_number=run_number)
        dists_misc_i, dists_misc_w_i, dists_KS_i, dists_KS_w_i, dists_AD_i, dists_AD_w_i = compute_distances_sim_Kepler(load_dir, None, None, P_min, P_max, radii_min, radii_max, AD_mod=AD_mod, dists_exclude=dists_exclude, run_number=run_number)

        dist_tot_w_KS_all[i,m] = dists_KS_i['tot_weighted']
        dist_tot_w_AD_all[i,m] = dists_AD_i['tot_weighted']
        dist_delta_f_all[i,m] = dists_misc_i['delta_f']
        dist_CRPDr_all[i,m] = dists_misc_i['d_mult_CRPD_switched']
        
        dist_P_KS_all[i,m] = dists_KS_i['P']
        dist_Rm_KS_all[i,m] = dists_KS_i['Rm']
        dist_tdur_KS_all[i,m] = dists_KS_i['tdur']
        dist_logxi_res_KS_all[i,m] = dists_KS_i['logxi_res']
        dist_logxi_nonres_KS_all[i,m] = dists_KS_i['logxi_nonres']
        dist_D_KS_all[i,m] = dists_KS_i['D']
        dist_D_ratio_KS_all[i,m] = dists_KS_i['D_ratio']

        dist_P_AD_all[i,m] = dists_AD_i['P']
        dist_Rm_AD_all[i,m] = dists_AD_i['Rm']
        dist_tdur_AD_all[i,m] = dists_AD_i['tdur']
        dist_logxi_res_AD_all[i,m] = dists_AD_i['logxi_res']
        dist_logxi_nonres_AD_all[i,m] = dists_AD_i['logxi_nonres']
        dist_D_AD_all[i,m] = dists_AD_i['D']
        dist_D_ratio_AD_all[i,m] = dists_AD_i['D_ratio']





##### To plot histograms of the individual distances:

subdirectory = 'Paper_Figures/' #'Paper_Figures/'; 'Talk_Figures/'

fig_size = (8,3) #size of each panel (figure)
fig_lbrt = [0.15, 0.3, 0.95, 0.925]

n_bins = 20
lw = 3 #linewidth
alpha = 0.2 #transparency of histograms

afs = 20 #axes labels font size
tfs = 20 #text labels font size
lfs = 16 #legend labels font size

# Total weighted distance (KS):
plot_fig_pdf_simple(fig_size, [dist_tot_w_KS_all[:,m] for m in range(n_models)], [], n_bins=n_bins, log_x=False, c_sim=model_colors, ls_sim=model_linestyles, lw=lw, labels_sim=model_names, xlabel_text=r'$\mathcal{D}_W (\rm KS)$', afs=afs, tfs=tfs, lfs=lfs, legend=True, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + save_name + '_tot_weighted_KS.pdf', save_fig=savefigures)

# Total weighted distance (AD):
plot_fig_pdf_simple(fig_size, [dist_tot_w_AD_all[:,m] for m in range(n_models)], [], n_bins=n_bins, log_x=False, c_sim=model_colors, ls_sim=model_linestyles, lw=lw, labels_sim=model_names, xlabel_text=r'$\mathcal{D}_W (\rm AD^\prime)$', afs=afs, tfs=tfs, lfs=lfs, legend=True, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + save_name + '_tot_weighted_AD.pdf', save_fig=savefigures)

# Rate of planets:
plot_fig_pdf_simple(fig_size, [dist_delta_f_all[:,m] for m in range(n_models)], [], n_bins=n_bins, log_x=False, c_sim=model_colors, ls_sim=model_linestyles, lw=lw, labels_sim=model_names, xlabel_text=r'$D_f$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + save_name + '_delta_f.pdf', save_fig=savefigures)

# Multiplicity CRPDr:
plot_fig_pdf_simple(fig_size, [dist_CRPDr_all[:,m] for m in range(n_models)], [], n_bins=n_bins, log_x=False, c_sim=model_colors, ls_sim=model_linestyles, lw=lw, labels_sim=model_names, xlabel_text=r'$\rho_{\rm CRPD}$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + save_name + '_CRPDr.pdf', save_fig=savefigures)

#'''
# Periods (KS):
plot_fig_pdf_simple(fig_size, [dist_P_KS_all[:,m] for m in range(n_models)], [], n_bins=n_bins, log_x=False, c_sim=model_colors, ls_sim=model_linestyles, lw=lw, labels_sim=model_names, xlabel_text=r'$\mathcal{D}_{\rm KS}$ for $\{P\}$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + save_name + '_periods_KS.pdf', save_fig=savefigures)

# Period ratios (KS):
plot_fig_pdf_simple(fig_size, [dist_Rm_KS_all[:,m] for m in range(n_models)], [], n_bins=n_bins, log_x=False, c_sim=model_colors, ls_sim=model_linestyles, lw=lw, labels_sim=model_names, xlabel_text=r'$\mathcal{D}_{\rm KS}$ for $\{\mathcal{P}\}$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + save_name + '_periodratios_KS.pdf', save_fig=savefigures)

# Durations (KS):
plot_fig_pdf_simple(fig_size, [dist_tdur_KS_all[:,m] for m in range(n_models)], [], n_bins=n_bins, log_x=False, c_sim=model_colors, ls_sim=model_linestyles, lw=lw, labels_sim=model_names, xlabel_text=r'$\mathcal{D}_{\rm KS}$ for $\{t_{\rm dur}\}$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + save_name + '_durations_KS.pdf', save_fig=savefigures)

# Logxi-res (KS):
plot_fig_pdf_simple(fig_size, [dist_logxi_res_KS_all[:,m] for m in range(n_models)], [], n_bins=n_bins, log_x=False, c_sim=model_colors, ls_sim=model_linestyles, lw=lw, labels_sim=model_names, xlabel_text=r'$\mathcal{D}_{\rm KS}$ for $\{\xi_{\rm res}\}$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + save_name + '_logxi_res_KS.pdf', save_fig=savefigures)

# Logxi_nonres (KS):
plot_fig_pdf_simple(fig_size, [dist_logxi_nonres_KS_all[:,m] for m in range(n_models)], [], n_bins=n_bins, log_x=False, c_sim=model_colors, ls_sim=model_linestyles, lw=lw, labels_sim=model_names, xlabel_text=r'$\mathcal{D}_{\rm KS}$ for $\{\xi_{\rm non-res}\}$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + save_name + '_logxi_nonres_KS.pdf', save_fig=savefigures)

# Depths (KS):
plot_fig_pdf_simple(fig_size, [dist_D_KS_all[:,m] for m in range(n_models)], [], n_bins=n_bins, log_x=False, c_sim=model_colors, ls_sim=model_linestyles, lw=lw, labels_sim=model_names, xlabel_text=r'$\mathcal{D}_{\rm KS}$ for $\{\delta\}$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + save_name + '_depths_KS.pdf', save_fig=savefigures)

# Depth ratios (KS):
plot_fig_pdf_simple(fig_size, [dist_D_ratio_KS_all[:,m] for m in range(n_models)], [], n_bins=n_bins, log_x=False, c_sim=model_colors, ls_sim=model_linestyles, lw=lw, labels_sim=model_names, xlabel_text=r'$\mathcal{D}_{\rm KS}$ for $\{\delta_{i+1}/\delta_i\}$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + save_name + '_depthratios_KS.pdf', save_fig=savefigures)
#'''

#'''
# Periods (AD):
plot_fig_pdf_simple(fig_size, [dist_P_AD_all[:,m] for m in range(n_models)], [], n_bins=n_bins, log_x=False, c_sim=model_colors, ls_sim=model_linestyles, lw=lw, labels_sim=model_names, xlabel_text=r'$\mathcal{D}_{\rm AD^\prime}$ for $\{P\}$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + save_name + '_periods_AD.pdf', save_fig=savefigures)

# Period ratios (AD):
plot_fig_pdf_simple(fig_size, [dist_Rm_AD_all[:,m] for m in range(n_models)], [], n_bins=n_bins, log_x=False, c_sim=model_colors, ls_sim=model_linestyles, lw=lw, labels_sim=model_names, xlabel_text=r'$\mathcal{D}_{\rm AD^\prime}$ for $\{\mathcal{P}\}$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + save_name + '_periodratios_AD.pdf', save_fig=savefigures)

# Durations (AD):
plot_fig_pdf_simple(fig_size, [dist_tdur_AD_all[:,m] for m in range(n_models)], [], n_bins=n_bins, log_x=False, c_sim=model_colors, ls_sim=model_linestyles, lw=lw, labels_sim=model_names, xlabel_text=r'$\mathcal{D}_{\rm AD^\prime}$ for $\{t_{\rm dur}\}$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + save_name + '_durations_AD.pdf', save_fig=savefigures)

# Logxi-res (AD):
plot_fig_pdf_simple(fig_size, [dist_logxi_res_AD_all[:,m] for m in range(n_models)], [], n_bins=n_bins, log_x=False, c_sim=model_colors, ls_sim=model_linestyles, lw=lw, labels_sim=model_names, xlabel_text=r'$\mathcal{D}_{\rm AD^\prime}$ for $\{\xi_{\rm res}\}$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + save_name + '_logxi_res_AD.pdf', save_fig=savefigures)

# Logxi_nonres (AD):
plot_fig_pdf_simple(fig_size, [dist_logxi_nonres_AD_all[:,m] for m in range(n_models)], [], n_bins=n_bins, log_x=False, c_sim=model_colors, ls_sim=model_linestyles, lw=lw, labels_sim=model_names, xlabel_text=r'$\mathcal{D}_{\rm AD^\prime}$ for $\{\xi_{\rm non-res}\}$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + save_name + '_logxi_nonres_AD.pdf', save_fig=savefigures)

# Depths (AD):
plot_fig_pdf_simple(fig_size, [dist_D_AD_all[:,m] for m in range(n_models)], [], n_bins=n_bins, log_x=False, c_sim=model_colors, ls_sim=model_linestyles, lw=lw, labels_sim=model_names, xlabel_text=r'$\mathcal{D}_{\rm AD^\prime}$ for $\{\delta\}$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + save_name + '_depths_AD.pdf', save_fig=savefigures)

# Depth ratios (AD):
plot_fig_pdf_simple(fig_size, [dist_D_ratio_AD_all[:,m] for m in range(n_models)], [], n_bins=n_bins, log_x=False, c_sim=model_colors, ls_sim=model_linestyles, lw=lw, labels_sim=model_names, xlabel_text=r'$\mathcal{D}_{\rm AD^\prime}$ for $\{\delta_{i+1}/\delta_i\}$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + save_name + '_depthratios_AD.pdf', save_fig=savefigures)
#'''

plt.show()
plt.close()
