from ExoplanetsSysSim_functions import *





savefigures = False

savefigures_directory = 'ExoplanetsSysSim_Clusters/Figures/'
save_name = 'Models_Compare'





##### To load the underlying populations:

# Model 1:
loadfiles_directory1 = 'ACI/Simulated_Data/Julia_v0.7/Kepler_catalog_optimization/q1q17_dr25_gaia_fgk_stars79935/Clustered_P_R/f_high_incl_low_incl_mmr/Fit_rate_mult_P_Pratios_D_Dratios_dur_durratios_mmr/Some11_params_CRPDr_KS/Fixed_Rbreak3_Ncrit8/lc_0p2_5_lp_0p5_10_alphaP_-2_2_alphaR1_-4_2_alphaR2_-6_0_ecc_0_0p1_incl_inclmmr_0_90_sigmaR_0_0p5_sigmaP_0_0p3/targs79935_maxincl0_maxiters5000/sigma_i_greater_sigma_i_mmr/GP_med/' #'ExoplanetsSysSim_Clusters/clusters_v0.7/'
run_number1 = ''

N_sim, cos_factor, P_min, P_max, radii_min, radii_max = read_targets_period_radius_bounds(loadfiles_directory1 + 'periods%s.out' % run_number1)

param_vals_all1 = read_sim_params(loadfiles_directory1 + 'periods%s.out' % run_number1)
sssp_per_sys1, sssp1 = compute_summary_stats_from_cat_phys(file_name_path=loadfiles_directory1, run_number=run_number1)

# Model 2:
loadfiles_directory2 = 'ACI/Simulated_Data/Julia_v0.7/Kepler_catalog_optimization/q1q17_dr25_gaia_fgk_stars79935/Clustered_P/f_high_incl_low_incl_mmr/Fit_rate_mult_P_Pratios_D_Dratios_dur_durratios_mmr/Some10_params_CRPDr_KS/Fixed_Rbreak3_Ncrit8/lc_0p2_5_lp_0p5_10_alphaP_-2_2_alphaR1_-4_2_alphaR2_-6_0_ecc_0_0p1_incl_inclmmr_0_90_sigmaP_0_0p3/targs79935_maxincl0_maxiters5000/sigma_i_greater_sigma_i_mmr/GP_med/'
run_number2 = ''

param_vals_all2 = read_sim_params(loadfiles_directory2 + 'periods%s.out' % run_number2)
sssp_per_sys2, sssp2 = compute_summary_stats_from_cat_phys(file_name_path=loadfiles_directory2, run_number=run_number2)

# Model 3:
loadfiles_directory3 = 'ACI/Simulated_Data/Julia_v0.7/Kepler_catalog_optimization/q1q17_dr25_gaia_fgk_stars79935/Non_Clustered/f_high_incl_low_incl_mmr/Fit_rate_mult_P_Pratios_D_Dratios_dur_durratios_mmr/Some8_params_CRPDr_KS/Fixed_Rbreak3_Ncrit8/lc_1_8_alphaP_-2_2_alphaR1_-4_2_alphaR2_-6_0_ecc_0_0p1_incl_inclmmr_0_90/targs79935_maxincl0_maxiters5000/sigma_i_greater_sigma_i_mmr/GP_med/'
run_number3 = ''

param_vals_all3 = read_sim_params(loadfiles_directory3 + 'periods%s.out' % run_number3)
sssp_per_sys3, sssp3 = compute_summary_stats_from_cat_phys(file_name_path=loadfiles_directory3, run_number=run_number3)

model_names = ['Clustered P+R', 'Clustered P', 'Non-clustered'] # Make sure this matches the models loaded!
model_linestyles = ['-','-','-'] #['-', '-.', '--']
model_colors = ['b', 'g', 'r']





##### To plot the simulated catalog as marginal distributions:

subdirectory = 'Paper_Figures/Models/Underlying/' #'Paper_Figures/'; 'Talk_Figures/'

fig_size = (8,3) #size of each panel (figure)
fig_lbrt = [0.15, 0.3, 0.95, 0.925]

n_bins = 100
lw = 1 #linewidth

afs = 20 #axes labels font size
tfs = 20 #text labels font size
lfs = 16 #legend labels font size

#'''
# Multiplicities:
ax = setup_fig_single(fig_size, fig_lbrt[0], fig_lbrt[1], fig_lbrt[2], fig_lbrt[3])
Mtot_list = [sssp1['Mtot_all'], sssp2['Mtot_all'], sssp3['Mtot_all']]
for i,x in enumerate(Mtot_list):
    counts, bins = np.histogram(x, bins=np.max(x)+1, range=(-0.5, np.max(x)+0.5))
    counts[0] = N_sim - len(x) #to compute the number of systems with no planets
    bins_mid = (bins[:-1] + bins[1:])/2.
    #plt.plot(bins_mid, counts/float(np.sum(counts)), 'o-', color=model_colors[i], lw=lw, label=model_names[i])
    plt.plot(bins_mid, counts/float(np.sum(counts)), drawstyle='steps-mid', color=model_colors[i], lw=lw, label=model_names[i])
ax.tick_params(axis='both', labelsize=afs)
plt.xlim([0, 12]) #max([np.max(x) for x in Mtot_list])
plt.xlabel('Planets per system', fontsize=tfs)
plt.ylabel('Fraction', fontsize=tfs)
plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, fontsize=lfs) #show the legend
if savefigures:
    plt.savefig(savefigures_directory + subdirectory + save_name + '_underlying_multiplicities.pdf')
    plt.close()

# Numbers of clusters:
Nc_max = np.max([np.max(sssp['clustertot_all']) for sssp in [sssp1, sssp2]])
plot_fig_pdf_simple(fig_size, [sssp1['clustertot_all'], sssp2['clustertot_all']], [], x_min=0.5, x_max=Nc_max+0.5, n_bins=Nc_max, log_y=True, c_sim=model_colors[:2], lw=lw, ls_sim=model_linestyles[:2], labels_sim=model_names[:2], xlabel_text=r'Clusters per system $N_c$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + save_name + '_underlying_clusters.pdf', save_fig=savefigures)

# Numbers of planets per cluster:
Np_max = np.max([np.max(sssp['pl_per_cluster_all']) for sssp in [sssp1, sssp2]])
plot_fig_pdf_simple(fig_size, [sssp1['pl_per_cluster_all'], sssp2['pl_per_cluster_all']], [], x_min=0.5, x_max=Np_max+0.5, n_bins=Np_max, c_sim=model_colors[:2], lw=lw, ls_sim=model_linestyles[:2], labels_sim=model_names[:2], xlabel_text=r'Planets per cluster $N_p$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + save_name + '_underlying_planets_per_cluster.pdf', save_fig=savefigures)

# Periods:
plot_fig_pdf_simple(fig_size, [sssp1['P_all'], sssp2['P_all'], sssp3['P_all']], [], x_min=P_min, x_max=P_max, n_bins=n_bins, log_x=True, c_sim=model_colors, lw=lw, ls_sim=model_linestyles, labels_sim=model_names, xticks_custom=[3,10,30,100,300], xlabel_text=r'$P$ (days)', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + save_name + '_underlying_periods.pdf', save_fig=savefigures)

# Period ratios (all):
plot_fig_pdf_simple(fig_size, [sssp1['Rm_all'], sssp2['Rm_all'], sssp3['Rm_all']], [], x_min=1., x_max=10., y_max=0.05, n_bins=n_bins, log_x=True, c_sim=model_colors, lw=lw, ls_sim=model_linestyles, labels_sim=model_names, xticks_custom=[1,2,3,4,5,10], xlabel_text=r'$P_{i+1}/P_i$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt)
plt.minorticks_off()
if savefigures:
    plt.savefig(savefigures_directory + subdirectory + save_name + '_underlying_periodratios.pdf')
    plt.close()

# Period ratios (< 5):
plot_fig_pdf_simple(fig_size, [sssp1['Rm_all'][sssp1['Rm_all'] < 5], sssp2['Rm_all'][sssp2['Rm_all'] < 5], sssp3['Rm_all'][sssp3['Rm_all'] < 5]], [], x_min=1., x_max=5., y_max=0.04, n_bins=n_bins, log_x=True, c_sim=model_colors, lw=lw, ls_sim=model_linestyles, labels_sim=model_names, xticks_custom=[1,2,3,4,5], xlabel_text=r'$P_{i+1}/P_i$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt)
plt.minorticks_off()
if savefigures:
    plt.savefig(savefigures_directory + subdirectory + save_name + '_underlying_periodratios_less5.pdf')
    plt.close()

# Eccentricities:
plot_fig_pdf_simple(fig_size, [sssp1['e_all'], sssp2['e_all'], sssp3['e_all']], [], x_min=0., x_max=0.1, n_bins=n_bins, c_sim=model_colors, lw=lw, ls_sim=model_linestyles, labels_sim=model_names, xlabel_text=r'$e$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + save_name + '_underlying_eccentricities.pdf', save_fig=savefigures)

# Planet masses:
plot_fig_pdf_simple(fig_size, [sssp1['mass_all'], sssp2['mass_all'], sssp3['mass_all']], [], x_min=0.07, x_max=1e3, y_max=0.05, n_bins=n_bins, log_x=True, c_sim=model_colors, lw=lw, ls_sim=model_linestyles, labels_sim=model_names, xlabel_text=r'$M_p$ ($M_\oplus$)', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + save_name + '_underlying_masses.pdf', save_fig=savefigures)

# Planet radii:
plot_fig_pdf_simple(fig_size, [sssp1['radii_all'], sssp2['radii_all'], sssp3['radii_all']], [], x_min=radii_min, x_max=radii_max, y_max=0.03, n_bins=n_bins, log_x=True, c_sim=model_colors, lw=lw, ls_sim=model_linestyles, labels_sim=model_names, xticks_custom=[0.5,1,2,4,10], xlabel_text=r'$R_p$ ($R_\oplus$)', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + save_name + '_underlying_radii.pdf', save_fig=savefigures)

# Planet radii ratios:
plot_fig_pdf_simple(fig_size, [sssp1['radii_ratio_all'], sssp2['radii_ratio_all'], sssp3['radii_ratio_all']], [], x_min=1e-1, x_max=10., n_bins=n_bins, log_x=True, c_sim=model_colors, lw=lw, ls_sim=model_linestyles, labels_sim=model_names, xlabel_text=r'$R_{p,i+1}/R_{p,i}$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + save_name + '_underlying_radii_ratios.pdf', save_fig=savefigures)

# Separations in mutual Hill radii:
plot_fig_pdf_simple(fig_size, [sssp1['N_mH_all'], sssp2['N_mH_all'], sssp3['N_mH_all']], [], x_min=8., x_max=200., y_max=0.03, n_bins=n_bins, log_x=True, c_sim=model_colors, lw=lw, ls_sim=model_linestyles, labels_sim=model_names, xlabel_text=r'$\Delta$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + save_name + '_underlying_deltas.pdf', save_fig=savefigures)

# Stellar radii:
plot_fig_pdf_simple(fig_size, [sssp1['Rstar_all'], sssp2['Rstar_all'], sssp3['Rstar_all']], [], x_min=0.5, x_max=2.5, n_bins=n_bins, c_sim=model_colors, lw=lw, ls_sim=model_linestyles, labels_sim=model_names, xlabel_text=r'$R_\star (R_\odot)$', afs=afs, tfs=tfs, lfs=lfs, fig_lbrt=fig_lbrt, save_name=savefigures_directory + subdirectory + save_name + '_underlying_stellar_radii.pdf', save_fig=savefigures)

plt.show()
plt.close()
#'''





'''
##### For ICS 2019 Poster:

lw = 3 #linewidth

fig = plt.figure(figsize=(12,9))
plot = GridSpec(3,2, left=0.15, bottom=0.1, right=0.98, top=0.98, wspace=0.2, hspace=0.35)

# Multiplicities:
ax = plt.subplot(plot[0,0])
Mtot_list = [sssp1['Mtot_all'], sssp2['Mtot_all'], sssp3['Mtot_all']]
for i,x in enumerate(Mtot_list):
    counts, bins = np.histogram(x, bins=np.max(x)+1, range=(-0.5, np.max(x)+0.5))
    counts[0] = N_sim - len(x) #to compute the number of systems with no planets
    bins_mid = (bins[:-1] + bins[1:])/2.
    plt.plot(bins_mid, counts, 'o-', color=model_colors[i], lw=lw, label=model_names[i])
ax.tick_params(axis='both', labelsize=afs)
plt.xlim([0, max([np.max(x) for x in Mtot_list])])
plt.xlabel('Planets per system', fontsize=tfs)
plt.ylabel('Number', fontsize=tfs)
plt.legend(loc='upper right', bbox_to_anchor=(0.95,0.95), ncol=1, frameon=False, fontsize=lfs) #show the legend

# Separations in mutual Hill radii:
ax = plt.subplot(plot[0,1])
plot_panel_pdf_simple(ax, [sssp1['N_mH_all'], sssp2['N_mH_all'], sssp3['N_mH_all']], [], y_max=0.05, n_bins=n_bins, log_x=True, c_sim=model_colors, ls_sim=model_linestyles, lw=lw, labels_sim=model_names, xlabel_text=r'$\Delta$', ylabel_text='', afs=afs, tfs=tfs, lfs=lfs)

# Periods:
ax = plt.subplot(plot[1,0])
plot_panel_pdf_simple(ax, [sssp1['P_all'], sssp2['P_all'], sssp3['P_all']], [], x_min=P_min, x_max=P_max, n_bins=n_bins, log_x=True, c_sim=model_colors, ls_sim=model_linestyles, lw=lw, labels_sim=model_names, xticks_custom=[3,10,30,100,300], xlabel_text=r'$P$ (days)', afs=afs, tfs=tfs, lfs=lfs)

# Period ratios (all):
ax = plt.subplot(plot[1,1])
plot_panel_pdf_simple(ax, [sssp1['Rm_all'], sssp2['Rm_all'], sssp3['Rm_all']], [], x_min=1., y_max=0.1, n_bins=n_bins, log_x=True, c_sim=model_colors, ls_sim=model_linestyles, lw=lw, labels_sim=model_names, xticks_custom=[1,2,3,4,5,10,20,50,100], xlabel_text=r'$P_{i+1}/P_i$', ylabel_text='', afs=afs, tfs=tfs, lfs=lfs)

# Planet radii:
ax = plt.subplot(plot[2,0])
plot_panel_pdf_simple(ax, [sssp1['radii_all'], sssp2['radii_all'], sssp3['radii_all']], [], x_min=0.5, x_max=10., y_max=0.03, n_bins=n_bins, log_x=True, c_sim=model_colors, ls_sim=model_linestyles, lw=lw, labels_sim=model_names, xticks_custom=[0.5,1,2,4,10], xlabel_text=r'$R_p$ ($R_\oplus$)', afs=afs, tfs=tfs, lfs=lfs)

# Planet radii ratios:
ax = plt.subplot(plot[2,1])
plot_panel_pdf_simple(ax, [sssp1['radii_ratio_all'], sssp2['radii_ratio_all'], sssp3['radii_ratio_all']], [], n_bins=n_bins, log_x=True, c_sim=model_colors, ls_sim=model_linestyles, lw=lw, labels_sim=model_names, xlabel_text=r'$R_{p,i+1}/R_{p,i}$', ylabel_text='', afs=afs, tfs=tfs, lfs=lfs)

if savefigures:
    plt.savefig(savefigures_directory + save_name + '_underlying.pdf')
    plt.close()
else:
    plt.show()
'''






'''
# Model 1:
loadfiles_directory1 = 'ACI/Simulated_Data/Julia_v0.7/Kepler_catalog_optimization/q1q17_dr25_gaia_fgk_stars80006/Clustered_P_R/f_high_incl_low_incl_mmr/Fit_rate_mult_P_Pratios_D_Dratios_dur_durratios_mmr/Some11_params_KSweightedrms/lc_lp_0p5_5_alphaP_-2_1_alphaR1_R2_-6_0_ecc_0_0p1_incl_inclmmr_0_90_sigmaR_0_0p5_sigmaP_0_0p3/Fixed_Rbreak3_Ncrit8/targs400030_maxincl0_maxiters5000/sigma_i_greater_sigma_i_mmr/best_N/' #'ExoplanetsSysSim_Clusters/clusters_v0.7/'

n1_1plus = np.zeros(10)
ptot1 = np.zeros(10)
for i in range(10):
    run_number1 = str(i)

    N_sim, cos_factor, P_min, P_max, radii_min, radii_max = read_targets_period_radius_bounds(loadfiles_directory1 + 'periods%s.out' % run_number1)

    param_vals_all1 = read_sim_params(loadfiles_directory1 + 'periods%s.out' % run_number1)
    sssp_per_sys1, sssp1 = compute_summary_stats_from_cat_phys(file_name_path=loadfiles_directory1, run_number=run_number1)

    n1_1plus[i] = len(sssp1['Mtot_all'])
    ptot1[i] = np.sum(sssp1['Mtot_all'])
n1_none = N_sim - n1_1plus
f1_none = n1_none/N_sim
fp1 = ptot1/float(N_sim)
print 'Clustered_P_R'
print 'Number of stars with no planets: ', n1_none
print 'Fraction of stars with no planets: ', f1_none
print 'Total number of planets: ', ptot1
print 'Ratio of planets to stars: ', fp1

# Model 2:
loadfiles_directory2 = 'ACI/Simulated_Data/Julia_v0.7/Kepler_catalog_optimization/q1q17_dr25_gaia_fgk_stars80006/Clustered_P/f_high_incl_low_incl_mmr/Fit_rate_mult_P_Pratios_D_Dratios_dur_durratios_mmr/Some10_params_KSweightedrms/Fixed_Rbreak3_Ncrit8/lc_lp_0p5_5_alphaP_-2_1_alphaR1_R2_-6_0_ecc_0_0p1_incl_inclmmr_0_90_sigmaP_0_0p3/targs400030_maxincl0_maxiters5000/sigma_i_greater_sigma_i_mmr/best_N/'

n2_1plus = np.zeros(10)
ptot2 = np.zeros(10)
for i in range(10):
    run_number2 = str(i)

    param_vals_all2 = read_sim_params(loadfiles_directory2 + 'periods%s.out' % run_number2)
    sssp_per_sys2, sssp2 = compute_summary_stats_from_cat_phys(file_name_path=loadfiles_directory2, run_number=run_number2)

    n2_1plus[i] = len(sssp2['Mtot_all'])
    ptot2[i] = np.sum(sssp2['Mtot_all'])
n2_none = N_sim - n2_1plus
f2_none = n2_none/N_sim
fp2 = ptot2/float(N_sim)
print 'Clustered_P'
print 'Number of stars with no planets: ', n2_none
print 'Fraction of stars with no planets: ', f2_none
print 'Total number of planets: ', ptot2
print 'Ratio of planets to stars: ', fp2

# Model 3:
loadfiles_directory3 = 'ACI/Simulated_Data/Julia_v0.7/Kepler_catalog_optimization/q1q17_dr25_gaia_fgk_stars80006/Non_Clustered/f_high_incl_low_incl_mmr/Fit_rate_mult_P_Pratios_D_Dratios_dur_durratios_mmr/Some8_params_KSweightedrms/Fixed_Rbreak3_Ncrit8/lc_1_10_alphaP_-2_1_alphaR1_R2_-6_0_ecc_0_0p1_incl_inclmmr_0_90/targs400030_maxincl0_maxiters5000/sigma_i_greater_sigma_i_mmr/best_N/'

n3_1plus = np.zeros(10)
ptot3 = np.zeros(10)
for i in range(10):
    run_number3 = str(i)

    param_vals_all3 = read_sim_params(loadfiles_directory3 + 'periods%s.out' % run_number3)
    sssp_per_sys3, sssp3 = compute_summary_stats_from_cat_phys(file_name_path=loadfiles_directory3, run_number=run_number3)

    n3_1plus[i] = len(sssp3['Mtot_all'])
    ptot3[i] = np.sum(sssp3['Mtot_all'])
n3_none = N_sim - n3_1plus
f3_none = n3_none/N_sim
fp3 = ptot3/float(N_sim)
print 'Non_clustered:'
print 'Number of stars with no planets: ', n3_none
print 'Fraction of stars with no planets: ', f3_none
print 'Total number of planets: ', ptot3
print 'Ratio of planets to stars: ', fp3
'''




