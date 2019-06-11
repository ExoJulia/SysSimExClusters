from ExoplanetsSysSim_functions import *
import csv





##### To load the files with the GP evaluated points:

savefigures = True
run_directory = 'Model_Optimization/Julia_v0.7/Kepler_catalog_optimization/q1q17_dr25_gaia_fgk_stars79935/Non_Clustered/f_high_incl_low_incl_mmr/Fit_rate_mult_P_Pratios_D_Dratios_dur_durratios_mmr/Some8_params_CRPDr_KS/Fixed_Rbreak3_Ncrit8/lc_1_8_alphaP_-2_2_alphaR1_-4_2_alphaR2_-6_0_ecc_0_0p1_incl_inclmmr_0_90/targs79935_maxincl0_maxiters5000/sigma_i_greater_sigma_i_mmr/GP_files/'
loadfiles_directory = 'ACI/' + run_directory #'ExoplanetsSysSim_Clusters/clusters_v0.7/GP_files/'
sub_directory = '' #'Transformed_rates/'
savefigures_directory = 'ExoplanetsSysSim_Clusters/Figures/' + run_directory + sub_directory #'ExoplanetsSysSim_Clusters/clusters_v0.7/GP_files/'
model_name = 'Non_Clustered_Model' #'Non_Clustered_Model', 'Clustered_P_Model', 'Clustered_P_R_Model'

active_params_symbols = [r'$f_{\sigma_{i,\rm high}}$',
                         #r'$\ln{(\lambda_c)}$',
                         r'$\ln{(\lambda_p)}$',
                         r'$\alpha_P$',
                         r'$\alpha_{R1}$',
                         r'$\alpha_{R2}$',
                         r'$\sigma_e$',
                         r'$\sigma_{i,\rm high}$ ($^\circ$)', # '\n',
                         r'$\sigma_{i,\rm low}$ ($^\circ$)', # '\n',
                         #r'$\sigma_R$',
                         #r'$\sigma_P$'
                         ] #this list of parameter symbols must match the order of parameters in the loaded table!
dims = len(active_params_symbols)

active_params_transformed_symbols = np.copy(active_params_symbols)
active_params_transformed_symbols[1] = r'$\ln{(\lambda_c \lambda_p)}$' #'\n'
active_params_transformed_symbols[2] = r'$\ln{(\lambda_p /\lambda_c)}$' #'\n'

# To load the training points:
best, every = 100000, 10
data_train = load_training_points(dims, file_name_path=loadfiles_directory, best=best, every=every)
active_params_names = np.array(data_train['active_params_names'])

# To load the tables of points drawn from the prior based on the GP model:
n_train, mean_f, sigma_f, lscales, vol = 2000, 75.0, 1.0, 32.31, 1.09
n_points, max_mean, max_std, max_post = 10000, 'Inf', 'Inf', -22.0 #100000, 'Inf', 'Inf', 'Inf'
file_name = 'GP_train%s_meanf%s_sigmaf%s_lscales%s_vol%s_points%s_mean%s_std%s_post%s.csv' % (n_train, mean_f, sigma_f, lscales, vol, n_points, max_mean, max_std, max_post)
xprior_accepted_table = load_GP_table_prior_draws(file_name, file_name_path=loadfiles_directory + sub_directory)
active_params_transformed_names = np.array(xprior_accepted_table.dtype.names[:dims])

GP_model_name = '_GP_train%s_meanf%s_sigmaf%s_lscales%s' % (n_train, mean_f, sigma_f, lscales)
model_name = model_name + GP_model_name





##### To plot the mean, std, and posterior draws as histograms:

#plot_fig_hists_GP_draws((16,8), xprior_accepted_table, save_name=savefigures_directory + model_name + '_vol%s_prior%s_GP_mean%s_std%s_post%s_hists.pdf' % (vol, len(xprior_accepted_table), max_mean, max_std, max_post), save_fig=savefigures)
plt.show()

##### To make corner plots for the GP draws:

mean_cut, std_cut, post_cut = np.inf, np.inf, -22.
xprior_accepts = make_cuts_GP_mean_std_post(active_params_transformed_names, xprior_accepted_table, max_mean=mean_cut, max_std=std_cut, max_post=post_cut)

#plot_cornerpy_wrapper(active_params_symbols, xprior_accepts, save_name=savefigures_directory + model_name + '_vol%s_prior%s_GP_mean%s_std%s_post%s_corner.pdf' % (vol, len(xprior_accepts), mean_cut, std_cut, post_cut), save_fig=savefigures)

#plot_cornerpy_wrapper(active_params_symbols, transform_sum_diff_params_inverse(xprior_accepts, 1, 2), save_name=savefigures_directory + model_name + '_vol%s_prior%s_GP_mean%s_std%s_post%s_corner.pdf' % (vol, len(xprior_accepts), mean_cut, std_cut, post_cut), save_fig=savefigures)
#plot_cornerpy_wrapper(active_params_transformed_symbols, xprior_accepts, save_name=savefigures_directory + model_name + '_vol%s_prior%s_GP_mean%s_std%s_post%s_transformed_corner.pdf' % (vol, len(xprior_accepts), mean_cut, std_cut, post_cut), save_fig=savefigures)
plt.close()




##### To make a custom plotting function for making 'corner' plots with contours based on an array of values instead of the density of points:
'''
grid_dims = 50
GP_grids = load_GP_2d_grids(dims, n_train, mean_f, sigma_f, lscales, file_name_path=loadfiles_directory, grid_dims=grid_dims)

dist_cut = -15. # after subtracting the mean function
GP_prob_below_dist_cut_2d_grids = cdf_normal(dist_cut, mu=GP_grids['mean_grids'], std=GP_grids['std_grids'])

xtrain_sample = data_train['xtrain'][np.random.choice(np.arange(len(data_train['xtrain'])), 1000)]

plot_contours_and_points_corner(active_params_symbols, GP_grids['xlower'], GP_grids['xupper'], GP_grids['mean_grids'], xpoints=xtrain_sample, points_alpha=0.1, save_name=savefigures_directory + model_name + '_grids2d_mean_corner.pdf', save_fig=savefigures)

plot_contours_and_points_corner(active_params_symbols, GP_grids['xlower'], GP_grids['xupper'], GP_grids['std_grids'], xpoints=xtrain_sample, points_alpha=0.1, save_name=savefigures_directory + model_name + '_grids2d_std_corner.pdf', save_fig=savefigures)

plot_contours_and_points_corner(active_params_symbols, GP_grids['xlower'], GP_grids['xupper'], GP_prob_below_dist_cut_2d_grids, xpoints=xtrain_sample, points_alpha=0.1, save_name=savefigures_directory + model_name + '_grids2d_frac_mean%s_corner.pdf' % dist_cut, save_fig=savefigures)

plt.close()
'''
