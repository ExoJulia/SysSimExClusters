from ExoplanetsSysSim_functions import *
import csv





##### To load the files with the GP evaluated points:

savefigures = True
run_directory = 'Model_Optimization/Julia_v0.7/Kepler_catalog_optimization/q1q17_dr25_gaia_fgk_stars79935/Clustered_P_R/f_high_incl_low_incl_mmr/Fit_rate_mult_P_Pratios_D_Dratios_dur_durratios_mmr/Some11_params_CRPDr_AD/Fixed_Rbreak3_Ncrit8/lc_0p2_5_lp_0p5_10_alphaP_-2_2_alphaR1_-4_2_alphaR2_-6_0_ecc_0_0p1_incl_inclmmr_0_90_sigmaR_0_0p5_sigmaP_0_0p3/targs79935_maxincl0_maxiters5000/sigma_i_greater_sigma_i_mmr/AD_mod/GP_files/'
loadfiles_directory = 'ACI/' + run_directory
sub_directory = 'Transformed_rates/'
savefigures_directory = 'ExoplanetsSysSim_Clusters/Figures/' + run_directory + sub_directory
model_name = 'Clustered_P_R_Model' #'Non_Clustered_Model', 'Clustered_P_Model', 'Clustered_P_R_Model'

active_params_symbols = [r'$f_{\sigma_{i,\rm high}}$',
                         r'$\ln{(\lambda_c)}$',
                         r'$\ln{(\lambda_p)}$',
                         r'$\alpha_P$',
                         r'$\alpha_{R1}$',
                         r'$\alpha_{R2}$',
                         r'$\sigma_e$',
                         r'$\sigma_{i,\rm high}$ ($^\circ$)',
                         r'$\sigma_{i,\rm low}$ ($^\circ$)',
                         r'$\sigma_R$',
                         r'$\sigma_N$'
                         ] #this list of parameter symbols must match the order of parameters in the loaded table!
dims = len(active_params_symbols)

active_params_transformed_symbols = np.copy(active_params_symbols)
active_params_transformed_symbols[1] = r'$\ln{(\lambda_c \lambda_p)}$'
active_params_transformed_symbols[2] = r'$\ln{(\lambda_p /\lambda_c)}$'

# To load the training points:
best, every = 100000, 10
data_train = load_training_points(dims, file_name_path=loadfiles_directory, best=best, every=every)
active_params_names = np.array(data_train['active_params_names'])

##### If we want to compute and plot the un-logged rates (i.e. lambda_c, lambda_p):
#data_train['xtrain'][:,[1,2]] = np.exp(data_train['xtrain'][:,[1,2]])
#active_params_symbols[1], active_params_symbols[2] = r'$\lambda_c$', r'$\lambda_p$'
#active_params_transformed_symbols[1], active_params_transformed_symbols[2] = r'$\lambda_c \lambda_p$', r'$\lambda_p /\lambda_c$'

# To make corner plots for the GP training points:

plot_cornerpy_wrapper(active_params_symbols, data_train['xtrain'], save_name=savefigures_directory + model_name + '_training_corner.pdf', save_fig=savefigures)
plot_cornerpy_wrapper(active_params_transformed_symbols, transform_sum_diff_params(data_train['xtrain'], 1, 2), save_name=savefigures_directory + model_name + '_training_transformed_corner.pdf', save_fig=savefigures)
plt.close()





##### To load the table of points minimizing the GP mean and overplot them:
'''
n_points_min = 100
file_name = 'GP_train%s_meanf%s_sigmaf%s_lscales%s_minimize_mean_iterations%s.csv' % (n_train, mean_f, sigma_f, lscales, n_points_min)
xmin_table = load_table_points_min_GP(file_name, file_name_path=loadfiles_directory + sub_directory)
xmins = xmin_table[active_params_names].view((float, dims))

plot_cornerpy_wrapper(active_params_symbols, data_train['xtrain'], xpoints_extra=xmins, save_name=savefigures_directory + model_name + '_training_corner.pdf', save_fig=savefigures)
plot_cornerpy_wrapper(active_params_transformed_symbols, transform_sum_diff_params(data_train['xtrain'], 1, 2), xpoints_extra=transform_sum_diff_params(xmins, 1, 2), save_name=savefigures_directory + model_name + '_training_transformed_corner.pdf', save_fig=savefigures)
plt.close()
'''
