##### GP regression to fit to the results of the model optimization (fit GP model to the distance of the model compared to the Kepler data as a function of the model parameters)

include("GP_functions.jl")

using Plots
using Optim





function load_data(dims::Int64, data_path::String)
    data_table_original = CSV.read(joinpath(data_path,"Active_params_distances_table_best100000_every10.txt"), delim=" ", allowmissing=:none)
    data_table_recomputed = CSV.read(joinpath(data_path,"Active_params_distances_recomputed_table_best100000_every10.txt"), delim=" ", allowmissing=:none)
    #data_table_recomputed = CSV.read("GP_files/Active_params_distances_table_best388529_every1.txt", delim=" ", allowmissing=:none)[1:388000,:]

    params_names = names(data_table_recomputed)[1:dims]
    dists_names = names(data_table_recomputed)[dims+1:end]

    params_array_original = convert(Matrix, data_table_original[1:end, params_names])
    dist_array_original = convert(Array, data_table_original[1:end, :dist_tot_weighted])

    params_array = convert(Matrix, data_table_recomputed[1:end, params_names])
    dist_array = convert(Array, data_table_recomputed[1:end, :dist_tot_weighted])

    data = Dict()
    data[:params_names] = params_names
    data[:dists_names] = dists_names
    data[:params_array_original] = params_array_original
    data[:dist_array_original] = dist_array_original
    data[:params_array] = params_array
    data[:dist_array] = dist_array

    return data
end

"""
Loads a dataset and trains a GP model on it by returning the training points, a mean function (constant), and a set of hyperparameters for the GP kernel.
Also gives the option to make a bunch of plots showing how well the GP model performs.
"""
function train_GP_emulator(; dims::Int64, data_path::String, f_err::Float64=0.8, n_train::Int64, n_cv::Int64, mean_f::Float64, kernel::Function, hparams_best::Vector{Float64}=zeros(dims+1), make_plots::Bool=false)
    @assert n_train > 0
    @assert n_cv > 0

    # To load and plot the data:

    data = load_data(dims, data_path)

    if make_plots
        fig = histogram([data[:dist_array_original] data[:dist_array]], fillalpha=0.5, xlabel="Distance", ylabel="Number of points", label=["Best distances during optimization", "Recomputed distances"])
        display(fig)

        fig1 = histogram(reshape(data[:dist_array], (div(length(data[:dist_array]), 5), 5)), fillalpha=0.5, xlabel="Distance", ylabel="Number of points")
        display(fig1)
    end



    # To choose a subset of the data for training and cross-validation sets:

    n_data = n_train + n_cv

    i_data = Random.randperm(sum(data[:dist_array] .< Inf))[1:n_data]
    #i_data = 1:n_data
    xdata = data[:params_array][i_data, 1:end]
    ydata = data[:dist_array][i_data] .- mean_f
    ydata_err = f_err .*ones(n_data)

    xtrain, xcheck, ytrain, ycheck, ytrain_err, ycheck_err = split_data_training_cv(xdata, ydata, ydata_err; n_train=n_train)

    if make_plots
        fig2 = histogram(ydata, fillalpha=0.5, xlabel="Distance - mean_f", ylabel="Number of points")
        display(fig2)
    end



    # To optimize the hyperparameters, plot the resulting GP model, and assess the fit of the model:

    hparams_guess = [1.; ones(dims)]
    if sum(hparams_best) == 0
        println("Optimizing hyperparameters with MLE...")
        hparams_best, log_p_best = optimize_hparams_with_MLE(hparams_guess, xtrain, ytrain, kernel; ydata_err=ytrain_err)
    end



    # To predict at the training points:
    mu_train, stdv_train, f_posterior_train = draw_from_posterior_given_kernel_and_data(xtrain, xtrain, ytrain, kernel, hparams_best; ydata_err=ytrain_err)

    # To predict at the checking points (cross-validation):
    mu_cv, stdv_cv, f_posterior_cv = draw_from_posterior_given_kernel_and_data(xcheck, xtrain, ytrain, kernel, hparams_best; ydata_err=ytrain_err)

    # To plot histograms of the residuals of the mean prediction compared to the data:
    ydiff_train = mu_train - ytrain
    ydiff_cv = mu_cv - ycheck

    if make_plots
        fig3 = histogram([ydiff_train ydiff_cv], fillalpha=0.5, xlabel="Mean prediction - Data", ylabel="Number of points", label=["Training", "Cross-validation"])

        fig4 = scatter([ytrain ycheck], [mu_train mu_cv], markersize=1, xlabel="Data", ylabel="Mean prediction", label=["Training", "Cross-validation"])
        plot!(range(0, stop=maximum(ycheck), length=100), range(0, stop=maximum(ycheck), length=100), label="Perfect prediction")

        fig5 = scatter([ydiff_train ydiff_cv], [stdv_train stdv_cv], markersize=1, xlabel="Mean prediction - Data", ylabel="Uncertainty of prediction", label=["Training", "Cross-validation"])

        fig6 = scatter([mu_train mu_cv], [stdv_train stdv_cv], markersize=1, xlabel="Mean prediction", ylabel="Uncertainty of prediction", label=["Training", "Cross-validation"])

        fig3_6 = plot(fig3,fig4,fig5,fig6, layout=(2,2), guidefontsize=8, legend=true, legendfontsize=4)
        display(fig3_6)
    end

    GP_model = Dict()
    GP_model[:params_names] = data[:params_names]
    GP_model[:xtrain] = xtrain
    GP_model[:ytrain] = ytrain
    GP_model[:ytrain_err] = ytrain_err
    GP_model[:mean_f] = mean_f
    GP_model[:kernel] = kernel
    GP_model[:hparams_best] = hparams_best
    if make_plots
        GP_model[:plots] = [fig, fig1, fig2, fig3, fig4, fig5, fig6, fig3_6]
    end

    return GP_model
end





data_path = "GP_files"

# Clustered_P_R:
#data_path = "/Users/Matthias/Documents/GradSchool/Eric_Ford_Research/ACI/Model_Optimization/Julia_v0.7/Kepler_catalog_optimization/q1q17_dr25_gaia_fgk_stars80006/Clustered_P_R/f_high_incl_low_incl_mmr/Fit_rate_mult_P_Pratios_D_Dratios_dur_durratios_mmr/Some11_params_KSweightedrms/lc_lp_0p5_5_alphaP_-2_1_alphaR1_R2_-6_0_ecc_0_0p1_incl_inclmmr_0_90_sigmaR_0_0p5_sigmaP_0_0p3/Fixed_Rbreak3_Ncrit8/targs400030_maxincl0_maxiters5000/sigma_i_greater_sigma_i_mmr/GP_files"

#hparams_best = [8.03106, 0.2036, 0.344798, 0.415796, 0.484306, 0.757346, 1.81563, 950.991, 36.9642, 0.851163, 0.159388, 395.641] # Clustered_P_R with 'n_data = 2000', 'mean_f = 30.'
#hparams_best = [1., 0.3, 1.2, 1.5, 1.2, 1.8, 3., 0.03, 60., 1., 0.15, 0.1]

# Clustered_P:
#data_path = "/Users/Matthias/Documents/GradSchool/Eric_Ford_Research/ACI/Model_Optimization/Julia_v0.7/Kepler_catalog_optimization/q1q17_dr25_gaia_fgk_stars80006/Clustered_P/f_high_incl_low_incl_mmr/Fit_rate_mult_P_Pratios_D_Dratios_dur_durratios_mmr/Some10_params_KSweightedrms/Fixed_Rbreak3_Ncrit8/lc_lp_0p5_5_alphaP_-2_1_alphaR1_R2_-6_0_ecc_0_0p1_incl_inclmmr_0_90_sigmaP_0_0p3/targs400030_maxincl0_maxiters5000/sigma_i_greater_sigma_i_mmr/GP_files"
#hparams_best = [1., 0.3, 0.8, 1., 0.5, 0.8, 2., 0.03, 60., 1., 0.1]

# Non_Clustered:
#data_path = "/Users/Matthias/Documents/GradSchool/Eric_Ford_Research/ACI/Model_Optimization/Julia_v0.7/Kepler_catalog_optimization/q1q17_dr25_gaia_fgk_stars80006/Non_Clustered/f_high_incl_low_incl_mmr/Fit_rate_mult_P_Pratios_D_Dratios_dur_durratios_mmr/Some8_params_KSweightedrms/Fixed_Rbreak3_Ncrit8/lc_1_10_alphaP_-2_1_alphaR1_R2_-6_0_ecc_0_0p1_incl_inclmmr_0_90/targs400030_maxincl0_maxiters5000/sigma_i_greater_sigma_i_mmr/GP_files"
#hparams_best = [1., 0.05, 0.25, 0.1, 0.8, 1.5, 0.01, 60., 0.5]

dims = 11
GP_model = train_GP_emulator(; dims=dims, data_path=data_path, n_train=2000, n_cv=2000, mean_f=50., kernel=kernel_SE_ndims, hparams_best=hparams_best, make_plots=true)
