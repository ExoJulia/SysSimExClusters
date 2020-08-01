##### GP regression to fit to the results of the model optimization (fit GP model to the distance of the model compared to the Kepler data as a function of the model parameters)

include("GP_functions.jl")

using PyPlot
using Optim





function load_data(dims::Int64, data_path::String)
    data_table_original =  CSV.read(joinpath(data_path,"Active_params_distances_table_best100000_every10.txt"), comment="#")
    data_table_recomputed = CSV.read(joinpath(data_path,"Active_params_recomputed_distances_table_best100000_every10.txt"), comment="#")

    params_names = names(data_table_recomputed)[1:dims]
    dists_names = names(data_table_recomputed)[dims+1:end]

    params_array_original = convert(Matrix, data_table_original[1:end, params_names])
    dist_array_original = convert(Array, data_table_original[1:end, :dist_tot_weighted])

    params_array = convert(Matrix, data_table_recomputed[1:end, params_names])
    dist_array = convert(Array, data_table_recomputed[1:end, :dist_tot_weighted])

    # If we want to transform some parameters by sum and difference:
    transform_sum_diff_params!(params_names, params_array, 4, 5) #####

    data = Dict()
    data[:params_names] = params_names
    data[:dists_names] = dists_names
    #data[:params_array_original] = params_array_original
    data[:dist_array_original] = dist_array_original
    data[:params_array] = params_array
    data[:dist_array] = dist_array

    return data
end

"""
Transforms a dataset ('params_names' and 'params_array') in place, by performing a sum and difference transformation of two parameters.
"""
function transform_sum_diff_params!(params_names::Vector{String}, params_array::Array{Float64,2}, i::Int64, j::Int64)
    @assert length(params_names) == size(params_array, 2)
    dims = length(params_names)
    @assert i != j
    @assert (i <= dims) && (j <= dims)

    sum_params_name = "sum_"*params_names[i]*"_"*params_names[j]
    diff_params_name = "diff_"*params_names[j]*"_"*params_names[i]

    println("Transforming ($(params_names[i]), $(params_names[j])) to ($sum_params_name, $diff_params_name).")

    params_names[i] = sum_params_name
    params_names[j] = diff_params_name
    params_array[:,i], params_array[:,j] = params_array[:,i] .+ params_array[:,j], params_array[:,j] .- params_array[:,i]
end

"""
Loads a dataset and trains a GP model on it by returning the training points, a mean function (constant), and a set of hyperparameters for the GP kernel.
Also gives the option to make a bunch of plots showing how well the GP model performs.
"""
function train_GP_emulator(; dims::Int64, data_path::String, f_err::Float64=0.8, n_train::Int64, n_cv::Int64, mean_f::Float64, kernel::Function, hparams_best::Vector{Float64}=zeros(dims+1), optimize_hparams::Bool=false, make_plots::Bool=false)
    @assert n_train > 0
    @assert n_cv > 0

    # To load and plot the data:

    data = load_data(dims, data_path)

    if make_plots
        #=
        fig = histogram([data[:dist_array_original] data[:dist_array]], fillalpha=0.5, xlabel="Distance", ylabel="Number of points", label=["Best distances during optimization", "Recomputed distances"])
        display(fig)

        fig1 = histogram(reshape(data[:dist_array], (div(length(data[:dist_array]), 5), 5)), fillalpha=0.5, xlabel="Distance", ylabel="Number of points")
        display(fig1)
        =#

        fig = figure()
        hist([data[:dist_array_original], data[:dist_array]], bins=100, histtype="step", stacked=false, fill=false, label=["Best distances during optimization", "Recomputed distances"])
        xlabel("Distance")
        ylabel("Number of points")
        legend()

        fig1 = figure()
        hist(reshape(data[:dist_array], (div(length(data[:dist_array]), 5), 5)), bins=100, histtype="step", stacked=false, fill=false)
        xlabel("Distance")
        ylabel("Number of points")
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
        #=
        fig2 = histogram(ydata, fillalpha=0.5, xlabel="Distance - mean_f", ylabel="Number of points")
        display(fig2)
        =#

        fig2 = figure()
        hist(ydata, bins=100, histtype="step", stacked=false, fill=false)
        xlabel("Distance - mean_f")
        ylabel("Number of points")
    end



    # To optimize the hyperparameters, plot the resulting GP model, and assess the fit of the model:

    hparams_guess = [1.; ones(dims)]
    if optimize_hparams
        #println("Optimizing all hyperparameters with MLE...")
        #hparams_best, log_p_best = optimize_hparams_with_MLE(hparams_guess, xtrain, ytrain, kernel; ydata_err=ytrain_err)

        #println("Optimizing 'sigma_f' and 'scale_l' hyperparameters with MLE, fixing 'lscales_rel' = $(hparams_best[2:end])...")
        #hparams_best, log_p_best = optimize_hparams_SE_ndims_fixed_relative_lscales_with_MLE([1., 1.], xtrain, ytrain; ydata_err=ytrain_err, lscales_rel=hparams_best[2:end])

        println("Optimizing 'lscales' hyperparameters with MLE, fixing 'sigma_f' = $(hparams_best[1])...")
        hparams_best, log_p_best = optimize_hparams_SE_ndims_fixed_sigmaf_with_MLE(hparams_guess[2:end], xtrain, ytrain; ydata_err=ytrain_err, sigma_f=hparams_best[1])
    end



    # To predict at the training points:
    mu_train, stdv_train, f_posterior_train = draw_from_posterior_given_kernel_and_data(xtrain, xtrain, ytrain, kernel, hparams_best; ydata_err=ytrain_err)

    # To predict at the checking points (cross-validation):
    mu_cv, stdv_cv, f_posterior_cv = draw_from_posterior_given_kernel_and_data(xcheck, xtrain, ytrain, kernel, hparams_best; ydata_err=ytrain_err)

    # To plot histograms of the residuals of the mean prediction compared to the data:
    ydiff_train = mu_train - ytrain
    ydiff_cv = mu_cv - ycheck

    if make_plots
        #=
        fig3 = histogram([ydiff_train ydiff_cv], fillalpha=0.5, xlabel="Mean prediction - Data", ylabel="Number of points", label=["Training", "Cross-validation"])

        fig4 = scatter([ytrain ycheck], [mu_train mu_cv], markersize=1, xlabel="Data", ylabel="Mean prediction", label=["Training", "Cross-validation"])
        plot!(range(0, stop=maximum(ycheck), length=100), range(0, stop=maximum(ycheck), length=100), label="Perfect prediction")

        fig5 = scatter([ydiff_train ydiff_cv], [stdv_train stdv_cv], markersize=1, xlabel="Mean prediction - Data", ylabel="Uncertainty of prediction", label=["Training", "Cross-validation"])

        fig6 = scatter([mu_train mu_cv], [stdv_train stdv_cv], markersize=1, xlabel="Mean prediction", ylabel="Uncertainty of prediction", label=["Training", "Cross-validation"])

        fig3_6 = plot(fig3,fig4,fig5,fig6, layout=(2,2), guidefontsize=8, legend=true, legendfontsize=4)
        display(fig3_6)
        =#

        fig3 = figure()

        subplot(2,2,1)
        hist([ydiff_train, ydiff_cv], bins=100, histtype="step", stacked=false, fill=false, label=["Training", "Cross-validation"])
        xlabel("Mean prediction - Data")
        ylabel("Number of points")
        legend()

        subplot(2,2,2)
        scatter(ytrain, mu_train, label="Training", s=1)
        scatter(ycheck, mu_cv, label="Cross-validation", s=1)
        plot(range(0, stop=maximum(ycheck), length=100), range(0, stop=maximum(ycheck), length=100), label="Perfect prediction")
        xlabel("Data")
        ylabel("Mean prediction")
        legend()

        subplot(2,2,3)
        scatter(ydiff_train, stdv_train, label="Training", s=1)
        scatter(ydiff_cv, stdv_cv, label="Cross-validation", s=1)
        xlabel("Mean prediction - Data")
        ylabel("Uncertainty of prediction")
        legend()

        subplot(2,2,4)
        scatter(mu_train, stdv_train, label="Training", s=1)
        scatter(mu_cv, stdv_cv, label="Cross-validation", s=1)
        xlabel("Mean prediction")
        ylabel("Uncertainty of prediction")
        legend()
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
        #GP_model[:plots] = [fig, fig1, fig2, fig3, fig4, fig5, fig6, fig3_6]
        GP_model[:plots] = [fig, fig1, fig2, fig3]
    end

    return GP_model
end





data_path = "GP_files"
prior_bounds = nothing

# Transformed:
hparams_best = [2.7, 0.2, 1., 1., 1.5, 3., 1., 0.5, 1., 0.2, 0.15, 0.15]
prior_bounds = [(0.6, 1.), (-0.6, 2.), (-1., 2.), (-1.5, 3.), (6., 15.), (-0.8, 1.6), (-2., -0.5), (-6., -3.), (0., 0.5), (0.1, 0.5), (0.1, 0.5)]





dims = 11
mean_f = 75. # KS
#mean_f = 150. # AD
GP_model = train_GP_emulator(; dims=dims, data_path=data_path, f_err=2.7, n_train=2000, n_cv=2000, mean_f=mean_f, kernel=kernel_SE_ndims, hparams_best=hparams_best, optimize_hparams=false, make_plots=false)
