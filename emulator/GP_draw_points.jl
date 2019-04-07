#include("GP_test.jl")
include("GP_emulator.jl")

"""
Evaluates 'predict_model_from_uniform_prior_until_accept_n_points' given a dataset and GP model (i.e. kernel and hyperparameters) to draw 'n_points' with GP mean and std better than 'max_mean' and 'max_std' and save the results as a table in a CSV file.
"""
function draw_points_with_GP_and_save(n_points::Int64; params_names::Array{Symbol,1}, xdata::Array{Float64,2}, mean_f::Float64, ydata::Vector{Float64}, ydata_err::Vector{Float64}, kernel::Function, hparams::Vector{Float64}, max_mean::Float64=Inf, max_std::Float64=Inf, max_post::Float64=Inf, save_path::String="")
    @assert size(xdata, 2) == length(params_names)
    @assert size(xdata, 1) == length(ydata) == length(ydata_err)
    @assert n_points > 0
    dims = length(params_names)
    n_train = length(ydata)

    println("Drawing $n_points points satisfying GP mean < $max_mean, GP std < $max_std, GP posterior draw < $max_post...")
    @time prior_draws_GP_table = predict_model_from_uniform_prior_until_accept_n_points(params_names, xdata, ydata, kernel, hparams, ydata_err, n_points; max_mean=max_mean, max_std=max_std, max_post=max_post)

    file_name = joinpath(save_path, "GP_train$(n_train)_meanf$(mean_f)_sigmaf$(hparams[1])_lscales$(round(sum(hparams[2:end]), digits=2))_points$(n_points)_mean$(max_mean)_std$(max_std)_post$(max_post).csv")
    f = open(file_name, "w")
    println(f, "# n_data = $n_train, mean_f = $mean_f")
    println(f, "# hparams = $hparams")
    println(f, "# n_accept = $n_points, max_mean = $max_mean, max_std = $max_std, max_post = $max_post")
    CSV.write(f, prior_draws_GP_table; append=true, writeheader=true)
    close(f)
end





# To use the optimized GP model to predict at a large number of points drawn from the prior:

n_points = 100000
max_mean, max_std, max_post = Inf, Inf, Inf
save_path = data_path
draw_points_with_GP_and_save(n_points; params_names=GP_model[:params_names], xdata=GP_model[:xtrain], mean_f=GP_model[:mean_f], ydata=GP_model[:ytrain], ydata_err=GP_model[:ytrain_err], kernel=kernel_SE_ndims, hparams=GP_model[:hparams_best], max_mean=max_mean, max_std=max_std, max_post=max_post, save_path=save_path)
