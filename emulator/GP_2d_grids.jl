include("GP_emulator.jl")





"""
Evaluates 'predict_mean_std_on_2d_grids_all_combinations' which computes the GP mean and std on each 2d grid of parameters given a dataset and GP model (i.e. kernel and hyperparameters), and saves the results as a reshaped 2d array to file (one for GP mean and one for GP std).
"""
function predict_GP_2d_grids_and_save(grid_dims::Int64; params_names::Array{Symbol,1}, xmin_guess::Array{Float64,2}, xdata::Array{Float64,2}, mean_f::Float64, ydata::Vector{Float64}, ydata_err::Vector{Float64}, kernel::Function, hparams::Vector{Float64}, save_path::String="")
    @assert size(xdata, 2) == length(params_names)
    @assert size(xdata, 1) == length(ydata) == length(ydata_err)
    dims = length(params_names)
    n_train = length(ydata)

    xdata_lower, xdata_upper = findmin(xdata, dims=1)[1], findmax(xdata, dims=1)[1]
    @time stacked_grids_dict = predict_mean_std_on_2d_grids_all_combinations(params_names, reshape(xmin_guess, dims), xdata, ydata, kernel, hparams, ydata_err; grid_dims=grid_dims)

    file_name = joinpath(save_path, "GP_train$(n_train)_meanf$(mean_f)_sigmaf$(hparams[1])_lscales$(round(sum(hparams[2:end]), digits=2))_grids2d_$(grid_dims)x$(grid_dims)_mean.csv")
    f = open(file_name, "w")
    println(f, "# n_data = $n_train, mean_f = $mean_f")
    println(f, "# hparams = $hparams")
    println(f, "# xlower = $xdata_lower")
    println(f, "# xupper = $xdata_upper")
    println(f, "# xmin_guess = $xmin_guess")
    CSV.write(f, stacked_grids_dict[:mean_grids]; append=true)
    close(f)

    file_name = joinpath(save_path, "GP_train$(n_train)_meanf$(mean_f)_sigmaf$(hparams[1])_lscales$(round(sum(hparams[2:end]), digits=2))_grids2d_$(grid_dims)x$(grid_dims)_std.csv")
    f = open(file_name, "w")
    println(f, "# n_data = $n_train, mean_f = $mean_f")
    println(f, "# hparams = $hparams")
    println(f, "# xlower = $xdata_lower")
    println(f, "# xupper = $xdata_upper")
    println(f, "# xmin_guess = $xmin_guess")
    CSV.write(f, stacked_grids_dict[:std_grids]; append=true)
    close(f)
end





# To use the optimized GP model to predict on a series of 2d grids of points, with the other dimensions (model parameters) set to best-fit values:

data = load_data(dims, data_path)
xmin_guess = Statistics.median(data[:params_array], dims=1)
println("Median data point: ", xmin_guess)

grid_dims = 50
save_path = data_path
predict_GP_2d_grids_and_save(grid_dims; params_names=GP_model[:params_names], xmin_guess=xmin_guess, xdata=GP_model[:xtrain], mean_f=GP_model[:mean_f], ydata=GP_model[:ytrain], ydata_err=GP_model[:ytrain_err], kernel=kernel_SE_ndims, hparams=GP_model[:hparams_best], save_path=save_path)
