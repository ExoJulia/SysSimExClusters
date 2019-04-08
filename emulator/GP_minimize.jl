include("GP_emulator.jl")





"""
Optimize the parameters of the model (i.e. find the point that minimizes the function) by using a dataset and a GP model (i.e. a given kernel and set of hyperparameters) to minimize the GP mean prediction.
"""
function optimize_model_GP_mean(xmin_guess::Array{Float64,2}, xdata::Array{Float64,2}, ydata::Vector{Float64}, kernel::Function, hparams::Vector{Float64}; ydata_err::Vector{Float64}=zeros(length(ydata)))
    @assert size(xdata)[1] == length(ydata) == length(ydata_err)
    dims = size(xmin_guess,2)

    xdata_lower = minimum(xdata, dims=1)
    xdata_upper = maximum(xdata, dims=1)

    L = compute_kernel_given_data(xdata, kernel, hparams; ydata_err=ydata_err)

    @time result1 = optimize(xpoint -> draw_from_posterior_given_precomputed_kernel_from_data(reshape(xpoint, (1,dims)), xdata, ydata, L, kernel, hparams)[1][1], xmin_guess, BFGS()) # unconstrained, no gradient

    #@time result2 = optimize(xpoint -> draw_from_posterior_given_precomputed_kernel_from_data(reshape(xpoint, (1,dims)), xdata, ydata, L, kernel, hparams)[1][1], xdata_lower, xdata_upper, xmin_guess, Fminbox(GradientDescent())) # constrained, no gradient

    println("Best (xpoint, GP_mean): ", (Optim.minimizer(result1), Optim.minimum(result1)))
    #println("Best (xpoint, GP_mean): ", (Optim.minimizer(result2), Optim.minimum(result2)))
    result = result1
    xmin, ymin = Optim.minimizer(result), Optim.minimum(result)
    return xmin, ymin
end

"""
Optimize the parameters of the model repeatedly (i.e. find the point that minimizes the function) by using a GP model (i.e. a given kernel and set of hyperparameters) and different iterations of the training data to minimize the GP mean prediction.
"""
function optimize_model_GP_mean_multiple_datasets(iterations::Int64, n_train::Int64, params_names::Array{Symbol,1}, xmin_guess::Array{Float64,2}, xdata_all::Array{Float64,2}, ydata_all::Vector{Float64}, mean_f::Float64, kernel::Function, hparams::Vector{Float64}; ydata_err_all::Vector{Float64}=zeros(length(ydata_all)), save_path::String="")
    @assert size(xdata_all, 1) == length(ydata_all) == length(ydata_err_all)
    @assert n_train < length(ydata_all)
    dims = length(xmin_guess)

    min_points_table = Array{Float64,2}(undef, iterations, dims+1)
    @time for i in 1:iterations
        i_train = Random.randperm(length(ydata_all))[1:n_train]
        xtrain = xdata_all[i_train,:]
        ytrain = ydata_all[i_train] .- mean_f
        ytrain_err = ydata_err_all[i_train]
        min_points_table[i,1:dims], min_points_table[i,dims+1] = optimize_model_GP_mean(xmin_guess, xtrain, ytrain, kernel, hparams; ydata_err=ytrain_err)
    end

    min_points_table = DataFrame(min_points_table, [params_names; :ymin])

    if save_path != ""
        file_name = joinpath(save_path, "GP_train$(n_train)_meanf$(mean_f)_sigmaf$(hparams[1])_lscales$(round(sum(hparams[2:end]), digits=2))_minimize_mean_iterations$(iterations).csv")
        f = open(file_name, "w")
        println(f, "# n_data = $n_train, mean_f = $mean_f")
        println(f, "# hparams = $hparams")
        CSV.write(f, min_points_table; append=true, writeheader=true)
        close(f)
    end

    return min_points_table
end





data = load_data(dims, data_path)

# To run an optimizer using the best GP model in order to find the point that minimizes the function (i.e. GP mean prediction):

xmin_guess = Statistics.median(data[:params_array], dims=1)
println("Median data point: ", xmin_guess)

# To optimize for the minimum once using a dataset:
xmin, ymin = optimize_model_GP_mean(xmin_guess, GP_model[:xtrain], GP_model[:ytrain], GP_model[:kernel], GP_model[:hparams_best]; ydata_err=GP_model[:ytrain_err])

# To optimize for the minimum multiple times using random iterations of the dataset:

iterations = 100
n_train = 2000
save_path = data_path
min_table = optimize_model_GP_mean_multiple_datasets(iterations, n_train, data[:params_names], xmin_guess, data[:params_array], data[:dist_array], GP_model[:mean_f], GP_model[:kernel], GP_model[:hparams_best]; ydata_err_all=0.8 .*ones(length(data[:dist_array])), save_path=save_path)
