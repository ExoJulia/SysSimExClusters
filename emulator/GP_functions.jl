using LinearAlgebra
using Random
using DataFrames
using CSV
using Sobol
using Statistics
using Distributed





##### To define various GP kernels:

"""
Computes the squared exponential kernel in one dimension.
Assumes the first element in 'hparams' is the amplitude and the second element is the length scale.
"""
function kernel_SE_1d(xrows::Vector{Float64}, xcols::Vector{Float64}, hparams::Vector{Float64})
    sigma_f, lscale = hparams
    sqdist = xrows.^2 .+ transpose(xcols.^2) .- 2*(xrows .* transpose(xcols))
    return (sigma_f^2).*exp.((-1/(2*lscale^2)) .* sqdist)
end

"""
Computes the squared exponential kernel in any number of dimensions, allowing for a different length scale for each dimension (i.e. it is a product of multiple 1d SE kernels).
Assumes the first element in 'hparams' is the amplitude and the rest are the length scales.
"""
function kernel_SE_ndims(xpoints1::Array{Float64,2}, xpoints2::Array{Float64,2}, hparams::Vector{Float64})
    sigma_f, lscales = hparams[1], hparams[2:end]
    @assert length(lscales) == size(xpoints1, 2) == size(xpoints2, 2)

    sqdist_normed = zeros(size(xpoints1, 1), size(xpoints2, 1))
    for i in 1:length(lscales)
        sqdist_normed += (xpoints1[:,i].^2 .+ transpose(xpoints2[:,i].^2) .- 2*(xpoints1[:,i] .* transpose(xpoints2[:,i])))/(lscales[i]^2)
    end
    return (sigma_f^2).*exp.(-0.5 .* sqdist_normed)
end





##### To define functions that make random draws from a given kernel (with and without data):

"""
Returns a number of draws from a GP with a given kernel as a 2-d array (where each column is a single draw), at the points provided.
"""
function draw_from_prior_given_kernel(xpoints::Array{Float64,2}, kernel::Function, hparams::Vector{Float64}; diag_noise::Float64=1e-5, draws::Integer=1)

    K_ss = kernel(xpoints, xpoints, hparams)
    L = transpose(cholesky(K_ss + diag_noise*I).U)
    f_prior = L * randn(size(xpoints, 1), draws)

    return f_prior
end

"""
Computes the GP model given a kernel and a set of data, by computing the mean and standard deviation of the prediction at the points provided, and also making a number of draws from the GP posterior.
"""
function draw_from_posterior_given_kernel_and_data(xpoints::Array{Float64,2}, xdata::Array{Float64,2}, ydata::Vector{Float64}, kernel::Function, hparams::Vector{Float64}; ydata_err::Vector{Float64}=zeros(length(ydata)), diag_noise::Float64=1e-5, draws::Integer=1)
    @assert size(xdata, 1) == length(ydata) == length(ydata_err)

    K_ss = kernel(xpoints, xpoints, hparams)
    K_s = kernel(xdata, xpoints, hparams)
    K = kernel(xdata, xdata, hparams)
    var_I = zeros(size(K))
    var_I[diagind(var_I)] = ydata_err.^2
    L = transpose(cholesky(K + var_I).U)

    # To compute the mean at our test points 'xpoints':
    Lk = L \ K_s
    mu = transpose(Lk) * (L \ ydata)

    # To compute the std at our test points 'xpoints':
    s2 = diag(K_ss) .- transpose(sum(Lk.^2, dims=1))
    stdv = sqrt.(s2)

    # To draw samples from the posterior at our test points 'xpoints':
    L = transpose(cholesky(K_ss - (transpose(Lk) * Lk) + diag_noise*I).U) # NOTE: is it OK to add a diagonal term here to help avoid PosDefException error?
    f_posterior = mu .+ (L * randn(size(xpoints, 1), draws))

    # Returns the mean, std, and draws from the posterior at 'xpoints':
    return mu, stdv, f_posterior
end

"""
Computes the cholesky decomposition (triangular) matrix of a kernel given a dataset.
"""
function compute_kernel_given_data(xdata::Array{Float64,2}, kernel::Function, hparams::Vector{Float64}; ydata_err::Vector{Float64}=zeros(length(ydata)))
    @assert size(xdata, 1) == length(ydata_err)

    K = kernel(xdata, xdata, hparams)
    var_I = zeros(size(K))
    var_I[diagind(var_I)] = ydata_err.^2
    L = transpose(cholesky(K + var_I).U)
    return L
end

"""
Does the same thing as the function 'draw_from_posterior_given_kernel_and_data', but takes in a precomputed cholesky decomposition matrix 'L' to save time.
"""
function draw_from_posterior_given_precomputed_kernel_from_data(xpoints::Array{Float64,2}, xdata::Array{Float64,2}, ydata::Vector{Float64}, L::Transpose{Float64,UpperTriangular{Float64,Array{Float64,2}}}, kernel::Function, hparams::Vector{Float64}; diag_noise::Float64=1e-5, draws::Integer=1)
    @assert size(xdata, 1) == length(ydata)

    K_ss = kernel(xpoints, xpoints, hparams)
    K_s = kernel(xdata, xpoints, hparams)

    # To compute the mean at our test points 'xpoints':
    Lk = L \ K_s
    mu = transpose(Lk) * (L \ ydata)

    # To compute the std at our test points 'xpoints':
    s2 = diag(K_ss) .- transpose(sum(Lk.^2, dims=1))
    stdv = sqrt.(s2)

    # To draw samples from the posterior at our test points 'xpoints':
    L = transpose(cholesky(K_ss - (transpose(Lk) * Lk) + diag_noise*I).U) # NOTE: is it OK to add a diagonal term here to help avoid PosDefException error?
    f_posterior = mu .+ (L * randn(size(xpoints, 1), draws))

    # Returns the mean, std, and draws from the posterior at 'xpoints':
    return mu, stdv, f_posterior
end

"""
Uses a GP model (i.e. a given kernel and set of hyperparameters), an estimate for the best-fit point 'xmin', and a set of data points, to compute the mean and standard deviation along a 2-d grid of parameters for the i-th and j-th parameters (holding all the other parameters fixed at 'xmin').
"""
function predict_mean_std_on_2d_grid(i::Int64, j::Int64, xi_axis::AbstractRange, xj_axis::AbstractRange, xmin::Vector{Float64}, xdata::Array{Float64,2}, ydata::Vector{Float64}, L::Transpose{Float64,UpperTriangular{Float64,Array{Float64,2}}}, kernel::Function, hparams::Vector{Float64})
    @assert size(xdata, 2) == length(xmin)
    @assert size(xdata, 1) == length(ydata)
    dims = length(xmin)
    @assert 1 <= i <= dims
    @assert 1 <= j <= dims

    mean_2d_grid = Array{Float64,2}(undef, length(xi_axis), length(xj_axis))
    std_2d_grid = Array{Float64,2}(undef, length(xi_axis), length(xj_axis))
    for (k,xi) in enumerate(xi_axis)
        for (l,xj) in enumerate(xj_axis)
            xpoint = copy(xmin)
            xpoint[[i,j]] = [xi, xj]
            GP_result = draw_from_posterior_given_precomputed_kernel_from_data(reshape(xpoint, (1,dims)), xdata, ydata, L, kernel, hparams)
            mean_2d_grid[k,l] = GP_result[1][1]
            std_2d_grid[k,l] = GP_result[2][1]
        end
    end

    return mean_2d_grid, std_2d_grid
end

"""
Uses a GP model (i.e. a given kernel and set of hyperparameters), an estimate for the best-fit point 'xmin', and a set of data points, to compute the mean and standard deviation along 2-d grids of parameters for all possible combinations of the parameters (holding all the other parameters fixed at 'xmin').
Also computes the function 'f_true' on the same 2-d grids, if provided.
"""
function predict_mean_std_on_2d_grids_all_combinations(params_names, xmin::Vector{Float64}, xdata::Array{Float64,2}, ydata::Vector{Float64}, kernel::Function, hparams::Vector{Float64}, ydata_err::Vector{Float64}; grid_dims::Int64=20, f_x::Union{Function, Nothing}=nothing)
    @assert size(xdata, 2) == length(xmin) == length(params_names)
    @assert size(xdata, 1) == length(ydata) == length(ydata_err)
    dims = length(params_names)
    pairs = sum(1:dims-1) # number of pairs of parameters; equivalent to dims choose 2

    xtrain_mins, xtrain_maxs = findmin(xdata, dims=1)[1], findmax(xdata, dims=1)[1]

    L = compute_kernel_given_data(xdata, kernel, hparams; ydata_err=ydata_err)

    if f_x !== nothing
        f_true_2d_grids = Array{Float64,3}(undef, grid_dims, grid_dims, pairs)
    end
    mean_2d_grids = Array{Float64,3}(undef, grid_dims, grid_dims, pairs)
    std_2d_grids = Array{Float64,3}(undef, grid_dims, grid_dims, pairs)
    grid_count = 1 # To index how many 2d grids are computed
    for i in 1:dims
        for j in 1:i-1
            xi_axis = range(xtrain_mins[i], stop=xtrain_maxs[i], length=grid_dims)
            xj_axis = range(xtrain_mins[j], stop=xtrain_maxs[j], length=grid_dims)
            if f_x !== nothing
                f_true_2d_grid = Array{Float64,2}(undef, grid_dims, grid_dims)
                for (k,xi) in enumerate(xi_axis)
                    for (l,xj) in enumerate(xj_axis)
                        xpoint = copy(xmin)
                        xpoint[[i,j]] = [xi, xj]
                        f_true_2d_grid[k,l] = f_x(reshape(xpoint, (1,dims)))[1]
                    end
                end
                f_true_2d_grids[:,:,grid_count] = f_true_2d_grid
            end
            mean_2d_grids[:,:,grid_count], std_2d_grids[:,:,grid_count] = predict_mean_std_on_2d_grid(i, j, xi_axis, xj_axis, xmin, xdata, ydata, L, kernel, hparams)
            grid_count += 1
        end
    end

    grids_dict = Dict{Symbol, DataFrame}()
    grids_dict[:mean_grids] = DataFrame(vcat([mean_2d_grids[:,:,i] for i in 1:pairs]...))
    grids_dict[:std_grids] = DataFrame(vcat([std_2d_grids[:,:,i] for i in 1:pairs]...))
    if f_x !== nothing
        grids_dict[:f_x_grids] = DataFrame(vcat([f_true_2d_grids[:,:,i] for i in 1:pairs]...))
    end

    return grids_dict
end

"""
Uses a GP model (i.e. a given kernel and set of hyperparameters) and a set of data points, to compute the mean, standard deviation, and a draw from the posterior at each point drawn from a uniform prior for 'n_draws' points.
"""
function predict_model_at_n_points_from_uniform_prior(params_names::Array{Symbol,1}, xdata::Array{Float64,2}, ydata::Vector{Float64}, kernel::Function, hparams::Vector{Float64}, ydata_err::Vector{Float64}, n_draws::Int64)
    @assert size(xdata, 2) == length(params_names)
    @assert size(xdata, 1) == length(ydata) == length(ydata_err)
    dims = length(params_names)

    xtrain_mins, xtrain_maxs = findmin(xdata, dims=1)[1], findmax(xdata, dims=1)[1]
    prior_bounds = [(xtrain_mins[i], xtrain_maxs[i]) for i in 1:dims] # To prevent predicting "outside" of the n-dim box of training points

    L = compute_kernel_given_data(xdata, kernel, hparams; ydata_err=ydata_err)

    prior_draws_GP_table = Array{Float64,2}(undef, n_draws, dims+3)
    for i in 1:n_draws
        prior_draws_GP_table[i, 1:dims] = map(j -> prior_bounds[j][1] + (prior_bounds[j][2] - prior_bounds[j][1])*rand(), 1:dims)
        GP_result = draw_from_posterior_given_precomputed_kernel_from_data(reshape(prior_draws_GP_table[i, 1:dims], (1,dims)), xdata, ydata, L, kernel, hparams)
        prior_draws_GP_table[i, dims+1:dims+3] = [GP_result[j][1] for j in 1:3]
    end

    return DataFrame(prior_draws_GP_table, [params_names; [:GP_mean, :GP_std, :GP_posterior_draw]])
end

"""
Uses a GP model (i.e. a given kernel and set of hyperparameters) with a pre-computed cholesky decomposition matrix 'L' and a set of data points, to compute the mean, standard deviation, and a draw from the posterior at points drawn from a uniform prior until a point passing the criteria for the mean and std is accepted.
"""
function predict_model_from_uniform_prior_until_accept_point(prior_bounds::Array{Tuple{Float64,Float64},1}, xdata::Array{Float64,2}, ydata::Vector{Float64}, kernel::Function, hparams::Vector{Float64}, L::Transpose{Float64,UpperTriangular{Float64,Array{Float64,2}}}, ydata_err::Vector{Float64}; n_accept::Int64=1, max_mean::Float64=Inf, max_std::Float64=Inf, max_post::Float64=Inf)
    @assert size(xdata, 1) == length(ydata) == length(ydata_err)
    dims = size(xdata, 2)

    points_accepted_GP = Array{Float64,2}(undef, n_accept, dims+3)
    count_accepted = 0
    count_draws = 0
    while count_accepted < n_accept
        prior_draw = map(j -> prior_bounds[j][1] + (prior_bounds[j][2] - prior_bounds[j][1])*rand(), 1:dims)
        GP_result = draw_from_posterior_given_precomputed_kernel_from_data(reshape(prior_draw, (1,dims)), xdata, ydata, L, kernel, hparams)
        count_draws += 1
        if GP_result[1][1] < max_mean && GP_result[2][1] < max_std && GP_result[3][1] < max_post
            count_accepted += 1
            points_accepted_GP[count_accepted,:] = [prior_draw; [GP_result[j][1] for j in 1:3]]
            #prior_draws_GP_table[count_accepted, 1:end] = [prior_draw; [GP_result[j][1] for j in 1:3]]
        end
    end

    return points_accepted_GP, count_draws
end

"""
Evaluates 'predict_model_from_uniform_prior_until_accept_point' until 'n_accept' total points are accepted (in series).
"""
function predict_model_from_uniform_prior_until_accept_n_points(params_names::Array{Symbol,1}, xdata::Array{Float64,2}, ydata::Vector{Float64}, kernel::Function, hparams::Vector{Float64}, ydata_err::Vector{Float64}, n_accept::Int64; max_mean::Float64=Inf, max_std::Float64=Inf, max_post::Float64=Inf)
    @assert size(xdata, 2) == length(params_names)
    @assert size(xdata, 1) == length(ydata) == length(ydata_err)
    @assert n_accept > 0
    dims = length(params_names)

    xtrain_mins, xtrain_maxs = findmin(xdata, dims=1)[1], findmax(xdata, dims=1)[1]
    prior_bounds = [(xtrain_mins[i], xtrain_maxs[i]) for i in 1:dims] # To prevent predicting "outside" of the n-dim box of training points

    L = compute_kernel_given_data(xdata, kernel, hparams; ydata_err=ydata_err)

    prior_draws_GP_table = Array{Float64,2}(undef, n_accept, dims+3)
    count_draws = 0
    for i in 1:n_accept
        prior_draws_GP_table[i, 1:end], counts = predict_model_from_uniform_prior_until_accept_point(prior_bounds, xdata, ydata, kernel, hparams, L, ydata_err; max_mean=max_mean, max_std=max_std, max_post=max_post)
        count_draws += counts

        if i == 10
            println("First ", i, " points accepted after ", count_draws, " draws...")
        end
    end

    println("Accepted ", n_accept, " points after ", count_draws, " draws.")

    return DataFrame(prior_draws_GP_table, [params_names; [:GP_mean, :GP_std, :GP_posterior_draw]])
end

"""
Evaluates 'predict_model_from_uniform_prior_until_accept_point' until 'n_accept' total points are accepted (in parallel).
"""
function predict_model_from_uniform_prior_until_accept_n_points_parallel(params_names::Array{Symbol,1}, xdata::Array{Float64,2}, ydata::Vector{Float64}, kernel::Function, hparams::Vector{Float64}, ydata_err::Vector{Float64}, n_accept::Int64; n_batch::Int64=n_accept, max_mean::Float64=Inf, max_std::Float64=Inf, max_post::Float64=Inf)
    @assert size(xdata, 2) == length(params_names)
    @assert size(xdata, 1) == length(ydata) == length(ydata_err)
    @assert rem(n_accept, n_batch) == 0 # 'n_batch' is the number of batches, NOT the number of points per batch
    dims = length(params_names)
    n_per_batch = div(n_accept, n_batch)

    xtrain_mins, xtrain_maxs = findmin(xdata, dims=1)[1], findmax(xdata, dims=1)[1]
    prior_bounds = [(xtrain_mins[i], xtrain_maxs[i]) for i in 1:dims] # To prevent predicting "outside" of the n-dim box of training points

    L = compute_kernel_given_data(xdata, kernel, hparams; ydata_err=ydata_err)

    println("Number of procs: ", nprocs())
    prior_draws_GP_table = SharedArray{Float64,2}(n_accept, dims+3)
    draws_per_accept = SharedArray{Int64,1}(n_batch)
    @sync @distributed for i in 1:n_batch
        prior_draws_GP_table[1+(i-1)*n_per_batch:i*n_per_batch,:], draws_per_accept[i] = predict_model_from_uniform_prior_until_accept_point(prior_bounds, xdata, ydata, kernel, hparams, L, ydata_err; n_accept=n_per_batch, max_mean=max_mean, max_std=max_std, max_post=max_post)
    end
    count_draws = sum(draws_per_accept)

    println("Accepted ", n_accept, " points after ", count_draws, " draws.")

    return DataFrame(prior_draws_GP_table, [params_names; [:GP_mean, :GP_std, :GP_posterior_draw]])
end





##### To define a function for calculating the log-marginal likelihood given the data, a kernel and its hyperparameters:

"""
Computes the log-marginal likelihood of the GP model given the data, a kernel, and a set of hyperparameters.
"""
function log_marginal_likelihood(xdata::Array{Float64,2}, ydata::Vector{Float64}, kernel::Function, hparams::Vector{Float64}; ydata_err::Vector{Float64}=zeros(length(ydata)))
    @assert size(xdata, 1) == length(ydata) == length(ydata_err)

    K_f = kernel(xdata, xdata, hparams)
    var_I = zeros(size(K_f))
    var_I[diagind(var_I)] = ydata_err.^2
    K_y = K_f + var_I
    L = transpose(cholesky(K_y).U)

    return -(1/2)*transpose(L \ ydata)*(L \ ydata) - logdet(transpose(L)) - (length(ydata)/2)*log(2*pi) # Note: det(L) == det(transpose(L)) in general
end

"""
Optimize the hyperparameters of a given kernel using maximum likelihood estimation, given a training dataset and an initial guess for the hyperparameters.
"""
function optimize_hparams_with_MLE(hparams_guess::Vector{Float64}, xdata::Array{Float64,2}, ydata::Vector{Float64}, kernel::Function; ydata_err::Vector{Float64}=zeros(length(ydata)), hparams_lower::Vector{Float64}=ones(length(hparams_guess)).*1e-2, hparams_upper::Vector{Float64}=ones(length(hparams_guess)).*1e2)
    @assert size(xdata, 1) == length(ydata) == length(ydata_err)
    @assert length(hparams_guess) == length(hparams_lower) == length(hparams_upper)

    #@time result1 = optimize(hparams -> -log_marginal_likelihood(xdata, ydata, kernel, hparams; ydata_err=ydata_err), hparams_guess, BFGS()) #unconstrained, no gradient

    @time result2 = optimize(hparams -> -log_marginal_likelihood(xdata, ydata, kernel, hparams; ydata_err=ydata_err), hparams -> -gradient_log_marginal_likelihood(xdata, ydata, kernel, hparams; ydata_err=ydata_err), hparams_guess, BFGS(); inplace = false) #unconstrained, with gradient

    #@time result3 = optimize(hparams -> -log_marginal_likelihood(xdata, ydata, kernel, hparams; ydata_err=ydata_err), hparams_lower, hparams_upper, hparams_guess, Fminbox(GradientDescent())) #constrained, no gradient

    #@time result4 = optimize(hparams -> -log_marginal_likelihood(xdata, ydata, kernel, hparams; ydata_err=ydata_err), hparams -> -gradient_log_marginal_likelihood(xdata, ydata, kernel, hparams; ydata_err=ydata_err), hparams_lower, hparams_upper, hparams_guess, Fminbox(GradientDescent()); inplace = false) #constrained, with gradient

    #println("Best (hparams, log_p): ", (Optim.minimizer(result1), -Optim.minimum(result1)))
    println("Best (hparams, log_p): ", (Optim.minimizer(result2), -Optim.minimum(result2)))
    #println("Best (hparams, log_p): ", (Optim.minimizer(result3), -Optim.minimum(result3)))
    #println("Best (hparams, log_p): ", (Optim.minimizer(result4), -Optim.minimum(result4)))

    result = result2
    hparams_best = Optim.minimizer(result)
    log_p_best = -Optim.minimum(result)
    return hparams_best, log_p_best
end

"""
Ad hoc metric for finding the best hyperparameters for a kernel: summing the std's of the GP predictions at the cross-validation points (idea is to minimize this sum so that the prediction is as accurate as possible).
Computes the sum of the std's of the GP model at a given set of cross-validation points, given the data, a kernel, and a set of hyperparameters.
TODO: maybe should sum the differences in the GP mean predictions as well as the std's  at the cross-validation points?
"""
function sum_var_prediction_at_CV_points(xdata::Array{Float64,2}, xcheck::Array{Float64,2}, ydata::Vector{Float64}, kernel::Function, hparams::Vector{Float64}; ydata_err::Vector{Float64}=zeros(length(ydata)))
    @assert size(xdata, 1) == length(ydata) == length(ydata_err)

    K_ss = kernel(xcheck, xcheck, hparams)
    K_s = kernel(xdata, xcheck, hparams)
    K = kernel(xdata, xdata, hparams)
    var_I = zeros(size(K))
    var_I[diagind(var_I)] = ydata_err.^2
    L = transpose(cholesky(K + var_I).U)

    # To compute the mean at our test points 'xpoints':
    Lk = L \ K_s
    mu = transpose(Lk) * (L \ ydata)

    # To compute the std at our test points 'xpoints':
    s2 = diag(K_ss) .- transpose(sum(Lk.^2, dims=1))
    stdv = sqrt.(s2)

    return sum(s2)
end

"""
Optimize the hyperparameters of a given kernel using our ad hoc metric of minimizing the total variance at the cross-validiation points, given a training dataset and an initial guess for the hyperparameters.
NOTE: this does not really work.
"""
function optimize_hparams_with_minimize_GP_var(xdata::Array{Float64,2}, xcheck::Array{Float64,2}, ydata::Vector{Float64}, kernel::Function; ydata_err::Vector{Float64}=zeros(length(ydata)))

    @time result = optimize(hparams -> sum_var_prediction_at_CV_points(xdata, xcheck, ydata, kernel, [1.; hparams]; ydata_err=ydata_err), ones(dims), BFGS()) #unconstrained, no gradient
    hparams_best = [1; Optim.minimizer(result)]
    sum_var_best = Optim.minimum(result)
    println("Best (hparams, sum_std): ", (hparams_best, sum_var_best))
    return hparams_best, sum_var_best
end





##### To define functions for calculating the partial derivatives of various kernels with respect to the hyperparameters, as well as the gradient for the log-marginal likelihood as a function of the hyperparameters:

"""
Computes the matrix of partial derivatives of the squared exponential kernel with respect to the amplitude ('sigma_f') hyperparameter.
Assumes the first element in 'hparams' is the amplitude and the rest are the length scales.
"""
function dK_dsigma_kernel_SE_ndims(xpoints1::Array{Float64,2}, xpoints2::Array{Float64,2}, hparams::Vector{Float64})
    sigma_f, lscales = hparams[1], hparams[2:end]
    @assert length(lscales) == size(xpoints1, 2) == size(xpoints2, 2)

    return (2/sigma_f).*kernel_SE_ndims(xpoints1, xpoints2, hparams)
end

"""
Computes the matrix of partial derivatives of the squared exponential kernel with respect to the j-th dimension length-scale ('lscales[j]') hyperparameter.
Assumes the first element in 'hparams' is the amplitude and the rest are the length scales.
"""
function dK_dlj_kernel_SE_ndims(xpoints1::Array{Float64,2}, xpoints2::Array{Float64,2}, hparams::Vector{Float64}, j::Integer)
    sigma_f, lscales = hparams[1], hparams[2:end]
    @assert length(lscales) == size(xpoints1, 2) == size(xpoints2, 2)
    @assert 1 <= j <= length(lscales)

    K = kernel_SE_ndims(xpoints1, xpoints2, hparams)
    r2_lj3 = (xpoints1[:,j].^2 .+ transpose(xpoints2[:,j].^2) .- 2*(xpoints1[:,j] .* transpose(xpoints2[:,j])))/(lscales[j]^3)
    return r2_lj3 .* K # this is an element-wise product of the two matrices!
end

"""
Computes the gradient of the log-marginal likelihood as a function of the hyperparameters, given the data, a kernel, its partial derivatives with respect to the hyperparameters, and a set of hyperparameters.
"""
function gradient_log_marginal_likelihood(xdata::Array{Float64,2}, ydata::Vector{Float64}, kernel::Function, hparams::Vector{Float64}; ydata_err::Vector{Float64}=zeros(length(ydata)))
    @assert size(xdata, 1) == length(ydata) == length(ydata_err)

    K_f = kernel(xdata, xdata, hparams)
    var_I = zeros(size(K_f))
    var_I[diagind(var_I)] = ydata_err.^2
    K_y = K_f + var_I
    L = transpose(cholesky(K_y).U)
    K_y_inv = transpose(L) \ (L \ I) # inverse of K_y using cholesky decomposition

    matrix = (K_y_inv*ydata)*transpose(K_y_inv*ydata) - K_y_inv
    if kernel == kernel_SE_ndims
        # Assumes the first element in 'hparams' is the amplitude and the rest are the length scales
        sigma_f, lscales = hparams[1], hparams[2:end]
        return [0.5 .* tr(matrix * dK_dsigma_kernel_SE_ndims(xdata, xdata, hparams)); map(i -> 0.5 .* tr(matrix * dK_dlj_kernel_SE_ndims(xdata, xdata, hparams, i)), 1:length(lscales))]
    else
        println("The partial derivatives for this kernel are unknown.")
        return NaN
    end
end





##### To define various functions for generating and splitting data:

"""
Sample a number of points 'n_points' from a Sobol (quasi-random) sequence in n-dimensions, in an n-dim hypercube with each dimension bounded by [a,b].
"""
function sample_points_a_to_b_ndims_Sobol(dims::Int64, a::Float64, b::Float64, n_points::Int64)
    @assert dims >= 1

    s = SobolSeq(dims)
    skip(s, n_points) # skip the first n_points to get better uniformity
    xdata = transpose(hcat([next!(s) for i = 1:n_data]...))
    xdata = a .+ (b - a).*xdata # Note: can only do this because we assume the same range for each dimension; need to re-scale each dimension separately if different
    return xdata
end

"""
Splits a given dataset randomly into two samples, a 'training' set with 'n_train' points and a 'cross-validation' set with the remaining number of points.
"""
function split_data_training_cv(xdata, ydata, ydata_err; n_train::Int64=div(size(xdata, 1), 2))
    @assert n_train < size(xdata, 1)
    dims = size(xdata, 2)

    randperm_data = Random.randperm(size(xdata, 1))
    i_train, i_check = randperm_data[1:n_train], randperm_data[n_train+1:end]
    xtrain, xcheck = xdata[i_train,1:dims], xdata[i_check,1:dims]
    ytrain, ycheck = ydata[i_train], ydata[i_check]
    ytrain_err, ycheck_err = ydata_err[i_train], ydata_err[i_check]
    return xtrain, xcheck, ytrain, ycheck, ytrain_err, ycheck_err
end

"""
Computes a simple quadratic function (centered at 'x0') in one dimension.
"""
function f_quadratic_1d(x, f0::Float64, x0::Float64)
    return f0 .+ (x .- x0).^2
end

"""
Computes a quadratic surface (centered at 'x0') in n-dimensions.
"""
function f_quadratic_gaussian_ndims(x, f0::Float64, x0::Vector{Float64}, cov)
    @assert length(x) == length(x0) == size(cov, 1) == size(cov, 2)
    return f0 + (x .- x0)'*cov*(x .- x0)
end
