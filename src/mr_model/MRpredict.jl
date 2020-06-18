if !@isdefined ExoplanetsSysSim
    using ExoplanetsSysSim
end
using CSV
using Distributions
using Roots #if we want to use the 'find_zero()' function
using Optim

dir_path_MR = dirname(@__FILE__)

##### This code implements the M-R prediction model of Ning, Wolfgang, and Ghosh (2018) (hereafter NWG-2018)
##### It is a translation of 'MRpredict.R' at 'https://github.com/Bo-Ning/Predicting-exoplanet-mass-and-radius-relationship/blob/master/MR-predict/MRpredict.R'

struct MR_param_Ning2018
    Mass_min::Float64 #log10 scale
    Mass_max::Float64 #log10 scale
    Radius_min::Float64 #log10 scale
    Radius_max::Float64 #log10 scale
    deg::Int64
    weights_mle::Vector{Float64}
    indices_keep::Vector{Int64}
    indices_div::Vector{Int64}
    indices_rem::Vector{Int64}
end

function pdfnorm_beta(x::Float64, x_obs::Float64, x_sd::Float64, x_max::Float64, x_min::Float64, shape1::Int64, shape2::Int64 ; log::Bool=true)
    if log == true
        norm_dist = Normal(10^x, x_sd)
    else
        norm_dist = Normal(x, x_sd)
    end
    #beta_dist = Beta(shape1, shape2)
    return pdf.(norm_dist, x_obs) .* pdf.(Beta(shape1, shape2), (x-x_min)/(x_max-x_min)) / (x_max-x_min)
end

function fn_for_integrate(x_obs::Float64, x_sd::Float64, deg::Int, degree::Int, x_max::Float64, x_min::Float64 ; log::Bool=false, abs_tol::Float64=1e-10)
    shape1, shape2 = degree, deg-degree+1
    return quadgk(x -> pdfnorm_beta(x, x_obs, x_sd, x_max, x_min, shape1, shape2; log=log), x_min, x_max, abstol=abs_tol)[1]
end

function conditional_density(y::Float64, y_max::Float64, y_min::Float64, x_max::Float64, x_min::Float64, deg::Int64, w_hat::Vector{Float64}, indices_keep::Vector{Int64}, indices_div::Vector{Int64}, indices_rem::Vector{Int64}; y_sd::Union{Nothing,Float64}=nothing, qtl::Vector{Float64}=[0.16, 0.84])
    deg_vec = 1:deg

    #To compute conditional mean, variance, quantile, distribution:
    y_beta_indv = zeros(length(indices_rem))
    if y_sd == nothing
        map!(degree -> pdf.(Beta(degree, deg-degree+1), (y-y_min)/(y_max-y_min)), y_beta_indv, deg_vec[indices_rem])
        y_beta_indv /= (y_max-y_min)
    else
        map!(degree -> fn_for_integrate(y, y_sd, deg, degree, y_max, y_min), y_beta_indv, deg_vec[indices_rem])
    end
    denominator = sum(w_hat .* y_beta_indv)

    #Mean:
    mean = sum(w_hat::Vector{Float64} .* (((deg_vec::UnitRange{Int64})/(deg+1)*(x_max-x_min) .+ x_min)[indices_div] .* y_beta_indv::Vector{Float64})::Vector{Float64}) / denominator::Float64

    #Variance:
    var = sum(w_hat::Vector{Float64} .* (((deg_vec::UnitRange{Int64}) .* (deg .- deg_vec .+ 1) / ((deg+1)^2*(deg+2))*(x_max-x_min)^2)[indices_div] .* y_beta_indv::Vector{Float64})::Vector{Float64}) / denominator::Float64

    #Quantile:
    function pbeta_conditional_density(x::Float64)::Float64
        function mix_density(j::Float64)::Float64
            @fastmath @inbounds @views return sum(w_hat::Vector{Float64} .* (map(degree -> cdf.(Beta(degree, deg::Int64-degree+1), (j-x_min)/(x_max-x_min)), (deg_vec::UnitRange{Int64})[indices_div]) .* y_beta_indv::Vector{Float64})) / denominator::Float64
        end

        return mix_density(x)
    end

    function conditional_quantile(q::Float64)::Float64
        function g(x::Float64)::Float64
            return pbeta_conditional_density(x) - q
        end
        function g2(x::Float64)::Float64
            return (pbeta_conditional_density(x) - q)^2
        end

        #return find_zero(g, (x_min, x_max), FalsePosition())
        #return find_zero(g, (x_min, x_max))::Float64 ###'find_zero' is from 'using Roots' (is this a good function to use?); Eric suggested brent's method in the Optim package

        result = optimize(g2, x_min, x_max, Brent(), rel_tol=1e-3, abs_tol=1e-3) #Brent(); GoldenSection()
        if Optim.converged(result)
            return Optim.minimizer(result)::Float64
        else
            return find_zero(g, (x_min, x_max))::Float64
        end
    end

    quantile = map(conditional_quantile, qtl)

    return (mean, var, quantile, denominator, y_beta_indv)
end

function conditional_density(y::Float64, y_max::Float64, y_min::Float64, x_max::Float64, x_min::Float64, deg::Int64, w_hat::Vector{Float64}, indices_keep::Vector{Int64}; y_sd::Union{Nothing,Float64}=nothing, qtl::Vector{Float64}=[0.16, 0.84])
    deg_vec = 1:deg

    #To compute conditional mean, variance, quantile, distribution:
    y_beta_indv = zeros(deg)
    if y_sd == nothing
        map!(degree -> pdf.(Beta(degree, deg-degree+1), (y-y_min)/(y_max-y_min)), y_beta_indv, deg_vec)
        y_beta_indv /= (y_max-y_min)
    else
        map!(degree -> fn_for_integrate(y, y_sd, deg, degree, y_max, y_min), y_beta_indv, deg_vec)
    end
    y_beta_pdf = kron(ones(deg), y_beta_indv)[indices_keep]
    denominator = sum(w_hat .* y_beta_pdf)

    #Mean:
    mean = sum(w_hat::Vector{Float64} .* kron(deg_vec::UnitRange{Int64}/(deg+1)*(x_max-x_min)+x_min, y_beta_indv::Vector{Float64})[indices_keep]) / denominator::Float64

    #Variance:
    var = sum(w_hat::Vector{Float64} .* kron(deg_vec::UnitRange{Int64} .* (deg-deg_vec+1) / ((deg+1)^2*(deg+2))*(x_max-x_min)^2, y_beta_indv::Vector{Float64})[indices_keep]) / denominator::Float64

    #Quantile:
    function pbeta_conditional_density(x::Float64)::Float64
        function mix_density(j::Float64)::Float64
            @fastmath @inbounds @views return sum(w_hat::Vector{Float64} .* kron(map(degree -> cdf.(Beta(degree, deg-degree+1), (j-x_min)/(x_max-x_min)), deg_vec::UnitRange{Int64}), y_beta_indv::Vector{Float64})[indices_keep]) / denominator::Float64
            #@fastmath @inbounds @views return sum(w_hat::Vector{Float64} .* (map(degree -> cdf.(Beta(degree, deg::Int64-degree+1), (j-x_min)/(x_max-x_min)), (deg_vec::UnitRange{Int64})[indices_div]) .* y_beta_indv::Vector{Float64})) / denominator::Float64
        end

        return mix_density(x)
    end

    function conditional_quantile(q::Float64)::Float64
        function g(x::Float64)::Float64
            return pbeta_conditional_density(x) - q
        end
        function g2(x::Float64)::Float64
            return (pbeta_conditional_density(x) - q)^2
        end

        #return find_zero(g, (x_min, x_max), FalsePosition())
        #return find_zero(g, (x_min, x_max))::Float64 ###'find_zero' is from 'using Roots' (is this a good function to use?); Eric suggested brent's method in the Optim package

        result = optimize(g2, x_min, x_max, Brent(), rel_tol=1e-3, abs_tol=1e-3) #Brent(); GoldenSection()
        if Optim.converged(result)
            return Optim.minimizer(result)::Float64
        else
            return find_zero(g, (x_min, x_max))::Float64
        end
    end

    quantile = map(conditional_quantile, qtl)

    return (mean, var, quantile, denominator, y_beta_indv)
end

function predict_mass_given_radius(Radius::Float64, param::MR_param_Ning2018 ; R_sigma::Union{Nothing,Float64}=nothing, posterior_sample::Bool=false, qtl::Vector{Float64}=[0.16, 0.84])
    #Convert data to log scale:
    l_radius = log10.(Radius)

    #The 'posterior_sample == false' condition can deal with two cases:
    #Case I: if input data do not have measurement errors
    #Case II: if input data have measurement error
    if posterior_sample == false
        predicted_value = conditional_density(l_radius, param.Radius_max, param.Radius_min, param.Mass_max, param.Mass_min, param.deg, param.weights_mle, param.indices_keep, param.indices_div, param.indices_rem; y_sd=R_sigma, qtl=qtl)
        ###predicted_value = conditional_density(l_radius, param.Radius_max, param.Radius_min, param.Mass_max, param.Mass_min, param.deg, param.weights_mle, param.indices_keep; y_sd=R_sigma, qtl=qtl)
        predicted_mean = predicted_value[1]
        predicted_quantiles = predicted_value[3]
    elseif posterior_sample == true
        #Case III: if the input are posterior samples
        radius_sample = log10(Radius)
        k = length(radius_sample)
        denominator_sample = zeros(k)
        mean_sample = zeros(k)
        y_beta_indv_sample = zeros(k, param.deg)

        #Calculate mean:
        for i in 1:k
            results = cond_density_estimation(y=radius_sample[i], y_max=param.Radius_max, y_min=param.Radius_min, x_max=param.Mass_max, x_min=param.Mass_min, deg=param.deg, w_hat=param.weights_mle, qtl=quantile, only_output_mean=true) ###what function is this??? Is it part of some CRAN package...
            mean_sample[i] = (results[1]) ###why are there parentheses here? It makes no difference...
            denominator_sample[i] = results[2]
            y_beta_indv_sample[i,:] = results[3:57]
        end
        predicted_mean = mean(mean_sample)

        #Calculate 16% and 84% quantiles:
        #Mixture of the CDF of k conditional densities
        function pbeta_conditional_density(x)
            function mix_density(j)
                deg_vec = 1:param.deg
                x_indv_cdf = map(degree -> cdf.(Beta(degree, param.deg-degree+1), (j-x_min)/(x_max-x_min)), deg_vec)
                quantile_numerator = zeros(k)
                p_beta_sample = zeros(k)
                for ii in 1:k
                    quantile_numerator[ii] = sum(param.weights_mle .* kron(x_indv_cdf, y_beta_indv_sample[ii,:]))
                    p_beta_sample[ii] = quantile_numerator[ii]/denominator_sample[ii]
                end
                return p_beta = mean(p_beta_sample)
            end

            return map(mix_density, x)
        end

        function mixture_conditional_quantile(q, x_min, x_max)
            function g(x)
                return pbeta_conditional_density(x) - q
            end
            return find_zero(g, (x_min, x_max)) ###see note above
            #return Optim.minimum(optimize(g, x_min, x_max))
        end

        predicted_quantiles = map(q -> mixture_conditional_quantile(q, param.Mass_min, param.Mass_max), qtl)
    end

    #Return the output:
    return (predicted_mean, predicted_quantiles)::Tuple{Float64, Vector{Float64}}
end

function draw_planet_mass_from_radius_Ning2018(Radius::Float64, param::MR_param_Ning2018)
    #This function takes in a Radius (in solar radii) and draws a mass (returning in solar masses) probabilistically
    Radius_in_earths = Radius/ExoplanetsSysSim.earth_radius
    log_Radius_in_earths = log10.(Radius_in_earths)
    @assert param.Radius_min < log_Radius_in_earths < param.Radius_max
    q::Float64 = rand()
    l_mass = predict_mass_given_radius(Radius_in_earths, param; qtl=[q])[2][1]
    return (10^l_mass)*ExoplanetsSysSim.earth_mass
end

function generate_planet_mass_from_radius_Ning2018(Radius::Float64, sim_param::SimParam)
    global MR_param
    return draw_planet_mass_from_radius_Ning2018(Radius, MR_param)
end

function interpolate_planet_mass_from_radius_quantile_Ning2018_table(log_Radius::Float64, quantile::Float64)
    #If using "GridInterpolations" package (very slow for repeated use):
    #=
    global grid, gridData
    x = [log_Radius, quantile]
    return interpolate(grid, gridData, x)
    =#

    #If using "Interpolations" package (fastest):
    global scaled_itp
    return scaled_itp(log_Radius, quantile)
end

function generate_planet_mass_from_radius_Ning2018_table(R::Float64, sim_param::SimParam)
    R_in_earths = R/ExoplanetsSysSim.earth_radius
    log_R_in_earths = log10(R_in_earths)
    q::Float64 = rand()
    log_M = interpolate_planet_mass_from_radius_quantile_Ning2018_table(log_R_in_earths, q)
    return (10^log_M)*ExoplanetsSysSim.earth_mass
end





##### To initialize the M-R model parameters:

weights_mle = CSV.read(joinpath(dir_path_MR, "weights.mle.csv"))[!,:x]
degrees = 55

N_keep = 25
indices_keep = sortperm(weights_mle, rev=true)[1:N_keep]
indices_div = 1 .+ div.(indices_keep .- 1, degrees)
indices_rem = 1 .+ rem.(indices_keep .- 1, degrees)
weights = weights_mle[indices_keep]
MR_param = MR_param_Ning2018(-1., 3.809597, -0.302, 1.357509, degrees, weights, indices_keep, indices_div, indices_rem) #(-1., 3.809597, -0.3, 1.357509, 55, weights_mle)





##### Examples:

#=
#Observation without measurement errors:

Radius = 5. #original scale, not log scale
predict_result = predict_mass_given_radius(Radius, MR_param)
println(predict_result)

#Observation with a measurement error:

Radius = 5. #original scale, not log scale
R_sigma = 0.1
predict_result = predict_mass_given_radius(Radius, MR_param; R_sigma=0.1)
println(predict_result)

#If want to change 16% and 84% quantiles to 5% and 95% quantiles:

Radius = 5. #original scale, not log scale
predict_result = predict_mass_given_radius(Radius, MR_param; qtl=[0.05, 0.95])
println(predict_result)
=#

#Input are posterior samples: ###currently broken because the function 'cond_density_estimation' is undefined
#=
Radius_dist = Normal(5, 0.5) #original scale, not log scale
Radius_samples = rand(Radius_dist, 100)
predict_result = predict_mass_given_radius(Radius_samples, MR_param; R_sigma=0.1, posterior_sample=true)
println(predict_result)
=#





##### For plotting and timing the model:
#=
#Pkg.add("PyPlot")
using PyPlot

N_keep_array = [3025, 19, 18]
colors = ["blue", "green", "red"]
Radii = 10 .^(range(MR_param.Radius_min+0.01, stop=MR_param.Radius_max-0.01, length=1000))

fig, ax = subplots()
for (i,N_keep) in enumerate(N_keep_array)
    indices_keep = sortperm(weights_mle, rev=true)[1:N_keep]
    indices_div = 1 .+ div.(indices_keep .- 1, degrees)
    indices_rem = 1 .+ rem.(indices_keep .- 1, degrees)
    weights = weights_mle[indices_keep]
    MR_param = MR_param_Ning2018(-1., 3.809597, -0.302, 1.357509, degrees, weights, indices_keep, indices_div, indices_rem)

    @time Masses_predict = map(r -> predict_mass_given_radius(r, MR_param; qtl=[0.16, 0.84]), Radii)

    Masses_mean = 10 .^[x[1] for x in Masses_predict]
    Masses_lower = 10 .^[x[2][1] for x in Masses_predict]
    Masses_upper = 10 .^[x[2][2] for x in Masses_predict]

    loglog(Radii, Masses_mean, color=colors[i], linewidth=2.0, linestyle="-", label="M = "*string(N_keep))
    loglog(Radii, Masses_lower, color=colors[i], linewidth=1.0, linestyle="--")
    loglog(Radii, Masses_upper, color=colors[i], linewidth=1.0, linestyle="--")
end
xlabel(L"$R_p (R_\oplus)$")
ylabel(L"$M_p (M_\oplus)$")
legend(loc="upper left")
=#





##### If we want to generate and save a pre-computed table of masses on a radius vs. quantile grid:
#=
N_keep = 3025
indices_keep = sortperm(weights_mle, rev=true)[1:N_keep]
indices_div = 1 .+ div.(indices_keep .- 1, degrees)
indices_rem = 1 .+ rem.(indices_keep .- 1, degrees)
weights = weights_mle[indices_keep]
MR_param = MR_param_Ning2018(-1., 3.809597, -0.302, 1.357509, degrees, weights, indices_keep, indices_div, indices_rem)

N_radii = 1001
N_quantiles = 1001
log_Radii = collect(range(MR_param.Radius_min, stop=MR_param.Radius_max, length=N_radii))
Radii = 10 .^log_Radii
quantiles = collect(range(0., stop=1.0, length=N_quantiles))

file_name = "MRpredict_table_weights"*string(N_keep)*"_R"*string(N_radii)*"_Q"*string(N_quantiles)*".txt"
f = open(file_name, "w")
println(f, "# All masses are in log10; weights used = ", N_keep)
println(f, "# First uncommented line is the header with the column labels (first column contains the radii in log10 while the remaining columns correspond to the quantiles)")
writedlm(f, reshape(append!(["log_R"], "q".*string.(quantiles)), (1,:)), ", ")

log_Mass_table = zeros(N_radii, N_quantiles + 1)
for (i,r) in enumerate(Radii)
    log_Mass_quantiles = predict_mass_given_radius(r, MR_param; qtl=quantiles)[2]
    log_Mass_table[i,:] = append!([log_Radii[i]], log_Mass_quantiles)
    println(i, ", Radius = ", r)
end
writedlm(f, log_Mass_table, ", ")
close(f)
=#





##### To load a pre-computed table for interpolating the mass radius relation:

log_Mass_table = CSV.read(joinpath(dir_path_MR, "MRpredict_table_weights3025_R1001_Q1001.txt"), header=3)

#One way to do the interpolation is to use "GridInterpolations" but this is very slow:
#=
#Pkg.add("GridInterpolations")
using GridInterpolations

N_quantiles = 1001
quantiles = collect(range(0., stop=1.0, length=N_quantiles))
grid = RectangleGrid(log_Mass_table[1:end,1], quantiles)
gridData = convert(Matrix, log_Mass_table[1:end,2:end])

#x = [0., 0.5] #an example data point to interpolate at
#interpolate(grid, gridData, x) #interpolate at the example data point
=#

#A faster way to do the interpolation is to first construct an interpolation object:

#Pkg.add("Interpolations") #NOTE: the README is somewhat out-dated so some of the examples/syntax do not actually work as the way they are shown
using Interpolations

N_radii = 1001
N_quantiles = 1001
log_Radii = range(MR_param.Radius_min, stop=MR_param.Radius_max, length=N_radii)
quantiles = range(0., stop=1.0, length=N_quantiles)

table_data = convert(Matrix, log_Mass_table[1:end,2:end])
itp = interpolate(table_data, BSpline(Cubic(Line(OnGrid())))) #interpolation object where the x and y axes are indices
scaled_itp = Interpolations.scale(itp, log_Radii, quantiles) #scaled interpolation object where the x and y axes are scaled to their physical units

#scaled_itp(0., 0.5) #calls the interpolation object to perform an interpolation at a given point





##### To draw from an alternate M-R model for small planets (e.g. Li Zeng's model for Earth-like rocky):
##### https://www.cfa.harvard.edu/~lzeng/planetmodels.html
##### Will assume a mean density drawn from a normal distribution centered the Earth-like rocky model, with a std that increases towards larger sizes (a linear function for log(σ_ρ) vs. R between the min and max radius for this M-R model)

# Define some useful constants and functions:
earth_radius_in_cm = ExoplanetsSysSim.earth_radius_eq_in_m_IAU2015 * 100.
earth_mass_in_g = (ExoplanetsSysSim.G_mass_earth_in_mks / ExoplanetsSysSim.G_in_mks_IAU2015) * 1e3

# Compute mass M (Earth masses) given radius R (Earth radii) and mean density ρ (g/cm^3)
mass_given_radius_density(R::Float64, ρ::Float64) = ρ * (4π/3) * (R*earth_radius_in_cm)^3 * (1/earth_mass_in_g)

# Compute mean density ρ (g/cm^3) given mass M (Earth masses) and radius R (Earth radii)
density_given_mass_radius(M::Float64, R::Float64) = (M*earth_mass_in_g) / ((4π/3)*(R*earth_radius_in_cm)^3)



# To load Li Zeng's M-R tables:
MR_earthlike_rocky = CSV.read(joinpath(dir_path_MR, "MR_earthlike_rocky.txt"), delim="\t", comment="#", header=["mass","radius"]) # Earth units

# To construct an interpolation function for planet mass as a function of radius:
using Dierckx
M_earthlike_rocky_spl = Spline1D(MR_earthlike_rocky[!,"radius"], MR_earthlike_rocky[!,"mass"])



# To define a function for computing the std in mean density as a function of radius:
# NOTE: while the following constants could be defined in sim_param, we define them here to prevent users from easily changing them
radius_switch = 1.472 # Earth radii; where the Earth-like rocky model intersects the mean prediction from NWG-2018
radius_min = 0.5 # Earth radii; hardcoded here in addition to "min_radius" in sim_param so changing that does not change the M-R model here
σ_ρ_at_radius_switch = 3. # g/cm^3; std of mean density at radius_switch
σ_ρ_at_radius_min = 1. # g/cm^3; std of mean density at radius_min

"""
    σ_ρ_loglinear_given_radius(R; σ_ρ_r1=σ_ρ_at_radius_min, σ_ρ_r2=σ_ρ_at_radius_switch, r1=radius_min, r2=radius_switch)

Compute the standard deviation in mean density (σ_ρ) at a given planet radius (R), assuming a linear relation between log(σ_ρ) and R defined by two points.

# Arguments:
- `R::Float64`: planet radius (Earth radii).
- `σ_ρ_r1::Float64=σ_ρ_at_radius_min`: standard deviation of mean density (g/cm^3) at radius `r1` (Earth radii).
- `σ_ρ_r2::Float64=σ_ρ_at_radius_switch`: standard deviation of mean density (g/cm^3) at radius `r2` (Earth radii).
- `r1::Float64=radius_min`: radius at one point (Earth radii).
- `r2::Float64=radius_switch`: radius at another point (Earth radii).

# Returns:
The standard deviation in mean density, σ_ρ (g/cm^3), at radius `R`.
"""
function σ_ρ_loglinear_given_radius(R::Float64; σ_ρ_r1::Float64=σ_ρ_at_radius_min, σ_ρ_r2::Float64=σ_ρ_at_radius_switch, r1::Float64=radius_min, r2::Float64=radius_switch)
    log_σ_ρ_r1, log_σ_ρ_r2 = log10(σ_ρ_r1), log10(σ_ρ_r2)
    slope = (log_σ_ρ_r2 - log_σ_ρ_r1) / (r2 - r1)
    log_σ_ρ = slope*(R - r1) + log_σ_ρ_r1
    return 10^log_σ_ρ
end



# To define functions for drawing from the new M-R relation for planets below "radius_switch":

"""
    generate_planet_mass_from_radius_normal_density_around_earthlike_rocky(R, sim_param; ρ_min=1., ρ_max=100.)

Draw a planet mass given a planet radius, from a normal distribution of mean density centered around the Earth-like rocky model.

# Arguments:
- `R::Float64`: planet radius (solar radii).
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.
- `ρ_min::Float64=1.`: minimum density (g/cm^3) to truncate at.
- `ρ_max::Float64=100.`: maximum density (g/cm^3) to truncate at.
NOTE: requires an interpolation object `M_earthlike_rocky_spl` to be defined as a global.

# Returns:
A planet mass (in solar masses).
"""
function generate_planet_mass_from_radius_normal_density_around_earthlike_rocky(R::Float64, sim_param::SimParam; ρ_min::Float64=1., ρ_max::Float64=100.)
    global M_earthlike_rocky_spl
    R_in_earths = R/ExoplanetsSysSim.earth_radius
    M_in_earths = M_earthlike_rocky_spl(R_in_earths)
    @assert(M_in_earths > 0)

    μ_ρ = density_given_mass_radius(M_in_earths, R_in_earths)
    σ_ρ = σ_ρ_loglinear_given_radius(R_in_earths)
    ρ_dist = Truncated(Normal(μ_ρ, σ_ρ), ρ_min, ρ_max) # truncated normal distribution for mean density (g/cm^3)

    ρ = rand(ρ_dist) # draw a mean density
    @assert(ρ > 0)
    return mass_given_radius_density(R_in_earths, ρ)*ExoplanetsSysSim.earth_mass
end

"""
    generate_planet_mass_from_radius_Ning2018_table_above_earthlike_rocky_below(R, sim_param; R_switch=radius_switch, ρ_min=1., ρ_max=100.)

Draw a planet mass given a planet radius, from interpolating the NWG-2018 table if R > R_switch, or from the normal distribution of mean density centered around the Earth-like rocky model if R <= R_switch.

# Arguments:
- `R::Float64`: planet radius (solar radii).
- `sim_param::SimParam`: a SimParam object containing various simulation parameters.
- `R_switch::Float64=radius_switch`: planet radius (Earth radii) where the M-R relation switches.
- `ρ_min::Float64=1.`: minimum density (g/cm^3) to truncate at (only for below R_switch).
- `ρ_max::Float64=100.`: maximum density (g/cm^3) to truncate at (only for below R_switch).
NOTE: requires an interpolation object `M_earthlike_rocky_spl` to be defined as a global.

# Returns:
A planet mass (in solar masses).
"""
function generate_planet_mass_from_radius_Ning2018_table_above_earthlike_rocky_below(R::Float64, sim_param::SimParam; R_switch::Float64=radius_switch, ρ_min::Float64=1., ρ_max::Float64=100.)
    @assert(0 < R)
    R_in_earths = R/ExoplanetsSysSim.earth_radius

    if R_in_earths > R_switch
        return generate_planet_mass_from_radius_Ning2018_table(R, sim_param::SimParam)
    else
        return generate_planet_mass_from_radius_normal_density_around_earthlike_rocky(R, sim_param; ρ_min=ρ_min, ρ_max=ρ_max)
    end
end





##### For timing the functions:

sim_param = setup_sim_param_model()
#Radii = ones(10000)*ExoplanetsSysSim.earth_radius
Radii = (10 .^(range(MR_param.Radius_min+0.01, stop=MR_param.Radius_max-0.01, length=10000)))*ExoplanetsSysSim.earth_radius
@time Masses = map(r -> generate_planet_mass_from_radius_Ning2018(r, sim_param), Radii)
@time Masses = map(r -> generate_planet_mass_from_radius_Ning2018_table(r, sim_param), Radii)
@time Masses = map(r -> generate_planet_mass_from_radius_normal_density_around_earthlike_rocky(r, sim_param), Radii)
@time Masses = map(r -> generate_planet_mass_from_radius_Ning2018_table_above_earthlike_rocky_below(r, sim_param), Radii)
