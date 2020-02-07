##### Functions for information theory metrics taken from or inspired by Gilbert & Fabrycky (2020):

"""
Compute the Shannon entropy.
"""
function Shannon_entropy(p::Vector{Float64})
    @assert all(p .>= 0)
    @assert all(p .<= 1)
    H = -sum(p .* log.(p))
    return H
end



"""
Compute the disequilibrium.
"""
function disequilibrium(p::Vector{Float64})
    @assert all(p .>= 0)
    @assert all(p .<= 1)
    D = sum((p .- (1/length(p))).^2)
    return D
end



"""
    LMC_complexity(K, p)

Compute the Lopez-Ruiz, Mancini, & Calbet (LMC) complexity (1995); product of Shannon entropy and disequilibrium.

# Arguments:
- `K::Float64}`: a normalization (positive, real) constant.
- `p::Vector{Float64}`: list of occupancy probabilities.

# Returns:
- `C::Float64`: the LMC complexity.
"""
function LMC_complexity(K::Float64, p::Vector{Float64})
    @assert K > 0
    H = Shannon_entropy(p)
    D = disequilibrium(p)
    C = K*H*D
    return C
end



"""
Compute the Pearson correlation coefficient between two variables.
"""
function Pearson_correlation_coefficient(x::Vector{Float64}, y::Vector{Float64})
    xmean, ymean = mean(x), mean(y)
    r_xy = sum((x .- xmean).*(y .- ymean)) / sqrt(sum((x .- xmean).^2)*sum((y .- ymean).^2))
    return r_xy
end



"""
Compute the Spearman correlation coefficient between two variables. This is the Pearson correlation of the ranks of the variables.
"""
function Spearman_correlation_coefficient(x::Vector{Float64}, y::Vector{Float64})
    xsort, ysort = sortperm(x), sortperm(y)
    xranks, yranks = zeros(length(x)), zeros(length(y))
    xranks[xsort], yranks[ysort] = 1:length(x), 1:length(y)
    rho_S = Pearson_correlation_coefficient(xranks, yranks)
    return rho_S
end





"""
    partitioning(x)

Compute the "partitioning" of quantity x. Analogous to the "mass partitioning" as defined in Gilbert & Fabrycky (2020).

# Arguments:
- `x::Vector{Float64}`: quantities (e.g. masses) in a given system.

# Returns:
- `Q::Float64`: the "partitioning" of quantity x, which is the disequilibrium normalized to (0,1).
"""
function partitioning(x::Vector{Float64})
    @assert all(x .>= 0)
    xnorm = x./sum(x)
    Q = (length(x)/(length(x)-1)) * disequilibrium(xnorm)
    return Q
end



"""
    monotonicity_GF2020(x)

Compute the "monotonicity" of quantity x. Analogous to the "monotonicity" as defined in Gilbert & Fabrycky (2020) for mass ordering.

# Arguments:
- `x::Vector{Float64}`: quantities (e.g. masses) in a given system.

# Returns:
- `M::Float64`: the "monotonicity" of quantity x, which is the Spearman correlation coefficient times a factor involving the "partitioning" of the same quantity (in order to capture the magnitude of the ordering trends).
"""
function monotonicity_GF2020(x::Vector{Float64})
    rho_S = Spearman_correlation_coefficient(convert(Vector{Float64},1:length(x)), x)
    Q = partitioning(x)
    M = rho_S*(Q^(1/length(x)))
    return M
end



"""
    gap_complexity_GF2020(P)

Compute the "gap complexity" given a set of periods as defined in Gilbert & Fabrycky (2020).

# Arguments:
- `P::Vector{Float64}`: periods in a system (do not have to be sorted).

# Returns:
- `C::Float64`: the gap complexity.
"""
function gap_complexity_GF2020(P::Vector{Float64})
    @assert length(P) >= 3
    n = length(P)-1

    P = sort(P)
    Pmin, Pmax = minimum(P), maximum(P)
    Pratios = P[2:end]./P[1:end-1]

    pnorm = log.(Pratios) ./ log(Pmax/Pmin) # NOTE: assuming GF2020 also used natural log
    Cmax = n < 10 ? Cmax_table_GF2020[n] : Cmax_approx_GF2020(n)
    K = 1/Cmax
    C = LMC_complexity(K, pnorm)
    return C
end


# Normalization constant (C_max = 1/K) for the gap complexity function:

Cmax_table_GF2020 = Dict{Int64,Float64}(2=>0.106, 3=>0.212, 4=>0.291, 5=>0.350, 6=>0.398, 7=>0.437, 8=>0.469, 9=>0.497) # Table 2 of Gilbert & Fabrycky (2020)
Cmax_approx_GF2020(n) = 0.262*log(0.766*n) # Eq. 15 of Gilbert & Fabrycky (2020)
