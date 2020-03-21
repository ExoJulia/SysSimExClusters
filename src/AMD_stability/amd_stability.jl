#=
Date:   June 15, 2019
Author: Daniel Carrera (dcarrera@gmail.com)
Modified by: Matthias Yang He

Based on:
Laskar & Petit (2017): https://arxiv.org/pdf/1703.07125.pdf
Petit, Laskar, & Boue (2017): https://www.aanda.org/articles/aa/pdf/2017/11/aa31196-17.pdf
=#

using Random

##################################################
#
#     AMD STABILITY
#
##################################################
function version_1(γ::Real,α::Real; tol::Real = 1e-15)
	@assert(α <= 1)
	#
	# F(e,γ,α) from Equation 35 of Laskar & Petit (2017).
	#
	F(e) = α*e + γ*e / sqrt(α*(1-e)*(1+e) + γ^2*e^2) - 1.0 + α
	#
	# dF/de from Equation 36 of Laskar & Petit (2017). dF/de > 0 always.
	#
	dF(e) = α + α*γ / (α*(1-e)*(1+e) + γ^2*e^2)^(3//2)
	#
	# Find the root of F() by Newton's method + bisection.
	#
	eL = 0.0
	eR = α < 0.5 ? 1.0 : 1/α - 1 # a1*(1 + e) = a2
	e0 = (eL + eR)/2 # Initial guess.
	#@info("Bisect: e0 = $e0")
	F0 = F(e0)
	if F0 >= 0
		eR = e0
		FR = F0
		ex = eL
		FL = F(eL)
	else
		eL = e0
		FL = F0
		ex = eR
		FR = F(eR)
	end
	@assert(FL <= 0 <= FR)
	
	while abs(ex - e0) > tol
		#
		# Try a Newton step first.
		#
		ex = e0 - F0 / dF(e0)
		if (ex <= eL) || (ex >= eR)
			#
			# Newton's method jumped out of the interval.
			# Switch to linear interpolation.
			#
			m = (FR - FL) / (eR - eL)
			b = F0 - m*e0
			ex = -b/m
			#@info("Linear: ex = $ex")
		else
			#@info("Newton: ex = $ex")
		end
		#
		# Swap and update.
		#
		e0,ex = ex,e0
		F0 = F(e0)
		if F0 >= 0
			eR = e0
			FR = F0
		else
			eL = e0
			FL = F0
		end
	end
	return ex
end

function version_2(γ::Real,α::Real; tol::Real = 1e-15)
	@assert(α <= 1)
	
	function foo(e)
		e2  = e*e
		tmp = 1 / sqrt(α*(1 - e2) + γ^2*e2)
		#
		# F(e,γ,α) from Equation 35 of Laskar & Petit (2017).
		#
		F = α*e + γ*e * tmp - 1.0 + α
		#
		# dF/de from Equation 36 of Laskar & Petit (2017). dF/de > 0 always.
		#
		dF = α + α*γ * tmp^3
		
		return F, dF
	end
	
	#
	# Initial box.
	#
	eL = 0.0
	eR = min(1.0, 1/α - 1) # a1*(1 + e) = a2
	yL, dyL = foo(eL)
	yR, dyR = foo(eR)
	@assert(yL <= 0 <= yR)
	#
	# Try Newton steps from the edges. Chances are, one of them is
	# way off and the other is a good initial guess.
	#
	e_l = eL - yL / dyL # Newton's method starting from eL.
	e_r = eR - yR / dyR # Newton's method starting from eR.
	
	#@info("Newton: e0 = $e_l")
	#@info("Newton: e0 = $e_r")
	
	if (eL < e_l < eR)
		ex = e_l
	elseif (eL < e_r < eR)
		ex = e_r
	else
		# That's weird... well. Use bisection then.
		ex = (eL + eR)/2
	end
	e0 = 0.0
	
	while abs(ex - e0) > tol
		#
		# Swap and update.
		#
		e0,ex = ex,e0
		y0, dy0 = foo(e0)
		if y0 >= 0
			eR = e0
			yR = y0
		else
			eL = e0
			yL = y0
		end
		#
		# Try a Newton step first.
		#
		ex = e0 - y0 / dy0
		if (ex < eL) || (ex > eR)
			#
			# Switch to linear interpolation.
			#
			m = (yR - yL) / (eR - eL)
			b = y0 - m*e0
			ex = -b/m
			#@info("Linear: ex = $ex")
		else
			#@info("Newton: ex = $ex")
		end
	end
	return ex
end

"""
Finds the root of F(e,γ,α) from Equation 35 of Laskar & Petit (2017).

INPUT:
		γ		m_in / m_out ; planet-planet mass ratio.
		α		a_in / a_out ; semimajor axis ratio
		[tol]	Tolenrance (default: 1e-15)

OUTPUT: Critical eccentricity of the inner planet.
"""
function critical_eccentricity(γ::Real,α::Real; tol::Real = 1e-15)
	@assert(α <= 1)
	#
	# F(e,γ,α) from Equation 35 of Laskar & Petit (2017).
	#
	F(e) = α*e + γ*e / sqrt(α*(1-e)*(1+e) + γ^2*e^2) - 1.0 + α
	#
	# dF/de from Equation 36 of Laskar & Petit (2017). dF/de > 0 always.
	#
	dF(e) = α + α*γ / (α*(1-e)*(1+e) + γ^2*e^2)^(3//2)
	#
	# Find the root of F() by Newton's method + bisection.
	#
	e0 = 0.5 # Initial guess.
	F0 = F(e0)
	if F0 >= 0
		eR = e0
		FR = F0
		eL = ex = 0.0
		FL = F(eL)
		@assert(FL <= 0)
	else
		eL = e0
		FL = F0
		eR = ex = 1.0
		FR = F(eR)
		@assert(FR >= 0)
	end
	
	while abs(ex - e0) > tol
		#
		# Try a Newton step first.
		#
		ex = e0 - F0 / dF(e0)
		if (ex <= eL) || (ex >= eR)
			#
			# Newton's method jumped out of the interval.
			# Switch to linear interpolation.
			#
			m = (FR - FL) / (eR - eL)
			b = F0 - m*e0
			ex = -b/m
		end
		#
		# Swap and update.
		#
		e0,ex = ex,e0
		F0 = F(e0)
		if F0 >= 0
			eR = e0
			FR = F0
		else
			eL = e0
			FL = F0
		end
	end
	return ex
end

"""
Determines the AMD stability of a pair of planets using the
collision condition and the MMR overlap conditions of Laskar
& Petit (2017) and Petit & Laskar (2017).

INPUT:
		μ1		Planet/star mass ratio of inner planet.
		μ2		Planet/star mass ratio of outer planet.
		a1		Semimajor axis of inner planet.
		a2		Semimajor axis of outer planet.
		Cx		Relative AMD (Eq. 29 of Laskar & Petit 2017).

OUTPUT 1:
		:Stable			AMD stable.
		:Collison		Fails collision condition
		:MMR_circular	Fails MMR overlap condition for circular orbits.
		:MMR_eccentric	Fails MMR overlap condition for eccentric orbits.

OUTPUT 2:
		Ratio of relative AMD vs relative AMD needed for instability.
		This is ∞ for :MMR_circular and Cx/min(C_coll,C_mmr) otherwise.
"""
function AMD_stability(μ1::Real,μ2::Real,a1::Real,a2::Real,Cx::Real)
	# Convention: 1 == inner planet; 2 == outer planet.
	@assert(a1 <= a2)

	γ = μ1/μ2
	α = a1/a2
	ϵ = μ1+μ2 # Equation 4 of Petit, Laskar, & Boue 2017

    # Equation 76 and Section 3.6 of Petit, Laskar, & Boue 2017:
	α_crit = 1 - 1.46*ϵ^(2/7)
	if α > α_crit
		return :MMR_circular, Inf
	end

    # Relative AMD from collision condition (Laskar & Petit 2017):
    C_coll = relative_AMD_collision(μ1, μ2, a1, a2)

    # Relative AMD from MMR overlap condition (Petit, Laskar, & Boue 2017):
	C_mmr = relative_AMD_MMR_overlap(μ1, μ2, a1, a2)

    # Final result:
    C_crit = min(C_coll, C_mmr)
	ratio = Cx/C_crit
	if Cx < C_crit
		return :Stable, ratio
	elseif C_coll < C_mmr
		return :Collision, ratio
	else
		return :MMR_eccentric, ratio
	end
end

"""
    relative_AMD_collision(μ1, μ2, a1, a2)

Compute the critical (minimum) relative AMD for collision, based on Equations 29 & 39 in Laskar & Petit (2017).

# Arguments:
- `μ1::Real`: planet/star mass ratio of inner planet.
- `μ2::Real`: planet/star mass ratio of outer planet.
- `a1::Real`: semimajor axis of inner planet.
- `a2::Real`: semimajor axis of outer planet.

# Returns:
- `C_coll::Float64`: critical relative AMD for collision.
"""
function relative_AMD_collision(μ1::Real, μ2::Real, a1::Real, a2::Real)
    @assert(a1 <= a2)

    γ = μ1/μ2
    α = a1/a2

    e1 = critical_eccentricity(γ,α)
    e2 = 1 - α - α*e1
    @assert(a1*(1 + e1) ≈ a2*(1 - e2))
    C_coll = γ*sqrt(α)*(1 - sqrt((1-e1)*(1+e1)) + 1 - sqrt((1-e2)*(1+e2)))
    return C_coll
end

"""
    relative_AMD_MMR_overlap(μ1, μ2, a1, a2)

Compute the critical (minimum) relative AMD for MMR overlap, based on Equation 74 in Petit, Laskar, & Boue (2017).

# Arguments:
- `μ1::Real`: planet/star mass ratio of inner planet.
- `μ2::Real`: planet/star mass ratio of outer planet.
- `a1::Real`: semimajor axis of inner planet.
- `a2::Real`: semimajor axis of outer planet.

# Returns:
- `C_mmr::Float64`: critical relative AMD for MMR overlap.
"""
function relative_AMD_MMR_overlap(μ1::Real, μ2::Real, a1::Real, a2::Real)
    @assert(a1 <= a2)

    γ = μ1/μ2
    α = a1/a2
    ϵ = μ1+μ2
    r = 0.80199 # Equation 28 of Petit, Laskar, & Boue (2017)

    g = ((3^4*(1-α)^5) / (2^9*r*ϵ)) - ((32*r*ϵ) / (9*(1-α)^2))
    C_mmr = (g^2*γ*sqrt(α)) / (2+2γ*sqrt(α))
    return C_mmr
end

"""
    critical_relative_AMD(μ1, μ2, a1, a2)

Compute the critical (minimum) relative AMD for collision or MMR overlap.

# Arguments:
- `μ1::Real`: planet/star mass ratio of inner planet.
- `μ2::Real`: planet/star mass ratio of outer planet.
- `a1::Real`: semimajor axis of inner planet.
- `a2::Real`: semimajor axis of outer planet.

# Returns:
The minimum of the relative AMD for collision and for MMR overlap.
"""
function critical_relative_AMD(μ1::Real, μ2::Real, a1::Real, a2::Real)
    @assert(a1 <= a2)

    C_coll = relative_AMD_collision(μ1, μ2, a1, a2)
    C_mmr = relative_AMD_MMR_overlap(μ1, μ2, a1, a2)
    return min(C_coll, C_mmr)
end



"""
Determines the AMD stability of a planetary system using the
collision condition and the MMR overlap conditions of Laskar
& Petit (2017) and Petit & Laskar (2017).

INPUT:
		μ		Planet/star mass ratios.
		a		Semimajor axes.
		e		Eccentricities.
		I		Inclinations (in radians).

NOTE:	The planets must be sorted from inner to outer.

OUTPUT:
		stat	Return status:
				:Stable			AMD stable.
				:Collison		Collision condition
				:MMR_eccentric	MMR overlap condition for eccentric orbits.
				:MMR_circular	MMR overlap condition for circular orbits.

		k		Index of inner planet that may be unstable,
				or 0 if the system is AMD stable.
"""
function AMD_stability(μ::Array{T}, a::Array{T}, e::Array{T},
						I::Array{T}) where T <: Real
	#
	# WARNING:	Laskar & Petit define μ = G*Mstar
	#			But I define μ = m/Mstar and set G*Mstar = 1.
	#
	N = length(μ)
	@assert(length(a) == N)
	@assert(length(e) == N)
	@assert(length(I) == N)
	#
	# We'll work in units where G*Mstar = 1.
	#
	Λ = μ .* sqrt.(a) # Ang. momentum of a circular orbit.
	AMD = 0.0
	for k in 1:N
		AMD += Λ[k] * (1 - sqrt((1-e[k])*(1+e[k])) * cos(I[k]))
	end
	#
	# Return a list of possibly unstable pairs, where:
	#
	#   1 -->  star + innermost planet
	#   k -->  kth planet + (k+1)th planet.
	#
	stats = Symbol[]
	pairs = Int[]
	ratios = Float64[]
	#=
	NOTE:	1 -->	ratio 	= AMD / Λ[1]
			k -->	ratio	= Cx / C_crit
							= (AMD/Λ[k]) / C_crit
							= AMD / (Λ[k] * C_crit)
							= AMD / Λ_crit[k]
	
	If we let AMD_crit be the critical ang.momentum deficit that
	we can put on the kth body to make it unstable w.r.t. the
	(k-1)st body, then it is sensible to define
	
			AMD_crit[1] = Λ[1]
	
	In other words, for the 1st planet to hit the star, the AMD
	you need to put in the planet is equal to the entire ang.
	momentum of a circular orbit. And for the other planets, it
	is set by the orbit-crossing and MMR criteria.
	=#
	#
	# Report on whether planet 1 could collide with the star.
	#
	push!(pairs,1)
	push!(stats,AMD > Λ[1] ? :Collision : :Stable)
	push!(ratios,AMD/Λ[1])
	
	for k in 2:N
		Cx = AMD / Λ[k]
        #@info("k = $k; Cx = $Cx")
		stat,ratio = AMD_stability(μ[k-1],μ[k],a[k-1],a[k],Cx)
		#
		# Track the results.
		#
		push!(pairs,k)
		push!(stats,stat)
		push!(ratios,ratio)
	end
	return stats,pairs,ratios
end



"""
    critical_AMD_system(μ, a)

Compute the critical AMD for a system to be stable.

# Arguments:
- `μ::Vector{T}`: list of planet/star mass ratios.
- `a::Vector{T}`: list of semimajor axes.
Note: μ and a must be sorted from innermost to outermost.

# Returns:
The minimum of the critical AMD for each pair of planets (including the star-planet1 pair). Any greater AMD value for the system would render it AMD-unstable.
"""
function critical_AMD_system(μ::Vector{T}, a::Vector{T}) where T <: Real
    # As above, we define μ = m/Mstar and set G*Mstar = 1

    N = length(μ)
    @assert(length(a) == N)

    Λ = μ .* sqrt.(a) # angular momentum of circular orbits

    # Compute the critical AMD for each pair of planets:
    # AMD_crit_pairs[1] = Λ[1] is simply the star-planet 1 pair
    # AMD_crit_pairs[k] = Λ[k]*min(C_coll, C_mmr) for k=2,...,N (each (k-1,k) planet pair)
    AMD_crit_pairs = deepcopy(Λ)
    for k in 2:N
        AMD_crit_pairs[k] *= critical_relative_AMD(μ[k-1], μ[k], a[k-1], a[k])
    end
    return minimum(AMD_crit_pairs)
end





"""
    distribute_AMD_planets_equal(AMD_tot, N)

Distribute a total AMD amount equally amongst the planets in the system.

# Arguments:
- `AMD_tot::Real`: total AMD of the system.
- `N::Integer`: number of planets in the system.

# Returns:
- `AMD::Vector{Real}`: list of AMD per planet.
"""
function distribute_AMD_planets_equal(AMD_tot::Real, N::Integer)
    AMD = (AMD_tot/N) .* ones(N)
    return AMD
end

"""
    distribute_AMD_planets_per_mass(AMD_tot, μ)

Distribute a total AMD amount amongst the planets in the system, equally per unit mass.

# Arguments:
- `AMD_tot::Real`: total AMD of the system.
- `μ::Vector{T}`: planet/star mass ratios (can also be planet masses).

# Returns:
- `AMD::Vector{Real}`: list of AMD per planet.
"""
function distribute_AMD_planets_per_mass(AMD_tot::Real, μ::Vector{T}) where T <: Real
    N = length(μ)
    AMD = (AMD_tot/sum(μ)) .* μ
    return AMD
end

# Distribute a total AMD amount amongst the planets in the system:
function distribute_AMD_planets(AMD_tot::Real, μ::Vector{T}) where T <: Real
    #AMD = distribute_AMD_planets_equal(AMD_tot, length(μ))
    AMD = distribute_AMD_planets_per_mass(AMD_tot, μ)
    return AMD
end



"""
    distribute_AMD_planet_ecc_incl_random(AMD, μ, a)

Distribute the AMD of a planet amongst its eccentricity and inclination components randomly (sin(i), e*sin(w), and e*cos(w)).

# Arguments:
- `AMD::Real`: AMD assigned to the planet.
- `μ::Real`: planet/star mass ratio.
- `a::Real`: semimajor axis of the planet.

# Returns:
- `e::Real`: orbital eccentricity.
- `ω::Real`: argument of pericenter (rad).
- `i::Real`: inclination (rad) relative to invariant plane.
"""
function distribute_AMD_planet_ecc_incl_random(AMD::Real, μ::Real, a::Real)
    # Let x = e*sin(w), y = e*cos(w), z = sin(i)
    # Thus e = sqrt(x^2 + y^2) and i = asin(z)

    Λ = μ*sqrt(a) # angular momentum of circular orbit
    @assert(AMD <= Λ) # too much AMD if AMD >= Λ (planet can collide with star)
    # NOTE: allowing AMD = Λ for case of a single planet, where the critical AMD is Λ
    sumsq_xyz = (AMD/Λ)*(2 - AMD/Λ) # x^2+y^2+z^2, related to the total AMD

    # Randomly assign x^2, y^2, z^2 such that their sum equals sumsq_xyz:
    split = sort(Random.rand(2))
    xsq, ysq, zsq = [split[1], split[2]-split[1], 1-split[2]] .* sumsq_xyz
    x, y = sign(randn())*sqrt(xsq), sign(randn())*sqrt(ysq)

    e = sqrt(xsq + ysq) # eccentricity
    ω = atan(x, y) # argument of pericenter
    i = asin(sqrt(zsq)) # inclination relative to invariant plane (rad)
    #@assert(sqrt(1 - e^2)*cos(i) >= 1 - AMD/Λ) # NOTE: I expected this to be equal given how we assigned x^2+y^2+z^2 = (AMD/Λ)*(2 - AMD/Λ), but in practice it is >= (which does not break AMD stability since this implies the true AMD given (e,i) is less than the AMD provided)
    #@info("e = $e, ω = $(ω*180/π) deg, i = $(i*180/π) deg")

    return e, ω, i
end

# Distribute the AMD of a planet amongst its eccentricity and inclination:
function distribute_AMD_planet_ecc_incl(AMD::Real, μ::Real, a::Real)
    e, ω, i = distribute_AMD_planet_ecc_incl_random(AMD, μ, a)
    return e, ω, i
end



"""
    draw_ecc_incl_system_critical_AMD(μ, a; check_stability=false)

Draw eccentricities and inclinations for the planets in a system by distributing the critical AMD of the system.

# Arguments:
- `μ::Vector{T}`: list of planet/star mass ratios.
- `a::Vector{T}`: list of semimajor axes.
- `check_stability::Bool=false`: whether to double check (if true) or not (if false) that the system is AMD-stable.
NOTE 1: μ and a must be sorted from innermost to outermost.
NOTE 2: by definition, the system should ALWAYS be AMD-stable if drawn using this function. Thus the `check_stability` flag is more of a debugging/testing tool.

# Returns:
- `AMD::Vector{Real}`: list of AMD assigned to the planets.
- `e::Vector{Float64}`: list of eccentricities drawn for the planets.
- `ω::Vector{Float64}`: list of arguments of pericenter (rad) drawn for the planets.
- `i::Vector{Float64}`: list of inclinations (rad) relative to the invariant plane drawn for the planets.
"""
function draw_ecc_incl_system_critical_AMD(μ::Vector{T}, a::Vector{T}; check_stability::Bool=false) where T <: Real
    N = length(μ)
    @assert(length(a) == N)

    AMD_tot = critical_AMD_system(μ, a) # total (critical) AMD of system
    AMD = distribute_AMD_planets(AMD_tot, μ) # list of AMD distributed to each planet

    e = Vector{Float64}(undef, N)
    ω = Vector{Float64}(undef, N)
    i = Vector{Float64}(undef, N)
    for n in 1:N
        e[n], ω[n], i[n] = distribute_AMD_planet_ecc_incl(AMD[n], μ[n], a[n])
    end

    # To double check that the system is AMD stable given the drawn e and i:
    if check_stability
        stats, pairs, ratios = AMD_stability(μ, a, e, i)
        @assert(all(stats .== :Stable))
    end

    return AMD, e, ω, i
end
