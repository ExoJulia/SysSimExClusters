#=
Date:   June 15, 2019
Author: Daniel Carrera (dcarrera@gmail.com)
=#
using Printf

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
	#
	# Convention: 1 == inner planet; 2 == outer planet.
	#
	@assert(a1 <= a2)
	γ = μ1/μ2
	α = a1/a2
	ϵ = μ1+μ2   # Equation  4 of Petit & Laskar (2017).
	r = 0.80199 # Equation 28 of Petit & Laskar (2017).
	#
	# Equantion 76 and Section 3.6 of Petit & Laskar (2017).
	#
	α_crit = 1 - 1.46 * ϵ^(2/7)
	if α > α_crit
		return :MMR_circular, Inf
	end
	#
	# Relative AMD fron collision condition (Laskar & Petit 2017).
	#
	e1 = critical_eccentricity(γ,α)
	e2 = 1 - α - α*e1
	@assert(a1*(1 + e1) ≈ a2*(1 - e2))
	sqrtα = sqrt(α)
	C_coll = γ*sqrtα*(1 - sqrt((1-e1)*(1+e1)) + 1 - sqrt((1-e2)*(1+e2)))
	#
	# Relative AMD fron MMR overlap condition (Petit & Laskar 2017).
	#
	g = 3^4 * (1 - α)^5 / (2^9 * r * ϵ) - 32r * ϵ / (9*(1 - α)^2)
	C_mmr = g^2 * γ*sqrtα / (2 + 2γ*sqrtα)
	#
	# Final result.
	#
	ratio = Cx/min(C_coll,C_mmr)
	if Cx < min(C_coll,C_mmr)
		return :Stable,ratio
	elseif C_coll < C_mmr
		return :Collision,ratio
	else
		return :MMR_eccentric,ratio
	end
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
