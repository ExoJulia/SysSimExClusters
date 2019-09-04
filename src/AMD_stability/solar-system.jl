
include("amd_stability.jl")

using Printf

##################################################
#
#     EXAMPLE -- SOLAR SYSTEM
#
##################################################
μ = [
	1.66e-7, 2.448e-6, 3e-6, 3.227e-7, # Terrestrial planets.
	9.546e-4, 2.858e-4, 4.366e-5, 5.151e-5 # Jovian planets.
]
a = [
	0.387098, 0.723332, 1.00, 1.523679,
	5.2044, 9.5826, 19.2184, 30.11
]
e = [
	0.205630, 0.006772, 0.0167086, 0.0934,
	0.0489, 0.0565, 0.046381, 0.009456
]
I = [ # Relative to the invariant plane.
	6.34, 2.19, 1.57869, 1.67,
	0.32, 0.93, 1.02, 0.72
] .* π/180


function report(name,stats,pairs)
	if all(stats .== :Stable)
		@info("$name: AMD stable.")
	else
		@info("$name: AMD unstable")
		for j in 1:length(stats)
			k = pairs[j]
			@info("    Pair ($(k-1),$k) -- $(stats[j])")
		end
	end
	println("")
end

#
# TEST 1
#
stats, pairs, ratios = AMD_stability(μ, a, e, I)
report("Solar System", stats, pairs)
#
# TEST 2
#
#stats, pairs, ratios = AMD_stability(μ[1:4], a[1:4], e[1:4], I[1:4])
#report("Terrestrial Planets", stats, pairs)
#
# TEST 3
#
#stats, pairs, ratios = AMD_stability(μ[5:8], a[5:8], e[5:8], I[5:8])
#report("Jovian Planets", stats, pairs)
