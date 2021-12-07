# SysSimExClusters

[He, Ford, & Ragozzine (2021b), AJ, 162, 216 (22pp)](https://ui.adsabs.harvard.edu/abs/2021AJ....162..216H/abstract) [[arXiv](https://arxiv.org/abs/2105.04703)]


This paper utilizes the same model ("Maximum AMD model") from "Paper III"; refer to that paper for a full description of the model.

Here, we provide a framework for users to simulate planetary systems *conditioned* on a given planet (i.e. within a given period, radius, and mass range for the properties of the chosen planet). For example, you may have a system with a known transiting planet of measured period and radius, and wish to compute the conditional distribution of planets in such systems.


## To simulate systems conditioned on a given planet of your choosing:

The procedure for generating simulated catalogs (physical and observed pairs) with a user defined set of model parameters is the same as before. Refer to the Paper III branch ("He_et_al_2020b") for a description of the procedure and the simulation outputs, and refer to the paper itself for a description of all the model parameters.

To generate a simulated catalog only containing systems with a conditioned planet:

1. Move into the "src/" directory and edit the file "param_common.jl". Make sure the `generate_planetary_system_clustered_conditional` option is selected:
```julia
#add_param_fixed(sim_param,"generate_planetary_system", generate_planetary_system_clustered) # this is for simulating a full catalog without conditioning (e.g. for previous papers)
add_param_fixed(sim_param,"generate_planetary_system", generate_planetary_system_clustered_conditional) # this is for simulating a catalog conditioned on a given planet
```

2. Set a range for the period, radius, and mass of the conditioned planet:
```julia
add_param_fixed(sim_param,"cond_period_min", 8.) # minimum period for the conditioned planets
add_param_fixed(sim_param,"cond_period_min", 12.) # maximum period for the conditioned planets
add_param_fixed(sim_param,"cond_radius_min", 0.9*ExoplanetsSysSim.earth_radius) # minimum radius for the conditioned planets
add_param_fixed(sim_param,"cond_radius_max", 1.1*ExoplanetsSysSim.earth_radius) # maximum radius for the conditioned planets
add_param_fixed(sim_param,"cond_mass_min", 0.5*ExoplanetsSysSim.earth_mass) # minimum mass for the conditioned planets
add_param_fixed(sim_param,"cond_mass_max", 2.0*ExoplanetsSysSim.earth_mass) # maximum mass for the conditioned planets
```
**Note**: the ranges for the planet radii and masses are in solar units, and thus one should multiply by `ExoplanetsSysSim.earth_radius` and `ExoplanetsSysSim.earth_mass` respectively if entering values in Earth units, as shown in this example.

If you do not want to condition on one or more properties (e.g. planet mass, if you only know the period and radius of the conditioned planet), simply input `0` and `Inf` for the min and max values.

3. Choose whether or not the conditioned planet must also transit:
```julia
add_param_fixed(sim_param,"cond_also_transits", true) # true for transiting, false for isotropic distribution (may or may not transit)
```

4. Due to the nature of the rejection sampling, the narrower the conditioned ranges are (and especially if the conditioned planets must also transit), the longer the simulations will take! We also advise reducing the total number of simulated systems:
```julia
add_param_fixed(sim_param,"num_targets_sim_pass_one", 1000) # how many systems with conditioned planets to collect in total
```

5. Move into a directory where you want your simulated catalogs to be saved and run the script "generate_catalogs.jl" in "examples/" in Julia. For example, you can navigate to "examples/", start Julia, and run:
```julia
include("generate_catalogs.jl")
```

As before, the physical catalog will contain all the simulated systems (you can check that each system contains at least one planet satisfying the conditioning bounds), while the observed catalog will only contain the measured properties of Kepler-detected planets.
