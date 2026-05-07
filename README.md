# SysSimExClusters

[He \& Ford (2026), ApJ, 1001, 78 (34pp)](https://ui.adsabs.harvard.edu/abs/2026ApJ..1001...78H/abstract) ("Paper IV") \[[arXiv](https://arxiv.org/abs/2601.13480)\]

<center><img src="/cartoons/Hybrid_model_cartoon_no_steps.png" alt="Cartoon illustration of the hybrid models" width="600"/></center>  



## To download simulated catalogs from our models:

* Go to the [SysSim Catalogs Google Drive](https://drive.google.com/drive/folders/1hUsNwjtF0Y8Y8PCxORrQfcWPlFKh-rFl).
* Navigate into "HM-C/" and pick a folder ("100 catalogs, depth_min 0.29" is recommended for the best models). 
* Download the "physical_catalogs.zip" file, which contains many "physical_catalogX.csv" tables (X = an index/number) including all the physical planets in each simulated catalog.

| target_id | star_id | planet_mass    | planet_radius | clusterid | period     | ecc      | incl      | omega     | asc_node   | mean_anom | incl_invariable | asc_node_invariable | star_mass      | star_radius |
|-----------|---------|----------------|---------------|-----------|------------|----------|-----------|-----------|------------|-----------|-----------------|---------------------|----------------|-------------|
|           |         | (solar masses) | (solar radii) |           | (days)     |          | (radians) | (radians) | (radians)  | (radians) | (radians)       | (radians)           | (solar masses) | (solar radii) |
| 1         | 60465   | 3.8166e-7      | 0.0048        | 1         | 21.1215    | 0.2292   | 1.3397    | 2.9827    | 5.8627     | 2.5189    | 0.1540          | 2.7316              | 1.098          | 1.254       |
| ...       | ...     | ...            | ...           | ...       | ...        | ...      | ...       | ...       | ...        | ...       | ...             | ...                 | ...            | ...         |

* Download the "observed_catalogs.zip" file, which contains many "observed_catalogX.csv" tables (X = an index/number) including all the observed planets that a simulated Kepler mission would detect given the true planetary systems listed in "physical_catalogX.csv".

| target_id | star_id | period    | period_err | depth   | depth_err | duration | duration_err    | star_mass      | star_radius |
|-----------|---------|-----------|------------|---------|-----------|----------|-----------------|----------------|-------------|
|           |         | (days)    | (days)     |         |           | (days)   | (days)          | (solar masses) | (solar radii) |
| 193       | 1132    | 5.9486    | 0.0009     | 0.0007  | 1.9547e-5 | 0.1209   | 0.00182         | 0.983          | 0.973       |
| ...       | ...     | ...       | ...        | ...     | ...       | ...      | ...             | ...            | ...         |

In each of these files, the header contains all the parameters of the model used to generate that catalog.

:mega: **NOTE: The formats of these tables have NOT changed compared to the simulated catalogs from the previous models, so any functions used to load the previous tables should also work on these new files.**

For each planet (row):
* **target_id**: the index of the star in the simulation (e.g. 1 for the first star in the simulation) the planet orbits
* **star_id**: the index of the star based on where it is in the input stellar catalog (i.e. the row number in the "q1_q17_dr25_gaia_fgk_cleaned.csv" catalog, which can be found in the "plotting/" directory)
* **clusterid**: a cluster identifier (i.e., which "cluster" in the system the planet belongs to)
* **incl**: inclination of the orbit relative to the sky plane
* **omega**: argument of periapse relative to the sky plane
* **asc_node**: argument of ascending node relative to the sky plane
* **mean_anom**: mean anomaly relative to the sky plane
* **incl_invariable**: inclination relative to the system invariable plane
* **asc_node_invariable**: argument of ascending node relative to the system invariable plane

All other fields should be self explanatory.
Note that indexing starts at 1 in Julia. Stars without any planets are not included in these tables.



## To simulate new catalogs from the "Hybrid models" (He \& Ford 2026) on your own:

The procedure for simulating catalogs from the new model is the same as before. To generate one simulated catalog (physical and observed) with a user defined set of model parameters:

1. Move into the "src/" directory and edit the file "param_common.jl". Set a value for each of the model parameters. The relevant parameters of the hybrid models are:
```julia
add_param_active(sim_param, "log_rate_clusters", log(0.40)) # the (log) mean number of clusters per system
add_param_active(sim_param, "log_rate_planets_per_cluster", log(2.12)) # the (log) mean number of planets per cluster
add_param_active(sim_param, "power_law_P", -0.1) # the period power-law index
add_param_fixed(sim_param, "sigma_logperiod_per_pl_in_cluster", 0.25) # the cluster width (std) in log-period, per planet in the cluster, of each cluster's period distribution
add_param_active(sim_param, "mean_ln_mass", 0.35) # the mean log-mass of the initial planet mass distribution, in ln(Earth masses)
add_param_active(sim_param, "sigma_ln_mass", 1.6) # the standard deviation in log-mass of the initial planet mass distribution, in ln(Earth masses)
add_param_active(sim_param, "norm_radius", 2.1) # the normalization radius (corresponding to an initial mass of 1 Earth mass), in Earth radii
add_param_fixed(sim_param, "break1_mass", 20.) # the break mass, in Earth masses
add_param_active(sim_param, "power_law_γ0", 0.07) # the radius-mass power-law index below the break mass
add_param_fixed(sim_param, "power_law_γ1", 0.5) # the radius-mass power-law index above the break mass
add_param_active(sim_param, "power_law_σ0", 0.23) # the spread (std) around the mean radius, below the break mass
add_param_fixed(sim_param, "power_law_σ1", 0.3) # the spread (std) around the mean radius, above the break mass
add_param_active(sim_param, "log_α_pret", log(8.)) # normalization factor for the envelope retention probability
add_param_active(sim_param, "sigma_hk", 0.25) # the eccentricity (Rayleigh) scale for true singles
```
:twisted_rightwards_arrows: **There are two types of Hybrid Models: "Unclustered" (HM-U) vs. "Clustered" (HM-C).** To simulate catalogs from the **HM-C** (the preferred model), also set the cluster width (std) in log-mass of the initial planet mass distribution:
```julia
add_param_active(sim_param, "sigma_ln_mass_in_cluster", 0.48) # for clustered initial masses, in ln(Earth masses)
```
If you want to simulate from the **HM-U** instead, do NOT set this parameter (i.e., this line should be commented out).

:bulb: **NOTE:** as before, there is no difference between setting a parameter as "active" vs. "fixed" for the purposes of simulating a catalog. This is only used for choosing which parameters to vary if one wants to run optimization algorithms on the models.

2. Move into a directory where you want your simulated catalogs to be saved and run the script "generate_catalogs.jl" in "examples/" in Julia. For example, you can navigate to "examples/", start Julia, and run:
```julia
include("generate_catalogs.jl")
```
This will generate the following files:
* Physical and observed catalogs of planets (and stars) in table format:
  * "physical_catalog.csv"
  * "physical_catalog_stars.csv"
  * "observed_catalog.csv"
  * "observed_catalog_stars.csv"

These files are analogous to the simulated catalogs we provide as described above.

In addition, the following files will also be generated:
* Individual files with the true cluster id's, periods, orbital eccentricities, mutual inclinations, sky inclinations, planet radii, planet masses, stellar radii, and stellar masses, of all the planets per system (and stars with planets) in the physical catalog:
  * "clusterids_all.out"
  * "periods_all.out"
  * "eccentricities_all.out"
  * "mutualinclinations_all.out"
  * "inclinations_all.out"
  * "radii_all.out"
  * "masses_all.out"
  * "stellar_radii_with_planets.out"
  * "stellar_masses_with_planets.out"

The data in these files are the same as those in "physical_catalog.csv", just organized in a different format.

* Individual files with the planets' *initial* radii, masses, and envelope masses (i.e., before envelope mass-loss due to photoevaporation), as well as the mass-loss timescales, envelope retention probabilities, and envelope retention booleans, of all the planets per system in the physical catalog:
  * **"initial_radii_all.out"**
  * **"initial_masses_all.out"**
  * **"envelope_masses_all.out"**
  * **"mass_loss_timescales_all.out"**
  * **"prob_retained_all.out"**
  * **"envelope_retained_all.out"**

**NOTE:** these files (bolded) are only generated for the new models (i.e. **HM-U** and **HM-C**). This data is NOT provided in the tables/csv files.

* Individual files with the observed periods, transit depths, transit durations, stellar radii, and stellar masses of all the planets (and stars with observed planets) in the observed catalog:
  * "periods.out"
  * "depths.out"
  * "durations.out"
  * "stellar_radii_obs.out"
  * "stellar_masses_obs.out"

The data in these files are the same as those in "observed_catalog.csv", just organized in a different format (sorted into systems with 1, 2, 3, ..., and 8 observed planets). For the simulated catalogs in the SysSim Catalogs folder, these files are also provided in the "catalogs_separate.zip" file.
