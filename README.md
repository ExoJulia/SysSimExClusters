# SysSimExClusters

[He et al (2020), AJ, 160, 276 (38pp)](https://ui.adsabs.harvard.edu/abs/2020AJ....160..276H/abstract) ("Paper III") \[[arXiv](https://arxiv.org/abs/2007.14473)\]



## To download simulated catalogs from our models:

* Go to the [SysSimExClusters Simulated Catalogs](https://pennstateoffice365-my.sharepoint.com/:f:/g/personal/myh7_psu_edu/Ei7QJqnmaCBGipPM4uMzrusBjw_hUwo0KfIDBe-0UTYyMw) folder.
* Navigate into "He_et_al_2020b/" and through the subdirectories. 
* Download the "physical_catalogs.zip" file, which contains many "physical_catalogX.csv" tables (X = an index/number) including all the physical planets in each simulated catalog.

| target_id | star_id | planet_mass    | planet_radius | clusterid | period     | ecc      | incl      | omega     | asc_node   | mean_anom | incl_invariable | asc_node_invariable | star_mass      | star_radius |
|-----------|---------|----------------|---------------|-----------|------------|----------|-----------|-----------|------------|-----------|-----------------|---------------------|----------------|-------------|
|           |         | (solar masses) | (solar radii) |           | (days)     |          | (radians) | (radians) | (radians)  | (radians) | (radians)       | (radians)           | (solar masses) | (solar radii) |
| 1         | 32722   | 2.3983e-5      | 0.0366        | 1         | 13.0340    | 0.0124   | 1.1409    | -2.6147   | 5.5608     | 1.1570    | 0.0298          | 4.6866              | 1.031          | 1.32        |
| ...       | ...     | ...            | ...           | ...       | ...        | ...      | ...       | ...       | ...        | ...       | ...             | ...                 | ...            | ...         |

* Download the "observed_catalogs.zip" file, which contains many "observed_catalogX.csv" tables (X = an index/number) including all the observed planets that a simulated Kepler mission would detect given the true planetary systems listed in "physical_catalogX.csv".

| target_id | star_id | period    | period_err | depth   | depth_err | duration | duration_err    | star_mass      | star_radius |
|-----------|---------|-----------|------------|---------|-----------|----------|-----------------|----------------|-------------|
|           |         | (days)    | (days)     |         |           | (days)   | (days)          | (solar masses) | (solar radii) |
| 65        | 77706   | 4.3019    | 0.0004     | 0.0015  | 2.6409e-5 | 0.0835   | 0.00097         | 0.778          | 0.745       |
| ...       | ...     | ...       | ...        | ...     | ...       | ...      | ...             | ...            | ...         |

In each of these files, the header contains all the parameters of the model used to generate that catalog.

For each planet (row),
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



## To simulate new catalogs from our Maximum AMD model on your own:

To generate one simulated catalog (physical and observed) with a user defined set of model parameters:

1. Move into the "src/" directory and edit the file "param_common.jl". Set a value for each of the model parameters. For example, to change the eccentricity scale of intrinsic single planets:
```julia
add_param_active(sim_param,"sigma_hk", 0.25) # Set the eccentricity (Rayleigh) scale for true singles
```
**NOTE:** the line:
```julia
add_param_fixed(sim_param,"osd_file","dr25fgk_relaxcut_osds.jld2")
```
will load a file that requires about 8gb of memory to read. If your machine cannot handle that, you can comment out this line and uncomment the previous line to use a smaller osd file ("dr25fgk_small_osds.jld2"). However, the results will not be as accurate.

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
* An individual file for the true cluster id's, periods, orbital eccentricities, mutual inclinations, sky inclinations, planet radii, planet masses, stellar radii, and stellar masses, of all the planets per system (and stars with planets) in the physical catalog:
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
* An individual file for the observed periods, transit depths, transit durations, stellar radii, and stellar masses of all the planets (and stars with observed planets) in the observed catalog:
  * "periods.out"
  * "depths.out"
  * "durations.out"
  * "stellar_radii_obs.out"
  * "stellar_masses_obs.out"

The data in these files are the same as those in "observed_catalog.csv", just organized in a different format (sorted into systems with 1, 2, 3, ..., and 8 observed planets). For the simulated catalogs in the SysSimExClusters Simulated Catalogs folder, these files are also provided in the "catalogs_separate.zip" file.
