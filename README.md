# SysSimExClusters

He, Ford, and Ragozzine (2020) (Paper II): \[submitted to MNRAS\] \[[arXiv](https://arxiv.org/abs/2003.04348)\]



## To download simulated catalogs from our models:

* Go to the [SysSimExClusters Simulated Catalogs](https://psu.box.com/s/v09s9fhbmyele911drej29apijlxsbp3) folder.
* Navigate into "He_Ford_Ragozzine_2020/", then into the "Clustered_P_R_fswp_bprp/"\*, and then into either the KS or AD subdirectories. (Use KS if you don't know the difference or have trouble deciding.) 
* Download a "physical_catalogX.csv" file (X = an index/number) for a table including all the physical planets in a simulated catalog.

| target_id | star_id | planet_mass    | planet_radius | clusterid | period     | ecc      | incl_mut  | star_mass      | star_radius |
|-----------|---------|----------------|---------------|-----------|------------|----------|-----------|----------------|-------------|
|           |         | (solar masses) | (solar radii) |           | (days)     |          | (radians) | (solar masses) | (solar radii) |
| 1.0       | 5441.0  | 2.5670e-5      | 0.0182        | 1.0       | 105.5113   | 0.0199   | 0.5802    | 0.896          | 0.818       |
| ...       | ...     | ...            | ...           | ...       | ...        | ...      | ...       | ...            | ...         |

Note that compared to the physical catalogs from Paper I, these files have an extra column (the mutual inclination, "incl_mut").

* Download an "observed_catalogX.csv" file (X = an index/number) for a table including all the observed planets that a simulated Kepler mission would detect given the true planetary systems listed in "physical_catalogX.csv".

| target_id | star_id | period    | depth   | duration | star_mass      | star_radius |
|-----------|---------|-----------|---------|----------|----------------|-------------|
|           |         | (days)    |         | (days)   | (solar masses) | (solar radii) |
| 1.0       | 5441.0  | 256.4506  | 0.00076 | 0.2215   | 0.896          | 0.818       |
| ...       | ...     | ...       | ...     | ...      | ...            | ...         |

In each of these files, the header contains all the parameters of the model used to generate that catalog.

For each planet (row),
* **target_id** refers to the index of the star in the simulation (e.g. 1 for the first star in the simulation) the planet orbits,
* **star_id** refers to the index of the star based on where it is in the input stellar catalog (i.e. the row number in the "q1_q17_dr25_gaia_fgk_cleaned.csv" catalog, which can be found in the "plotting/" directory), and
* **clusterid** is a cluster identifier (i.e., which "cluster" in the system the planet belongs to).
Note that indexing starts at 1 in Julia. Stars without any planets are not included in these tables.

\*Note: We also provide simulated catalogs from a "Clustered_P_R_fswp" model, which assumes the same fraction of stars with planets across all stars. This model is also described in the paper linked above.



## To simulate new catalogs from our Linear f<sub>swpa</sub>(b<sub>p</sub> - r<sub>p</sub>) model on your own:

To generate one simulated catalog (physical and observed) with a user defined set of model parameters:

1. Move into the "src/" directory and edit the file "param_common.jl". Set a value for each of the model parameters. For example, to change the fraction of stars with planets at the median colour and the slope for the linear f<sub>swpa</sub>(b<sub>p</sub> - r<sub>p</sub>) relation:
```julia
add_param_active(sim_param,"f_stars_with_planets_attempted_color_slope", 0.6)  # Set the slope
add_param_active(sim_param,"f_stars_with_planets_attempted_at_med_color", 0.6) # Set the normalization
```
**NOTE:** for most accurate results, you should also add/uncomment the line:
```julia
add_param_fixed(sim_param,"osd_file","dr25fgk_relaxcut_osds.jld2")
```
However, this file requires about 8gb of memory to read.

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
* An individual file for the true cluster id's, periods, orbital eccentricities, planet radii, planet masses, stellar radii, and stellar masses, of all the planets per system (and stars with planets) in the physical catalog:
  * "clusterids_all.out"
  * "periods_all.out"
  * "eccentricities_all.out"
  * "mutualinclinations_all.out"
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

The data in these files are the same as those in "observed_catalog.csv", just organized in a different format (sorted into systems with 1, 2, 3, ..., and 8 observed planets).
