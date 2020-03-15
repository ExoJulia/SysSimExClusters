# SysSimExClusters

He, Ford, and Ragozzine (2019) (Paper I): \[published in [MNRAS](https://doi.org/10.1093/mnras/stz2869)\] \[[arXiv](https://arxiv.org/abs/1907.07773)\]



## To download simulated catalogs from our models:

* Go to the [SysSimExClusters Simulated Catalogs](https://psu.box.com/s/v09s9fhbmyele911drej29apijlxsbp3) folder.
* Navigate into "He_Ford_Ragozzine_2019/", then into "Clustered_P_R/"\*, and then into either the KS or AD subdirectories. (Use KS if you don't know the difference or have trouble deciding.) 
* Download a "physical_catalogX.csv" file (X = an index/number) for a table including all the physical planets in a simulated catalog.

| target_id | star_id | planet_mass    | planet_radius | clusterid | period     | ecc      | star_mass      | star_radius |
|-----------|---------|----------------|---------------|-----------|------------|----------|----------------|-------------|
|           |         | (solar masses) | (solar radii) |           | (days)     |          | (solar masses) | (solar radii) |
| 1.0       | 17593.0 | 3.3830e-7      | 0.0061        | 1.0       | 159.0465   | 0.0107   | 0.992          | 0.985       |
| ...       | ...     | ...            | ...           | ...       | ...        | ...      | ...            | ...         |

* Download an "observed_catalogX.csv" file (X = an index/number) for a table including all the observed planets that a simulated Kepler mission would detect given the true planetary systems listed in "physical_catalogX.csv".

| target_id | star_id | period    | depth   | duration | star_mass      | star_radius |
|-----------|---------|-----------|---------|----------|----------------|-------------|
|           |         | (days)    |         | (days)   | (solar masses) | (solar radii) |
| 75.0      | 66728.0 | 52.7541   | 0.00083 | 0.1090   | 0.928          | 0.7874      |
| ...       | ...     | ...       | ...     | ...      | ...            | ...         |

In each of these files, the header contains all the parameters of the model used to generate that catalog.

For each planet (row),
* **target_id** refers to the index of the star in the simulation (e.g. 1 for the first star in the simulation) the planet orbits,
* **star_id** refers to the index of the star based on where it is in the input stellar catalog (i.e. the row number in the "q1_q17_dr25_gaia_fgk_cleaned.csv" catalog, which can be found in the "plotting/" directory), and
* **clusterid** is a cluster identifier (i.e., which "cluster" in the system the planet belongs to).
Note that indexing starts at 1 in Julia. Stars without any planets are not included in these tables.

\*Note: We also provide simulated catalogs from two other models, "Clustered_P" and "Non_Clustered", in their respective subdirectories. These models are also described in the paper linked above but we do not recommend using these models for science as they are not as good of a description of the data as our "Clustered_P_R" model.



## To simulate new catalogs from our Clustered Periods and Sizes model:

To generate one simulated catalog (physical and observed) with a user defined set of model parameters:

1. Move into the "src/" directory and edit the file "param_common.jl". Set a value for each of the model parameters. For example, to change the mean number of attempted clusters (&lambda;<sub>c</sub>), planets per cluster (&lambda;<sub>p</sub>), and period power-law (&alpha;<sub>P</sub>):
```julia
add_param_active(sim_param,"log_rate_clusters", log(1.6))            # Set the (log) mean number of clusters
add_param_active(sim_param,"log_rate_planets_per_cluster", log(1.6)) # Set the (log) mean number of planets per cluster
add_param_active(sim_param,"power_law_P", 0.)                        # Set the period power-law index (-1 = flat in log-period)
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
