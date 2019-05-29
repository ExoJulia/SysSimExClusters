# SysSimExClusters

This repository provides a comprehensive model for Clustered Planetary Systems for SysSim. The details of the specific models for exoplanetary systems are described in our paper (link to be provided).



## How to use our Clustered model for studying exoplanetary system populations:

We provide a large set of simulated catalogs from our models directly in this repository. We refer to a *physical catalog* as a set of intrinsic, physical planetary systems, and an *observed catalog* as a set of transiting and detected planet candidates that a Kepler-like mission might produce. If you simply wish to use these simulated catalogs as examples of our models, then no installation is required! Simply download any of these tables and use them for your own science.
* Navigate to the "examples/best_models/Clustered_P_R/" directory. (*Note: you will find that we also provide simulated catalogs from two other models, "Clustered_P" and "Non_Clustered", in their respective subdirectories. These models are also described in the paper linked above but we do not recommend using these models for science as they are not as good of a description of the data as our "Clustered_P_R" model.*)
* Download a "physical_catalogX.csv" file (X = an index/number) for a table including all the physical planets in a simulated catalog.
* Download an "observed_catalogX.csv" file (X = an index/number) for a table including all the observed planets that a simulated Kepler mission would detect given the true planetary systems listed in "physical_catalogX.csv".



## How to simulate catalogs (physical and observed) from our Clustered model on your own:

You will need to first install the [ExoplanetsSysSim](https://github.com/ExoJulia/ExoplanetsSysSim.jl) package and set up some additional repositories; follow the instructions listed in the README of that page.



## If you wish to make similar plots as those included in our paper:

While the core ExoplanetsSysSim and SysSimExClusters code is written in Julia, almost all of the figures produced for the paper are generated from Python (2.7) code that was written by Matthias He. We provide these Python scripts but do not fully maintain or document them.



## If you wish to reproduce the results and methodology described in our paper:

First read the paper, especially Section 2 (Methods), in detail.
