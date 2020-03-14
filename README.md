# SysSimExClusters

This repository provides a comprehensive forward modelling framework for studying planetary systems using the Kepler planet candidates.

<center><img src="/best_models/Clustered_P_R_observed.gif" alt="Best fit clustered models" width="800"/></center>  

We develop and provide several statistical models for describing the intrinsic planetary systems, their architectures, and the correlations within multi-planet systems. Our specific models are described in the following papers:

* He, Ford, and Ragozzine (2019) (Paper I): \[published in [MNRAS](https://doi.org/10.1093/mnras/stz2869)\] \[[arXiv](https://arxiv.org/abs/1907.07773)\]
* He, Ford, and Ragozzine (2020) (Paper II): \[submitted to MNRAS\] \[[arXiv](https://arxiv.org/abs/2003.04348)\]

**Important:** We have a separate code branch for each paper (e.g. see "He_Ford_Ragozzine_2019", "He_Ford_Ragozzine_2020"); these should be used if you want to run our code instead of the master branch, which is actively updated. In addition, the README file is different for each branch, and we provide more details for the models and code usage specific to each paper/branch.



## How to use our models for studying exoplanetary system populations:

We provide a large set of simulated catalogs from our models in the [SysSimExClusters Simulated Catalogs](https://psu.box.com/s/v09s9fhbmyele911drej29apijlxsbp3) folder. If you simply wish to use these simulated catalogs as examples of our models, then no installation is required! Simply download any of these tables and use them for your own science. To be able to use them, you must understand that we provide two types of catalogs:

* *Physical catalog:* a set of intrinsic, physical planetary systems (before any observations; contains properties like the true orbital periods, planet radii, etc.)
* *Observed catalog:* a set of transiting and detected planet candidates derived from a *physical catalog* (after a Kepler-like mission; contains properties like the measured orbital periods, transit depths, etc.)

Refer to the README of the branch specific to each paper for complete details on what each set of catalogs contains.



## How to simulate catalogs (physical and observed) on your own:

### Installation:

* You will need to first install the [ExoplanetsSysSim](https://github.com/ExoJulia/ExoplanetsSysSim.jl) package and set up some additional repositories; follow the instructions listed in the README of that page.
* Clone this repository.
```
git clone git@github.com:ExoJulia/SysSimExClusters.git
```
* Switch to the branch of this repository containing the model you want to simulate from. For example, to simulate models from Paper II, do:
```
git checkout He_Ford_Ragozzine_2020
```

### Usage:

Refer to the README of the branch containing the model you want to simulate from for steps.



## If you wish to make similar plots as those included in our paper:

While the core ExoplanetsSysSim and SysSimExClusters code is written in Julia, almost all of the figures produced for the paper are generated from Python (3.7) code that was written by Matthias He. We provide these Python scripts in the "plotting/" directory but do not fully maintain or document them yet.



## What if I need help?

Feel free to email Matthias He at myh7@psu.edu!
