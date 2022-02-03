# SysSimExClusters

This repository provides a comprehensive forward modelling framework for studying planetary systems.

<center><img src="/best_models/Clustered_P_R_observed.gif" alt="Best fit clustered models" width="800"/></center>  

We develop and provide several *statistical* models for describing the intrinsic planetary systems, their architectures, and the correlations within multi-planet systems, using the Kepler population of exoplanet candidates. Our specific models are described in the following papers:

* [He, Ford, and Ragozzine (2019), MNRAS, 490, 4575 (30pp)](https://ui.adsabs.harvard.edu/abs/2019MNRAS.490.4575H/abstract) ("Paper I") \[[arXiv](https://arxiv.org/abs/1907.07773)\]
* [He, Ford, and Ragozzine (2021a)\*, AJ, 161, 16 (24pp)](https://ui.adsabs.harvard.edu/abs/2021AJ....161...16H/abstract) ("Paper II") \[[arXiv](https://arxiv.org/abs/2003.04348)\]
* [He et al. (2020), AJ, 160, 276 (38pp)](https://ui.adsabs.harvard.edu/abs/2020AJ....160..276H/abstract) ("Paper III") \[[arXiv](https://arxiv.org/abs/2007.14473)\]

\*Paper II was published shortly after Paper III due to an extended referee process, but the most updated model is the one described in Paper III.

In addition to these papers describing the new models, the simulated catalogs from these models have been directly used for several other publications:

* [Gilbert and Fabrycky (2020), AJ, 159, 281 (17pp)](https://ui.adsabs.harvard.edu/abs/2020AJ....159..281G/abstract) \[[arXiv](https://arxiv.org/abs/2003.11098)\]
* [Millholland, He, Ford, et al. (2021), AJ, 162, 166 (17pp)](https://ui.adsabs.harvard.edu/abs/2021arXiv210615589M/abstract), \[[arXiv](https://arxiv.org/abs/2106.15589)\]
* [He, Ford, and Ragozzine (2021b), AJ, 162, 216 (22pp)](https://ui.adsabs.harvard.edu/abs/2021arXiv210504703H/abstract), \[[arXiv](https://arxiv.org/abs/2105.04703)\]

**Important:** We have a separate code branch for each paper that provides new models or functionality:

* "He_Ford_Ragozzine_2019" branch for He, Ford, and Ragozzine (2019)
* "He_Ford_Ragozzine_2020" branch for He, Ford, and Ragozzine (2021a)
* "He_et_al_2020" branch for He et al. (2020)
* "He_Ford_Ragozzine_2021b" branch for He, Ford, and Ragozzine (2021b)

These should be used if you want to run our code instead of the master branch, which is actively being updated. In addition, the README file is different for each branch, and we provide more details for the models and code usage specific to each paper/branch.



## How do I use these models?

We provide a large set of simulated catalogs from our models in the [SysSimExClusters Simulated Catalogs](https://pennstateoffice365-my.sharepoint.com/:f:/g/personal/myh7_psu_edu/Ei7QJqnmaCBGipPM4uMzrusBjw_hUwo0KfIDBe-0UTYyMw) folder. If you simply wish to use these simulated catalogs as examples of our models, then no installation is required! Simply download any of these tables and use them for your own science. To be able to use them, you must understand that we provide two types of catalogs:

* *Physical catalog:* a set of intrinsic, physical planetary systems (before any observations; contains properties like the true orbital periods, planet radii, etc.)
* *Observed catalog:* a set of transiting and detected planet candidates derived from a *physical catalog* (after a Kepler-like mission; contains properties like the measured orbital periods, transit depths, etc.)

Refer to the README of the branch specific to each paper for complete details on what each set of catalogs contains.



## How do I simulate my own (physical and observed) catalogs?

### Installation:

* You will need to first install the [ExoplanetsSysSim](https://github.com/ExoJulia/ExoplanetsSysSim.jl) package and set up some additional repositories; follow the instructions listed in the README of that page.
* Clone this repository.
```
git clone git@github.com:ExoJulia/SysSimExClusters.git
```
* Switch to the branch of this repository containing the model you want to simulate from. For example, to simulate models from the most recent paper, do:
```
git checkout He_et_al_2020b
```

### Usage:

Refer to the README of the branch containing the model you want to simulate from for steps.



## How do I make plots similar to those in the papers?

While the core ExoplanetsSysSim and SysSimExClusters code is written in Julia, almost all of the figures produced for the paper are generated from Python (3.7) code written by Matthias He. These Python scripts are accessible in the [SysSim_Plotting](https://github.com/hematthi/SysSim_Plotting) package.



## What if I need help?

Feel free to email Matthias He at myh7@psu.edu!
