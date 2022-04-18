
<!-- README.md is generated from README.Rmd. Please edit that file -->

# methanogenlab 
<!-- badges: start -->

[![Git\_Hub\_Version](https://img.shields.io/github/r-package/v/mankeldy/Methanogen_package?label=Github&logo=Github)](/commits)
<!-- badges: end -->

## About
**Note: This model is still work in progress, so some functions may become depreciated as new updates are released.**
## Installation

You can install **[methanogenlab]()** from github with the
devtools package.

``` r
# install.packages("devtools")
devtools::install_github("hoehler-labgroup/methanogenlab")
```
## Usage
There are currently two available models within **[methanogenlab]()**, `methanogenesis.time` and `methanogenesis.DIC`. The underlying computations for each model are quite similar and will produce similar results; the main distinction is whether the model steps in time or through DIC consumed. For the latter, each step consumes a constant amount DIC (and hydrogen) leading to the production of methane, changing the overall chemistry of the closed system and the amount of biomass generated. For `methanogenesis.time`, however, Monod kinetics dictate the number of new cells generated at each step and can be translated to an amount of DIC necessary for that growth (i.e. a DIC step that will vary as the reaction progresses).

Below is an example using `methanogenesis.time` with an initial concentrations [CH4], [H2], and [DIC] at 1 &#956;M, 500 &#956;M, and 3.2 mM, respectively. 
``` r
library(methanogenlab)
library(CHNOSZ) #some parts of CHNOSZ are not imported with methanogenlab, causing certain parts of the model to break. 
                #This will hopefully be corrected for in future updates. You can install CHNOSZ from CRAN.

model <- methanogensis.time(
CH4.initial= 1e-6,
H2.initial= 5e-4,
is.H2.limiting=FALSE,
Ks.H2=10e-6,
DIC.initial= 3.2e-3,
pH.initial= 7.5,
is.CO2.limiting=TRUE,
Ks.CO2=20e-6,
umax=0.6,
temperature= 313.15,
VolumeSolution = 80e-3,
VolumeHeadspace = 20e-3,
time.step=0.1,
total.time=30,
inoculum.cell.number = 1e6,
biomass.yield=2.4,
carbon.fraction=0.44,
cell.weight=30e-15)
)
```

Documentation for all functions within methanogenlab can be found in the [manual](methanogenlab_0.1.2.pdf).

------------------------------------------------------------------------

## Interactive plot

You can access an online application for methanogenlab [here](https://hoehler-labgroup.shinyapps.io/methanogenlab) **(CURRENTLY DEPRECIATED, NOT UPDATED TO v0.1.2)**. The repository is located in [hoehler-labgroup/shinyApps](https://github.com/hoehler-labgroup/shinyApps)

