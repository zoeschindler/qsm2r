# qsm2r: Import and analyze QSMs from TreeQSM in R 

## Description <img src="https://github.com/zoeschindler/qsm2r/blob/main/inst/figures/logo.png" align="right" width = 290/>

`qsm2r` is an R package to read in and analyze QSMs (quantitative structure models) created in Matlab using <a href = "https://github.com/InverseTampere/TreeQSM">TreeQSM</a>. QSMs stored in `*.mat` files can be read in using `readQSM()`. After changing the cylinder structure of QSMs, the `updateQSM()` function can be used to update the remaining data in the QSMs (e.g. branch data, crown statistics, ...) based on the edited data. Various summaries of the cylinder and branch data can be derived using a set of `summary_cylinder_*()` and `summary_branch_*()` functions. Different types of pruning treatments can be simulated using several `pruning_*()` functions.

## Installation from source

To install the package from github, the package `remotes` is required.

```R
# install package from github
remotes::install_github("zoeschindler/qsm2r")
```

## Usage

```R
# load package
library(qsm2r)

# load qsm
file_path <- system.file("extdata", "QSM_Juglans_regia_M.mat", package="qsm2r")
qsm <- readQSM(file_path)

# get summary
print(qsm)

# show qsm
plot(qsm)

# pruning
qsm <- pruning_conventional(qsm, threshold_m = 3, method = "height", remove = TRUE)

# show pruned qsm
plot(qsm)

# check if branch & overview are up-to-date with cylinder
checkQSM(qsm)

# update branch & overview from cylinder
qsm <- updateQSM(qsm)

# check if branch & overview are up-to-date with cylinder
checkQSM(qsm)

# get summary
print(qsm)
```

## About

Author: Zoe Schindler, <a href = "https://www.iww.uni-freiburg.de/">Chair of Forest Growth and Dendroecology</a>, <a href = "https://uni-freiburg.de/">University of Freiburg</a>
