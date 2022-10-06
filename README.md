# qsm2r: Import and analyze QSMs from TreeQSM in R

## Description

`qsm2r` is an R package to read in and analyze QSMs (quantitative structure models) created in Matlab using <a href = "https://github.com/InverseTampere/TreeQSM">TreeQSM</a>.

## Installation from source

```R
# install dependencies
install.packages(c("data.table", "R.matlab", "rgl"))

# install package from github
remotes::install_github("zoeschindler/qsm2r")
```

## Usage

```R
# load qsm
file_path <- system.file("extdata", "medium.mat", package="qsm2r")
qsm <- readQSM(file_path)

# get summary
print(qsm)

# check if branch & overview are up-to-date with cylinder
checkQSM(qsm)

# update branch & overview from cylinder
updateQSM(qsm)

# plot qsm
plot(qsm)
```

## About

Author: Zoe Schindler, <a href = "https://www.iww.uni-freiburg.de/">Chair of Forest Growth and Dendroecology</a>, <a href = "https://uni-freiburg.de/">University of Freiburg</a>