# Load required libraries
library("MASS")
suppressPackageStartupMessages(library("Rfast"))
# library("fastmatrix")
library("parallel")
# library("progress")
suppressPackageStartupMessages(library("loo"))
suppressWarnings(library(CVXR, warn.conflicts=FALSE))
if(require(Rmosek)){
  suppressPackageStartupMessages(library(Rmosek))
}else{
  message("MOSEK not installed. Use ECOS/SCS solvers. To install R-interface for MOSEK,
  visit https://docs.mosek.com/latest/rmosek/install-interface.html.")
}
suppressPackageStartupMessages(library("geoR"))

suppressPackageStartupMessages(library("dplyr"))
library("ggplot2")
library("MBA")
library("latex2exp")
library("knitr")

# Source required scripts

source("../src/distributions.R")
source("../src/functions.R")
source("../src/projection.R")
source("../src/cholesky.R")
source("../src/posterior.R")
source("../src/predictive.R")
source("../src/sim_EFdata.R")
source("../src/stacking_weights.R")
source("../src/cm_stacking.R")
source("../src/pointrefplot.R")
