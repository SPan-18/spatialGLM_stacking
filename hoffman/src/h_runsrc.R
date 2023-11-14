# Load required libraries
library("MASS")
suppressPackageStartupMessages(library("Rfast"))
# library("fastmatrix")
# library("parallel")
# library("progress")
suppressPackageStartupMessages(library("loo"))
suppressPackageStartupMessages(library("geoR"))

suppressPackageStartupMessages(library("dplyr"))
library("ggplot2")
library("MBA")
library("latex2exp")
library("knitr")

# Source required scripts

source("../src/h_distributions.R")
source("../src/h_functions.R")
source("../src/h_projection.R")
source("../src/h_cholesky.R")
source("../src/h_posterior.R")
source("../src/h_predictive.R")
source("../src/h_sim_EFdata.R")
source("../src/h_cm_stacking.R")
source("../src/h_pointrefplot.R")
