# Load required libraries
library(MASS)
library(fastmatrix)
library(progress)
library(dplyr)
library(ggplot2)
library(akima)

# Source required scripts

source("../src/distributions.R")
source("../src/functions.R")
source("../src/projection.R")
source("../src/posterior.R")
source("../src/predictive.R")
source("../src/sim_countdata.R")
source("../src/cm_stacking.R")