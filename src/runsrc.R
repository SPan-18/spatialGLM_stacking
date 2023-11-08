# Load required libraries
library("MASS")
suppressPackageStartupMessages(library("Rfast"))
library("fastmatrix")
library("parallel")
library("progress")
suppressPackageStartupMessages(library("loo"))
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
source("../src/cm_stacking.R")
source("../src/pointrefplot.R")