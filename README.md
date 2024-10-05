<a name="readme-top"></a>

[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]

# spatialGLM_stacking: Bayesian Inference For Geostatistical Count Data Using Predictive Stacking <img src="plots/logo4_ubuntu.png" align="right" height="120" alt="" />

>[!IMPORTANT]
>The R package [**spStack 1.0.0**](https://span-18.github.io/spStack-dev/) has been
recently released on [CRAN](https://CRAN.R-project.org/package=spStack). See the
website [spStack: Byaesian Geostatistics Using Predictive Stacking](https://span-18.github.io/spStack-dev/articles/spStack.html)
for a brief tutorial consisting of illustrative examples demonstrating the main 
functions of the package. Currently, the package implements predictive stacking 
for Bayesian hierarchical spatial models for Gaussian, Poisson, binomial and 
binary data. The next minor release of spStack is aimed at including more 
functions to implement different spatial-temporal models.

***

This repository contains code to implement different analyses, as it appears in the manuscript "Bayesian Inference for Spatial-temporal Non-Gaussian Data Using Predictive Stacking".

Soumyakanti Pan, Lu Zhang, Jonathan R. Bradley, and Sudipto Banerjee. 2024. _Bayesian Inference for Spatial-temporal Non-Gaussian Data Using Predictive Stacking._ https://arxiv.org/abs/2406.04655.

Implemented models include Bayesian spatial and spatial-temporal (discrete as well as continuous time domains) regressions on Poisson and binomial point-referenced count data with application to modeling avian counts recorded in the [North American Breeding Bird Survey](https://www.usgs.gov/data/2022-release-north-american-breeding-bird-survey-dataset-1966-2021).

## Authors

| Name   | Email       |              |
|:------ |:----------- | :----------- |
| Soumyakanti Pan | span18@ucla.edu | PhD student, UCLA Biostatistics |

## Instructions

To reproduce results or to try out the example code, issue
```bash
git clone https://github.com/SPan-18/spatialGLM_stacking.git
```
on your personal computer or, remote server to clone this repository.

The following instructions will assume installation and use of the R statistical environment and installation of some standard packages. Details on the required packages are given in the document `vignette/vignette.html`. Install them along with the packages `rmarkdown` and `bookdown`. To run the example code in `vignette.Rmd`, run the following commands. 
```bash
cd PATH-TO-CLONED-RESPOSITORY
cd vignette
make
```
If you intend to run on a server, you may need to load some necessary module files before issuing *make*. For example, if you are working on UCLA's Hoffman2 Cluster, then run the following first.
```bash
module load R/4.2.2
module load pandoc/2.17.1.1
module load gcc/11.3.0
```
If you are installing the package `ggpubr` on Hoffman2, then you may also need to include the module files *cmake/3.19.5* and *nlopt/2.7.1*.

The Makefile runs the vignette.Rmd based on cached output. For a fresh run (which may take a very long time!) of the Rmd file, the user is required to toggle the cache options suitably in the respective code chunk headers.

The directory `data` contains synthetic datasets used in different simulation experiments as well as a dataset on migrant bird sightings from the [North American Breeding Bird Survey](https://www.usgs.gov/data/2022-release-north-american-breeding-bird-survey-dataset-1966-2021) (2010-19). The directory `src` contains all R functions required to implement our algorithm. The directory `test` contains executable `.R` scripts, which the user needs to source to run them. For example, to run the script *test_stack_poisson.R*, on Mac, run the following in terminal, or "source" the script from RStudio.
```bash
cd PATH-TO-CLONED-RESPOSITORY
cd test
Rscript test_stack_poisson.R
```

### Directory structure

```bash
.
├── ... (license, README, etc.)
├── data
│   └── ... (contains different datasets)
├── hoffman2
│   └── ... (scripts to run on cluster)
├── plots
│   └── ... (contains useful plots)
├── src
│   └── ... (contains all necessary functions)
├── test
│   └── ... (sourceable R scripts)
└── vignette
    ├── Makefile
    ├── cache
    │   └── ... (Rmd outputs cache)
    ├── refs.bib
    ├── vignette.Rmd
    └── vignette.html

```

Licensing
---------
* Code &copy; 2024, Soumyakanti Pan, licensed under [CC0 1.0 Universal](http://creativecommons.org/ns#).

<!-- MARKDOWN LINKS & IMAGES -->
[contributors-shield]: https://img.shields.io/github/contributors/SPan-18/spatialGLM_stacking.svg?style=for-the-badge
[contributors-url]: https://github.com/SPan-18/spatialGLM_stacking/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/SPan-18/spatialGLM_stacking.svg?style=for-the-badge
[forks-url]: https://github.com/SPan-18/spatialGLM_stacking/network/members
[license-shield]: https://img.shields.io/github/license/SPan-18/spatialGLM_stacking?label=license&style=for-the-badge
[license-url]: https://github.com/SPan-18/spatialGLM_stacking/blob/master/LICENSE
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://www.linkedin.com/in/soumyakanti-pan-9b660b145/