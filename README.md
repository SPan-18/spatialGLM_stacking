<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a name="readme-top"></a>
<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Don't forget to give the project a star!
*** Thanks again! Now go create something AMAZING! :D
-->



<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]



<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/SPan-18/spatialGLM_stacking">
    <img src="images/logo.png" alt="Logo" width="80" height="80">
  </a>

<h3 align="center">project_title</h3>

  <p align="center">
    project_description
    <br />
    <a href="https://github.com/SPan-18/spatialGLM_stacking"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/SPan-18/spatialGLM_stacking">View Demo</a>
    ·
    <a href="https://github.com/SPan-18/spatialGLM_stacking/issues/new?labels=bug&template=bug-report---.md">Report Bug</a>
    ·
    <a href="https://github.com/SPan-18/spatialGLM_stacking/issues/new?labels=enhancement&template=feature-request---.md">Request Feature</a>
  </p>
</div>


We develop Bayesian predictive stacking for geostatistical models involving non-Gaussian responses. We carefully build upon the distribution theory proposed by Bradley et. al. (2023) which assumes the data to be distributed as some member of the exponential family and provides a straightforward strategy to sample exactly from the posterior distribution of the latent processes and the model parameters conditional on certain hyperparameters (including spatial process parameters) hence overcoming the computational burden of traditional Markov Chain Monte Carlo (MCMC) algorithms. Following the above strategy, the posterior samples are often found to be sensitive towards the choice of these hyperparameters which often remains elusive in addition to being weakly identifiable. Under this setting, stacking of predictive densities, as a model averaging procedure can prove to be effective as we combine the inference by stacking these individual models on a grid of model hyperparameters. 
