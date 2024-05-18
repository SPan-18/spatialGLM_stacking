<a name="readme-top"></a>
[![CC0 1.0 Universal][license-shield]][https://opensource.org/licenses/BSD-3-Clause]
[![LinkedIn][linkedin-shield]][https://www.linkedin.com/in/soumyakanti-pan-9b660b145/]

# spatialGLM_stacking: Predictive Stacking in Bayesian Hierarchical Models for Dependent Data from the Natural Exponential Family

We develop Bayesian predictive stacking for geostatistical models involving non-Gaussian responses. We carefully build upon the distribution theory proposed by Bradley et. al. (2023) which assumes the data to be distributed as some member of the exponential family and provides a straightforward strategy to sample exactly from the posterior distribution of the latent processes and the model parameters conditional on certain hyperparameters (including spatial process parameters) hence overcoming the computational burden of traditional Markov Chain Monte Carlo (MCMC) algorithms. Following the above strategy, the posterior samples are often found to be sensitive towards the choice of these hyperparameters which often remains elusive in addition to being weakly identifiable. Under this setting, stacking of predictive densities, as a model averaging procedure can prove to be effective as we combine the inference by stacking these individual models on a grid of model hyperparameters. 
