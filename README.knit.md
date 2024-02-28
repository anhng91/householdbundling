---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# familyenrollment

<!-- badges: start -->
<!-- badges: end -->

This package includes all functions used in data cleaning and estimation of "Household Bundling in Health Insurance". The latest version of the paper can be found here <https://www.andrew.cmu.edu/user/anhnguye/Household-Bundling.pdf>

## Installation
This package is not available on CRAN. To install this package, first navigate the directory to `familyenrollment`, then run `devtools::install()`. 

## List of functions 

### CARA\_fun\_Gauss  
This function computes the following integral using Gauss-Laguerre approximation: 
    \[
      \mathbb{E}\left[\frac{R^{1 - \omega} - 1}{1 - \omega}\right]
    \]
when $\omega \sim \mathcal{N}(\mu, \sigma^2)$ truncated at 0. 


### second_taylor_CARA
This function computes the following expectation w.r.t $U$ and $r$, where $r$ is a truncated normal distribution with mean $mu$ and standard deviation $\sigma$, and $U$ is uncorrelated with $r$. Let $\phi$ denote the standard normal distribution's pdf and $\Phi$ denote its CDF. 
\[
  \begin{align*}
    \mathbb{E}\left[ - \exp(-r U) \right] & \approx  \mathbb{E}\left[ - \exp(-r \mathbb{E}\left[U\right]) \right] +  \mathbb{E}\left[ - \frac{r^2}{2} \exp(-r \mathbb{E}\left[U\right]) \left(U - \mathbb{E}\left[U\right] \right)^2  \right] \\
    & = \left[ - \mathbb{E}\exp(-r \mathbb{E}\left[U\right]) \left(1 + \frac{r^2}{2} Var(U) \right) \right] \\
    & = -\exp \left(\frac{\sigma^2 (\mathbb{E} U)^2}{2} + \mu \mathbb{E} U  \right) \left[ 1 + \frac{Var(U)}{2} \sigma^2  \left[1 + \frac{-\mu/\sigma \phi(-\mu/\sigma)}{1 - \Phi(-\mu/\sigma)} - \left(\phi(-\mu/\sigma)/(1 - \Phi(-\mu/\sigma))\right)^2 \right] \right]
  \end{align*} 
\]
The arguments of `second_taylor_CARA` is `a`, which is $\mathbb{E}(U)$, `b`, which is $Var(U)$, and `mu` and `sigma`. 

