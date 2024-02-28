
<!-- README.md is generated from README.Rmd. Please edit that file -->

# familyenrollment

<!-- badges: start -->

<!-- badges: end -->

This package includes all functions used in data cleaning and estimation
of “Household Bundling in Health Insurance”. The latest version of the
paper can be found here
<https://www.andrew.cmu.edu/user/anhnguye/Household-Bundling.pdf>

## Installation

This package is not available on CRAN. To install this package, first
navigate the directory to `familyenrollment`, then run
`devtools::install()`.

## List of functions

### CARA\_fun\_Gauss

This function computes the following integral using Gauss-Laguerre
approximation:   
![&#10; \\mathbb{E}\\left\[\\frac{R^{1 - \\omega} - 1}{1 -
\\omega}\\right\]&#10;
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%20%20%20%20%20%20%5Cmathbb%7BE%7D%5Cleft%5B%5Cfrac%7BR%5E%7B1%20-%20%5Comega%7D%20-%201%7D%7B1%20-%20%5Comega%7D%5Cright%5D%0A%20%20%20%20
"
      \\mathbb{E}\\left[\\frac{R^{1 - \\omega} - 1}{1 - \\omega}\\right]
    ")  
when ![\\omega \\sim \\mathcal{N}(\\mu,
\\sigma^2)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Comega%20%5Csim%20%5Cmathcal%7BN%7D%28%5Cmu%2C%20%5Csigma%5E2%29
"\\omega \\sim \\mathcal{N}(\\mu, \\sigma^2)") truncated at 0.

### second\_taylor\_CARA

This function computes the following expectation w.r.t
![U](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;U
"U") and
![r](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;r
"r"), where
![r](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;r
"r") is a truncated normal distribution with mean
![mu](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;mu
"mu") and standard deviation
![\\sigma](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Csigma
"\\sigma"), and
![U](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;U
"U") is uncorrelated with
![r](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;r
"r"). Let
![\\phi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cphi
"\\phi") denote the standard normal distribution’s pdf and
![\\Phi](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CPhi
"\\Phi") denote its CDF.   
![&#10; \\begin{align\*}&#10; \\mathbb{E}\\left\[ - \\exp(-r U)
\\right\] & \\approx \\mathbb{E}\\left\[ - \\exp(-r
\\mathbb{E}\\left\[U\\right\]) \\right\] + \\mathbb{E}\\left\[ -
\\frac{r^2}{2} \\exp(-r \\mathbb{E}\\left\[U\\right\]) \\left(U -
\\mathbb{E}\\left\[U\\right\] \\right)^2 \\right\] \\\\&#10; & =
\\left\[ - \\mathbb{E}\\exp(-r \\mathbb{E}\\left\[U\\right\]) \\left(1 +
\\frac{r^2}{2} Var(U) \\right) \\right\] \\\\&#10; & \\propto -
\\int\_0^{\\infty} \\exp(-r \\mathbb{E}\\left\[U\\right\]) \\left(1 +
\\frac{r^2}{2} Var(U) \\right) \\exp\\left\[-0.5 \\frac{(r -
\\mu)^2}{\\sigma^2} \\right\] \\\\&#10; & = -\\exp
\\left(\\frac{\\sigma^2 (\\mathbb{E} U)^2}{2} - \\mu \\mathbb{E} U
\\right) \\frac{1 - \\Phi\\left( -(\\mu - a \\sigma^2)/\\sigma
\\right)}{1 - \\Phi(-\\mu/\\sigma)} \\\\ &#10; & \\left\[ 1 +
\\frac{Var(U)}{2} \\left( \\sigma^2 \\left\[1 - \\frac{-\\mu/\\sigma
\\phi(-\\mu/\\sigma)}{1 - \\Phi(-\\mu/\\sigma)} -
\\left(\\phi(-\\mu/\\sigma)/(1 - \\Phi(-\\mu/\\sigma))\\right)^2
\\right\] + \\left(\\mu + \\phi(-\\mu/\\sigma)/(1 -
\\Phi(-\\mu/\\sigma))\\sigma \\right)^2 \\right) \\right\]&#10;
\\end{align\*}
&#10;](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%20%20%5Cbegin%7Balign%2A%7D%0A%20%20%20%20%5Cmathbb%7BE%7D%5Cleft%5B%20-%20%5Cexp%28-r%20U%29%20%5Cright%5D%20%26%20%5Capprox%20%20%5Cmathbb%7BE%7D%5Cleft%5B%20-%20%5Cexp%28-r%20%5Cmathbb%7BE%7D%5Cleft%5BU%5Cright%5D%29%20%5Cright%5D%20%2B%20%20%5Cmathbb%7BE%7D%5Cleft%5B%20-%20%5Cfrac%7Br%5E2%7D%7B2%7D%20%5Cexp%28-r%20%5Cmathbb%7BE%7D%5Cleft%5BU%5Cright%5D%29%20%5Cleft%28U%20-%20%5Cmathbb%7BE%7D%5Cleft%5BU%5Cright%5D%20%5Cright%29%5E2%20%20%5Cright%5D%20%5C%5C%0A%20%20%20%20%26%20%3D%20%5Cleft%5B%20-%20%5Cmathbb%7BE%7D%5Cexp%28-r%20%5Cmathbb%7BE%7D%5Cleft%5BU%5Cright%5D%29%20%5Cleft%281%20%2B%20%5Cfrac%7Br%5E2%7D%7B2%7D%20Var%28U%29%20%5Cright%29%20%5Cright%5D%20%5C%5C%0A%20%20%20%20%26%20%5Cpropto%20-%20%5Cint_0%5E%7B%5Cinfty%7D%20%5Cexp%28-r%20%5Cmathbb%7BE%7D%5Cleft%5BU%5Cright%5D%29%20%5Cleft%281%20%2B%20%5Cfrac%7Br%5E2%7D%7B2%7D%20Var%28U%29%20%5Cright%29%20%5Cexp%5Cleft%5B-0.5%20%5Cfrac%7B%28r%20-%20%5Cmu%29%5E2%7D%7B%5Csigma%5E2%7D%20%5Cright%5D%20%5C%5C%0A%20%20%20%20%26%20%3D%20-%5Cexp%20%5Cleft%28%5Cfrac%7B%5Csigma%5E2%20%28%5Cmathbb%7BE%7D%20U%29%5E2%7D%7B2%7D%20-%20%5Cmu%20%5Cmathbb%7BE%7D%20U%20%20%5Cright%29%20%5Cfrac%7B1%20-%20%5CPhi%5Cleft%28%20-%28%5Cmu%20-%20a%20%5Csigma%5E2%29%2F%5Csigma%20%5Cright%29%7D%7B1%20-%20%5CPhi%28-%5Cmu%2F%5Csigma%29%7D%20%5C%5C%20%0A%20%20%20%20%26%20%5Cleft%5B%201%20%2B%20%5Cfrac%7BVar%28U%29%7D%7B2%7D%20%5Cleft%28%20%5Csigma%5E2%20%20%5Cleft%5B1%20-%20%5Cfrac%7B-%5Cmu%2F%5Csigma%20%5Cphi%28-%5Cmu%2F%5Csigma%29%7D%7B1%20-%20%5CPhi%28-%5Cmu%2F%5Csigma%29%7D%20-%20%5Cleft%28%5Cphi%28-%5Cmu%2F%5Csigma%29%2F%281%20-%20%5CPhi%28-%5Cmu%2F%5Csigma%29%29%5Cright%29%5E2%20%5Cright%5D%20%2B%20%5Cleft%28%5Cmu%20%2B%20%5Cphi%28-%5Cmu%2F%5Csigma%29%2F%281%20-%20%5CPhi%28-%5Cmu%2F%5Csigma%29%29%5Csigma%20%5Cright%29%5E2%20%5Cright%29%20%5Cright%5D%0A%20%20%5Cend%7Balign%2A%7D%20%0A
"
  \\begin{align*}
    \\mathbb{E}\\left[ - \\exp(-r U) \\right] & \\approx  \\mathbb{E}\\left[ - \\exp(-r \\mathbb{E}\\left[U\\right]) \\right] +  \\mathbb{E}\\left[ - \\frac{r^2}{2} \\exp(-r \\mathbb{E}\\left[U\\right]) \\left(U - \\mathbb{E}\\left[U\\right] \\right)^2  \\right] \\\\
    & = \\left[ - \\mathbb{E}\\exp(-r \\mathbb{E}\\left[U\\right]) \\left(1 + \\frac{r^2}{2} Var(U) \\right) \\right] \\\\
    & \\propto - \\int_0^{\\infty} \\exp(-r \\mathbb{E}\\left[U\\right]) \\left(1 + \\frac{r^2}{2} Var(U) \\right) \\exp\\left[-0.5 \\frac{(r - \\mu)^2}{\\sigma^2} \\right] \\\\
    & = -\\exp \\left(\\frac{\\sigma^2 (\\mathbb{E} U)^2}{2} - \\mu \\mathbb{E} U  \\right) \\frac{1 - \\Phi\\left( -(\\mu - a \\sigma^2)/\\sigma \\right)}{1 - \\Phi(-\\mu/\\sigma)} \\\\ 
    & \\left[ 1 + \\frac{Var(U)}{2} \\left( \\sigma^2  \\left[1 - \\frac{-\\mu/\\sigma \\phi(-\\mu/\\sigma)}{1 - \\Phi(-\\mu/\\sigma)} - \\left(\\phi(-\\mu/\\sigma)/(1 - \\Phi(-\\mu/\\sigma))\\right)^2 \\right] + \\left(\\mu + \\phi(-\\mu/\\sigma)/(1 - \\Phi(-\\mu/\\sigma))\\sigma \\right)^2 \\right) \\right]
  \\end{align*} 
")  
The arguments of `second_taylor_CARA` is `a`, which is
![\\mathbb{E}(U)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbb%7BE%7D%28U%29
"\\mathbb{E}(U)"), `b`, which is
![Var(U)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Var%28U%29
"Var(U)"), and `mu` and `sigma`.
