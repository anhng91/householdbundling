library(knitr)
library(tidyverse)
library(lfe)
library(sandwich)
library(plm)
library(stargazer)
library(parallel)
library(randomForest)
library(randtoolbox)
library(familyenrollment)

# Estimate the probability of getting sick 
estimation = function(numcores, sample_index) {
  data_here = do.call('cbind', data_hh_list[sample_index])
  sick_parameters = optim(rep(0, ncol(data_here)), fn = function(x) llh_sick(x, data_here), gr = function(x) grad_llh_sick(x, data_here), control=list(maxit = 1e4), method='BFGS')
  xi_parameters = optim(rep(0, 2 * ncol(data_here)), fn = function(x) llh_xi(x, data_here), gr = function(x) grad_llh_xi(x, data_here), control=list(maxit = 1e4), method='BFGS')
  
}
