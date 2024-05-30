args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) { 
  numcores = 24;
  job_index = 1;  
} else {
  job_index = as.numeric(args[1]);
  numcores = as.numeric(args[2]); 
}
if (dir.exists('work/teckyongtan/tecktan/Oil/Data/R_lib')) {
  .libPaths('work/teckyongtan/tecktan/Oil/Data/R_lib')
}
options("install.lock"=FALSE)
library(knitr)
library(tidyverse)
library(lfe)
library(sandwich)
library(plm)
library(stargazer)
library(parallel)
library(randomForest)
library(randtoolbox)

# setwd('./familyenrollment')
devtools::install(upgrade='never')
library(familyenrollment)

message('constructing the list of Compulsory households')

Com_HH_list_index = lapply(1:length(data_hh_list), function(hh_index) {
  data = data_hh_list[[hh_index]]; 
  if (0 == nrow(data %>% filter(Bef_sts + Com_sts + Std_w_ins == 0))) {
    return(hh_index);
  }
  else {
    return(NA); 
  }
}) %>% unlist(); 

Com_HH_list_index = Com_HH_list_index[!(is.na(Com_HH_list_index))]

message('constructing the list of Voluntary households')
Vol_HH_list_index = lapply(1:length(data_hh_list), function(hh_index) {
  data = data_hh_list[[hh_index]]; 
  if (nrow(data) > nrow(data %>% filter(Bef_sts + Com_sts + Std_w_ins == 1))) {
    return(hh_index);
  }
  else {
    return(NA); 
  }
}) %>% unlist(); 

Vol_HH_list_index = Vol_HH_list_index[!(is.na(Vol_HH_list_index))] 

full_index = 1:length(data_hh_list); 

message('Identify the sample for theta')
sample_identify_theta = full_index[which((lapply(full_index, function(sample_index_i) {
      data_mini = data_hh_list[[sample_index_i]]
      if ((data_mini$Year[1] == 2008) & (data_mini$Income[1] < 0) & (data_mini$HHsize_s[1] == 0)) {
        return(1)
      } else {
        return(0)
      }
    }) %>% unlist()) == 1)]

sample_no_sick = full_index[which((lapply(full_index, function(sample_index_i) {
      data_mini = data_hh_list[[sample_index_i]]
      if (((data_mini$HHsize_s[1] == 0) & (sum(data_mini$sick_dummy == 0)))) {
        return(1)
      } else {
        return(0)
      }
    }) %>% unlist()) == 1)]

message('Identify the sample for preference')
sample_identify_pref = lapply(Com_HH_list_index, function(x) ifelse(x %in% c(sample_identify_theta, sample_no_sick), NA, x)) %>% unlist()
sample_identify_pref = sample_identify_pref[!(is.na(sample_identify_pref))]

# Bootstrapping indices 
message('bootstrapping indices')
set.seed(job_index);
sample_index = sample(1:length(data_hh_list), length(data_hh_list), replace=TRUE)
sample_r_theta = Vol_HH_list_index[which(!(is.na(lapply(Vol_HH_list_index, function(x) ifelse(nrow(data_hh_list[[x]]) <= 4, x, NA)))))]
sample_r_theta = sample(sample_r_theta, 10, replace=TRUE)
sample_identify_pref = sample(sample_identify_pref, 10,  replace=TRUE)
sample_identify_theta = sample(sample_identify_theta, 50,  replace=TRUE)
# sample_r_theta = sample(sample_r_theta, length(sample_r_theta), replace=TRUE)
# sample_identify_pref = sample(sample_identify_pref, length(sample_identify_pref), replace=TRUE)
# sample_identify_theta = sample(sample_identify_theta, length(sample_identify_theta), replace=TRUE)


message('merging household frames')
data = do.call('rbind', data_hh_list[sample_index])

# Estimate the probability of getting sick 
message('computing the sick parameters')
sick_parameters = optim(rep(0, ncol(var_ind(data_hh_list[[1]]))), fn = function(x) llh_sick(x, data), gr = function(x) grad_llh_sick(x, data), control=list(maxit = 1e4), method='BFGS')

message('computing the coverage parameters')
xi_parameters = optim(rep(0, 2 * ncol(var_ind(data_hh_list[[1]]))), fn = function(x) llh_xi(x, data), gr = function(x) grad_llh_xi(x, data), control=list(maxit = 1e4), method='BFGS')

param_trial = transform_param(return_index=TRUE, init=TRUE);
transform_param_trial = transform_param(param_trial, return_index=TRUE)
x_transform = transform_param_trial

var_list = c('sigma_thetabar', 'beta_theta', 'beta_theta_ind');
active_index = unlist(transform_param_trial[[2]][var_list]);
active_index = setdiff(unlist(transform_param_trial[[2]][var_list]),tail(transform_param_trial[[2]][['beta_theta']], n=4));

message('Preparing data for estimation')
data_hh_list_theta = list()
for (index in 1:length(sample_identify_theta)) {
  message(paste0('constructing data for index', index))
  data_hh_list_theta[[index]] = household_draw_theta_kappa_Rdraw(hh_index=sample_identify_theta[index], param=transform_param_trial[[1]], n_draw_halton = 1000, n_draw_gauss = 10, sick_parameters, xi_parameters, short=FALSE)
}


message('Estimating theta')

if (Sys.info()[['sysname']] == 'Windows') {
  cl = makeCluster(numcores);
  clusterEvalQ(cl, library('tidyverse'))
  clusterEvalQ(cl, library('familyenrollment'))
  clusterExport(cl, c('active_index', 'param_trial', 'data_hh_list_theta'))
}

n_draw_halton = 100; 
if (Sys.info()[['sysname']] == 'Windows') {
  clusterExport(cl, c('active_index', 'param_trial', 'transform_param_trial', 'sample_identify_pref', 'xi_parameters', 'sick_parameters'))
  clusterExport(cl,'n_draw_halton')
}


message('compute covariates before estimation')

if (Sys.info()[['sysname']] == 'Windows') {
  data_hh_list_pref = parLapply(cl, sample_identify_pref,function(index) {
    output = tryCatch(household_draw_theta_kappa_Rdraw(hh_index=index, param=transform_param_trial[[1]], n_draw_halton = n_draw_halton, n_draw_gauss = 10, sick_parameters, xi_parameters, short=FALSE),error=function(e) e)
    return(output)
  })
} else {
  data_hh_list_pref = mclapply(sample_identify_pref, function(index) tryCatch(household_draw_theta_kappa_Rdraw(hh_index=index, param=transform_param_trial[[1]], n_draw_halton = n_draw_halton, n_draw_gauss = 10, sick_parameters, xi_parameters, short=FALSE), error=function(e) e), mc.cores=numcores)
}

mat_M = do.call('rbind', lapply(data_hh_list_pref, function(x) {
      return(cbind(x$data$M_expense, x$data$M_expense^2))
    }))
mat_Year = do.call('rbind', lapply(data_hh_list_pref, function(x) {
        return(cbind(x$data$Year == 2004, x$data$Year == 2006, x$data$Year == 2010, x$data$Year == 2012))
      }))
X_ind_pref = do.call('rbind', lapply(data_hh_list_pref, function(x) x$X_ind))
X_hh_pref = do.call('rbind', lapply(data_hh_list_pref, function(x) x$X_hh))
X_ind_pref_with_year = cbind(X_ind_pref, mat_Year)
mat_M_rtheta = do.call('rbind', lapply(data_hh_list[sample_r_theta], function(x) {
    return(cbind(x$M_expense, x$M_expense^2))
  }))

mat_Y_rtheta = do.call('c', lapply(data_hh_list[sample_r_theta], function(x) x$Income[1]))

mat_Y_rtheta = pnorm(mat_Y_rtheta, mean = mean(mat_Y_rtheta), sd = sd(mat_Y_rtheta))

full_insurance_indicator = do.call('c', lapply(data_hh_list[sample_r_theta], function(x) {
    return(x$HHsize_s[1] == x$N_vol[1])
  }))

full_insurance_indicator_ind_level = do.call('c', lapply(data_hh_list[sample_r_theta], function(x) {
    return(rep(x$HHsize_s[1] == x$N_vol[1], x$HHsize[1]))
  }))

X_r_all = do.call('rbind', lapply(data_hh_list[sample_r_theta], function(x) var_hh(x)))

X_hh_theta_r = do.call('rbind',lapply(sample_r_theta, function(output_hh_index) var_hh(data_hh_list[[output_hh_index]])))

n_involuntary = do.call('c', lapply(sample_r_theta, function(output_hh_index) data_hh_list[[output_hh_index]]$N_com[1] + data_hh_list[[output_hh_index]]$N_bef[1] + data_hh_list[[output_hh_index]]$N_std_w_ins[1]))

n_halton_at_r = 100; 

param_trial[x_transform[[2]]$beta_theta[1]] = 0

initial_param_trial = param_trial

compute_inner_loop = function(x_stheta, return_result=FALSE, estimate_theta=TRUE, estimate_pref=TRUE) {
  param_trial_here = initial_param_trial
  param_trial_here[transform_param_trial[[2]]$sigma_theta] = x_stheta; 
  transform_param_trial = transform_param(param_trial_here, return_index=TRUE);
  var_list = c('beta_delta', 'beta_omega', 'beta_gamma', 'sigma_delta', 'sigma_gamma', 'sigma_omega')
  x_transform = transform_param(param_trial_here, return_index=TRUE)
  print(x_transform[[1]])

  aggregate_moment_pref = function(x_transform) {
    if (Sys.info()[['sysname']] == 'Windows') {
      clusterExport(cl, c('x_transform', 'n_halton_at_r'),envir=environment())
      data_hh_list_pref = parLapply(cl, sample_identify_pref,function(index) {
        output = tryCatch(household_draw_theta_kappa_Rdraw(hh_index=index, param=x_transform[[1]], n_draw_halton = n_draw_halton, n_draw_gauss = 10, sick_parameters, xi_parameters, short=FALSE),error=function(e) e)
        return(output)
      })
      mat_YK = do.call('cbind', parLapply(cl, data_hh_list_pref, function(x) {
            output = rbind(colMeans(matrix(x$kappa_draw[[1]], nrow=n_draw_halton)), colMeans(matrix(x$kappa_draw[[1]]^2, nrow=n_draw_halton)), x$income[1], x$income[1]^2, colMeans(matrix(x$kappa_draw[[1]], nrow=n_draw_halton))*x$income[1])
            return(output)
          }))
    } else {
      data_hh_list_pref = mclapply(sample_identify_pref, function(index) tryCatch(household_draw_theta_kappa_Rdraw(hh_index=index, param=x_transform[[1]], n_draw_halton = n_draw_halton, n_draw_gauss = 10, sick_parameters, xi_parameters, short=FALSE), error=function(e) e), mc.cores=numcores)
      mat_YK = do.call('cbind', mclapply(data_hh_list_pref, function(x) {
            output = rbind(colMeans(matrix(x$kappa_draw[[1]], nrow=n_draw_halton)), colMeans(matrix(x$kappa_draw[[1]]^2, nrow=n_draw_halton)), x$income[1], x$income[1]^2, colMeans(matrix(x$kappa_draw[[1]], nrow=n_draw_halton))*x$income[1])
            return(output)
          }, mc.cores=numcores))
    }

    if (Sys.info()[['sysname']] == 'Windows') {
      clusterExport(cl, 'x_transform',envir=environment())
      moment_ineligible_hh_output = parLapply(cl, data_hh_list_pref,function(mini_data) {
        output = tryCatch(moment_ineligible_hh(mini_data, x_transform[[1]]),error=function(e) e)
        return(output)
      })
    } else {
      moment_ineligible_hh_output = mclapply(data_hh_list_pref, function(mini_data) tryCatch(moment_ineligible_hh(mini_data, x_transform[[1]]), error=function(e) e), mc.cores=numcores)
    }

    output_1 = do.call('c', lapply(moment_ineligible_hh_output, function(x) x[[1]]))
    print('predicted');print(summary(output_1))
    print('actual'); print(summary(mat_M[,1]))
    output_2 = do.call('c', lapply(moment_ineligible_hh_output, function(x) x[[2]]))
    d_output_1 = list();
    d_output_2 = list();
    for (varname in var_list) {
      if (grepl('sigma', varname)) {
        d_output_1[[varname]] = do.call('c', lapply(moment_ineligible_hh_output, function(x) x[[3]][[varname]]))
        d_output_2[[varname]] = do.call('c', lapply(moment_ineligible_hh_output, function(x) x[[4]][[varname]]))
      } else {
        d_output_1[[varname]] = do.call('cbind', lapply(moment_ineligible_hh_output, function(x) x[[3]][[paste0('mean_',varname)]]))
        d_output_2[[varname]] = do.call('cbind', lapply(moment_ineligible_hh_output, function(x) x[[4]][[paste0('mean_',varname)]]))
      }
    }

    moment = NULL
    d_moment = list()
    output = list(); 
    output[[3]] = list();
    for (moment_index in c(1:5)) {
      output[[3]][[moment_index]] = ((output_1 - mat_M[,1]) * mat_YK[moment_index,])^2
      moment[moment_index] = mean(output[[3]][[moment_index]]); 
      d_moment[[moment_index]] = 2*((output_1 - mat_M[,1]) * mat_YK[moment_index,]);
    }
    for (moment_index in c(1:5)) {
      output[[3]][[moment_index + 5]] = ((output_2 - mat_M[,2]) * mat_YK[moment_index,])^2
      moment[moment_index + 5] = mean(output[[3]][[moment_index + 5]]); 
      d_moment[[moment_index + 5]] = 2*((output_2 - mat_M[,2]) * mat_YK[moment_index,]);
    }

    d_moment = do.call('cbind', d_moment)

    output[[1]] = sum(moment); 
    active_index_pref = x_transform[[2]][var_list] %>% unlist()
    output[[2]] = rep(0, length(active_index_pref)); 

    output[[2]][x_transform[[2]][['beta_delta']]] = apply(d_output_1[['beta_delta']], 1, function(x) sum((x %*% d_moment[,1:5])/length(x))) + apply(d_output_2[['beta_delta']], 1, function(x) sum((x %*% d_moment[,6:10])/length(x)))
    output[[2]][x_transform[[2]][['beta_omega']]] = apply(d_output_1[['beta_omega']], 1, function(x) sum((x %*% d_moment[,1:5])/length(x))) + apply(d_output_2[['beta_omega']], 1, function(x) sum((x %*% d_moment[,6:10])/length(x)))
    output[[2]][x_transform[[2]][['beta_gamma']]] = apply(d_output_1[['beta_gamma']], 1, function(x) sum((x %*% d_moment[,1:5])/length(x))) + apply(d_output_2[['beta_gamma']], 1, function(x) sum((x %*% d_moment[,6:10])/length(x)))

    output[[2]][x_transform[[2]][['sigma_delta']]] =  sum((d_output_1[['sigma_delta']] %*% d_moment[,1:5])/length(active_index)) + sum((d_output_2[['sigma_delta']] %*% d_moment[,6:10])/length(active_index))
    output[[2]][x_transform[[2]][['sigma_gamma']]] =  sum((d_output_1[['sigma_gamma']] %*% d_moment[,1:5])/length(active_index)) + sum((d_output_2[['sigma_gamma']] %*% d_moment[,6:10])/length(active_index))
    output[[2]][x_transform[[2]][['sigma_omega']]] =  sum((d_output_1[['sigma_omega']] %*% d_moment[,1:5])/length(active_index)) + sum((d_output_2[['sigma_omega']] %*% d_moment[,6:10])/length(active_index))

    output[[2]] = output[[2]][active_index_pref]
    output[[3]] = do.call('cbind', output[[3]])
    return(output)
  }

  aggregate_moment_theta = function(x_transform) {
    # computing moment using realized expenses
    if (Sys.info()[['sysname']] == 'Windows') {
      clusterExport(cl, 'x_transform',envir=environment())
      moment_realized_expense = parLapply(cl, data_hh_list_theta,function(mini_data) {
        output = tryCatch(identify_theta(mini_data, x_transform[[1]], n_draw_halton = n_draw_halton),error=function(e) e)
        return(output)
      })
    } else {
      moment_realized_expense = mclapply(data_hh_list_theta, function(mini_data) tryCatch(identify_theta(mini_data, x_transform[[1]], n_draw_halton = n_draw_halton), error=function(e) e), mc.cores=numcores)
    }
    
    moment_realized_expense_val = do.call('c', lapply(moment_realized_expense, function(x) x[[1]])) %>% mean
    print(paste0('moment from theta = ', moment_realized_expense_val))
    moment_realized_expense_deriv = rep(0, length(param_trial));
    moment_realized_expense_deriv[x_transform[[2]]$beta_theta] = do.call('rbind', lapply(moment_realized_expense, function(x) x[[2]]$beta_theta)) %>% colMeans
    moment_realized_expense_deriv[x_transform[[2]]$beta_theta_ind] = do.call('rbind', lapply(moment_realized_expense, function(x) x[[2]]$beta_theta_ind)) %>% colMeans
    moment_realized_expense_deriv[x_transform[[2]]$sigma_thetabar] = do.call('c', lapply(moment_realized_expense, function(x) x[[2]]$sigma_thetabar)) %>% mean

    return(list(moment_realized_expense_val, moment_realized_expense_deriv))
  }

  active_index = x_transform[[2]][c(var_list, 'beta_theta', 'beta_theta_ind', 'sigma_thetabar')] %>% unlist()
  tol = 1e-4

  iteration = 1; 

  init_pref = aggregate_moment_pref(transform_param(param_trial_here, return_index=TRUE))
  init_theta = aggregate_moment_theta(transform_param(param_trial_here, return_index=TRUE))

  if (is.na(init_pref[[1]]) | is.nan(init_pref[[1]]) | is.na(init_theta[[1]]) | is.nan(init_theta[[1]])) {
    return(NA)
  }

  if (estimate_theta) {
    optim_pref_theta = splitfngr::optim_share(param_trial_here[active_index], function(x_pref_theta) {
      if (max(abs(x_pref_theta)) > 10) {
        return(list(NA, rep(NA, length(active_index))))
      }
      param_trial_inner_theta = param_trial_here
      param_trial_inner_theta[active_index] = x_pref_theta 
      x_transform = transform_param(param_trial_inner_theta, return_index = TRUE)
      message('computing preference moment')
      pref_moment = aggregate_moment_pref(transform_param(param_trial_inner_theta, return_index=TRUE))
      deriv_theta_index = c(x_transform[[2]]$beta_theta[1], x_transform[[2]]$beta_theta_ind[1], x_transform[[2]]$sigma_thetabar)
      if (is.na(pref_moment[[1]])) {
        return(list(NA, rep(NA, length(x_pref_theta))))
      }
      message('computing preference moment derivative')
      for (i in deriv_theta_index) {
        param_trial_i = param_trial_inner_theta; param_trial_i[i] = param_trial_inner_theta[i] + tol
        fi = aggregate_moment_pref(transform_param(param_trial_i, return_index=TRUE))
        if (i == deriv_theta_index[1]) {
          pref_moment[[2]] = c(pref_moment[[2]],(rowMeans((fi[[3]] - pref_moment[[3]])/tol) %*% X_ind_pref_with_year)/nrow(X_ind_pref))
        } else if (i == deriv_theta_index[2]) {
          pref_moment[[2]] = c(pref_moment[[2]],(rowMeans((fi[[3]] - pref_moment[[3]])/tol) %*% X_ind_pref)/nrow(X_ind_pref))
        } else {
          pref_moment[[2]] = c(pref_moment[[2]], (fi[[1]] - pref_moment[[1]])/tol)
        }
      }
      message('computing theta moment')
      output = list()
      output[[1]] =  pref_moment[[1]]
      output[[2]] =  pref_moment[[2]]
      print(paste0('output here = ',output[[1]]))
      iteration <<- iteration + 1; 
      print(paste0('at iteration = ', iteration))
      output[[2]][length(x_pref_theta)] = output[[2]][length(x_pref_theta)] * exp(x_pref_theta[length(x_pref_theta)])/(1 + exp(x_pref_theta[length(x_pref_theta)]))^2

      output_theta = aggregate_moment_theta(x_transform)
      output[[1]] = output[[1]] + output_theta[[1]] 
      output[[2]] = output[[2]] + output_theta[[2]][active_index] 
      return(output)
    }, control=list(maxit=1e1), method='BFGS')

    param_trial_here[active_index] = optim_pref_theta$par
  }
  
  active_index = c(x_transform[[2]]$sigma_r,  x_transform[[2]]$sigma_theta, x_transform[[2]]$beta_r); 

  x_transform = transform_param(param_trial_here, return_index=TRUE)
  
  min_theta_R = do.call('c', lapply(sample_r_theta, function(output_hh_index) min(cbind(var_ind(data_hh_list[[output_hh_index]]), data_hh_list[[output_hh_index]]$Year == 2004, data_hh_list[[output_hh_index]]$Year == 2006, data_hh_list[[output_hh_index]]$Year == 2010, data_hh_list[[output_hh_index]]$Year == 2012) %*% x_transform[[1]]$beta_theta)))

  max_theta_R = do.call('c', lapply(sample_r_theta, function(output_hh_index) max(cbind(var_ind(data_hh_list[[output_hh_index]]), data_hh_list[[output_hh_index]]$Year == 2004, data_hh_list[[output_hh_index]]$Year == 2006, data_hh_list[[output_hh_index]]$Year == 2010, data_hh_list[[output_hh_index]]$Year == 2012) %*% x_transform[[1]]$beta_theta)))

  mean_theta_R = do.call('c', lapply(sample_r_theta, function(output_hh_index) mean(cbind(var_ind(data_hh_list[[output_hh_index]]), data_hh_list[[output_hh_index]]$Year == 2004, data_hh_list[[output_hh_index]]$Year == 2006, data_hh_list[[output_hh_index]]$Year == 2010, data_hh_list[[output_hh_index]]$Year == 2012) %*% x_transform[[1]]$beta_theta)))

  message('start estimation of r')

  if (Sys.info()[['sysname']] == 'Windows') {
    clusterExport(cl, c('x_transform', 'sick_parameters', 'xi_parameters'),envir=environment())
    moment_eligible_hh_output = parLapply(cl, sample_r_theta, function(mini_data_index) {
      output = (household_draw_theta_kappa_Rdraw(mini_data_index, x_transform[[1]], n_halton_at_r, 10, sick_parameters, xi_parameters, u_lowerbar = -1))
      return(output)
    })
  } else {
    moment_eligible_hh_output = mclapply(sample_r_theta, function(mini_data_index) tryCatch(household_draw_theta_kappa_Rdraw(mini_data_index, x_transform[[1]], n_halton_at_r, 10, sick_parameters, xi_parameters, u_lowerbar = -1), error=function(e) e), mc.cores=numcores)
  }

  root_r = do.call('c',lapply(moment_eligible_hh_output, function(output_hh) output_hh$root_r))
  hh_theta = do.call('c',lapply(moment_eligible_hh_output, function(output_hh) output_hh$hh_theta))
  fx_r = function(x_r, derivative=FALSE, silent=TRUE) {
    if (max(abs(x_r)) > 10) {
      return(NA)
    }
    
    sd = exp(x_r[length(x_r)-1])
    correlation = x_r[length(x_r)]
    mean_vec = rep(X_hh_theta_r %*% x_r[1:(length(x_r)-2)], each = n_halton_at_r) + correlation * hh_theta
    prob = pnorm((root_r -mean_vec)/sd)
    prob[which(root_r == 4)] = 1
    prob[which(root_r < 0)] = 0
    prob[which(root_r >= 0)] = (prob[which(root_r >= 0)] -pnorm(- mean_vec[which(root_r >= 0)]/sd)) /(1 - pnorm(- mean_vec[which(root_r >= 0)]/sd))
    prob[which(is.na(prob))] = 1;

    # output = - sum(log(matrix(prob, nrow = n_halton_at_r) %>% colMeans + 1e-20) * full_insurance_indicator + log(matrix(1 - prob, nrow = n_halton_at_r) %>% colMeans + 1e-20) * (1 - full_insurance_indicator))
    output = sum((full_insurance_indicator - (matrix((1 - prob),nrow=n_halton_at_r) %>% colMeans))^2 * mat_Y_rtheta^2) + 
      sum((full_insurance_indicator - (matrix((1 - prob),nrow=n_halton_at_r) %>% colMeans))^2 * mean_theta_R^2) +
      sum((full_insurance_indicator - (matrix((1 - prob),nrow=n_halton_at_r) %>% colMeans))^2 * min_theta_R^2) + 
      sum((full_insurance_indicator - (matrix((1 - prob),nrow=n_halton_at_r) %>% colMeans))^2 * max_theta_R^2) + 
      sum((full_insurance_indicator - (matrix((1 - prob),nrow=n_halton_at_r) %>% colMeans))^2 * n_involuntary^2)

    if (!silent) {
      print('actual insurance = '); print(summary(full_insurance_indicator))
      print('predicted insurance '); print(summary(matrix((1 - prob),nrow=n_halton_at_r) %>% colMeans))
      print('x_r = '); print(x_r);
      print('output = '); print(output);
    }

    if (is.nan(output) | is.infinite(output)) {
      return(NA)
    }
    else {
      return(output)
    }
  }

  init_val = fx_r(rep(0,9)); 
  if (is.nan(init_val) | is.na(init_val)) {
    return(NA)
  }
  optim_r = optim(rep(0,9), function(x) fx_r(x, silent=TRUE), control=list(maxit=1000), method='Nelder-Mead') 

  fx_r(optim_r$par, silent=FALSE)

  param_trial_here[c(x_transform[[2]]$beta_r, x_transform[[2]]$sigma_r, x_transform[[2]]$correlation)] = optim_r$par
  x_transform = transform_param(param_trial_here,return_index=TRUE); 

  output_2 = do.call('c', lapply(moment_eligible_hh_output, function(output_hh) {
      sd = exp(optim_r$par[length(optim_r$par)-1])
      correlation = optim_r$par[length(optim_r$par)]
      mean_vec = rep(X_hh_theta_r %*% optim_r$par[1:(length(optim_r$par)-2)], each = n_halton_at_r) + correlation * hh_theta
      prob = pnorm((root_r -mean_vec)/sd)
      prob[which(root_r == 4)] = 1
      prob[which(root_r < 0)] = 0
      prob[which(root_r >= 0)] = (prob[which(root_r >= 0)] -pnorm(- mean_vec[which(root_r >= 0)]/sd)) /(1 - pnorm(- mean_vec[which(root_r >= 0)]/sd))
      prob[which(is.na(prob))] = 1;
      if (sum(prob) == 0) {
        Em = NA
      } else {
        Em = colSums(apply(output_hh$m, 2, function(x) x * prob/sum(prob)))
      }
      return(Em)
    }))

  print(output_2 %>% summary)
  output_2 =  mean(((output_2 - mat_M_rtheta[,1]) * full_insurance_indicator_ind_level)^2, na.rm=TRUE) * 1e9
  print(paste0('output_2 = ',output_2))
  
  print(mat_M_rtheta[which(full_insurance_indicator_ind_level == 1),1] %>% summary)
  print(paste0('optim_r')); print(optim_r$par)
  print('------')
  if (!(estimate_theta)) {
    optim_pref_theta = list()
    optim_pref_theta$value = 0 
  }
  if (return_result) {
    return(param_trial_here)
  } else {
    return(output_2 + optim_r$value + optim_pref_theta$value) 
  }
  
}

param_trial = compute_inner_loop(-3, return_result=TRUE)

initial_param_trial = param_trial; 

estimate_r_thetabar = optimize(function(xs) {
  output = try(compute_inner_loop(xs))
  if ('try-error' %in% class(output)) {
    return(NA)
  } else {
    return(output)
  }
}, c(-6,1)) 


param_trial = compute_inner_loop(estimate_r_thetabar$minimum, return_result=TRUE, estimate_theta=TRUE, estimate_pref = TRUE)

message('computing final param_trial')

param_final = list(); 
param_final$other = param_trial; 
param_final$xi = xi_parameters;
param_final$sick = sick_parameters

if (dir.exists('../../householdbundling_estimate')) {
  saveRDS(param_final, file=paste0('../../householdbundling_estimate/estimate_',job_index,'.rds'))
} else {
  dir.create('../../householdbundling_estimate') 
  saveRDS(param_final, file=paste0('../../householdbundling_estimate/estimate_',job_index,'.rds'))
}

