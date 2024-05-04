args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) { 
  numcores = 2;
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
sample_r_theta = sample(sample_r_theta, length(sample_r_theta), replace=TRUE)
# sample_r_theta = sample(Vol_HH_list_index[which(!(is.na(lapply(Vol_HH_list_index, function(x) ifelse(nrow(data_hh_list[[x]]) <= 4, x, NA)))))], 5000)
sample_identify_pref = sample(sample_identify_pref, length(sample_identify_pref), replace=TRUE)
sample_identify_theta = sample(sample_identify_theta, length(sample_identify_theta), replace=TRUE)


message('merging household frames')
data = do.call('rbind', data_hh_list[sample_index])

# Estimate the probability of getting sick 
message('computing the sick parameters')
sick_parameters = optim(rep(0, ncol(var_ind(data_hh_list[[1]]))), fn = function(x) llh_sick(x, data), gr = function(x) grad_llh_sick(x, data), control=list(maxit = 1e4), method='BFGS')

message('computing the coverage parameters')
xi_parameters = optim(rep(0, 2 * ncol(var_ind(data_hh_list[[1]]))), fn = function(x) llh_xi(x, data), gr = function(x) grad_llh_xi(x, data), control=list(maxit = 1e4), method='BFGS')

param_trial = transform_param(return_index=TRUE, init=TRUE);
transform_param_trial = transform_param(param_trial, return_index=TRUE)

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


compute_moment_theta = function(x_transform) { 
    var_list = c('sigma_thetabar', 'beta_theta', 'beta_theta_ind')
    active_index = setdiff(unlist(transform_param_trial[[2]][var_list]),tail(transform_param_trial[[2]][['beta_theta']], n=4))
    if (Sys.info()[['sysname']] == 'Windows') {
      clusterExport(cl, 'x_transform',envir=environment())
      theta_estimate_output = parLapply(cl, data_hh_list_theta,function(mini_data) {
        output = tryCatch(identify_theta(data_set = mini_data, param=x_transform[[1]]),error=function(e) e)
        return(output)
      })
    } else {
      theta_estimate_output = mclapply(data_hh_list_theta, function(mini_data) tryCatch(identify_theta(data_set = mini_data, param=x_transform[[1]]), error=function(e) e), mc.cores=numcores)
    }

    result=list();

    result[[1]] = lapply(theta_estimate_output, function(x) x[[1]]) %>% unlist; 
    if (any(!(is.numeric(result[[1]])))) {
      return(list(NA, rep(NA, length(x))))
    }

    na_index = which(is.na(result[[1]])) 
    if (length(na_index) > 10) {
      result[[1]] = NA;
      result[[2]] = rep(NA, length(x)); 
      return(result[1:2])
    } else {
      if (length(na_index) > 0) {
        result[[1]] = sum(result[[1]][-na_index], na.rm=TRUE)
      } else {
        result[[1]] = sum(result[[1]], na.rm=TRUE)
      }
      
    }
    derivative= list()
    result[[2]] = rep(0, length(param_trial)); 
    
    for (varname in var_list) {
      if (length(na_index) == 0) {
        derivative[[varname]] = Reduce('+', lapply(theta_estimate_output, function(x) x[[2]][[varname]]))
      } else {
        derivative[[varname]] = Reduce('+', lapply(theta_estimate_output[-na_index], function(x) x[[2]][[varname]]))
      }
      
      result[[2]][c(x_transform[[2]][[varname]])] = derivative[[varname]]
    }
  
    return(result[1:2])
  } 


n_draw_halton = 100; 
if (Sys.info()[['sysname']] == 'Windows') {
  clusterExport(cl, c('active_index', 'param_trial', 'transform_param_trial', 'sample_identify_pref', 'xi_parameters', 'sick_parameters'))
  clusterExport(cl,'n_draw_halton')
}

if (Sys.info()[['sysname']] == 'Windows') {
  data_hh_list_pref = parLapply(cl, sample_identify_pref,function(index) {
    output = tryCatch(household_draw_theta_kappa_Rdraw(hh_index=index, param=transform_param_trial[[1]], n_draw_halton = n_draw_halton, n_draw_gauss = 10, sick_parameters, xi_parameters, short=FALSE),error=function(e) e)
    return(output)
  })
} else {
  data_hh_list_pref = mclapply(sample_identify_pref, function(index) tryCatch(household_draw_theta_kappa_Rdraw(hh_index=index, param=transform_param_trial[[1]], n_draw_halton = n_draw_halton, n_draw_gauss = 10, sick_parameters, xi_parameters, short=FALSE), error=function(e) e), mc.cores=1)
}

mat_M = do.call('rbind', lapply(data_hh_list_pref, function(x) {
      return(cbind(x$data$M_expense, x$data$M_expense^2))
    }))

mat_Year = do.call('rbind', lapply(data_hh_list_pref, function(x) {
      return(cbind(x$data$Year == 2004, x$data$Year == 2006, x$data$Year == 2010, x$data$Year == 2012))
    }))
outer_output = 0; 

var_list = c('beta_delta', 'beta_omega', 'beta_gamma', 'sigma_delta', 'sigma_gamma', 'sigma_omega')

X_ind_pref = do.call('rbind', lapply(data_hh_list_pref, function(x) x$X_ind))
X_hh_pref = do.call('rbind', lapply(data_hh_list_pref, function(x) x$X_hh))
X_ind_pref_with_year = cbind(X_ind_pref, mat_Year)

max_sigma_thetabar = sd(data$M_expense) * 2; 

x_transform = transform_param(param_trial, return_index=TRUE)

aggregate_moment_pref = function(x_transform) {
  if (Sys.info()[['sysname']] == 'Windows') {
    clusterExport(cl, 'x_transform',envir=environment())
    data_hh_list_pref = parLapply(cl, sample_identify_pref,function(index) {
      output = tryCatch(household_draw_theta_kappa_Rdraw(hh_index=index, param=x_transform[[1]], n_draw_halton = n_draw_halton, n_draw_gauss = 10, sick_parameters, xi_parameters, short=FALSE),error=function(e) e)
      return(output)
    })
    mat_YK = do.call('cbind', parLapply(cl, data_hh_list_pref, function(x) {
          output = rbind(colMeans(matrix(x$kappa_draw[[1]], nrow=n_draw_halton)), colMeans(matrix(x$kappa_draw[[1]]^2, nrow=n_draw_halton)), x$income[1], x$income[1]^2, colMeans(matrix(x$kappa_draw[[1]], nrow=n_draw_halton))*x$income[1])
          return(output)
        }))
  } else {
    data_hh_list_pref = mclapply(sample_identify_pref, function(index) tryCatch(household_draw_theta_kappa_Rdraw(hh_index=index, param=x_transform[[1]], n_draw_halton = n_draw_halton, n_draw_gauss = 10, sick_parameters, xi_parameters, short=FALSE), error=function(e) e), mc.cores=1)
    mat_YK = do.call('cbind', mclapply(data_hh_list_pref, function(x) {
          output = rbind(colMeans(matrix(x$kappa_draw[[1]], nrow=n_draw_halton)), colMeans(matrix(x$kappa_draw[[1]]^2, nrow=n_draw_halton)), x$income[1], x$income[1]^2, colMeans(matrix(x$kappa_draw[[1]], nrow=n_draw_halton))*x$income[1])
          return(output)
        }, mc.cores=1))
  }

  if (Sys.info()[['sysname']] == 'Windows') {
    clusterExport(cl, 'x_transform',envir=environment())
    moment_ineligible_hh_output = parLapply(cl, data_hh_list_pref,function(mini_data) {
      output = tryCatch(moment_ineligible_hh(mini_data, x_transform[[1]]),error=function(e) e)
      return(output)
    })
  } else {
    moment_ineligible_hh_output = mclapply(data_hh_list_pref, function(mini_data) tryCatch(moment_ineligible_hh(mini_data, x_transform[[1]]), error=function(e) e), mc.cores=1)
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
active_index = x_transform[[2]][c(var_list, 'beta_theta', 'beta_theta_ind', 'sigma_thetabar')] %>% unlist()
tol = 1e-4

iteration = 1; 
optim_pref_theta = splitfngr::optim_share(param_trial[active_index], function(x) {
  if (max(abs(x)) > 10) {
    return(list(NA, rep(NA, length(active_index))))
  }
  print('x = '); print(x)
  param_trial[active_index] <<- x 
  param_trial[x_transform[[2]]$sigma_thetabar] = log(exp(x[length(x)])/(exp(x[length(x)]) + 1) * max_sigma_thetabar)
  message('computing preference moment')
  pref_moment = aggregate_moment_pref(transform_param(param_trial, return_index=TRUE))
  deriv_theta_index = c(x_transform[[2]]$beta_theta[1], x_transform[[2]]$beta_theta_ind[1], x_transform[[2]]$sigma_thetabar)
  if (is.na(pref_moment[[1]])) {
    return(list(NA, rep(NA, length(x))))
  }
  message('computing preference moment derivative')
  for (i in deriv_theta_index) {
    param_trial_i = param_trial; param_trial_i[i] = param_trial[i] + tol
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
  theta_moment = aggregate_moment_pref(transform_param(param_trial, return_index=TRUE))
  output = list()
  output[[1]] = theta_moment[[1]] + pref_moment[[1]]*1e4
  output[[2]] = theta_moment[[2]] + pref_moment[[2]]*1e4
  print(paste0('output here = ',output[[1]]))
  iteration <<- iteration + 1; 
  print(paste0('at iteration = ', iteration))
  output[[2]][length(x)] = output[[2]][length(x)] * exp(x[length(x)])/(1 + exp(x[length(x)]))^2
  return(output)
}, control=list(maxit=1e2), method='BFGS')

x_transform = transform_param(param_trial, return_index=TRUE)

active_index = c(x_transform[[2]]$sigma_r,  x_transform[[2]]$sigma_theta, x_transform[[2]]$beta_r); 

param_trial[active_index] = 0; 

mat_M = do.call('rbind', lapply(data_hh_list[sample_r_theta], function(x) {
    return(cbind(x$M_expense, x$M_expense^2))
  }))

mat_Y = do.call('c', lapply(data_hh_list[sample_r_theta], function(x) x$Income))

full_insurance_indicator = do.call('c', lapply(data_hh_list[sample_r_theta], function(x) {
    return(x$HHsize_s[1] == x$N_vol[1])
  }))

full_insurance_indicator_ind_level = do.call('c', lapply(data_hh_list[sample_r_theta], function(x) {
    return(rep(x$HHsize_s[1] == x$N_vol[1], x$HHsize[1]))
  }))

X_r_all = do.call('rbind', lapply(data_hh_list[sample_r_theta], function(x) var_hh(x)))

save_obj_outside = NULL;
save_param_outside = NULL; 


# some draws will be insensitive to values of sigma_thetabar, so we do not need to simulate over these draws multiple times. 
x_new = param_trial; 
x_new[x_transform[[2]]$sigma_theta] = 0;
x_transform = transform_param(x_new, return_index=TRUE)
n_halton_at_r = 100; 
if (Sys.info()[['sysname']] == 'Windows') {
  x_new = param_trial; 
  x_new[x_transform[[2]]$sigma_theta] = 0;
  x_transform = transform_param(x_new, return_index=TRUE)
  clusterExport(cl, c('x_transform', 'sick_parameters', 'xi_parameters', 'n_halton_at_r'))
  moment_eligible_hh_output_upper = parLapply(cl, sample_r_theta, function(mini_data_index) {
    output = (household_draw_theta_kappa_Rdraw(mini_data_index, x_transform[[1]], n_halton_at_r, 10, sick_parameters, xi_parameters, u_lowerbar = -1))
    return(output)
  })
} else {
  moment_eligible_hh_output_upper = mclapply(sample_r_theta, function(mini_data_index) tryCatch(household_draw_theta_kappa_Rdraw(mini_data_index, x_transform[[1]], n_halton_at_r, 10, sick_parameters, xi_parameters, u_lowerbar = -1), error=function(e) e), mc.cores=numcores)
}

x_new = param_trial; 
x_new[x_transform[[2]]$sigma_theta] = 0;
x_transform = transform_param(x_new, return_index=TRUE)
if (Sys.info()[['sysname']] == 'Windows') {
  clusterExport(cl, c('x_transform', 'sick_parameters', 'xi_parameters', 'n_halton_at_r'))
  moment_eligible_hh_output_lower = parLapply(cl, sample_r_theta, function(mini_data_index) {
    output = (household_draw_theta_kappa_Rdraw(mini_data_index, x_transform[[1]], n_halton_at_r, 10, sick_parameters, xi_parameters, u_lowerbar = -1))
    return(output)
  })
} else {
  moment_eligible_hh_output_lower = mclapply(sample_r_theta, function(mini_data_index) tryCatch(household_draw_theta_kappa_Rdraw(mini_data_index, x_transform[[1]], n_halton_at_r, 10, sick_parameters, xi_parameters, u_lowerbar = -1), error=function(e) e), mc.cores=numcores)
}


relevant_index = lapply(1:length(sample_r_theta), function(index) {
  output = list();
  output[[1]] = sample_r_theta[index]
  output[[2]] = which(!(is.infinite(moment_eligible_hh_output_lower[[index]]$root_r)) | !(is.infinite(moment_eligible_hh_output_upper[[index]]$root_r)))
  return(output)    
})

relevant_index_na = lapply(1:length(relevant_index), function(id) ifelse(length(relevant_index[[id]][[2]]) > 0, id, NA)) %>% unlist()
relevant_index_nona = relevant_index[which(!(is.na(relevant_index_na)))]
full_insurance_indicator_nona = full_insurance_indicator[which(!(is.na(relevant_index_na)))]

full_insurance_indicator_ind_level_nona = do.call('c', lapply(data_hh_list[sample_r_theta[which(!(is.na(relevant_index_na)))]], function(x) {
    return(rep(x$HHsize_s[1] == x$N_vol[1], x$HHsize[1]))
  }))

mat_M = do.call('rbind', lapply(data_hh_list[sample_r_theta[which(!(is.na(relevant_index_na)))]], function(x) {
    return(cbind(x$M_expense, x$M_expense^2))
  }))

X_hh_all = do.call('rbind',lapply(moment_eligible_hh_output_lower, function(output_hh) output_hh$X_hh))

vec_Y_nona = do.call('c', lapply(data_hh_list[sample_r_theta[which(!(is.na(relevant_index_na)))]], function(x) {
    return(x$Income[1])
  }))


estimate_r_thetabar = optimize(function(x) {
      # optim_rf_trial = optim(trial_paraam[-zero_index], function(x) {
        x_new = param_trial; 
        x_new[x_transform[[2]]$sigma_theta] = x; 
        x_transform = transform_param(x_new, return_index=TRUE)
        print('----------------');
        print('x = '); print(x_transform[[1]])


        if (Sys.info()[['sysname']] == 'Windows') {
          clusterExport(cl, c('x_transform', 'sick_parameters', 'xi_parameters'),envir=environment())
          moment_eligible_hh_output = parLapply(cl, relevant_index_nona, function(mini_data_index) {
            output = (household_draw_theta_kappa_Rdraw(mini_data_index[[1]], x_transform[[1]], n_halton_at_r, 10, sick_parameters, xi_parameters, u_lowerbar = -1, relevant_index = mini_data_index[[2]]))
            return(output)
          })
        } else {
          moment_eligible_hh_output = mclapply(relevant_index_nona, function(mini_data_index) tryCatch(household_draw_theta_kappa_Rdraw(mini_data_index, x_transform[[1]], n_halton_at_r, 10, sick_parameters, xi_parameters, u_lowerbar = -1), error=function(e) e), mc.cores=numcores)
        }

        save_param_outside <<- x_transform[[1]]

        save_obj_outside <<- moment_eligible_hh_output


        root_r = do.call('c',lapply(moment_eligible_hh_output, function(output_hh) output_hh$root_r))

        root_r_none_inf = do.call('c', lapply(relevant_index, function(index) c(1:n_halton_at_r) %in% index[[2]]))
        fx_r = function(x_r, derivative=FALSE) {
          mean_vec = rep(X_hh_all %*% x_r[-length(x_r)], each = n_halton_at_r)
          sd = exp(x_r[length(x_r)])
          prob = pnorm(root_r, mean = mean_vec, sd = sd)
          prob_0 = pnorm(0, mean = mean_vec, sd = sd)
          if (derivative) {
            dprob = dnorm((root_r-mean_vec)/sd);
            dprob_0 = dnorm((0-mean_vec)/sd);
            d_prob_mean = apply(X_hh_all,2, function(x) x * dprob*(-1)) ;
            d_prob_sd = dprob *(-1)/sd^2; 
            d_prob_0_mean = apply(X_hh_all,2, function(x) x * dprob_0)*(-1) ;
            d_prob_0_sd = dprob_0*(-1)/sd^2; 
            d_output_mean = do.call('c', lapply(1:ncol(d_prob_mean), function(id_col) mean(matrix(((-d_prob_mean[,id_col])*(1 - prob_0) + d_prob_0_mean[,id_col]*(1 - prob))/(1 - prob_0)^2, nrow=n_halton_at_r) %>% colMeans * (full_insurance_indicator_nona - (matrix((1 - prob)/(1 - prob_0),nrow=n_halton_at_r) %>% colMeans)))))
            d_output_sd = mean(matrix(((-d_prob_sd)*(1 - prob_0) + d_prob_0_sd*(1 - prob))/(1 - prob_0)^2, nrow=n_halton_at_r) %>% colMeans * (full_insurance_indicator_nona - (matrix((1 - prob)/(1 - prob_0),nrow=n_halton_at_r) %>% colMeans)))
          }
          
          # output = -mean(full_insurance_indicator_nona * log(matrix((1 - prob)/(1 - prob_0),nrow=n_halton_at_r) %>% colMeans + 1e-5) + (1 - full_insurance_indicator_nona) * log((1 - (matrix((1 - prob)/(1 - prob_0),nrow=n_halton_at_r) %>% colMeans)) + 1e-5))
          # output = sum((full_insurance_indicator_nona - (matrix((1 - prob)/(1 - prob_0),nrow=n_halton_at_r)) %>% colMeans))^2 + sum((full_insurance_indicator_nona - (matrix((1 - prob)/(1 - prob_0),nrow=n_halton_at_r)) %>% colMeans)^2 * vec_Y_nona^2)
          output = sum((full_insurance_indicator_nona - (matrix((1 - prob)/(1 - prob_0),nrow=n_halton_at_r)) %>% colMeans)^2 * vec_Y_nona^2)
          print('actual insurance = '); print(summary(full_insurance_indicator_nona))
          print('predicted insurance '); print(summary(matrix((1 - prob)/(1 - prob_0),nrow=n_halton_at_r)))
          print('x_r = '); print(x_r);
          print('output = '); print(output);
          if (derivative) {
            return(list(output, c(- d_output_mean, - d_output_sd * sd)))
          } else{
            if (is.nan(output) | is.infinite(output)) {
              return(NA)
            }
            else {
              return(output)
            }
          }
        }

        # optim_r = splitfngr::optim_share(rep(0, length(param_trial[c(x_transform[[2]]$beta_r, x_transform[[2]]$sigma_r)])), function(x) fx_r(x, derivative=TRUE), control=list(maxit=1000), method='BFGS') 

        optim_r = optim(c(rep(0,7),-2), fx_r, control=list(maxit=1000), method='Nelder-Mead') 

        param_trial[c(x_transform[[2]]$beta_r, x_transform[[2]]$sigma_r)] <<- optim_r$par
        param_trial[x_transform[[2]]$sigma_theta] <<- x; 
        x_transform = transform_param(param_trial,return_index=TRUE); 

        output_2 = do.call('c', lapply(moment_eligible_hh_output, function(output_hh) {
            prob_full_insured = (1 - pnorm(output_hh$root_r, mean = output_hh$X_hh %*% x_transform[[1]]$beta_r, sd = exp(x_transform[[1]]$sigma_r)))/(1 - pnorm(0, mean = output_hh$X_hh %*% x_transform[[1]]$beta_r, sd = exp(x_transform[[1]]$sigma_r)))
            Em = colMeans(apply(output_hh$m, 2, function(x) x * prob_full_insured))
            return(Em)
          }))

        print(output_2 %>% summary)
        output_2 = mean(((output_2 - mat_M[,1] * full_insurance_indicator_ind_level_nona))^2, na.rm=TRUE) * 1e4
        print(paste0('output_2 = ',output_2))
        
        print((mat_M[,1] * full_insurance_indicator_ind_level_nona)%>% summary)
        print(paste0('optim_r')); print(optim_r$par)
        print('------')
        return(output_2)
      }, c(0,exp(x_transform[[1]]$sigma_thetabar))) 

param_trial[x_transform[[2]]$sigma_theta] <- estimate_r_thetabar$minimum

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

