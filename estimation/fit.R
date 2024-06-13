args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) { 
  numcores = 24;
  job_index = 1;  
} else {
  job_index = as.numeric(args[1]);
  numcores = as.numeric(args[2]); 
}

library(knitr)
library(tidyverse)
library(lfe)
library(sandwich)
library(plm)
library(stargazer)
library(parallel)
library(randomForest)
library(randtoolbox)
library(Hmisc)

# setwd('./familyenrollment')
devtools::install(upgrade='never')
library(familyenrollment)

param_final <- readRDS(paste0('../../householdbundling_estimate/estimate_',job_index,'.rds'))
transform_param_final = transform_param(param_final$other)

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

if (Sys.info()[['sysname']] == 'Windows') {
  numcores = 10; 
  cl = makeCluster(numcores);
  clusterEvalQ(cl, library('tidyverse'))
  clusterEvalQ(cl, library('familyenrollment'))
}

if (Sys.info()[['sysname']] == 'Windows') {
  clusterExport(cl, c('transform_param_final', 'param_final','counterfactual_household_draw_theta_kappa_Rdraw'))
  fit_values = parLapply(cl, Vol_HH_list_index, function(id) {
	output = counterfactual_household_draw_theta_kappa_Rdraw(id, transform_param_final, 100, 10, param_final$sick, param_final$xi, u_lowerbar = -1, policy_mat_hh = policy_mat[[id]], seed_number = 1, constraint_function = function(x) x)
	output = as.data.frame(output)
	output$Y = data_hh_list[[id]]$Income; 
	output$m_observed = data_hh_list[[id]]$M_expense; 
	return(output)
	})
} else {
  fit_values = mclapply(Vol_HH_list_index, function(id) {
	output = counterfactual_household_draw_theta_kappa_Rdraw(id, transform_param_final, 100, 10, param_final$sick, param_final$xi, u_lowerbar = -1, policy_mat_hh = policy_mat[[id]], seed_number = 1, constraint_function = function(x) x)
	output = as.data.frame(output)
	output$Y = data_hh_list[[id]]$Income; 
	output$m_observed = data_hh_list[[id]]$M_expense; 
	return(output)}, mc.cores=numcores)
}

fit_values = do.call('rbind', fit_values)

observed_data_voluntary = do.call('rbind', data_hh_list[Vol_HH_list_index])

fit_values = as.data.frame(fit_values)
fit_values$Y2 <- as.numeric(Hmisc::cut2(fit_values$Y, g=5))

observed_data_voluntary = as.data.frame(observed_data_voluntary)
observed_data_voluntary$Y2 <- as.numeric(Hmisc::cut2(observed_data_voluntary$Income, g=5))

predicted_data_summary = fit_values %>% group_by(Y2) %>% summarise(mean_Vol_sts = mean(vol_sts_counterfactual), mean_m = mean(m))
predicted_data_summary$type = 'predicted'
actual_data_summary = observed_data_voluntary %>% group_by(Y2) %>% summarise(mean_Vol_sts = mean(Vol_sts), mean_m = mean(M_expense))
actual_data_summary$type = 'actual'

plot_1 = ggplot(data = rbind(predicted_data_summary, actual_data_summary), aes(x = Y2, y = mean_Vol_sts, color=type)) + geom_line() 
plot_2 = ggplot(data = rbind(predicted_data_summary, actual_data_summary), aes(x = Y2, y = mean_m, color=type)) + geom_line() 
gridExtra::grid.arrange(plot_1, plot_2, nrow=1)