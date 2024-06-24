args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) { 
  numcores = 4;  
} else {
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

Vol_HH_list_index = lapply(1:length(data_hh_list), function(hh_index) {
  data = data_hh_list[[hh_index]]; 
  if (nrow(data) > nrow(data %>% filter(Bef_sts + Com_sts + Std_w_ins == 1)) & !(data$Year[1] == 2006 & data$HHsize_s[1] > 1)) {
    return(hh_index);
  }
  else {
    return(NA); 
  }
}) %>% unlist(); 

Vol_HH_list_index = Vol_HH_list_index[!(is.na(Vol_HH_list_index))] 

Com_HH_list_index = lapply(1:length(data_hh_list), function(hh_index) {
  data = data_hh_list[[hh_index]]; 
  if (0 == nrow(data %>% filter(Bef_sts + Com_sts + Std_w_ins == 0))) {
    return(hh_index);
  }
  else {
    return(NA); 
  }
}) %>% unlist(); 

out_sample_index = lapply(1:length(data_hh_list), function(hh_index) {
  data = data_hh_list[[hh_index]]; 
  if (nrow(data) > nrow(data %>% filter(Bef_sts + Com_sts + Std_w_ins == 1)) & (data$Year[1] == 2006 & data$HHsize_s[1] > 1)) {
    return(hh_index);
  }
  else {
    return(NA); 
  }
}) %>% unlist(); 

out_sample_index = out_sample_index[!(is.na(out_sample_index))] 

Com_HH_list_index = Com_HH_list_index[!(is.na(Com_HH_list_index))]

if (Sys.info()[['sysname']] == 'Windows') {
  numcores = 10; 
  cl = makeCluster(numcores);
  clusterEvalQ(cl, library('tidyverse'))
  clusterEvalQ(cl, library('familyenrollment'))
}

for (job_index in 0:1) {
	if (file.exists(paste0('../../householdbundling_estimate/estimate_',job_index,'.rds'))) {
		param_final <- readRDS(paste0('../../householdbundling_estimate/estimate_',job_index,'.rds'))
		transform_param_final = transform_param(param_final$other)
		if (Sys.info()[['sysname']] == 'Windows') {
		  clusterExport(cl, c('transform_param_final', 'param_final','counterfactual_household_draw_theta_kappa_Rdraw'))
		  fit_values = parLapply(cl, c(Vol_HH_list_index, Com_HH_list_index), function(id) {
			output = counterfactual_household_draw_theta_kappa_Rdraw(id, transform_param_final, 100, 10, param_final$sick, param_final$xi, u_lowerbar = -1, policy_mat_hh = policy_mat[[id]], seed_number = 1, constraint_function = function(x) x)
			output = as.data.frame(output)
			output$Y = data_hh_list[[id]]$Income; 
			output$m_observed = data_hh_list[[id]]$M_expense; 
			output$fit_type = ifelse(id %in% Vol_HH_list_index, 2, ifelse(id %in% Com_HH_list_index, 1, 3))
			output$id = id; 
			return(output)
			})
		} else {
		  fit_values = mclapply(c(Vol_HH_list_index, Com_HH_list_index), function(id) {
			output = counterfactual_household_draw_theta_kappa_Rdraw(id, transform_param_final, 100, 10, param_final$sick, param_final$xi, u_lowerbar = -1, policy_mat_hh = policy_mat[[id]], seed_number = 1, constraint_function = function(x) x)
			output = as.data.frame(output)
			output$Y = data_hh_list[[id]]$Income; 
			output$m_observed = data_hh_list[[id]]$M_expense; 
			output$fit_type = ifelse(id %in% Vol_HH_list_index, 2, ifelse(id %in% Com_HH_list_index, 1, 3))
			output$id = id
			return(output)}, mc.cores=numcores)
		}
		fit_values = do.call('rbind', fit_values)
		fit_values = as.data.frame(fit_values)
		fit_values$job_index = job_index; 

		if (Sys.info()[['sysname']] == 'Windows') {
		  clusterExport(cl, c('transform_param_final', 'param_final','counterfactual_household_draw_theta_kappa_Rdraw'))
		  out_sample_values = parLapply(cl, out_sample_index, function(id) {
			output = counterfactual_household_draw_theta_kappa_Rdraw(id, transform_param_final, 100, 10, param_final$sick, param_final$xi, u_lowerbar = -1, policy_mat_hh = policy_mat[[id]], seed_number = 1, constraint_function = function(x) x[c(1,length(x))])
			output = as.data.frame(output)
			output$Y = data_hh_list[[id]]$Income; 
			output$m_observed = data_hh_list[[id]]$M_expense; 
			output$fit_type = ifelse(id %in% Vol_HH_list_index, 2, ifelse(id %in% Com_HH_list_index, 1, 3))
			output$id = id
			return(output)
			})
		} else {
		  out_sample_values = mclapply(out_sample_index, function(id) {
			output = counterfactual_household_draw_theta_kappa_Rdraw(id, transform_param_final, 100, 10, param_final$sick, param_final$xi, u_lowerbar = -1, policy_mat_hh = policy_mat[[id]], seed_number = 1, constraint_function = function(x) x[c(1,length(x))])
			output = as.data.frame(output)
			output$Y = data_hh_list[[id]]$Income; 
			output$m_observed = data_hh_list[[id]]$M_expense; 
			output$fit_type = ifelse(id %in% Vol_HH_list_index, 2, ifelse(id %in% Com_HH_list_index, 1, 3))
			output$id = id 
			return(output)}, mc.cores=numcores)
		}
		out_sample_values = do.call('rbind', out_sample_values)
		out_sample_values = as.data.frame(out_sample_values)
		out_sample_values$job_index = job_index; 

		fit_values = rbind(fit_values, out_sample_values); 
	}
}
fit_values$Y2 <- as.numeric(Hmisc::cut2(fit_values$Y, g=5))
observed_vol_data = do.call('rbind', data_hh_list[Vol_HH_list_index]); observed_vol_data$fit_type = 2; 
observed_com_data = do.call('rbind', data_hh_list[Com_HH_list_index]); observed_com_data$fit_type = 1; 
observed_out_sample_data = do.call('rbind', data_hh_list[out_sample_index]); observed_out_sample_data$fit_type = 3; 
observed_data = rbind(observed_vol_data, observed_com_data, observed_out_sample_data)
observed_data = as.data.frame(observed_data)
observed_data$Y2 <- as.numeric(Hmisc::cut2(observed_data$Income, g=5))

fit_values = fit_values %>% mutate_at(c('m_observed', 'average_theta', 'wtp', 'cost_to_insurance', 'Y', 'm', 'optional_care'), function(x) x * unit_inc)
observed_data = observed_data %>% mutate_at(c('M_expense', 'Income'), function(x) x * unit_inc)

saveRDS(fit_values, file='../../Obj_for_manuscript/fit_values.rds')
saveRDS(observed_data, file='../../Obj_for_manuscript/observed_data.rds')


all_param_final = list()
job_index_normalized = 1;
for (job_index in 0:1) {
	if (file.exists(paste0('../../householdbundling_estimate/estimate_',job_index,'.rds'))) {
		param_final <- readRDS(paste0('../../householdbundling_estimate/estimate_',job_index,'.rds'))
		all_param_final[[job_index_normalized]] = transform_param(param_final$other)
		all_param_final[[job_index_normalized]]$beta_sick = param_final$sick$par
		job_index_normalized = job_index_normalized + 1;
	}
}

name_summary_param = names(all_param_final[[1]]); 
all_param_final_2 = lapply(name_summary_param, function(name_x) do.call('rbind', lapply(all_param_final, function(x) x[[name_x]])))
names(all_param_final_2) = name_summary_param
all_param_final_2$beta_theta = as.data.frame(all_param_final_2$beta_theta); names(all_param_final_2$beta_theta) = c(names(var_ind(data_hh_list[[1]]) %>% as.data.frame), 'year_dummy_2004', 'year_dummy_2006', 'year_dummy_2010', 'year_dummy_2012')
all_param_final_2$beta_theta_ind = as.data.frame(all_param_final_2$beta_theta_ind); names(all_param_final_2$beta_theta_ind) = c(names(var_ind(data_hh_list[[1]]) %>% as.data.frame))
all_param_final_2$beta_gamma = as.data.frame(all_param_final_2$beta_gamma); names(all_param_final_2$beta_gamma) = c(names(var_ind(data_hh_list[[1]]) %>% as.data.frame))
all_param_final_2$beta_sick = as.data.frame(all_param_final_2$beta_sick); names(all_param_final_2$beta_sick) = c(names(var_ind(data_hh_list[[1]]) %>% as.data.frame))
all_param_final_2$beta_delta = as.data.frame(all_param_final_2$beta_delta); names(all_param_final_2$beta_delta) = c(names(var_ind(data_hh_list[[1]]) %>% as.data.frame))
all_param_final_2$beta_r = as.data.frame(all_param_final_2$beta_r); names(all_param_final_2$beta_r) = c(names(var_hh(data_hh_list[[1]]) %>% as.data.frame))
all_param_final_2$beta_omega = as.data.frame(all_param_final_2$beta_omega); names(all_param_final_2$beta_omega) = c(names(var_hh(data_hh_list[[1]]) %>% as.data.frame))

saveRDS(all_param_final_2, file='../../Obj_for_manuscript/all_param_final_2.rds')

