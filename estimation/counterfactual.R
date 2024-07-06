args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) { 
  numcores = 20;  
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

benchmark = readRDS('../../Obj_for_manuscript/fit_values.rds')

list_hh_2012 = unlist(lapply(1:length(data_hh_list), function(hh_index) ifelse(data_hh_list[[hh_index]]$Year[1] == 2012, hh_index, NA)))
list_hh_2012 = list_hh_2012[which(!(is.na(list_hh_2012)))]

data_2012 = benchmark %>% filter(id %in% list_hh_2012);

budget_2012 = data_2012 %>% group_by(id) %>% group_modify(function(data_hh_i, ...) {
	return(data.frame(net_budget = Income_net_premium[[data_hh_i$id[1]]][1 + sum(data_hh_i$vol_sts_counterfactual)] - sum(data_hh_i$cost_to_insurance)))
},.keep=TRUE) %>% ungroup() %>% pull(net_budget) %>% sum()

budget_2012 = budget_2012/unit_inc

if (Sys.info()[['sysname']] == 'Windows') {
  numcores = 10; 
  cl = makeCluster(numcores);
  clusterEvalQ(cl, library('tidyverse'))
  clusterEvalQ(cl, library('familyenrollment'))
  clusterExport(cl,c('Vol_HH_list_index', 'Com_HH_list_index', 'out_sample_index'))
}

for (job_index in 0:1) {
	if (file.exists(paste0('../../householdbundling_estimate/estimate_',job_index,'.rds'))) {
		param_final <- readRDS(paste0('../../householdbundling_estimate/estimate_',job_index,'.rds'))
		param_final$other = init_param
		transform_param_final = transform_param(param_final$other)

		counterfactual_premium = function(premium_param, type, id) {
			baseprice = Income_net_premium[[id]][1] - Income_net_premium[[id]][2]
			if (type == 'bundle discount') {
				hhsize_s = length(Income_net_premium[[id]]) - 1; 
				income_vec = Income_net_premium[[id]][1] - c(0, premium_param[[hhsize_s]][1:hhsize_s]) * baseprice;
			} else if (type == 'pure bundling') {
				income_vec = rep(Income_net_premium[[id]][1], hhsize_s); 
				income_vec[hhsize_s] = Income_net_premium[[id]][1] - premium_param[[hhsize_s]] * baseprice; 
			} else {
				income_vec =  Income_net_premium[[id]][1] - premium_param[[hhsize_s]] * c(0:hhsize_s) * baseprice;
			}
		}

		for (i in )


		if (Sys.info()[['sysname']] == 'Windows') {
		  clusterExport(cl, c('transform_param_final', 'param_final','counterfactual_household_draw_theta_kappa_Rdraw'))
		  fit_values = parLapply(cl, c(Vol_HH_list_index, Com_HH_list_index), function(id) {
			output = tryCatch(counterfactual_household_draw_theta_kappa_Rdraw(id, transform_param_final, 100, 10, param_final$sick, param_final$xi, u_lowerbar = -1, policy_mat_hh = policy_mat[[id]], seed_number = 1, constraint_function = function(x) x), error=function(x) x)
			output = as.data.frame(output)
			output$Y = data_hh_list[[id]]$Income; 
			output$m_observed = data_hh_list[[id]]$M_expense; 
			output$fit_type = ifelse(id %in% Vol_HH_list_index, 2, ifelse(id %in% Com_HH_list_index, 1, 3))
			output$id = id; 
			return(output)
			})

			no_heterogeneity_values = parLapply(cl, c(Vol_HH_list_index), function(id) {
			output = tryCatch(counterfactual_household_draw_theta_kappa_Rdraw(id, transform_param_final, 100, 10, param_final$sick, param_final$xi, u_lowerbar = -1, policy_mat_hh = policy_mat[[id]], seed_number = 1, constraint_function = function(x) x, within_hh_heterogeneity = list(omega=FALSE, gamma=FALSE, delta=FALSE, theta_bar=FALSE)), error=function(x) x)
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

			no_heterogeneity_values = mclapply(c(Vol_HH_list_index), function(id) {
			output = tryCatch(counterfactual_household_draw_theta_kappa_Rdraw(id, transform_param_final, 100, 10, param_final$sick, param_final$xi, u_lowerbar = -1, policy_mat_hh = policy_mat[[id]], seed_number = 1, constraint_function = function(x) x, within_hh_heterogeneity = list(omega=FALSE, gamma=FALSE, delta=FALSE, theta_bar=FALSE)), error=function(x) x)
			output = as.data.frame(output)
			output$Y = data_hh_list[[id]]$Income; 
			output$m_observed = data_hh_list[[id]]$M_expense; 
			output$fit_type = ifelse(id %in% Vol_HH_list_index, 2, ifelse(id %in% Com_HH_list_index, 1, 3))
			output$id = id; 
			return(output)
			}, mc.cores = numcores)
		}
	}
}



