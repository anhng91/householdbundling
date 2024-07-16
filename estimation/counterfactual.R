args = commandArgs(trailingOnly=TRUE)
if (length(args)<2) { 
  numcores = 2;  
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

list_hh_2012 = unlist(lapply(1:length(data_hh_list), function(hh_index) ifelse(data_hh_list[[hh_index]]$Year[1] == 2012 & data_hh_list[[hh_index]]$HHsize_s[1] == 2, hh_index, NA)))
list_hh_2012 = list_hh_2012[which(!(is.na(list_hh_2012)))]


list_hh_2012 = sample(list_hh_2012, 500)

data_2012 = benchmark %>% filter(id %in% list_hh_2012);

budget_2012 = data_2012 %>% filter(Com_sts + Bef_sts + Std_w_ins == 0) %>% group_by(id) %>% group_modify(function(data_hh_i, ...) {
	return(data.frame(net_budget = Income_net_premium[[data_hh_i$id[1]]][1] - Income_net_premium[[data_hh_i$id[1]]][1 + sum(data_hh_i$vol_sts_counterfactual)] - sum(data_hh_i$cost_to_insurance)/unit_inc))
},.keep=TRUE) %>% ungroup() %>% pull(net_budget) %>% sum()


if (Sys.info()[['sysname']] == 'Windows') {
  numcores = 10; 
  cl = makeCluster(numcores);
  clusterEvalQ(cl, library('tidyverse'))
  clusterEvalQ(cl, library('familyenrollment'))
  clusterExport(cl,c('Vol_HH_list_index', 'Com_HH_list_index', 'out_sample_index'))
}

for (job_index in 4:4) {
	if (file.exists(paste0('../../householdbundling_estimate/estimate_',job_index,'.rds'))) {
		param_final <- readRDS(paste0('../../householdbundling_estimate/estimate_',job_index,'.rds'))
		param_final$other = init_param
		transform_param_final = transform_param(param_final$other)

		counterfactual_premium = function(premium_param, type, id) {
			MW = (Income_net_premium[[id]][1] - Income_net_premium[[id]][2]) * 4/3 / 0.06 # 2012 premium is at 4.5\% MW for the first member
			if (type == 'bundle discount') {
				hhsize_s = length(Income_net_premium[[id]]) - 1; 
				income_vec = Income_net_premium[[id]][1] - c(0, premium_param[[hhsize_s]][1:hhsize_s]) * MW;
			} else if (type == 'pure bundling') {
				income_vec = rep(Income_net_premium[[id]][1], hhsize_s); 
				income_vec[hhsize_s] = Income_net_premium[[id]][1] - premium_param[[hhsize_s]] * MW; 
			} else {
				income_vec =  Income_net_premium[[id]][1] - premium_param[[hhsize_s]] * c(0:hhsize_s) * MW;
			}
			return(income_vec)
		}

		bd_prem = list()
		pb_prem = list()
		prem_id = 0; 
		for (i1 in seq(0, 0.06, by = 0.005)) {
			for (i2 in seq(0, 1, by = 0.05)) {
				prem_id = prem_id + 1; 
				bd_prem[[prem_id]] = list(); 
				bd_prem[[prem_id]][[2]] = c(i1, i1 * (1 + i2))
				pb_prem[[prem_id]] = list(); 
				pb_prem[[prem_id]][[2]] = c(0, i1 * (1 + i2))
				if (i1 == 0) {
					break
				}
			}
		} 

		bd_output = list(); 
		bd_output$budget = NULL; 
		bd_output$surplus = NULL

		pb_output = list(); 
		pb_output$budget = NULL; 
		pb_output$surplus = NULL

		print_index = 0; 
		for (prem in bd_prem) {
			print_index = print_index + 1; print(paste0('computing bundle discount premium index ', print_index));
			counterfactual_values = mclapply(c(list_hh_2012), function(id) {
				income_vec = counterfactual_premium(prem, 'bundle discount', id)
				output = counterfactual_household_draw_theta_kappa_Rdraw(id, transform_param_final, 100, 10, param_final$sick, param_final$xi, u_lowerbar = -1, policy_mat_hh = policy_mat[[id]], seed_number = 1, constraint_function = function(x) x, income_vec = income_vec)
				output = as.data.frame(output)
				output$Y = data_hh_list[[id]]$Income; 
				output$m_observed = data_hh_list[[id]]$M_expense; 
				output$fit_type = ifelse(id %in% Vol_HH_list_index, 2, ifelse(id %in% Com_HH_list_index, 1, 3))
				output$id = id
				output$premium_optimal = income_vec[1] - income_vec[sum(output$vol_sts_counterfactual) + 1]
				return(output)}, mc.cores=numcores)
			counterfactual_values = do.call('rbind', counterfactual_values) %>% filter(Com_sts + Bef_sts + Std_w_ins == 0)
			bd_output$surplus = c(bd_output$surplus, counterfactual_values %>% group_by(id) %>% slice(1) %>% pull(wtp) %>% sum - sum(counterfactual_values %>% group_by(id) %>% slice(1) %>% pull(premium_optimal))) 
			bd_output$budget = c(bd_output$budget, sum(counterfactual_values %>% group_by(id) %>% slice(1) %>% pull(premium_optimal)) - (counterfactual_values$cost_to_insurance %>% sum()))
		}

		print_index = 0; 
		for (prem in pb_prem) {
			print_index = print_index + 1; print(paste0('computing pure bundling premium index ', print_index));
			counterfactual_values = mclapply(c(list_hh_2012), function(id) {
				income_vec = counterfactual_premium(prem, 'bundle discount', id)
				output = counterfactual_household_draw_theta_kappa_Rdraw(id, transform_param_final, 100, 10, param_final$sick, param_final$xi, u_lowerbar = -1, policy_mat_hh = policy_mat[[id]], seed_number = 1, constraint_function = function(x) {x_new = x; x_new[-c(1, length(x))] = -Inf; return(x_new)}, income_vec = income_vec)
				output = as.data.frame(output)
				output$Y = data_hh_list[[id]]$Income; 
				output$m_observed = data_hh_list[[id]]$M_expense; 
				output$fit_type = ifelse(id %in% Vol_HH_list_index, 2, ifelse(id %in% Com_HH_list_index, 1, 3))
				output$id = id
				output$premium_optimal = income_vec[1] - income_vec[sum(output$vol_sts_counterfactual) + 1]
				return(output)}, mc.cores=numcores)
			counterfactual_values = do.call('rbind', counterfactual_values) %>% filter(Com_sts + Bef_sts + Std_w_ins == 0)
			pb_output$surplus = c(pb_output$surplus, counterfactual_values %>% group_by(id) %>% slice(1) %>% ungroup() %>% pull(wtp) %>% sum - sum(counterfactual_values %>% group_by(id) %>% slice(1) %>% ungroup() %>% pull(premium_optimal)))
			pb_output$budget = c(pb_output$budget, sum(counterfactual_values %>% group_by(id) %>% slice(1) %>% ungroup() %>% pull(premium_optimal)) - (counterfactual_values$cost_to_insurance %>% sum()))
		}	

		bd_output = as.data.frame(bd_output)
		pb_output = as.data.frame(pb_output)
		pb_premium_vector = as.data.frame(do.call('rbind',lapply(pb_prem, function(x) x[[2]]))); names(pb_premium_vector) = c('p1','p2')
		bd_premium_vector = as.data.frame(do.call('rbind',lapply(bd_prem, function(x) x[[2]]))); names(bd_premium_vector) = c('p1','p2')
		bd_output = cbind(bd_output, bd_premium_vector); names(bd_output)[c(3,4)] = c('p1','p2')
		bd_output = bd_output %>% mutate(p_ratio = p2/p1); 
		bd_output$p_ratio[is.nan(bd_output$p_ratio)] = 1; 
		pb_output = cbind(pb_output, pb_premium_vector)
		plot_1 = ggplot(data = bd_output , aes(x = p1, y = p_ratio, fill = surplus, color = budget >= budget_2012)) + geom_tile(linewidth=0.3) + geom_text(aes(label = round(surplus,3)), color = 'white',size=2) + theme_bw() + theme(legend.position = "none") + xlab(expression(p[1])) + ylab(expression(p[2]/p[1])) + scale_y_continuous(breaks=seq(1,2,length.out=6))

		graph_data = rbind(pb_output, bd_output %>% select(-p_ratio) %>% filter(p2 == 2 * p1), bd_output %>% select(-p_ratio) %>% filter(p2 == 1.9 * p1))
		graph_data$type = c(rep('PB', nrow(pb_output)), rep('IP', nrow(bd_output %>% filter(p2 == 2 * p1))), rep('BD', nrow(bd_output %>% filter(p2 == 1.9 * p1))))
		ggplot(data = graph_data, aes(x = p2, y = surplus, linetype= budget >= budget_2012, color=type)) + geom_line()


		stop()

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



