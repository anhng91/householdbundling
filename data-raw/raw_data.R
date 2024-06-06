## code to prepare `raw_data` dataset goes here
library(tidyverse)

data = read.csv('data-raw/raw_data.csv')

unit_inc = data %>% group_by(hhid,Year) %>% slice(1) %>% select(Income) %>% unlist() %>% mean(); 
upper_bound = data %>% group_by(hhid,Year) %>% slice(1) %>% select(Income) %>% unlist() %>% quantile(0.975); 
lower_bound = data %>% group_by(hhid,Year) %>% slice(1) %>% select(Income) %>% unlist() %>% quantile(0.025);

data = data %>% filter(Income < upper_bound & Income > lower_bound)

age_hh_mean = data %>% group_by(hhid, Year) %>% slice(1) %>% select(hhage) %>% unlist() %>% mean();
age_hh_sd = data %>% group_by(hhid, Year) %>% slice(1) %>% select(hhage) %>% unlist() %>% sd();
age_hh_max_mean = data %>% group_by(hhid, Year) %>% slice(1) %>% select(hhmaxage) %>% unlist() %>% mean();
age_hh_max_sd = data %>% group_by(hhid, Year) %>% slice(1) %>% select(hhmaxage) %>% unlist() %>% sd(); 
data = data %>% mutate(hhage = (hhage - age_hh_mean)/age_hh_sd) %>% mutate(hhmaxage = (hhmaxage - age_hh_max_mean)/age_hh_max_sd)



data$lnhhavginc = data$lnhhavginc - log(unit_inc); 
data$hhavgeduc = (data$hhavgeduc -  mean(data$hhavgedu))/sd(data$hhavgedu); 

# Subtract subsistence income level
data$Income[which(data$Year == 2004)] = (data$Income[which(data$Year == 2004)] - 1897.272 * data$HHsize[which(data$Year == 2004)]) / unit_inc; 
data$Income[which(data$Year == 2006)] = (data$Income[which(data$Year == 2006)] - 2328.836 * data$HHsize[which(data$Year == 2006)]) / unit_inc; 
data$Income[which(data$Year == 2008)] = (data$Income[which(data$Year == 2008)] - 3698.429 * data$HHsize[which(data$Year == 2008)]) / unit_inc; 
data$Income[which(data$Year == 2010)] = (data$Income[which(data$Year == 2010)] - 6067.262 * data$HHsize[which(data$Year == 2010)]) / unit_inc; 
data$Income[which(data$Year == 2012)] = (data$Income[which(data$Year == 2012)] - 9284.135 * data$HHsize[which(data$Year == 2012)]) / unit_inc; 

relationship = factor(data$relationship);
HHtype = factor(data$HHtype); 

relationship_dummy = model.matrix(~relationship);

HHtype_dummy = model.matrix(~HHtype);

relationship_dummy = relationship_dummy[,2:dim(relationship_dummy)[2]]

HHtype_dummy = HHtype_dummy[,2:dim(HHtype_dummy)[2]]

geographic_dummy = cbind(data$geo2, data$geo3, data$geo4);

age_dummy = cbind(data$age2, data$age3, data$age4, data$age5); 

data$indshare = data$indincome/unit_inc / data$Income;


time_dummy = cbind(data$year3, data$year4, data$year5, data$year6);

hh_char = cbind(data$hhmaxage, data$hhfemale, data$hhage, data$HHsize, data$hhavgeduc, HHtype_dummy);
ind_char = cbind(data$college, data$married, data$female, data$employed, data$student);

data$student = 0;
data$student[which(data$study == 1 | data$studying == 1)] = 1; 


# Generate coinsurance 
	data$Baseprice = data$Baseprice/unit_inc; 
	data$Baseprice_s = data$Baseprice_s/unit_inc; 
	data$M_expense = data$M_expense /unit_inc; 
	Premium_fun <- function( N_vol,  year,  baseprice,  hhtype,  hhsize) { 
	    if (N_vol == 0) {
	    	output = 0; 
	    }
	    else {
	        if (year == 2004) {
	          if (N_vol == 1) {
	            output = baseprice; 
	          }
	          else {
	            output = baseprice + baseprice * 0.95 * (N_vol - 1); 
	          }
	        }
	        else if ((year <= 2008) && (year > 2004)) {
	          if (N_vol <= 2){
	            output = baseprice*N_vol ;
	          }
	          else if (N_vol == 3){
	            output = baseprice*2 + baseprice * 0.9 ; 
	          }
	          else if (N_vol >= 4){
	            output = baseprice*2 + baseprice * 0.9 + baseprice * (N_vol - 3) * 0.8;
	          }
	        }
	        else {
	          # // Note that only if your family is in the agriculture sector or near-poor that you have
	          # // to have everyone in the household paying for health insurance in order to receive discount
	          # // If on the other hand, the household has eligible-for-compulsory members in the household, 
	          # // then no need to have everyone committing to having health insurance to receive discount. 
	          # // Self employed household does not receive any discount (Kinda odd)
	          if (hhtype == 1) { # Agricultural sector/near poor 
	            if (N_vol == hhsize) {
	              if (N_vol == 1) {
	                output = baseprice;
	              }
	              else if (N_vol == 2) {
	                output = baseprice  + baseprice * 0.9;
	              }
	              else if (N_vol == 3) {
	                output = baseprice + baseprice * 0.9 + baseprice * 0.8; 
	              }
	              else {
	                output = baseprice + baseprice * 0.9 + baseprice * 0.8 + baseprice * (N_vol - 3) * 0.7;
	              }
	            }
	            else {
	              output = baseprice * N_vol;
	            }
	          }
	          else if (hhtype == 2) { # having compulsory members in the household 
	              if (N_vol == 1) {
	                output = baseprice;
	              }
	              else if (N_vol == 2) {
	                output = baseprice  + baseprice * 0.9;
	              }
	              else if (N_vol == 3) {
	                output = baseprice + baseprice * 0.9 + baseprice * 0.8; 
	              }
	              else {
	                output = baseprice + baseprice * 0.9 + baseprice * 0.8 + baseprice * (N_vol - 3) * 0.7;
	              }
	          }
	          else if (hhtype == 3) { # Self employed
	            output = baseprice * N_vol; 
	          }
	        }
	    }
  		return(output);  
  	}

  	data$HHsize_s = data$HHsize - data$N_bef - data$N_com - data$N_std_w_ins; 
	Premium_fun_hh <- function(data_hh_i,...) { 
		output = lapply(0:data_hh_i$HHsize_s[1], function(n_i) data_hh_i$Income[1] - Premium_fun(n_i, data_hh_i$Year[1], data_hh_i$Baseprice[1], 
			data_hh_i$HHtype[1], data_hh_i$HHsize_s[1]) - data_hh_i$Baseprice_s[1] * data_hh_i$N_std_w_ins[1]);
		return(unlist(output));
	}      	

	Income_net_premium = lapply(data %>% group_by(hhid, Year) %>% group_split(), Premium_fun_hh) 
      	
  	data$tot_cost_normalized = data$tot_cost/unit_inc; #this is m

  	# p is effective price now 
  	data = data %>% mutate(tot_cost_normalized = ifelse(zeta_observed < 0, M_expense, tot_cost_normalized)) 

  	data = data %>% mutate(effective_price = M_expense/tot_cost * unit_inc) %>% mutate(effective_price = ifelse(sick_dummy == 1 & M_expense == 0, 0, effective_price)) %>% mutate(effective_price = ifelse(zeta_observed < 0, -1, effective_price)) %>% mutate(effective_price = ifelse(Bef_sts + Vol_sts + Com_sts + Std_w_ins == 0, 1, effective_price)) %>% mutate(effective_price = ifelse(is.nan(effective_price), -1, effective_price))


# Generate truncated spending 
	data$upper_index = data$M_expense == 0 & data$sick_dummy == 1
	
	Upper_bound <- function(Bef_sts, Vol_sts, Com_sts, Std_w_ins, Year, upper_index) {
		output = -1; 
		if (upper_index == 0) {

		}
		else {
			if ((Bef_sts == 1) && (Year <= 2008)) {
				output = (0.5); 
			}
			else if ((Bef_sts == 1) && (Year > 2008)) {
				output = (100/unit_inc); 
			}
			else if ((Com_sts == 1) || (Vol_sts ==1) || (Std_w_ins == 1)) {
				if (Year == 2004) {
					output = (20/unit_inc);
				}
				else if (Year == 2006) {
					output = (7000/unit_inc);
				}
				else if ((Year == 2008) && (Vol_sts == 1)) {
					output = (100/unit_inc);
				}
				else if ((Year == 2008) && (Com_sts == 1)) {
					output = (7000/unit_inc);
				}
				else {
					output = (100/unit_inc);
				}
			}
		}
		return(output);
	}

	data$upper_bound = mapply(Upper_bound, data$Bef_sts, data$Vol_sts, data$Com_sts, 
		data$Std_w_ins, data$Year, data$upper_index, SIMPLIFY = FALSE) %>% unlist(); 

# Generate price vector for insurance policy (3 X N rows, 3 cols)
	policy_fun <- function(Bef_sts, Com_sts, Year) {
		#initiate
		kappa_vec = rep(NA, 10)
		m_threshold_vec = rep(NA, 10)
		cost_threshold_vec = rep(NA, 10)

		if ((Bef_sts == 1) && (Year <= 2008)) {
			kappa_vec[1] = 0; m_threshold_vec[1] = 0; cost_threshold_vec[1]=0;
		}
		else if ((Bef_sts == 1) && (Year >= 2010)) {
			kappa_vec[1:2] = c(0.0, 0.05); 
			m_threshold_vec[1:2] = c(0, 100/unit_inc); 
			cost_threshold_vec[1:2] = c(0, 0); 
		}
		else if ((Bef_sts == 0) && (Year == 2004)) {
			kappa_vec[1:2] = c(0.2, 0.0); 
			m_threshold_vec[1:2] = c(0, 1500/unit_inc); 
			cost_threshold_vec[1:2] = c(0, 0.2*1500/unit_inc); 
		}
		else if ((Bef_sts == 1) && (Year >= 2010)) {
			kappa_vec[1:2] = c(0, 0.05); 
			m_threshold_vec[1:2] = c(0, 100/unit_inc); 
			cost_threshold_vec[1:2] = c(0, 0); 
		}
		else if ((Bef_sts == 0) && (Year == 2006)) {
			kappa_vec[1:3] = c(0, 1, 0.4); 
			m_threshold_vec[1:3] = c(0, 7000/unit_inc, 7000/0.6/unit_inc)
			cost_threshold_vec[1:3] = c(0, 0, 7000/unit_inc);
		}
		else if ((Com_sts == 1) && (Year == 2008)) {
			kappa_vec[1:3] = c(0, 1, 0.4); 
			m_threshold_vec[1:3] = c(0, 7000/unit_inc, 7000/0.6/unit_inc)
			cost_threshold_vec[1:3] = c(0, 0, 7000/unit_inc);
		}
		else if (((Com_sts + Bef_sts == 0) && (Year >= 2008)) || ((Com_sts == 1) && (Year > 2008))) {
			kappa_vec[1:2] = c(0, 0.2); 
			m_threshold_vec[1:2] = c(0, 100/unit_inc); 
			cost_threshold_vec[1:2] = c(0, 0);
		}
		return(list(kappa_vec, m_threshold_vec, cost_threshold_vec));
	}

	policy_fun_hh <- function(data_hh_i) {
		output = list()
		for (i in 1:nrow(data_hh_i)) {
  			output[[i]] = policy_fun(data_hh_i$Bef_sts[i], data_hh_i$Com_sts[i], data_hh_i$Year[i]);
		}
		kappa_mat = do.call('cbind',lapply(output, function(x) x[[1]]))
		m_threshold_mat = do.call('cbind',lapply(output, function(x) x[[2]]))
		cost_threshold_mat = do.call('cbind',lapply(output, function(x) x[[3]]))
		return(list(kappa_mat, m_threshold_mat, cost_threshold_mat))
	}
	policy_mat = lapply(data %>% group_by(hhid, Year) %>% group_split(), policy_fun_hh)

# Export household into small dataframes 
data_hh_list = data %>% group_by(hhid, Year) %>% group_split()

# Export data objects
usethis::use_data(policy_mat, overwrite = TRUE)
usethis::use_data(Income_net_premium, overwrite = TRUE)
usethis::use_data(data_hh_list, overwrite = TRUE)

