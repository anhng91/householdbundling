#' Sample data used to test moment of eligible households
#'
#' A subset of households' constructed data after estimation of the preference parameters and the unconditional health shocks
#'
#' @format ## `sample_data_and_parameter`
#' A list with two elements
#' \describe{
#'   \item{data}{A list with different vectors characterizing a household}
#'   \item{param}{A list with different vectors characterizing parameters}
#' }
#' @source sample_data_and_parameter
"sample_data_and_parameter"


#' Cleaned data on household enrollment and medical expenditure
#'
#' @format ## `data_hh_list`
#' A list, each list element is the data frame for a household. 
#' \describe{
#'   \item{M_expense}{medical expense out of pocket}
#' 	 \item{tot_cost}{Total medical expenditure (not equal to M_expense only if proportion of coverage is recorded)}
#' 	 \item{Vol_sts}{Whether a member is voluntarily insured}
#' 	 \item{Com_sts}{Whether a member is compulsorily insured}
#'   \item{Bef_sts}{Whether a member has free insurance}
#' }
#' @source raw_data
"data_hh_list"


#' Constructed coinsurance matrices 
#'
#' @format ## `policy_mat`
#' A list, with three elements, each is a matrix of HHsize columns and 10 rows . 
#' \describe{
#'   \item{List 1}{The coinsurance rate of each linear segment}
#' 	 \item{List 2}{The threshold at which coinsurance rates change}
#' 	 \item{List 3}{The out of pocket costs at the nearest threshold}
#' }
#' @source raw_data
"policy_mat"


#' Constructed premium vector for each household 
#' 
#' @format ## `Income_net_premium`
#' A vector
#' \describe{
#'   \item{Element 1}{Income if no member is insured}
#' 	 \item{Element 2}{Income if 1 member is insured }
#' 	 \item{Element 3}{Income if 2 members are insured}
#' }
#' @source raw_data
"Income_net_premium"


#' Example of estimates for the parameters characterizing whether a household member is going to be sick
#' 
#' @format ## `sick_parameters_sample`
#' A list produced by optim
#' 
#' @source raw_data
"sick_parameters_sample"


#' Example of estimates for the parameters characterizing whether a household member is going to be covered under SHI
#' 
#' @format ## `xi_parameters_sample`
#' A list produced by optim
#' 
#' @source raw_data
"xi_parameters_sample"


#' unit of income in estimation
#' 
#' @format ## `unit_inc`
#' a numerical value
#' 
#' @source raw_data
"unit_inc"



#' starting values (estimated values from a subsample)
#' 
#' @format ## `init_param`
#' a vector 
#' 
#' @source estimation.R 

