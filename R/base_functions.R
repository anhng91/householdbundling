#' Compute the expectation of X^y, where y is a truncated normal random variable
#'
#' @param X this is a constant
#' @param mu The mean of the truncated normal random variable 
#' @param sigma The standard deviation of the truncated normal random variable (truncated at 0) 
#'
#' @return A numerical value.
#' @export
#'
#' @examples
#' mean_E_XW(1, 1, 1)

mean_E_XW = function(mu, sigma, X) {
	output = exp(1/2 * (sigma^2 * (log(X))^2 + 2 * mu * log(X))) * (1 - pnorm((-mu - sigma^2 * log(X))/sigma))/(1 - pnorm(-mu/sigma));
	output[which(X == 0)] = 0 
	return(output)
}

#' Compute the derivative of the expectation of X^y, where y is a truncated normal random variable
#'
#' @param X this is a constant
#' @param mu The mean of the truncated normal random variable 
#' @param sigma The standard deviation of the truncated normal random variable (truncated at 0) 
#'
#' @return A list with two values, derivative for the mean and the derivative for the standard deviation
#' @export
#'
#' @examples
#' d_mean_E_XW(1, 1, 1)

d_mean_E_XW = function(mu, sigma, X) {
	d_output = list();

	d_output$mu = (exp((sigma^2*log(X)^2)/2 + mu*log(X))*log(X)*(pracma::erfc((2^(1/2)*(log(X)*sigma^2 + mu))/(2*sigma))/2 - 1))/(pracma::erfc((2^(1/2)*mu)/(2*sigma))/2 - 1) - (2^(1/2)*exp((sigma^2*log(X)^2)/2 + mu*log(X))*exp(-(log(X)*sigma^2 + mu)^2/(2*sigma^2)))/(2*sigma*pi^(1/2)*(pracma::erfc((2^(1/2)*mu)/(2*sigma))/2 - 1)) + (2^(1/2)*exp(-mu^2/(2*sigma^2))*exp((sigma^2*log(X)^2)/2 + mu*log(X))*(pracma::erfc((2^(1/2)*(log(X)*sigma^2 + mu))/(2*sigma))/2 - 1))/(2*sigma*pi^(1/2)*(pracma::erfc((2^(1/2)*mu)/(2*sigma))/2 - 1)^2)

	d_output$sigma = (sigma*exp((sigma^2*log(X)^2)/2 + mu*log(X))*log(X)^2*(pracma::erfc((2^(1/2)*(log(X)*sigma^2 + mu))/(2*sigma))/2 - 1))/(pracma::erfc((2^(1/2)*mu)/(2*sigma))/2 - 1) - (exp((sigma^2*log(X)^2)/2 + mu*log(X))*exp(-(log(X)*sigma^2 + mu)^2/(2*sigma^2))*(2^(1/2)*log(X) - (2^(1/2)*(log(X)*sigma^2 + mu))/(2*sigma^2)))/(pi^(1/2)*(pracma::erfc((2^(1/2)*mu)/(2*sigma))/2 - 1)) - (2^(1/2)*mu*exp(-mu^2/(2*sigma^2))*exp((sigma^2*log(X)^2)/2 + mu*log(X))*(pracma::erfc((2^(1/2)*(log(X)*sigma^2 + mu))/(2*sigma))/2 - 1))/(2*sigma^2*pi^(1/2)*(pracma::erfc((2^(1/2)*mu)/(2*sigma))/2 - 1)^2)

	d_output$X = (exp((sigma^2*log(X)^2)/2 + mu*log(X))*(pracma::erfc((2^(1/2)*(log(X)*sigma^2 + mu))/(2*sigma))/2 - 1)*((log(X)*sigma^2)/X + mu/X))/(pracma::erfc((2^(1/2)*mu)/(2*sigma))/2 - 1) - (2^(1/2)*sigma*exp((sigma^2*log(X)^2)/2 + mu*log(X))*exp(-(log(X)*sigma^2 + mu)^2/(2*sigma^2)))/(2*X*pi^(1/2)*(pracma::erfc((2^(1/2)*mu)/(2*sigma))/2 - 1))

	d_output$mu[which(X == 0)] = 0; 
	d_output$sigma[which(X == 0)] = 0; 
	d_output$X[which(X == 0)] = 0; 

	return(d_output)
}


#' Compute the mean of truncated normal random variable
#'
#' @param mu The mean of the truncated normal random variable 
#' @param sigma The standard deviation of the truncated normal random variable (truncated at 0) 
#'
#' @return mean of truncated normal distribution at 0
#' @export
#'
#' @examples
#' truncated_normal_mean(1, 1)

truncated_normal_mean = function(mu, sigma) {
	output = mu  + dnorm(0, mean=mu, sd=sigma)/(1 - pnorm(0, mean=mu, sd=sigma))*sigma
	if (is.infinite(output) | is.nan(output)) {
		return(NA)
	} else {
		return(output)
	}
}


#' Compute the derivative of the mean of truncated normal random variable with respect to its mean and variance 
#'
#' @param mu The mean of the truncated normal random variable 
#' @param sigma The standard deviation of the truncated normal random variable (truncated at 0) 
#'
#' @return a list of 2 derivatives, one w.r.t to the mean and one w.r.t to standard deviation
#' @export
#'
#' @examples
#' d_truncated_normal_mean(1, 1)

d_truncated_normal_mean = function(mu, sigma) {
	output = mu  + dnorm(0, mean=mu, sd=sigma)/(1 - pnorm(0, mean=mu, sd=sigma))*sigma
	d_output = list();

	# the derivatives are produced by matlab symbolic functions, the numerical values are related to pi.

	d_output$mu = (2251799813685248*mu*exp(-mu^2/(2*sigma^2)))/(5644425081792261*sigma^2*(pracma::erfc((2^(1/2)*mu)/(2*sigma))/2 - 1)) - (1125899906842624*2^(1/2)*exp(-mu^2/sigma^2))/(5644425081792261*sigma*pi^(1/2)*(pracma::erfc((2^(1/2)*mu)/(2*sigma))/2 - 1)^2) + 1

 	d_output$sigma = (1125899906842624*2^(1/2)*mu*exp(-mu^2/sigma^2))/(5644425081792261*sigma^2*pi^(1/2)*(pracma::erfc((2^(1/2)*mu)/(2*sigma))/2 - 1)^2) - (2251799813685248*mu^2*exp(-mu^2/(2*sigma^2)))/(5644425081792261*sigma^3*(pracma::erfc((2^(1/2)*mu)/(2*sigma))/2 - 1))
 
 
	if (is.infinite(output) | is.nan(output)) {
		return(list(NA, NA))
	} else {
		return(d_output)
	}
}

#' Compute the variance of truncated normal random variable
#'
#' @param mu The mean of the truncated normal random variable 
#' @param sigma The standard deviation of the truncated normal random variable (truncated at 0) 
#'
#' @return variance of truncated normal distribution at 0
#' @export
#'
#' @examples
#' truncated_normal_variance(1, 1)

truncated_normal_variance = function(mu, sigma) {
	Z = (1 - pnorm(0, mean=mu, sd=sigma))
	output = sigma^2 * (1 + (0 - mu)/sigma * dnorm(0, mean=mu, sd=sigma)/Z - (dnorm(0, mean=mu, sd=sigma))^2/Z^2)
	if (is.infinite(output) | is.nan(output)) {
		return(NA)
	} else {
		return(output)
	}
}

#' Compute the derivative of the variance of truncated normal random variable with respect to its mean and variance 
#'
#' @param mu The mean of the truncated normal random variable 
#' @param sigma The standard deviation of the truncated normal random variable (truncated at 0) 
#'
#' @return a list of 2 derivatives, one w.r.t to the mean and one w.r.t to standard deviation
#' @export
#'
#' @examples
#' d_truncated_normal_variance(1, 1)

d_truncated_normal_variance = function(mu, sigma) {
	Z = (1 - pnorm(0, mean=mu, sd=sigma))
	output = sigma^2 * (1 + (0 - mu)/sigma * dnorm(0, mean=mu, sd=sigma)/Z - (dnorm(0, mean=mu, sd=sigma))^2/Z^2)

	d_output = list()
	d_output$mu = sigma^2*((2251799813685248*exp(-mu^2/(2*sigma^2)))/(5644425081792261*sigma^2*(pracma::erfc((2^(1/2)*mu)/(2*sigma))/2 - 1)) + (10141204801825835211973625643008*mu*exp(-mu^2/sigma^2))/(31859534503965572279823959492121*sigma^4*(pracma::erfc((2^(1/2)*mu)/(2*sigma))/2 - 1)^2) - (2251799813685248*mu^2*exp(-mu^2/(2*sigma^2)))/(5644425081792261*sigma^4*(pracma::erfc((2^(1/2)*mu)/(2*sigma))/2 - 1)) - (5070602400912917605986812821504*2^(1/2)*exp(-mu^2/sigma^2)*exp(-mu^2/(2*sigma^2)))/(31859534503965572279823959492121*sigma^3*pi^(1/2)*(pracma::erfc((2^(1/2)*mu)/(2*sigma))/2 - 1)^3) + (1125899906842624*2^(1/2)*mu*exp(-mu^2/sigma^2))/(5644425081792261*sigma^3*pi^(1/2)*(pracma::erfc((2^(1/2)*mu)/(2*sigma))/2 - 1)^2));


	d_output$sigma = 2*sigma*((2251799813685248*mu*exp(-mu^2/(2*sigma^2)))/(5644425081792261*sigma^2*(pracma::erfc((2^(1/2)*mu)/(2*sigma))/2 - 1)) - (5070602400912917605986812821504*exp(-mu^2/sigma^2))/(31859534503965572279823959492121*sigma^2*(pracma::erfc((2^(1/2)*mu)/(2*sigma))/2 - 1)^2) + 1) + sigma^2*((10141204801825835211973625643008*exp(-mu^2/sigma^2))/(31859534503965572279823959492121*sigma^3*(pracma::erfc((2^(1/2)*mu)/(2*sigma))/2 - 1)^2) - (4503599627370496*mu*exp(-mu^2/(2*sigma^2)))/(5644425081792261*sigma^3*(pracma::erfc((2^(1/2)*mu)/(2*sigma))/2 - 1)) - (10141204801825835211973625643008*mu^2*exp(-mu^2/sigma^2))/(31859534503965572279823959492121*sigma^5*(pracma::erfc((2^(1/2)*mu)/(2*sigma))/2 - 1)^2) + (2251799813685248*mu^3*exp(-mu^2/(2*sigma^2)))/(5644425081792261*sigma^5*(pracma::erfc((2^(1/2)*mu)/(2*sigma))/2 - 1)) - (1125899906842624*2^(1/2)*mu^2*exp(-mu^2/sigma^2))/(5644425081792261*sigma^4*pi^(1/2)*(pracma::erfc((2^(1/2)*mu)/(2*sigma))/2 - 1)^2) + (5070602400912917605986812821504*2^(1/2)*mu*exp(-mu^2/sigma^2)*exp(-mu^2/(2*sigma^2)))/(31859534503965572279823959492121*sigma^4*pi^(1/2)*(pracma::erfc((2^(1/2)*mu)/(2*sigma))/2 - 1)^3))


	if (is.infinite(output) | is.nan(output)) {
		return(list(NA, NA))
	} else {
		return(d_output)
	}

}


#' Compute the expectation of a CRRA function with Gaussian random variable
#'
#' @param log_R A constant.
#' @param mu The mean of the truncated normal random variable 
#' @param sigma The standard deviation of the truncated normal random variable (truncated at 0) 
#' @param n_draw The number of points used to approximate the integral
#'
#' @return A numerical value.
#' @export
#'
#' @examples
#' log_R = 1; mu = 1; sigma = 1; n_draw = 50;
#' CRRA_fun_Gauss(log_R, mu, sigma, n_draw)

CRRA_fun_Gauss = function(log_R, mu, sigma, n_draw = 50) {
	Gauss_laguerre = pracma::gaussLaguerre(n_draw, a = 0); 
	Gauss_laguerre_integral = sum(Gauss_laguerre$w * (exp((1 - Gauss_laguerre$x) * log_R - ((Gauss_laguerre$x - mu)/sigma)^2/2 + Gauss_laguerre$x - log(1 - pnorm(-mu/sigma))) - exp(-((Gauss_laguerre$x - mu)/sigma)^2/2 + Gauss_laguerre$x - log(1 - pnorm(-mu/sigma))))/sigma/sqrt(2 * pi)/(1 - Gauss_laguerre$x));
	return(Gauss_laguerre_integral)
}

#' Compute the mean and variance of a truncated from below normal distribution
#'
#' @param mu mean of the (untruncated) normal variable 
#' @param sigma standard deviation of the (untruncated) normal variable
#' @param a is the truncation point. DEFAULT to 0.
#'
#' @return a list of two values, mean and variance.
#' @export
#'
#' @examples
#' trunc_norm(1,1,2)

trunc_norm = function(mu, sigma, a=0) {
	alpha = (a - mu)/sigma;
	Z = 1 - pnorm(alpha);
	output = list(); 
	output$mean = mu + exp(-1/2*alpha^2 - log(Z))/sqrt(2*pi)*sigma; 
	output$variance = sigma^2 * (1 + alpha * exp(-1/2 * alpha^2 - log(Z))/sqrt(2 * pi) - (exp(-1/2 * alpha^2 - log(Z))/sqrt(2*pi))^2)
	return(output)
}

#' Compute the expectation of a CARA expected utility function
#'
#' @param a A constant, which is the mean of (direct) utility.
#' @param b A constant, which is the variance of (indirect) utility. 
#' @param sigma The standard deviation of the truncated normal random variable (truncated at 0) 
#' @param mu The mean of the truncated normal random variable
#'
#' @return A numerical value.
#' @export
#'
#' @examples
#' second_taylor_CARA(1,1,1,1)
second_taylor_CARA = function(a, b, mu, sigma) {
	truncated_output = trunc_norm(mu = mu - a * sigma^2, sigma = sigma, a = 0);
	output = - exp(sigma^2 * a^2/2 - mu * a) * ((1 + b/2 * (truncated_output$mean^2 + truncated_output$variance))) * ((1 - pnorm(-(mu - a * sigma^2)/sigma))/((1 - pnorm(-mu/sigma))))
	return(output)
}


#' Compute the expectation of a CARA expected utility function using numerical simulations 
#'
#' @param a A constant, which is the mean of (direct) utility.
#' @param b A constant, which is the variance of (indirect) utility. 
#' @param mu The mean of the truncated normal random variable. 
#' @param sigma The standard deviation of the truncated normal random variable (truncated at 0) 
#'
#' @return A list. The first element is the numerical integration, the second element is the computation based on the second-order Taylor approximation.
#' @export
#'
#' @examples
#' second_taylor_CARA_test(1,1,1,1)
second_taylor_CARA_test = function(a, b, mu, sigma) {
	r_draw = qnorm(pnorm(-mu/sigma) +  (1 - pnorm(-mu/sigma)) * randtoolbox::halton(1000)) * sigma + mu; 
	# plot = ggplot2::ggplot(data = data.frame(x = r_draw), ggplot2::aes(x=x)) + ggplot2::geom_histogram()
	output = list()
	output$numerical = mean(-exp(-r_draw * a) * (1 + r_draw^2/2 * b))
	output$manual = second_taylor_CARA(a,b,mu,sigma); 
	# print(paste0('value 2 = ', mean(-exp(- r_draw * a))))
	return(list(output))
}


#' Compute the moments of ineligible households
#'
#' @param data_set This is a list that contains all pre-computed elements of each household .
#' @param param This is a list, each element is a coefficient (i.e., a vector) of a preference parameter (beta_gamma, beta_omega, beta_delta, beta_theta, beta_theta_ind) and the standard deviations of the unobserved heterogeneity (sigma_gamma, sigma_omega, sigma_delta, sigma_theta, sigma_thetabar)
#'
#' @return a list. The first element is the expected value of m|Y,k. The second element is the expected value of m^2|Y,k. The third and the fourth elements are their respective derivatives with respect to the preference parameters.
#' 
#' @export
#'
#' @examples
#' mini_data = household_draw_theta_kappa_Rdraw(3, sample_data_and_parameter$param, 1000, 10, sick_parameters_sample, xi_parameters_sample, short=FALSE);
#' moment_ineligible_hh(mini_data, sample_data_and_parameter$param)
#' 
#' 
moment_ineligible_hh = function(data_set, param) {
	tol = 1e-3;
	data_hh_i = data_set$data;
	policy_mat_hh_index = data_set$policy_mat;
	hh_index = data_set$index; 
	theta_draw = data_set$theta_draw;
	income_vec = data_set$income;
	X_ind = data_set$X_ind; 
	X_hh = data_set$X_hh;
	sick_p = data_set$sick_p;
	R_draw = data_set$R_draw[[1]];
	kappa_draw = data_set$kappa_draw[[1]]; 
	

	HHsize = data_hh_i$HHsize[1]

	mean_beta_omega = X_hh %*% param$beta_omega; 
	mean_beta_delta = X_ind %*% param$beta_delta;
	mean_beta_gamma = X_ind %*% param$beta_gamma; 


	m_draw = NULL; 
	m_draw2 = NULL; 

	d_m_draw = list()

	d_m_draw$mean_beta_delta = NULL; 
	d_m_draw$mean_beta_omega = NULL; 
	d_m_draw$mean_beta_gamma = NULL; 
	d_m_draw$sigma_delta = NULL; 
	d_m_draw$sigma_omega = NULL; 
	d_m_draw$sigma_gamma = NULL; 

	d_m_draw2 = list()

	d_m_draw2$mean_beta_delta = NULL; 
	d_m_draw2$mean_beta_omega = NULL; 
	d_m_draw2$mean_beta_gamma = NULL; 
	d_m_draw2$sigma_delta = NULL; 
	d_m_draw2$sigma_omega = NULL; 
	d_m_draw2$sigma_gamma = NULL; 

	# mat1 = c(Y, colMeans(kappa));
	# mat2 = c(Y, colMeans(kappa), Y * colMeans(kappa), Y^2, colMeans(kappa^2));

	for (mem_index in 1:HHsize) {
		obj1 =  mean_E_XW(mean_beta_gamma[mem_index], exp(param$sigma_gamma), 1/(1 + kappa_draw[,mem_index]))
		obj2 = truncated_normal_mean(mean_beta_delta[mem_index], exp(param$sigma_delta))
		obj3 = mean_E_XW(mean_beta_omega, exp(param$sigma_omega), R_draw)
		obj4 =  mean_E_XW(mean_beta_omega, exp(param$sigma_omega), R_draw);

		obj5 = d_truncated_normal_mean(mean_beta_delta[mem_index], exp(param$sigma_delta)); 
		obj6 = d_mean_E_XW(mean_beta_omega, exp(param$sigma_omega), R_draw);
		obj7 = d_mean_E_XW(mean_beta_gamma[mem_index], exp(param$sigma_gamma), 1/(1 + kappa_draw[,mem_index]))

		obj8 = truncated_normal_variance(mean_beta_delta[mem_index], exp(param$sigma_delta))
		obj9 = mean_E_XW(mean_beta_omega, exp(param$sigma_omega), R_draw^2)
		obj10 = d_truncated_normal_variance(mean_beta_delta[mem_index], exp(param$sigma_delta))
		obj11 =  mean_E_XW(mean_beta_gamma[mem_index], exp(param$sigma_gamma), 1/(1 + kappa_draw[,mem_index])^2)
		obj12 = d_mean_E_XW(mean_beta_gamma[mem_index], exp(param$sigma_gamma), 1/(1 + kappa_draw[,mem_index])^2)
		obj13 = d_mean_E_XW(mean_beta_omega, exp(param$sigma_omega), R_draw^2)

		m_draw[mem_index] = mean(theta_draw[,mem_index] * kappa_draw[,mem_index] + obj2 * theta_draw * obj4 * obj1 * kappa_draw[,mem_index])
		
		d_m_draw$mean_beta_delta[mem_index] = mean(obj5$mu * theta_draw * obj4 * obj1 * kappa_draw[,mem_index])

		d_m_draw$sigma_delta[mem_index] = mean(obj5$sigma * theta_draw * obj4 * obj1 * exp(param$sigma_delta) * kappa_draw[,mem_index])

		d_m_draw$sigma_omega[mem_index] = mean(obj2 * theta_draw * obj6$sigma * exp(param$sigma_omega) * obj1 * kappa_draw[,mem_index]);

		d_m_draw$mean_beta_omega[mem_index] = mean(obj2 * theta_draw * obj6$mu * obj1 * kappa_draw[,mem_index]);

		d_m_draw$sigma_gamma[mem_index] = mean(obj2 * theta_draw * obj4 * obj7$sigma * exp(param$sigma_gamma) * kappa_draw[,mem_index])

		d_m_draw$mean_beta_gamma[mem_index] = mean(obj2 * theta_draw * obj4 * obj7$mu * kappa_draw[,mem_index])


		m_draw2[mem_index] = mean(theta_draw[,mem_index]^2 * kappa_draw[,mem_index]^2 + 
			2 * theta_draw[,mem_index]^2 * obj2 * obj4 * obj1 * kappa_draw[,mem_index]^2 + 
			theta_draw[,mem_index]^2 * (obj2^2 + obj8) * obj9 * obj11 * kappa_draw[,mem_index]^2)

		d_m_draw2$mean_beta_gamma[mem_index] = mean(
			2 * theta_draw[,mem_index]^2 * obj2 * obj4 * obj7$mu * kappa_draw[,mem_index]^2 + 
			theta_draw[,mem_index]^2 * (obj2^2 + obj8) * obj9 * obj12$mu * kappa_draw[,mem_index]^2)

		d_m_draw2$sigma_gamma[mem_index] = mean(
			2 * theta_draw[,mem_index]^2 * obj2 * obj4 * obj7$sigma * exp(param$sigma_gamma) * kappa_draw[,mem_index]^2+ 
			theta_draw[,mem_index]^2 * (obj2^2 + obj8) * obj9 * obj12$sigma * exp(param$sigma_gamma) * kappa_draw[,mem_index]^2)

		d_m_draw2$mean_beta_omega[mem_index] = mean(
			2 * theta_draw[,mem_index]^2 * obj2 * obj6$mu * obj1 * kappa_draw[,mem_index]^2 + 
			theta_draw[,mem_index]^2 * (obj2^2 + obj8) * obj13$mu * obj11 * kappa_draw[,mem_index]^2)

		d_m_draw2$sigma_omega[mem_index] = mean(
			2 * theta_draw[,mem_index]^2 * obj2 * obj6$sigma * exp(param$sigma_omega) * obj1 * kappa_draw[,mem_index]^2 + 
			theta_draw[,mem_index]^2 * (obj2^2 + obj8) * obj13$sigma * exp(param$sigma_omega) * obj11 * kappa_draw[,mem_index]^2)

		d_m_draw2$mean_beta_delta[mem_index] = mean(
			2 * theta_draw[,mem_index]^2 * obj5$mu * obj4 * obj1 * kappa_draw[,mem_index]^2 + 
			theta_draw[,mem_index]^2 * (2 * obj2 * obj5$mu + obj10$mu) * obj9 * obj11 * kappa_draw[,mem_index]^2)

		d_m_draw2$sigma_delta[mem_index] = mean(
			2 * theta_draw[,mem_index]^2 * obj5$sigma * exp(param$sigma_delta) * obj4 * obj1 * kappa_draw[,mem_index]^2 + 
			theta_draw[,mem_index]^2 * (2 * obj2 * obj5$sigma * exp(param$sigma_delta) + obj10$sigma * exp(param$sigma_delta)) * obj9 * obj11 * kappa_draw[,mem_index]^2)
	}

	d_m_draw$mean_beta_delta = apply(X_ind, 1, function(x) x * d_m_draw$mean_beta_delta); 
	d_m_draw$mean_beta_omega = matrix(X_hh, ncol=1) %*% t(d_m_draw$mean_beta_omega); 
	d_m_draw$mean_beta_gamma = apply(X_ind, 1, function(x) x * d_m_draw$mean_beta_gamma); 

	d_m_draw2$mean_beta_delta = apply(X_ind, 1, function(x) x * d_m_draw2$mean_beta_delta); 
	d_m_draw2$mean_beta_omega = matrix(X_hh, ncol=1) %*% t(d_m_draw2$mean_beta_omega); 
	d_m_draw2$mean_beta_gamma = apply(X_ind, 1, function(x) x * d_m_draw2$mean_beta_gamma); 

	return(list(m_draw, m_draw2, d_m_draw, d_m_draw2))
}


#' Compute all permutations (with repeated draws) 
#'
#' @param vec_ A vector to draw from .
#' @param n default to 1, is the sample size.
#' 
#' @return a matrix with vec_^n rows and n columns. 
#' 
#' @export
#'
#' @examples
#' all_combn(c(1:3), 2)
#' 
all_combn = function(vec_, n = 1) {
	if (n == 1) {
		output = matrix(vec_, nrow = length(vec_)); 
		return(output)
	} else {
		output_1 = all_combn(vec_, n = n-1);
		output = do.call('rbind', lapply(1:length(vec_), function(x) cbind(output_1, vec_[x])))
		return(output)
	}
}


#' Compute all possible sets of insured individuals among eligible members
#'
#' @param vec_ A vector of indices of eligible members.
#' 
#' @return combinations of all possible insurance bundles. Each is a vector of indices of members with insurance
#' 
#' @export
#'
#' @examples
#' all_insurance(c(1,2,3))
#' 
all_insurance = function(vec_) {
	n = length(vec_)
	if (n == 1) {
		output = as.list(vec_)
		return(output)
	} else {
		output_1 = all_insurance(vec_[1:(n-1)]);
		output_2 = lapply(output_1, function(x) c(x, vec_[n]))
		return(c(output_1, output_2))
	}
}


#' Compute the draws (preference parameters independent).
#'
#' @param hh_index index (from object data_hh_list) of a particular household
#' @param param This is a list, each element is a coefficient (i.e., a vector) of a preference parameter (beta_gamma, beta_omega, beta_delta, beta_theta, beta_theta_ind) and the standard deviations of the unobserved heterogeneity (sigma_gamma, sigma_omega, sigma_delta, sigma_theta, sigma_thetabar)
#' @param n_draw_gauss this is number of Gauss  points.
#' @param n_draw_halton this is the number of halton draws, default to 1000
#' @param u_lowerbar value of utility when income constraint is violated
#' @param sick_parameters a list produced from optim 
#' @param xi_parameters a list produced from optim
#' @param short if TRUE, does not return dataframe of the original data
#'
#' @return a list that includes the draws of household-related objects,taking into account the sick parameters and the distribution of coverage. This does not take into account estimated preference parameters or unconditional distribution of health shocks. See `compute_expected_U_m` for the draws post-estimation of preference parameters and health shocks distribution.
#' 
#' @export
#'
#' @examples
#' household_draw_theta_kappa_Rdraw(1, sample_data_and_parameter$param, 1000, 10, sick_parameters_sample, xi_parameters_sample)
household_draw_theta_kappa_Rdraw = function(hh_index, param, n_draw_halton = 1000, n_draw_gauss = 10, sick_parameters, xi_parameters, u_lowerbar = -10, short=TRUE) {
	set.seed(1);
	data_hh_i = data_hh_list[[hh_index]]; 
	HHsize = nrow(data_hh_i);
	X_ind = var_ind(data_hh_i)
	X_hh = var_hh(data_hh_i)
	draw_p_xi_0 = X_ind %*% xi_parameters$par[1:ncol(X_ind)]
	draw_p_xi_1 = X_ind %*% xi_parameters$par[(ncol(X_ind) + 1):(2*ncol(X_ind))]
	p_0 = exp(draw_p_xi_0)/(1 + exp(draw_p_xi_0) + exp(draw_p_xi_1))
	p_1 = exp(draw_p_xi_1)/(1 + exp(draw_p_xi_1) + exp(draw_p_xi_0))
	sick_p = c(1/(1 + exp(X_ind %*% sick_parameters$par)))
	halton_mat = randtoolbox::halton(n_draw_halton, HHsize * 6 + 2)
	halton_mat_list = list()
	halton_mat_list$household_random_factor = qnorm(halton_mat[,HHsize * 6 + 1]) ; 
	halton_mat_list$omega = halton_mat[,HHsize * 6 + 2] ; 
	halton_mat_list$individual_factor = qnorm(halton_mat[,1:HHsize]) %>% matrix(ncol = HHsize);
	halton_mat_list$coverage = (halton_mat[,(1 * HHsize + 1):(2 * HHsize)]) %>% matrix(ncol = HHsize);
	halton_mat_list$sick = (halton_mat[,(2 * HHsize + 1):(3 * HHsize)]) %>% matrix(ncol = HHsize);
	upper_r = 5

	for (i in 1:HHsize) {
		halton_mat_list$sick[,i] = (halton_mat_list$sick[,i] >= sick_p[i])
	}
	halton_mat_list$gamma = (halton_mat[,(3 * HHsize + 1):(4 * HHsize)]) %>% matrix(ncol = HHsize);
	halton_mat_list$delta = (halton_mat[,(4 * HHsize + 1):(5 * HHsize)]) %>% matrix(ncol = HHsize);
	halton_mat_list$theta =  halton_mat[,(5 * HHsize + 1):(6 * HHsize)] %>% matrix(ncol=HHsize);
	
	kappa_draw = list(); 

	theta_draw = matrix(NA, nrow=n_draw_halton, ncol = HHsize) # observed insurance 
	kappa_draw[[1]] = matrix(NA, nrow=n_draw_halton, ncol = HHsize) # observed coinsurance 

	policy_mat_hh_index = list(); 
	policy_mat_hh_index[[1]] = policy_mat[[hh_index]];

	if (data_hh_i$HHsize_s[1] != 0) {
		elig_member_index = which(data_hh_i$Bef_sts + data_hh_i$Com_sts + data_hh_i$Std_w_ins == 0)
	}

	common_household_factor = matrix(NA, nrow=n_draw_halton, ncol = HHsize);

	income_vec = Income_net_premium[[hh_index]]
	
	R_draw = list()

	if (data_hh_i$HHsize_s[1] == 0) {
		s_thetabar = exp(param$sigma_thetabar); 
		s_theta = exp(param$sigma_theta); 

		theta_bar = matrix(NA, nrow = halton_mat %>% nrow, ncol = HHsize)
		theta_draw = halton_mat_list$theta * 0; 
		theta_draw = theta_draw %>% matrix(ncol=HHsize)
		for (i in 1:HHsize) {
			theta_bar[, i] = halton_mat_list$individual_factor[,i] * s_thetabar + halton_mat_list$household_random_factor * (X_ind[i,] %*% param$beta_theta_ind) + t(c(X_ind[i,], data_hh_i$Year[i] == 2004, data_hh_i$Year[i] == 2006, data_hh_i$Year[i] == 2010, data_hh_i$Year[i] == 2012)) %*% param$beta_theta; 
			
			theta_draw[, i] = qnorm(halton_mat_list$theta[, i] * pnorm(-theta_bar[, i]/s_theta) + (1 - halton_mat_list$theta[, i])) * s_theta + theta_bar[, i]

			theta_draw[, i] = lapply(theta_draw[, i], function(x) {
				x[which(is.na(x) | is.infinite(x) | x < 0)] = 0
				return(x)
			}) %>% unlist()

			theta_draw[, i] = theta_draw[, i] * halton_mat_list$sick[, i]

			random_xi_draws = lapply(halton_mat_list$coverage[,i], function(x) ifelse(x <= p_0[i], 0, ifelse(x <= p_0[i] + p_1[i], 1, (x - p_0[i] - p_1[i])/(1 - p_0[i] - p_1[i])))) %>% unlist()
			kappa_draw[[1]][,i] = (lapply(1:nrow(theta_draw), function(j) policy_mat_hh_index[[1]][[1]][max(which((theta_draw[j,i] * random_xi_draws[j]) >= policy_mat_hh_index[[1]][[2]][,i])),i]) %>% unlist()) * random_xi_draws + 1 - random_xi_draws
		}

		R_draw[[1]] =  lapply(income_vec[1] - rowSums(theta_draw * kappa_draw[[1]]), function(x) max(x,0)) %>% unlist()		
	}

	if (data_hh_i$HHsize_s[1] > 0) {
		beta_gamma = X_ind %*% param$beta_gamma; 
		s_gamma = exp(param$sigma_gamma); 
		gamma = matrix(NA, nrow = halton_mat %>% nrow, ncol = HHsize)
		for (i in 1:HHsize) {
			lower_threshold = pnorm(0, mean = beta_gamma[i], sd = s_gamma) 
			gamma[, i] = qnorm(lower_threshold * (1 - halton_mat_list$gamma[,i]) + 1 * halton_mat_list$gamma[,i]) * s_gamma + beta_gamma[i]; 
			which_infinite = which(is.infinite(gamma[,i]))
			if (length(which_infinite) > 0) {
				gamma[which_infinite,i] = 0; 
			}
			which_negative = which(gamma[,i] < 0)
			if (length(which_negative) > 0) {
				gamma[which_negative,i] = 0; 
			}
		}
		beta_delta = X_ind %*% param$beta_delta; 
		s_delta = exp(param$sigma_delta); 
		delta = matrix(NA, nrow = halton_mat %>% nrow, ncol = HHsize)
		for (i in 1:HHsize) {
			lower_threshold = pnorm(0, mean = beta_delta[i], sd = s_delta) 
			delta[, i] = qnorm(lower_threshold * (1 - halton_mat_list$delta[,i]) + 1 * halton_mat_list$delta[,i]) * s_delta + beta_delta[i]; 
			which_infinite = which(is.infinite(delta[,i]))
			if (length(which_infinite) > 0) {
				delta[which_infinite,i] = 0; 
			}
			which_negative = which(delta[,i] < 0)
			if (length(which_negative) > 0) {
				delta[which_negative,i] = 0; 
			}
		}
		 

		s_thetabar = exp(param$sigma_thetabar); 
		s_theta = exp(param$sigma_theta); 

		theta_bar = matrix(NA, nrow = halton_mat %>% nrow, ncol = HHsize)

		for (i in 1:HHsize) { 
			theta_bar[, i] = halton_mat_list$individual_factor[,i] * s_thetabar + halton_mat_list$household_random_factor * (X_ind[i,] %*% param$beta_theta_ind) + t(c(X_ind[i,], data_hh_i$Year[i] == 2004, data_hh_i$Year[i] == 2006, data_hh_i$Year[i] == 2010, data_hh_i$Year[i] == 2012)) %*% param$beta_theta; 
		}

		beta_omega = X_hh %*% param$beta_omega 
		s_omega = exp(param$sigma_omega)
		lower_threshold = pnorm(0, mean = beta_omega, sd = s_omega)
		omega = qnorm(lower_threshold * (1 - halton_mat_list$omega) + halton_mat_list$omega) * s_omega + beta_omega
		which_infinite = which(is.infinite(omega))
		if (length(which_infinite) > 0) {
			omega[which_infinite] = 0; 
		}
		which_negative = which(omega < 0)
		if (length(which_negative) > 0) {
			omega[which_negative] = 0; 
		}
		prob_full_insured = NULL

		root_r_vec = rep(Inf, n_draw_halton) 
		
		m = matrix(0, nrow = nrow(halton_mat), ncol = HHsize)

		random_xi_draws = matrix(NA, nrow = halton_mat %>% nrow, ncol = HHsize)

		for (i in 1:HHsize) {
			random_xi_draws[,i] = lapply(halton_mat_list$coverage[,i], function(x) ifelse(x <= p_0[i], 0, ifelse(x <= p_0[i] + p_1[i], 1, (x - p_0[i] - p_1[i])/(1 - p_0[i] - p_1[i])))) %>% unlist()
		}

		for (j in 1:n_draw_halton) {
			# Full insurance: 
			theta_draw = matrix(t(apply(halton_mat_list$theta, 1, function(x) {
				output = qnorm(x * pnorm(-theta_bar[j,]/s_theta) + (1 - x)) * s_theta + theta_bar[j,] 
				output = lapply(output, function(y) {y[which(y < 0)] = 0; return(y)}) %>% unlist()
				return(output)
			})), ncol=HHsize) * halton_mat_list$sick

			theta_draw = matrix(apply(theta_draw, 2, function(x) {
				x[which(is.na(x) | is.infinite(x))] = 0
				return(x)
			}), ncol=HHsize)

			for (i in 1:HHsize) {
				kappa_draw[[1]][,i] = (lapply(1:nrow(theta_draw), function(j) policy_mat_hh_index[[1]][[1]][max(which((theta_draw[j,i] * random_xi_draws[j,i]) >= policy_mat_hh_index[[1]][[2]][,i])),i]) %>% unlist()) * random_xi_draws[,i] + 1 - random_xi_draws[,i]
			}


			R_draw[[1]] =  lapply(income_vec[1 + length(elig_member_index)] - rowSums(theta_draw * kappa_draw[[1]]), function(x) max(x,0)) %>% unlist()	

			m[j,] = colMeans(theta_draw * kappa_draw[[1]] * (1 + matrix(apply(theta_draw * kappa_draw[[1]], 2, function(x) R_draw[[1]]^omega[j]), ncol=HHsize) * matrix(t(apply(1 + kappa_draw[[1]], 1, function(x) x^(-gamma[j,]) * delta[j,])), ncol=HHsize)), na.rm=TRUE)

			U_full_insurance_without_ulowerbar = lapply((R_draw[[1]]^(1 - omega[j]) - 1)/(1 - omega[j]), function(x) ifelse(is.infinite(x), 0, x)) %>% unlist()  - rowSums(matrix(t(apply(theta_draw, 1, function(x) x * delta[j,])), ncol=HHsize) * matrix(t(apply(kappa_draw[[1]], 1, function(x) ((x + 1)^(1 - gamma[j,]) - 1)/(1 - gamma[j,]))), ncol=HHsize))
			if (length(U_full_insurance_without_ulowerbar[which(R_draw[[1]] > 0)]) > 0) {
				u_lowerbar = min(U_full_insurance_without_ulowerbar[which(R_draw[[1]] > 0)],na.rm=TRUE)
			} else {
				u_lowerbar = 0
			}

			U_full_insurance = (U_full_insurance_without_ulowerbar * (R_draw[[1]] > 0) + (R_draw[[1]] == 0) * (income_vec[1 + length(elig_member_index)] - rowSums(theta_draw * kappa_draw[[1]]) + u_lowerbar))
			# print('utility with fully insured')
			# print(summary(U_full_insurance))
			U_drop = list()

			root_r = NULL

			for (i in elig_member_index) {
				kappa_draw[[i + 1]] = kappa_draw[[1]]; 
				kappa_draw[[i + 1]][,i] = 1;  

				R_draw[[i+1]] =  lapply(income_vec[length(elig_member_index)] - rowSums(theta_draw * kappa_draw[[i+1]]), function(x) max(x,0)) %>% unlist()	

				U_without_ulowerbar = (lapply((R_draw[[i+1]]^(1 - omega[j]) - 1)/(1 - omega[j]), function(x) ifelse(is.infinite(x), 0, x)) %>% unlist() - rowSums(matrix(t(apply(theta_draw, 1, function(x) x * delta[j,])), ncol=HHsize) * matrix(t(apply(kappa_draw[[i+1]], 1, function(x) ((x + 1)^(1 - gamma[j,]) - 1)/(1 - gamma[j,]))), ncol=HHsize)))

				if (length(U_without_ulowerbar[which(R_draw[[1 + i]] > 0)]) > 0) {
					u_lowerbar = min(U_without_ulowerbar[which(R_draw[[1 + i]] > 0)],na.rm=TRUE)
				} else {
					u_lowerbar = 0
				}

				U_drop[[i+1]] =  U_without_ulowerbar * (R_draw[[i + 1]] > 0) + (R_draw[[i + 1]] == 0) * (income_vec[length(elig_member_index)] - rowSums(theta_draw * kappa_draw[[i+1]]) + u_lowerbar) 
				# print(paste0('utility if drop member', i))
				# print(summary(U_drop[[i + 1]]))
				fr = function(r) {
					if (r != 0) {
						return((mean((1 - exp(-r * U_full_insurance))/r) - mean((1 - exp(-r * U_drop[[i+1]]))/r))^2)
					} else if (r == 0) {
						return((mean(U_full_insurance) - mean(U_drop[[i + 1]]))^2)
					}
				}


				root_r_optimize = optimize(fr, c(-1,1))
				# print(root_r_optimize)
				# if (is.na(root_r_optimize$objective)) {
				# 	theta_draw[which(is.na(U_full_insurance)),] %>% print
				# }
				root_r[i] = root_r_optimize$minimum
				# if (root_r_optimize$objective < 1e-10) {
				# 	root_r[i] = root_r_optimize$minimum
				# } else {
				# 	if ((mean(U_full_insurance) - mean(U_drop[[i + 1]])) > 0) {
				# 		root_r[i] = -1
				# 	} else {
				# 		root_r[i] = 1
				# 	}

				# }
			}
			root_r_vec[j] = max(root_r, na.rm=TRUE)
			prob_full_insured[j] = 1 - pnorm(root_r_vec[j], mean = X_hh %*% param$beta_r, sd = exp(param$sigma_r))
		}
	}


	if (data_hh_i$HHsize_s[1] > 0) {
		output = list(); 
		output$Em = colMeans(apply(m, 2, function(x) x * prob_full_insured/sum(prob_full_insured)))
		output$Prob_full = mean(prob_full_insured)
		output$root_r = root_r_vec
		output$hh_theta = halton_mat_list$household_random_factor 
		output$m = m
		output$X_hh = var_hh(data_hh_i)
	} else {
		output = list(); 
		output$R_draw = R_draw; 
		output$kappa_draw = kappa_draw; 
		output$theta_draw = theta_draw;
		output$index = hh_index
		output$income = Income_net_premium[[hh_index]]
		if (!(short)){
			output$data = data_hh_i
			output$policy_mat = policy_mat[[hh_index]]
			output$X_ind = var_ind(data_hh_i)
			output$X_hh = var_hh(data_hh_i)
			output$sick_p = c(1/(1 + exp(output$X_ind %*% sick_parameters$par)))
		}	
	}


	return(output)
}


#' Compute the matrix representing individual characteristics of a household.
#'
#' @param data_mini A dataframe of a household's characteristics
#'
#' @return a matrix with nrow being equal to the number of household members.
#' 
#' @export
#'
#' @examples
#' var_ind(data_hh_list[[1]])
#' 

var_ind = function(data_mini) {
	output = data_mini %>% mutate(relationship_2 = relationship == 2, relationship_3 = relationship == 3, relationship_4 = relationship == 4, relationship_5 = relationship == 5, age2_female = female * age2) %>% select(age2, age3, age4, age5, relationship_2, relationship_3, relationship_4, relationship_5, indshare, employed, female, married, college, age2_female)
	return(cbind(1,as.matrix(output)))
}

#' Compute the matrix representing household characteristics of a household.
#'
#' @param data_mini A dataframe of a household
#'
#' @return a vector of X_hh.
#' 
#' @export
#'
#' @examples
#' var_hh(data_hh_list[[1]])
#'
var_hh = function(data_mini) {
	output = data_mini %>% mutate(HHtype_2 = HHtype == 2, HHtype_3 = HHtype == 3, Year_2006 = Year == 2006, Year_2008 = Year == 2008, Year_2010 = Year == 2010, Year_2012= Year == 2012) %>% mutate(HHsize = HHsize/4) %>% select(HHtype_2, HHtype_3, HHsize, hhfemale, hhmaxage, hhavgeduc) %>% slice(1)
	return(cbind(1,as.matrix(output)))
}



#' Compute the LLH of a household member getting sick.
#'
#' @param data A dataframe of a household or multiple households
#' @param beta A vector of length var_ind(data)
#'
#' @return a numerical value for LLH.
#' 
#' @export
#'
#' @examples
#' llh_sick(runif(ncol(var_ind(sample_data_and_parameter[[2]]$data))),sample_data_and_parameter[[2]]$data)
#'
llh_sick = function(beta, data) {
	sick_dummy = data %>% pull(sick_dummy);
	part1 = exp(var_ind(data) %*% beta)
	p0 = 1/(1 + part1); 
	p1 = part1/(1 + part1); 
	llh = -sum((sick_dummy == 0)*log(p0)) - sum((sick_dummy == 1)*log(p1));
	return(llh)
}


#' Compute the gradient of LLH of a household member getting sick.
#'
#' @param data A dataframe of a household or multiple households
#' @param beta A vector of length var_ind(data)
#'
#' @return a vector, each is a partial derivative w.r.t beta.
#' 
#' @export
#'
#' @examples
#' grad_llh_sick(runif(ncol(var_ind(sample_data_and_parameter[[2]]$data))),sample_data_and_parameter[[2]]$data)
#'
#' 
grad_llh_sick = function(beta, data) {
	sick_dummy = data %>% pull(sick_dummy);
	var_ind_data = var_ind(data)
	part1 = exp(var_ind_data %*% beta)
	p0 = 1/(1 + part1); 
	p1 = part1/(1 + part1); 
	grad_part1 = apply(var_ind_data, 2, function(x) x * part1)
	grad_p0 = - apply(grad_part1, 2, function(x) x/(1 + part1)^2)
	grad_p1 = -grad_p0
	
	grad_llh = - (colSums(apply(grad_p0, 2, function(x) (sick_dummy == 0)*x/p0)) + 
								colSums(apply(grad_p1, 2, function(x) (sick_dummy == 1)*x/p1)));
	return(grad_llh)
}


#' Compute the LLH of a household member getting coverage
#'
#' @param data A dataframe of a household or multiple households
#' @param beta A vector of length var_ind(data)
#'
#' @return a numerical value for LLH.
#' 
#' @export
#'
#' @examples
#' llh_xi(runif(2 * ncol(var_ind(data_hh_list[[34592]]))),data_hh_list[[34592]])
#'
llh_xi = function(beta, data) {
	zeta_observed = data %>% filter(Year == 2008) %>% pull(zeta_observed);

	var_ind_data = var_ind(data %>% filter(Year == 2008))
	beta_x0 = beta[1:ncol(var_ind_data)]; 
	beta_x1 = beta[(ncol(var_ind_data) + 1):length(beta)]; 
	X_mat = var_ind_data
	

	part0 = exp(X_mat  %*% beta_x0) 
	part1 = exp(X_mat %*% beta_x1)
	p0 = part0/(1 + part0 + part1); 
	p1 = part1/(1 + part0 + part1); 
	p_mid = 1/(1 + part0 + part1); 
	llh = -sum((zeta_observed == 0)*log(p0)) - sum((zeta_observed == 1)*log(p1)) - 
		sum((zeta_observed < 1 & zeta_observed > 0)*log(p_mid));
	# message(paste0('value = ', llh))
	return(llh)
}


#' Compute the gradient of LLH of a household member getting coverage.
#'
#' @param data A dataframe of a household or multiple households
#' @param beta A vector of length var_ind(data)
#'
#' @return a vector, each is a partial derivative w.r.t beta.
#' 
#' @export
#'
#' @examples
#' grad_llh_xi(runif(2 * ncol(var_ind(data_hh_list[[34592]]))),data_hh_list[[34592]])
#'
#' 
grad_llh_xi = function(beta, data) {
	zeta_observed = data %>% filter(Year == 2008) %>% pull(zeta_observed);

	var_ind_data = var_ind(data %>% filter(Year == 2008))
	beta_x0 = beta[1:ncol(var_ind_data)]; 
	beta_x1 = beta[(ncol(var_ind_data) + 1):length(beta)]; 
	X_mat = var_ind_data

	part0 = exp(X_mat %*% beta_x0) 
	part1 = exp(X_mat %*% beta_x1)
	p0 = part0/(1 + part0 + part1); 
	p1 = part1/(1 + part0 + part1); 
	p_mid = 1/(1 + part0 + part1); 
	mat_0 = matrix(0, nrow=nrow(X_mat),ncol=ncol(X_mat))
	grad_part0 = cbind(apply(X_mat, 2, function(x) x * part0), mat_0)
	grad_part1 = cbind(mat_0, apply(X_mat, 2, function(x) x * part1))
	grad_p0 = matrix(apply(grad_part0, 2, function(x) x * (1 + part0 + part1)/(1 + part0 + part1)^2), ncol=ncol(X_mat) * 2) - matrix(apply(grad_part0 + grad_part1, 2, function(x) x * part0/(1 + part0 + part1)^2), ncol=ncol(X_mat) * 2)
	grad_p1 =  matrix(apply(grad_part1, 2, function(x) x * (1 + part0 + part1)/(1 + part0 + part1)^2), ncol=ncol(X_mat) * 2) - matrix(apply(grad_part0 + grad_part1, 2, function(x) x * part1/(1 + part0 + part1)^2), ncol=ncol(X_mat) * 2)
	grad_p_mid = -apply(grad_part0 + grad_part1,2,function(x) x /(1 + part0 + part1)^2)
	grad_llh = - (colSums(matrix(apply(grad_p0, 2, function(x) (zeta_observed == 0)*x/p0), ncol=ncol(X_mat) * 2)) + 
								colSums(matrix(apply(grad_p1, 2, function(x) (zeta_observed == 1)*x/p1), ncol=ncol(X_mat) * 2)) + 
								colSums(matrix(apply(grad_p_mid, 2, function(x) (zeta_observed < 1 & zeta_observed > 0)*x/p_mid), ncol=ncol(X_mat) * 2)));
	return(grad_llh)
}

#' Compute the LLH of medical expenditure.
#'
#' @param data_set A list produced by household_draw_theta_kappa_Rdraw
#' @param param This is a list, each element is a coefficient (i.e., a vector) of a preference parameter (beta_gamma, beta_omega, beta_delta, beta_theta, beta_theta_ind) and the standard deviations of the unobserved heterogeneity (sigma_gamma, sigma_omega, sigma_delta, sigma_theta, sigma_thetabar)
#' @param n_draw_halton number of halton draws
#' 
#' @return a list. The first element is the LLH, the second element is the derivative with respect to beta_theta, beta_theta_ind, and sigma_thetabar
#' 
#' @export
#'
#' @examples
#' identify_theta(sample_data_and_parameter[[2]], sample_data_and_parameter[[1]])
#'
#' 
identify_theta = function(data_set, param, n_draw_halton = 1000) {
	tol = 1e-3;
	oop_0_const = 1e-20;

	data_hh_i = data_set$data;
	policy_mat_hh_index = data_set$policy_mat;
	hh_index = data_set$index; 
	income_vec = data_set$income;
	X_ind = data_set$X_ind; 
	X_hh = data_set$X_hh;
	sick_p = data_set$sick_p;

	HHsize = data_hh_i$HHsize[1]
	mean_beta_theta = c(cbind(X_ind, data_hh_i$Year == 2004, data_hh_i$Year == 2006, data_hh_i$Year == 2010, data_hh_i$Year == 2012) %*% param$beta_theta);
	mean_beta_theta_ind = c(X_ind %*% param$beta_theta_ind);

	mean_beta_omega = X_hh %*% param$beta_omega; 
	mean_beta_delta = X_ind %*% param$beta_delta;
	mean_beta_gamma = X_ind %*% param$beta_gamma; 
	mean_beta_r = X_hh %*% param$beta_r;

	input_vec = c(param$sigma_thetabar, mean_beta_theta, mean_beta_theta_ind); 

	names_input_vec = c('sigma_thetabar','mean_beta_theta', 'mean_beta_theta_ind')
	size_input_vec = c(1, HHsize, HHsize)
	halton_mat_list = list()
	halton_mat = randtoolbox::halton(n_draw_halton, HHsize * 6 + 2)
	halton_mat_list$individual_factor = qnorm(halton_mat[,1:HHsize]) %>% matrix(ncol = HHsize);
	halton_mat_list$household_random_factor = qnorm(halton_mat[,HHsize * 6 + 1]) ; 

	transform_size_input_vec = function(input_vec) {
		cumsum_size_input_vec = cumsum(size_input_vec)
		input_vec_list = list()
		for (i in 1:length(size_input_vec)){
			if (i == 1) {
				input_vec_list[[names_input_vec[i]]] = input_vec[1:(cumsum_size_input_vec[i])]
			} else {
				input_vec_list[[names_input_vec[i]]] = input_vec[(cumsum_size_input_vec[i-1] + 1):cumsum_size_input_vec[i]]
			}
		}
		names(input_vec_list) = names_input_vec
		return(input_vec_list)
	}

	if (sum(size_input_vec) != length(input_vec)) {
		print(size_input_vec)
		print(length(input_vec))
		stop('Incorrect dimension')
	}

	ineligible_member = which(data_hh_i$Bef_sts + data_hh_i$Com_sts + data_hh_i$Std_w_ins == 1);
	eligible_member = which(data_hh_i$Bef_sts + data_hh_i$Com_sts + data_hh_i$Std_w_ins == 0); 

	observed_insurance = (data_hh_i$Vol_sts * (1 - (data_hh_i$Com_sts + data_hh_i$Bef_sts + data_hh_i$Std_w_ins)) + (-1) * (data_hh_i$Com_sts + data_hh_i$Bef_sts + data_hh_i$Std_w_ins)); 

	list_insurance = matrix(c(observed_insurance, income_vec[1]), nrow=1)
	observed_insurance_index = which(apply(list_insurance, 1, function(x) all(x[1:HHsize] == observed_insurance)))

	theta = data_hh_i$tot_cost_normalized;
	theta_pos_index = which(theta > 0);

	s_thetabar = exp(param$sigma_thetabar); 
	s_theta = exp(param$sigma_theta); 

	theta_bar = matrix(NA, nrow = halton_mat %>% nrow, ncol = HHsize)

	inner_f = function(input_vec) {
		input_vec_transformed = transform_size_input_vec(input_vec)
		for (i in 1:HHsize) {
			theta_bar[, i] = halton_mat_list$individual_factor[,i] * exp(input_vec_transformed$sigma_thetabar) + input_vec_transformed$mean_beta_theta_ind + input_vec_transformed$mean_beta_theta; 
		}
		likelihood = -log(apply(theta_bar, 1, function(x) prod(dnorm(theta[theta_pos_index], mean = x[theta_pos_index], sd = exp(input_vec_transformed$sigma_theta))/(1 - pnorm(0, mean=x[theta_pos_index], sd=exp(input_vec_transformed$sigma_theta))))) %>% mean + 1e-10)
		likelihood = ifelse(is.infinite(likelihood), -log(1e-10), likelihood) 
		return(likelihood)
		
	}

	if (length(theta_pos_index) > 0) {
		input_vec_all = list(); 
		input_vec_all[[1]] = input_vec; 

		
		for (i in 1:length(input_vec)) {
			input_vec_new = input_vec; 
			input_vec_new[i] = input_vec[i] + tol; 
			input_vec_all[[i + 1]] = input_vec_new; 
		}
		
		log_llh = lapply(input_vec_all, inner_f) %>% unlist()
	} else {
		log_llh = rep(0, length(input_vec) + 1)
	}
	
	f0 = log_llh[1]; 

	f1 = NULL; 
	for (i in 1:length(input_vec)) {
		f1[i] = log_llh[1 + i]; 
	}
	derivative = (f1 - f0)/tol; 
	
	derivative_param = list()
	start=1; end = start;derivative_param$sigma_thetabar = derivative[start:end];
	start = end + 1; end = start + HHsize - 1; derivative_param$beta_theta =apply(cbind(X_ind, data_hh_i$Year == 2004, data_hh_i$Year == 2006, data_hh_i$Year == 2010, data_hh_i$Year == 2012), 2, function(x) sum(x*derivative[start:end])); 
	start = end + 1; end = start + HHsize - 1; derivative_param$beta_theta_ind =apply(X_ind, 2, function(x) sum(x*derivative[start:end])); 

	return(list(f0, derivative_param))
	
}

#' Construct parameters list.
#'
#' @param param_trial A vector 
#' @param return_index LOGICAL, whether the position of the parameters is also returned. 
#' @param init LOGICAL, whether the initial value needs to supplied 
#'
#' @return if return_index is TRUE, return a list, where the first element is a list of parameters (beta_theta, beta_delta, etc..) and the second element is the position of the parameters from the original vector. If return_index is FALSE, only the list of parameters is returned. If init==TRUE, return a vector of 0 with dimension equal to the dimension required.
#' 
#' @export
#'
#' @examples
#' transform_param(return_index=TRUE, init=TRUE)
#'
#'
transform_param = function(param_trial, return_index=FALSE, init=FALSE) {
	param = list()
	index = list()
	
	if (init) {
		param_trial = rep(0, 200); 
	}

	start = 1; end = start + ncol(var_ind(data_hh_list[[1]])) + 4 -1; param$beta_theta = param_trial[start:end];  index$beta_theta = c(start:end)
	start = end + 1; end = start + ncol(var_hh(data_hh_list[[1]])) -1; param$beta_r = param_trial[start:end]; index$beta_r = c(start:end)
	start = end + 1; end = start + ncol(var_ind(data_hh_list[[1]])) -1; param$beta_gamma = param_trial[start:end]; index$beta_gamma = c(start:end)
	start = end + 1; end = start + ncol(var_ind(data_hh_list[[1]])) -1; param$beta_delta = param_trial[start:end]; index$beta_delta = c(start:end)
	# print(paste0('start index for theta = ', start, 'end index for theta = ', end))
	start = end + 1; end = start + ncol(var_ind(data_hh_list[[1]])) -1; param$beta_theta_ind = param_trial[start:end]; index$beta_theta_ind = c(start:end)
	start = end + 1; end = start; param$sigma_delta = param_trial[start] ; index$sigma_delta = c(start)
	start = end + 1; end = start; param$sigma_omega = param_trial[start] ;  index$sigma_omega = c(start)
	start = end + 1; end = start; param$sigma_theta = param_trial[start] ;  index$sigma_theta = c(start)
	start = end + 1; end = start + ncol(var_hh(data_hh_list[[1]])) -1; param$beta_omega = param_trial[start:end]; index$beta_omega = c(start:end)
	start = end + 1; end = start; param$sigma_r = (param_trial[start:end]); index$sigma_r = c(start:end)
	start = end + 1; end = start; param$sigma_gamma = (param_trial[start:end]); index$sigma_gamma = c(start:end)
	start = end + 1; end = start; param$sigma_thetabar=param_trial[start:end] ; index$sigma_thetabar = end; 
	start = end + 1; end = start; param$correlation=param_trial[end] ; index$correlation = end; 

	if (init == TRUE) {
		return(param_trial[1:end])
	}

	if (return_index) {
		return(list(param, index))
	} else {
		return(param)
	}	
}


#' Compute the counterfactual draws (preference parameters independent). This is different from household_draw_theta_kappa_Rdraw because the draws of thetabar
#'  is only conducted once, and there is no integration over thetabar. Also, the insurance bundles are not one-person deviation from the observed bundle, but
#'  all possible bundles. 
#'
#' @param hh_index index (from object data_hh_list) of a particular household
#' @param param This is a list, each element is a coefficient (i.e., a vector) of a preference parameter (beta_gamma, beta_omega, beta_delta, beta_theta, beta_theta_ind) and the standard deviations of the unobserved heterogeneity (sigma_gamma, sigma_omega, sigma_delta, sigma_theta, sigma_thetabar)
#' @param n_draw_gauss this is number of Gauss  points.
#' @param n_draw_halton this is the number of halton draws, default to 1000
#' @param sick_parameters estimated sick_parameters 
#' @param xi_parameters estimated xi_parameters 
#' @param u_lowerbar value of u when income is negative 
#' @param policy_mat_hh is the policy matrix (might be actual or a counterfactual), should follow the format of policy_mat[[hh_index]]
#' @param constraint_function is a function that can reflect whether we are imposing pure bundling or something else on the possible choice set. For example, for pure bundling, constraint_function = function(x) {x[-c(1,length(x))] = -Inf; return(x)}
#' @param within_hh_heterogeneity is a list with list names gamma, delta, and theta_bar. If a list element is true, within-hh-heterogeneity is allowed.
#' @param cap_theta_draw_normalized is an index to insure that the moment generating function is finite.
#'
#' @return a list that includes the draws of household-related objects,taking into account the sick parameters and the distribution of coverage, the optimal insurance choice, the amount of oop, and the cost to the insurance company
#' 
#' @export
#'
#' @examples
#' constant_f = function(x) x
#' counterfactual_household_draw_theta_kappa_Rdraw(3, sample_data_and_parameter$param, n_draw_halton = 1000, n_draw_gauss = 10, sick_parameters = sick_parameters_sample, xi_parameters = xi_parameters_sample, u_lowerbar = -10, policy_mat_hh=policy_mat[[3]], seed_number=1, constraint_function=constant_f, within_hh_heterogeneity = list(omega=TRUE, gamma=TRUE, delta=TRUE, theta_bar=TRUE), cap_theta_draw_normalized = 0.95)
#' 
#' 
counterfactual_household_draw_theta_kappa_Rdraw = function(hh_index, param, n_draw_halton = 1000, n_draw_gauss = 10, sick_parameters, xi_parameters, u_lowerbar = -10, policy_mat_hh, seed_number=1, constraint_function, within_hh_heterogeneity = list(omega=TRUE, gamma=TRUE, delta=TRUE, theta_bar=TRUE), cap_theta_draw_normalized = 0.95) {
	set.seed(hh_index);
	data_hh_i = data_hh_list[[hh_index]]; 
	HHsize = nrow(data_hh_i);
	X_ind = var_ind(data_hh_i)
	X_hh = var_hh(data_hh_i)
	draw_p_xi_0 = X_ind %*% xi_parameters$par[1:ncol(X_ind)]
	draw_p_xi_1 = X_ind %*% xi_parameters$par[(ncol(X_ind) + 1):(2*ncol(X_ind))]
	p_0 = exp(draw_p_xi_0)/(1 + exp(draw_p_xi_0) + exp(draw_p_xi_1))
	p_1 = exp(draw_p_xi_1)/(1 + exp(draw_p_xi_1) + exp(draw_p_xi_0))
	sick_p = c(1/(1 + exp(X_ind %*% sick_parameters$par)))
	halton_mat = randtoolbox::halton(n_draw_halton, HHsize * 6 + 2)
	halton_mat_list = list()
	halton_mat_list$household_random_factor = qnorm(halton_mat[,HHsize * 6 + 1]) ; 
	halton_mat_list$omega = halton_mat[,HHsize * 6 + 2] ; 
	halton_mat_list$individual_factor = qnorm(halton_mat[,1:HHsize]) %>% matrix(ncol = HHsize);
	halton_mat_list$coverage = (halton_mat[,(1 * HHsize + 1):(2 * HHsize)]) %>% matrix(ncol = HHsize);
	halton_mat_list$sick = (halton_mat[,(2 * HHsize + 1):(3 * HHsize)]) %>% matrix(ncol = HHsize);

	for (i in 1:HHsize) {
		halton_mat_list$sick[,i] = (halton_mat_list$sick[,i] >= sick_p[i])
	}
	halton_mat_list$gamma = (halton_mat[,(3 * HHsize + 1):(4 * HHsize)]) %>% matrix(ncol = HHsize);
	halton_mat_list$delta = (halton_mat[,(4 * HHsize + 1):(5 * HHsize)]) %>% matrix(ncol = HHsize);
	halton_mat_list$theta = qnorm(halton_mat[,(5 * HHsize + 1):(6 * HHsize)] %>% matrix(ncol=HHsize));
	beta_gamma = X_ind %*% param$beta_gamma; 
	s_gamma = exp(param$sigma_gamma); 
	gamma = NULL
	for (i in 1:HHsize) {
		lower_threshold = pnorm(0, mean = beta_gamma[i], sd = s_gamma) 
		random_draw_here = runif(1)
		gamma[i] = qnorm(lower_threshold * (1 - random_draw_here) + 1 * random_draw_here) * s_gamma + beta_gamma[i]; 
		if (is.infinite(gamma[i]) | (gamma[i] < 0)) {
			gamma[i] = 0
		}
	}

	beta_delta = X_ind %*% param$beta_delta; 
	s_delta = exp(param$sigma_delta); 
	delta = NULL
	for (i in 1:HHsize) {
		lower_threshold = pnorm(0, mean = beta_delta[i], sd = s_delta) 
		random_draw_here = runif(1)
		delta[i] = qnorm(lower_threshold * (1 - random_draw_here) + 1 * random_draw_here) * s_delta + beta_delta[i]; 
		if (is.infinite(delta[i]) | (delta[i] < 0)) {
			delta[i] = 0
		}
	}
	 
	s_thetabar = exp(param$sigma_thetabar); 
	s_theta = exp(param$sigma_theta); 
	theta_bar = NULL

	random_hh_factor = rnorm(1)

	for (i in 1:HHsize) { 
		random_draw_here = rnorm(1)
		theta_bar[i] = random_draw_here * s_thetabar + random_hh_factor * (X_ind[i,] %*% param$beta_theta_ind) + t(c(X_ind[i,], data_hh_i$Year[i] == 2004, data_hh_i$Year[i] == 2006, data_hh_i$Year[i] == 2010, data_hh_i$Year[i] == 2012)) %*% param$beta_theta; 
	}

	beta_omega = X_hh %*% param$beta_omega 
	s_omega = exp(param$sigma_omega)
	lower_threshold = pnorm(0, mean = beta_omega, sd = s_omega)
	random_draw_here = runif(1)
	omega = qnorm(lower_threshold * (1 - random_draw_here) + random_draw_here) * s_omega + beta_omega
	if (is.infinite(omega) | (omega < 0)) {
		omega = 0
	}

	if (data_hh_i$HHsize_s[1] != 0) {
		elig_member_index = which(data_hh_i$Bef_sts + data_hh_i$Com_sts + data_hh_i$Std_w_ins == 0)
	}

	# Should within-household heterogeneity be allowed: 
	if (data_hh_i$HHsize_s[1] > 0) {
		if (within_hh_heterogeneity$gamma) {
			gamma[elig_member_index] = rep(mean(gamma[elig_member_index]), data_hh_i$HHsize_s[1])
		} 
		if (within_hh_heterogeneity$theta_bar) {
			theta_bar[elig_member_index] = rep(mean(theta_bar[elig_member_index]), data_hh_i$HHsize_s[1])
		} 
		if (within_hh_heterogeneity$delta) {
			delta[elig_member_index] = rep(mean(delta[elig_member_index]), data_hh_i$HHsize_s[1])
		}
	}

	m = matrix(NA, nrow = nrow(halton_mat), ncol = HHsize)

	random_xi_draws = matrix(NA, nrow = halton_mat %>% nrow, ncol = HHsize)

	for (i in 1:HHsize) {
		random_xi_draws[,i] = lapply(halton_mat_list$coverage[,i], function(x) ifelse(x <= p_0[i], 0, ifelse(x <= p_0[i] + p_1[i], 1, (x - p_0[i] - p_1[i])/(1 - p_0[i] - p_1[i])))) %>% unlist()
	}

	theta_draw = matrix(t(apply(halton_mat_list$theta, 1, function(x) {
		output = qnorm(x * pnorm(-theta_bar/s_theta) + (1 - x)) * s_theta + theta_bar 
		output = lapply(output, function(y) {y[which(y < 0)] = 0; return(y)}) %>% unlist()
		return(output)
	})), ncol=HHsize) * halton_mat_list$sick

	theta_draw = matrix(apply(theta_draw, 2, function(x) {
		x[which(is.na(x) | is.infinite(x))] = 0
		return(x)
	}), ncol=HHsize)

	random_draw_here = runif(1)
	beta_r = X_hh %*% param$beta_r + random_hh_factor * param$correlation
	s_r = exp(param$sigma_r)
	lower_threshold = pnorm(0, mean = beta_r, sd = s_r)

	r = qnorm(lower_threshold * (1 - random_draw_here) + random_draw_here) * s_r + beta_r
	if (is.infinite(r) | (r < 0)) {
		r = 0
	}
	
	# r = qnorm(random_draw_here) * s_r + beta_r

	
	
	kappa_draw = list(); 

	kappa_draw[[1]] = matrix(NA, nrow=n_draw_halton, ncol = HHsize) # observed coinsurance 

	common_household_factor = matrix(NA, nrow=n_draw_halton, ncol = HHsize);

	income_vec = Income_net_premium[[hh_index]]
	
	random_xi_draws = matrix(NA, nrow = n_draw_halton, ncol=HHsize)
	for (i in 1:HHsize) {
		random_xi_draws[,i] = lapply(halton_mat_list$coverage[,i], function(x) ifelse(x <= p_0[i], 0, ifelse(x <= p_0[i] + p_1[i], 1, (x - p_0[i] - p_1[i])/(1 - p_0[i] - p_1[i])))) %>% unlist()
		kappa_draw[[1]][,i] = (lapply(1:nrow(theta_draw), function(j) policy_mat_hh[[1]][max(which((theta_draw[j,i] * random_xi_draws[j,i]) >= policy_mat_hh[[2]][,i])),i]) %>% unlist()) * random_xi_draws[,i] + 1 - random_xi_draws[,i]
	}

	R_draw = list()

	if (data_hh_i$HHsize_s[1] == 0) {
		R_draw[[1]] =  lapply(income_vec[1] - rowSums(theta_draw * kappa_draw[[1]]), function(x) max(x,0)) %>% unlist()		
		m = colMeans(theta_draw * kappa_draw[[1]] * (1 + matrix(apply(theta_draw * kappa_draw[[1]], 2, function(x) R_draw[[1]]^omega), ncol=HHsize) * matrix(t(apply(1 + kappa_draw[[1]], 1, function(x) x^(-gamma) * delta)), ncol=HHsize)), na.rm=TRUE)
		cost_to_insurance = colMeans(theta_draw * (1 - kappa_draw[[1]]) * (1 + matrix(apply(theta_draw * kappa_draw[[1]], 2, function(x) R_draw[[1]]^omega), ncol=HHsize) * matrix(t(apply(1 + kappa_draw[[1]], 1, function(x) x^(-gamma) * delta)), ncol=HHsize)), na.rm=TRUE)
		vol_sts_counterfactual = rep(0, HHsize);
	}


	if (data_hh_i$HHsize_s[1] > 0) {
		un_censored_R = NULL
		cara = function(x) {
			if (r != 0) {
				return(mean((1-exp(-r * x))/r, na.rm=TRUE))
			} else {
				return(mean(x, na.rm=TRUE))
			}
		}

		# compute full insurance coinsurance rates
		for (i in 1:HHsize) {
			kappa_draw[[1]][,i] = (lapply(1:nrow(theta_draw), function(j) policy_mat_hh[[1]][max(which((theta_draw[j,i] * random_xi_draws[j,i]) >= policy_mat_hh[[2]][,i])),i]) %>% unlist()) * random_xi_draws[,i] + 1 - random_xi_draws[,i] 
		}
		R_draw[[1]] =  lapply(income_vec[1 + length(elig_member_index)] - rowSums(theta_draw * kappa_draw[[1]]), function(x) max(x,0)) %>% unlist()	
		un_censored_R[data_hh_i$HHsize_s[1] + 1] = mean(income_vec[1 + length(elig_member_index)] - rowSums(theta_draw * kappa_draw[[1]]))

		U = list()

		U_without_ulowerbar = (lapply((R_draw[[1]]^(1 - omega) - 1)/(1 - omega), function(x) ifelse(is.infinite(x), 0, x)) %>% unlist()  - rowSums(matrix(t(apply(theta_draw, 1, function(x) x * delta)), ncol=HHsize) * matrix(t(apply(kappa_draw[[1]], 1, function(x) ((x + 1)^(1 - gamma) - 1)/(1 - gamma))), ncol=HHsize)))
		
		if (length(U_without_ulowerbar[which(R_draw[[1]] > 0)]) > 0) {
			u_lowerbar = min(U_without_ulowerbar[which(R_draw[[1]] > 0)],na.rm=TRUE)
		} else {
			u_lowerbar = 0
		}

		U[[data_hh_i$HHsize_s[1] + 1]] = cara(U_without_ulowerbar * (R_draw[[1]] > 0) + (R_draw[[1]] == 0) * (income_vec[1 + length(elig_member_index)] - rowSums(theta_draw * kappa_draw[[1]]) + u_lowerbar))

		# compute one-member insurance
		U_uni = NULL

		for (i in elig_member_index) {
			kappa_draw[[1 + i]] = kappa_draw[[1]]
			kappa_draw[[1 + i]][,-i] = 1
			R_draw[[1 + i]] =  lapply(income_vec[2] - rowSums(theta_draw * kappa_draw[[1 + i]]), function(x) max(x,0)) %>% unlist()	
			U_without_ulowerbar = (lapply((R_draw[[1 + i]]^(1 - omega) - 1)/(1 - omega), function(x) ifelse(is.infinite(x), 0, x)) %>% unlist()  - rowSums(matrix(t(apply(theta_draw, 1, function(x) x * delta)), ncol=HHsize) * matrix(t(apply(kappa_draw[[1 + i]], 1, function(x) ((x + 1)^(1 - gamma) - 1)/(1 - gamma))), ncol=HHsize)))
			if (length(U_without_ulowerbar[which(R_draw[[1 + i]] > 0)]) > 0) {
				u_lowerbar = min(U_without_ulowerbar[which(R_draw[[1 + i]] > 0)],na.rm=TRUE)
			} else {
				u_lowerbar = 0
			}
			U_uni[i] = cara(U_without_ulowerbar * (R_draw[[1 + i]] > 0) + (R_draw[[1 + i]] == 0) * (income_vec[2] - rowSums(theta_draw * kappa_draw[[1 + i]]) + u_lowerbar))
		}

		member_order = elig_member_index[sort(U_uni[elig_member_index], index.return=TRUE, decreasing=TRUE)$ix]


		kappa_draw_ordered = list()
		R_draw_ordered = list()
		for (n_insured in 0:(data_hh_i$HHsize_s[1] -1)) {
			if (n_insured == 0){
				kappa_draw_ordered[[n_insured + 1]] = kappa_draw[[1]]
				kappa_draw_ordered[[n_insured + 1]][,elig_member_index] = 1	
			} else {
				kappa_draw_ordered[[n_insured + 1]] = kappa_draw_ordered[[n_insured]];
				kappa_draw_ordered[[n_insured + 1]][,member_order[n_insured]] = kappa_draw[[1]][,member_order[n_insured]]
			}
			R_draw_ordered[[1 + n_insured]] =  lapply(income_vec[n_insured + 1] - rowSums(theta_draw * kappa_draw_ordered[[1 + n_insured]]), function(x) max(x,0)) %>% unlist()
			un_censored_R[1 + n_insured] = mean(income_vec[n_insured + 1] - rowSums(theta_draw * kappa_draw_ordered[[1 + n_insured]]))
			U_without_ulowerbar = (lapply((R_draw_ordered[[1 + n_insured]]^(1 - omega) - 1)/(1 - omega), function(x) ifelse(is.infinite(x), 0, x)) %>% unlist()  - rowSums(matrix(t(apply(theta_draw, 1, function(x) x * delta)), ncol=HHsize) * matrix(t(apply(kappa_draw_ordered[[1 + n_insured]], 1, function(x) ((x + 1)^(1 - gamma) - 1)/(1 - gamma))), ncol=HHsize)))
			if (length(U_without_ulowerbar[which(R_draw[[1 + n_insured]] > 0)]) > 0) {
				u_lowerbar = min(U_without_ulowerbar[which(R_draw[[1 + n_insured]] > 0)],na.rm=TRUE)
			} else {
				u_lowerbar = 0
			}
			U[[1 + n_insured]] = cara(U_without_ulowerbar * (R_draw[[1 + n_insured]] > 0) + (R_draw[[1 + n_insured]] == 0) * (income_vec[n_insured + 1] - rowSums(theta_draw * kappa_draw_ordered[[1 + n_insured]]) + u_lowerbar))
		}

		kappa_draw_ordered[[data_hh_i$HHsize_s[1] + 1]] = kappa_draw[[1]]
		R_draw_ordered[[data_hh_i$HHsize_s[1] + 1]] = R_draw[[1]]

		optimal_U_index = sort(constraint_function(unlist(U)), index.return=TRUE, decreasing=TRUE)$ix[1]

		vol_sts_counterfactual = rep(0, HHsize);

		if (optimal_U_index > 1) {
			vol_sts_counterfactual[member_order[1:(optimal_U_index - 1)]] = 1
		}

		f_wtp = function(wtp_) {
			non_censored_r = income_vec[optimal_U_index] + wtp_ - rowSums(theta_draw * kappa_draw_ordered[[optimal_U_index]])
			r_draw_here = lapply(non_censored_r, function(x) max(x,0)) %>% unlist() 

			U_without_ulowerbar = lapply(((r_draw_here)^(1 - omega) - 1)/(1 - omega), function(x) ifelse(is.infinite(x), 0, x)) %>% unlist()  - rowSums(matrix(t(apply(theta_draw, 1, function(x) x * delta)), ncol=HHsize) * matrix(t(apply(kappa_draw_ordered[[optimal_U_index]], 1, function(x) ((x + 1)^(1 - gamma) - 1)/(1 - gamma))), ncol=HHsize))
			if (length(U_without_ulowerbar[which(r_draw_here > 0)]) > 0) {
				u_lowerbar = min(U_without_ulowerbar[which(r_draw_here > 0)],na.rm=TRUE)
			} else {
				u_lowerbar = 0
			}
			output = cara(U_without_ulowerbar * (r_draw_here > 0) + (r_draw_here == 0) * (non_censored_r + u_lowerbar))
			return(output)
		}
		if (optimal_U_index == 1) {
			wtp = 0
		} else {
			init_val = -0.1

			if ((f_wtp(init_val) < unlist(U)[1]) & (f_wtp(0) >= unlist(U)[1])) {
				wtp = uniroot(function(x) {
				return(f_wtp(x) - unlist(U)[1])
				}, c(0, init_val))$root
			} else {
				wtp = optimize(function(x) {
				return((f_wtp(x) - unlist(U)[1])^2)
				}, c(0,-1))$minimum
				print(abs(f_wtp(wtp) - unlist(U)[1]))
			}
			
		}

		m = colMeans(theta_draw * kappa_draw_ordered[[optimal_U_index]] * (1 + matrix(apply(theta_draw * kappa_draw_ordered[[optimal_U_index]], 2, function(x) R_draw_ordered[[optimal_U_index]]^omega), ncol=HHsize) * matrix(t(apply(1 + kappa_draw_ordered[[optimal_U_index]], 1, function(x) x^(-gamma) * delta)), ncol=HHsize)), na.rm=TRUE)
		cost_to_insurance = colMeans(theta_draw * (1 - kappa_draw_ordered[[optimal_U_index]]) * (1 + matrix(apply(theta_draw * kappa_draw_ordered[[optimal_U_index]], 2, function(x) R_draw_ordered[[optimal_U_index]]^omega), ncol=HHsize) * matrix(t(apply(1 + kappa_draw_ordered[[optimal_U_index]], 1, function(x) x^(-gamma) * delta)), ncol=HHsize)), na.rm=TRUE)
	}

	output = list()
	output$m = m 
	output$delta = delta 
	output$omega = rep(omega, HHsize) 
	output$theta_bar = theta_bar
	output$r = rep(r, HHsize)
	output$gamma = gamma
	
	if (data_hh_i$HHsize_s[1] > 0) {
		output$vol_sts_counterfactual = vol_sts_counterfactual
		output$average_theta = colMeans(theta_draw)
		output$wtp = rep(wtp, HHsize)
		output$cost_to_insurance = cost_to_insurance
	} else {
		output$cost_to_insurance = cost_to_insurance
		output$average_theta = colMeans(theta_draw)
		output$wtp = rep(0, HHsize)
		output$vol_sts_counterfactual = 0
	}

	return(output)	
}

#' Function used in the constraint_function argument of counterfactual_household_draw_theta_kappa
#'
#' @param x a vector of size n + 1, where n is the size of eligible household members
#'
#' @return a vector of size 2 (all is insured or none is insured)
#' 
#' @export
#'
#' @examples
#' counterfactual_bundling(c(1:6))
#' 
counterfactual_bundling = function(x) {
	return(x[c(1,length(x))])
}


