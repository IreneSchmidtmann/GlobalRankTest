# Author: Irene Schmidtmann, IMBEI
# Status: 06.05.2017 00:10

# moments_untied: function to obtain mean and variance of the U statistic in the presence of
# death-censored observations
# based on the asumption that 
# - patients who die before the quantitative can be determined get ranks 
#   according to the time of death ("untied worst-rank scores")
# - time to death T_ij ~ exp(lambda_i)  (i = 0, 1)
# - the quantitative variable X_ij ~ N(mu_i, sigma_i) (i = 0, 1)

# Default values for all relevant parameters in function definition
moments_untied <- function(p_0 = 0.2, p_1 = 0.3, mu_0 = 1, mu_1 = 1, 
                           sigma_1 = 1, sigma_2=1, tau = 1, 
                           n_0 = 50, n_1 = 50, verbose = "off"){
  # Probability of survival
  q_0 <- 1 - p_0
  q_1 <- 1 - p_1
  # Non-null mortality risk in both groups
  if (p_0 > 0 && p_1 > 0){
    # Determine hazards and hazard ratio
    lambda_0 <- -log(q_0) / tau
    lambda_1 <- -log(q_1) / tau
    HR <- lambda_1 / lambda_0
    
    pi_t1 <- (p_0 * (p_1  - 1) * HR / (1 + HR) + p_1 / (1 + HR)) / (p_0 * p_1)
    pi_t2 <- (p_1 - 2 * HR / (1 + HR) * (1 - q_0 * q_1) + (1 - q_0^2 * q_1) * HR / (2 + HR) ) / (p_0^2 * p_1)
    pi_t3 <- (p_0 * q_1^2 + 2 * q_1 * (q_0 * q_1 - 1) / (1 + HR) - (q_0 * q_1^2 - 1) / (1 + 2*HR)) / (p_0 * p_1^2) 
    if (tolower(verbose) == "on") {
      print(paste0("HR = ", HR, "\n"))
      print(paste0("pi_t1 = ", pi_t1))
      print(paste0("pi_t2 = ", pi_t2))
      print(paste0("pi_t3 = ", pi_t3))
    }
  } else{
    # if there is no mortality risk in at least one group then the terms with 
    # pi_t1, pi_t2, and pi_t3 vanish, i. e. 
    pi_t1 <- 0
    pi_t2 <- 0
    pi_t3 <- 0
  }
  delta <- mu_1 - mu_0
  
  total_var <- sigma_1^2 + sigma_2^2
  delta_s <- delta / sqrt(total_var)
  rho_1 <- sigma_1^2 / total_var
  rho_2 <- sigma_2^2 / total_var

  integrand2 <- function(x){
    return(dnorm(x)*pnorm((delta_s - rho_2*x)/sqrt(1 - rho_2^2)))
  }
  
  integrand3 <- function(x){
    return(dnorm(x)*pnorm((delta_s - rho_1*x)/sqrt(1 - rho_1^2)))
  }
  
  pi_x1 <- pnorm(delta_s)
  pi_x2 <- integrate(integrand2, -Inf, delta_s)$value
  pi_x3 <- integrate(integrand3, -Inf, delta_s)$value

    if (tolower(verbose) == "on") {
    print(paste0("delta = ", delta))
    print(paste0("total_var = ", total_var))
    print(paste0("rho_1 = ", rho_1, ", rho_2 = ", rho_2))

    print(paste0("pi_x1 = ", pi_x1))
    print(paste0("pi_x2 = ", pi_x2))
    print(paste0("pi_x3 = ", pi_x3))
    }
  
  pi_u1 <- p_0 * p_1 * pi_t1 + p_0 * q_1 + q_0 * q_1 * pi_x1
  pi_u2 <- p_0^2 * q_1 + p_0^2 * p_1 * pi_t2 + 2 * p_0 * q_0 * q_1 * pi_x1 + q_0^2 * q_1 * pi_x2
  pi_u3 <- p_0 * q_1^2 + p_0 * p_1^2 * pi_t3 + 2 * p_0 * p_1 * q_1 * pi_t1 + q_0 * q_1^2 * pi_x3
  
  mu_u <- pi_u1
  sigma_u <- sqrt((pi_u1 * (1 - pi_u1) + (n_0 - 1)*(pi_u2 - pi_u1^2) + (n_1 - 1)*(pi_u3 - pi_u1^2)) / (n_0 * n_1) )
    
  list(mu_u, sigma_u)
  return(list(mean_u = mu_u, sigma_u = sigma_u))
}


# moments_tied: function to obtain mean and variance of the U statistic in the presence of
# death-censored observations
# based on the asumption that 
# - patients who die before the quantitative can be determined get identical ranks 
# - time to death T_ij ~ exp(lambda_i)  (i = 0, 1)
# - the quantitative variable X_ij ~ N(mu_i, sigma_i) (i = 0, 1)

# Default values for all relevant parameters in function definition
moments_tied <- function(p_0 = 0.2, p_1 = 0.3, mu_0 = 1, mu_1 = 1, 
                           sigma_1 = 1, sigma_2=1, tau = 1, 
                           n_0 = 50, n_1 = 50, verbose = "off"){
  # Probability of survival
  q_0 <- 1 - p_0
  q_1 <- 1 - p_1

  delta <- mu_1 - mu_0
  total_var <- sigma_1^2 + sigma_2^2
  delta_s <- delta / sqrt(total_var)
  rho_1 <- sigma_1^2 / total_var
  rho_2 <- sigma_2^2 / total_var
  
  integrand2 <- function(x){
    return(dnorm(x)*pnorm((delta_s - rho_2*x)/sqrt(1 - rho_2^2)))
  }
  
  integrand3 <- function(x){
    return(dnorm(x)*pnorm((delta_s - rho_1*x)/sqrt(1 - rho_1^2)))
  }
  
  pi_x1 <- pnorm(delta_s)
  pi_x2 <- integrate(integrand2, -Inf, delta_s)$value
  pi_x3 <- integrate(integrand3, -Inf, delta_s)$value
  
  if (tolower(verbose) == "on"){
    print(paste0("delta = ", delta))
    print(paste0("total_var = ", total_var))
    print(paste0("rho_1 = ", rho_1, ", rho_2 = ", rho_2))
    
    print(paste0("pi_x1 = ", pi_x1))
    print(paste0("pi_x2 = ", pi_x2))
    print(paste0("pi_x3 = ", pi_x3))
  }
  
  pi_u1 <- q_0 * q_1 * pi_x1 + p_0 * q_1 + 0.5 * p_0 * p_1
  pi_u2 <- q_0^2 * q_1 * pi_x2 + 2 * p_0 * q_0 * q_1 * pi_x1 + p_0^2 * q_1  + p_0^2 * p_1  / 3 
  pi_u3 <- q_0 * q_1^2 * pi_x3 + p_0 * q_1^2  + p_0 * p_1 * q_1 + p_0 * p_1^2 / 3
  
  mu_u <- pi_u1
  sigma_u <- sqrt((pi_u1 * (1 - pi_u1) + (n_0 - 1)*(pi_u2 - pi_u1^2 - p_0^2 * p_1 / 12) 
                   + (n_1 - 1)*(pi_u3 - pi_u1^2 - p_0 * p_1^2 / 12) - p_0 * p_1 / 4) / (n_0 * n_1) )
  
  list(mu_u, sigma_u)
  return(list(mean_u = mu_u, sigma_u = sigma_u))
}

# =============================
# below: testing function calls
moments_untied(p_0 = 0.1, p_1 = 0.4, mu_0 = 0.3, mu_1 = 0.2, sigma_1 = 0.1, sigma_2 = 0.1, tau = 1, 
               n_0 = 48, n_1 = 96, verbose = "on")
moments_tied(p_0 = 0.1, p_1 = 0.4, mu_0 = 0.3, mu_1 = 0.2, sigma_1 = 0.1, sigma_2 = 0.1, tau = 1, 
               n_0 = 48, n_1 = 96, verbose = "on")

# Moments under the hypothesis that both treatments are identical with zero mortality
mymoments_H_id <- moments_untied(p_0 = 0, p_1 = 0, mu_0 = 0.3, mu_1 = 0.3, sigma_1 = 0.1, sigma_2 = 0.1, tau = 1, 
                          n_0 = 48, n_1 = 96, verbose = "on")
print(mymoments_H_id)


z_1_a <- qnorm(1 - 0.0294)


# Moments under the hypothesis that both treatments have zero mortality but 
# differ with respect to quantitative endpoint
mymoments_H_diff1 <- moments_untied(p_0 = 0, p_1 = 0, mu_0 = 0.25, mu_1 = 0.3, sigma_1 = 0.1, sigma_2 = 0.1, tau = 1, 
                             n_0 = 48, n_1 = 96, verbose = "on")
print(mymoments_H_diff1)
# Power for ordinary one-sided Mann-Whitney test
power_H_diff1 <- pnorm((mymoments_H_diff1$mean_u - mymoments_H_id$mean_u) / mymoments_H_diff1$sigma_u - 
                         z_1_a * mymoments_H_id$sigma_u / mymoments_H_diff1$sigma_u )
print(power_H_diff1)

noether_power <- pnorm((mymoments_H_diff1$mean_u - mymoments_H_id$mean_u) / mymoments_H_id$sigma_u - 
                         z_1_a )
print(noether_power)

# Power for non-inferiority test with H0 as in H_diff_1


# Moments under the hypothesis that both treatments have identical non-zero mortality but 
# differ with respect to quantitative endpoint
mymoments_H_diff2 <- moments(p_0 = 0.1, p_1 = 0.1, mu_0 = 0.3, mu_1 = 0.25, sigma = 0.1, tau = 1, 
                             n_0 = 48, n_1 = 96)
print(mymoments_H_diff2)


# Moments under the hypothesis that both treatments differ with respect to mortality and 
# with respect to quantitative endpoint
mymoments_H_diff3 <- moments(p_0 = 0.05, p_1 = 0.5, mu_0 = 0.3, mu_1 = 0.25, sigma = 0.1, tau = 1, 
                             n_0 = 48, n_1 = 96)
print(mymoments_H_diff3)


# apply moments function to vectors of parameters
p_0 <- c(0, 0)
p_1 <- c(0, 0)
mu_0 <- c(0.25, 0.25)
mu_1 <- c(0.3, 0.3)
sigma <- c(0.1, 0.1)
tau <- c(1, 1)
n_0 <- c(48, 48)
n_1 <- c(96, 96)

result.matrix <- mapply(moments, p_0, p_1, mu_0, mu_1, sigma, tau, n_0, n_1)

