
# function moments computes mean and variance of global rank statistic
# it returns $\mu_u$ and $\sigma_u$

moments <- function(p_0 = 0.2, p_1 = 0.3, mu_0 = 1, mu_1 = 1, sigma = 1, tau = 1, 
                    n_0 = 50, n_1 = 50, verbose = FALSE){
  q_0 <- 1 - p_0
  q_1 <- 1 - p_1
  # Non-null mortality risk in both groups
  if (p_0 > 0 && p_1 > 0) {
    lambda_0 <- -log(q_0) / tau
    lambda_1 <- -log(q_1) / tau
    HR <- lambda_1 / lambda_0
    pi_t1 <- (p_0 * (p_1  - 1) * HR / (1 + HR) + p_1 / (1 + HR)) / (p_0 * p_1)
    pi_t2 <- (p_1 - 2 * HR / (1 + HR) * (1 - q_0 * q_1) + (1 - q_0^2 * q_1) * HR / (2 + HR) ) / (p_0^2 * p_1)
    pi_t3 <- (p_0 * q_1^2 + 2 * q_1 * (q_0 * q_1 - 1) / (1 + HR) - (q_0 * q_1^2 - 1) / (1 + 2*HR)) / (p_0 * p_1^2) 
    if (verbose) {
      print(paste0("HR = ", HR))
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
  shift <- delta / sigma
  if (verbose) {
    print(paste0("delta = ", delta))
    print(paste0("shift = ", shift))
  } 
  
  integrand1 <- function(x) {
    return(dnorm(x - shift)*pnorm(x))
  }  
  integrand2 <- function(x) {
    return(dnorm(x - shift)*(pnorm(x))^2)
  }  

  integrand3 <- function(x) {
    return(dnorm(x)*(pnorm(x - shift))^2)
  }  

  pi_x1 <- integrate(integrand1, -Inf, Inf)$value
  pi_x2 <- integrate(integrand2, -Inf, Inf)$value
  pi_x3 <- 2*pi_x1 - 1 + integrate(integrand3, -Inf, Inf)$value

  if (verbose) {
    print(paste0("pi_x1 = ", pi_x1))
    print(paste0("pi_x2 = ", pi_x2))
    print(paste0("pi_x3 = ", pi_x3))
  } 
  pi_u1 <- p_0 * p_1 * pi_t1 + p_0 * q_1 + q_0 * q_1 * pi_x1
  pi_u2 <- p_0^2 * q_1 + p_0^2 * p_1 * pi_t2 + 2 * p_0 * q_0 * q_1 * pi_x1 + q_0^2 * q_1 * pi_x2
  pi_u3 <- p_0 * q_1^2 + p_0 * p_1^2 * pi_t3 + 2 * p_0 * p_1 * q_1 * pi_t1 + q_0 * q_1^2 * pi_x3
  if (verbose) {
    print(paste0("pi_u1 = ", pi_u1))
    print(paste0("pi_u2 = ", pi_u2))
    print(paste0("pi_u3 = ", pi_u3))
  } 
  
  mu_u <- pi_u1
  sigma_u <- sqrt((pi_u1 * (1 - pi_u1) + (n_0 - 1)*(pi_u2 - pi_u1^2) + (n_1 - 1)*(pi_u3 - pi_u1^2)) / (n_0 * n_1) )
  if (verbose) {
    print(paste0("mu_u = ", mu_u))
    print(paste0("sigma_u = ", sigma_u))
  }
  list(mu_u, sigma_u)
  return(list(mu_u = mu_u, sigma_u = sigma_u))
}

# Determine power for global rank test
# - normal one-sided test for difference $\epsilon = 0$
# - non-inferiority margin $\epsilon > 0$
# epsilon is computed from \mu_0$ under $H_0$
# $H_0: \mu(U) = P(\tilde{X_{0k}} < \tilde{X_{1l}}) \leq \frac{1}{2} - \epsilon$ \\
# $H_1: \mu(U) = P(\tilde{X_{0k}} < \tilde{X_{1l}}) > \frac{1}{2} - \epsilon $ \\
# Parameters \\
# $p_0_0 = p_0$ under $H_0$ (probability of censoring by death in reference group under $H_0$) \\
# $p_1_0 = p_1$ under $H_0$ (probability of censoring by death in treatment group under $H_0$) \\
# $p_0_1 = p_0$ under $H_1$ (probability of censoring by death in reference group under $H_1$) \\
# $p_1_1 = p_1$ under $H_1$ (probability of censoring by death in treatment group under $H_1$) \\
# $\mu_0_0 = \mu_0$ under $H_0$ (mean in reference group under $H_0$) \\
# $\mu_1_0 = \mu_1$ under $H_0$ (mean in treatment group under $H_0$) \\
# $\mu_0_1 = \mu_0$ under $H_1$ (mean in reference group under $H_1$) \\
# $\mu_1_1 = \mu_1$ under $H_1$ (mean in treatment group under $H_1$) \\
# $\sigma =  var(X_i), i=0, 1$ under $H_0$ and $H_1 $ \\
# $\tau$ = time when quantitative endpoint is determined, this is needed to calculate 
# the hazards \\
# $n_0 = $ sample size in reference group \\
# $n_1 = $ sample size in treatment group \\

power_gr <- function(alpha, 
                     p_0_0 = 0.2, p_1_0 = 0.2, p_0_1 = 0.2, p_1_1 = 0.2, 
                     mu_0_0 = 1, mu_1_0 = 1, mu_0_1 = 1, mu_1_1 = 1, 
                     sigma = 1, tau = 1,
                     n_0 = 50, n_1 = 50, verbose = FALSE){
  # determine moments unter $H_0$
  moments_H_0 <- moments(p_0 = p_0_0, p_1 = p_1_0, mu_0 = mu_0_0, mu_1 = mu_1_0, 
                         sigma = sigma, tau = tau, n_0 = n_0, n_1 = n_1, verbose = verbose)
  # compute effective epsilon, given parameters under $H_0$
  epsilon <- 0.5 - moments_H_0$mu_u
  
  if (verbose) {
    print("moments unter H_0")
    print(paste0("mu_u_0 = ", moments_H_0$mu_u))
    print(paste0("sigma_u_0 = ", moments_H_0$sigma_u))
  }
  # determine moments unter $H_1$
  moments_H_1 <- moments(p_0 = p_0_1, p_1 = p_1_1, mu_0 = mu_0_1, mu_1 = mu_1_1, 
                         sigma = sigma, tau = tau, n_0 = n_0, n_1 = n_1, verbose = verbose)
  if (verbose) {
    print("moments unter H_1")
    print(paste0("mu_u_1 = ", moments_H_1$mu_u))
    print(paste0("sigma_u_1 = ", moments_H_1$sigma_u))
  }
  # compute $1 - \alpha$ percentile of the standard normal distribution
  z_alpha <- qnorm(alpha, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  power = pnorm((z_alpha * moments_H_0$sigma_u  + moments_H_1$mu_u - moments_H_0$mu_u) / moments_H_1$sigma_u)
  list(power, epsilon)
  return(list(power = power, epsilon = epsilon))
}
