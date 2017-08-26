moments <- function(p_0 = 0.2, p_1 = 0.3, mu_0 = 1, mu_1 = 1, sigma = 1, tau = 1, 
                    n_0 = 50, n_1 = 50){
  q_0 <- 1 - p_0
  q_1 <- 1 - p_1
  # Non-null mortality risk in both groups
  if (p_0 > 0 && p_1 > 0) {
    lambda_0 <- -log(q_0) / tau
    lambda_1 <- -log(q_1) / tau
    HR <- lambda_1 / lambda_0
    print(paste0("HR = ", HR))
    pi_t1 <- (p_0 * (p_1  - 1) * HR / (1 + HR) + p_1 / (1 + HR)) / (p_0 * p_1)
    print(paste0("pi_t1 = ", pi_t1))
    pi_t2 <- (p_1 - 2 * HR / (1 + HR) * (1 - q_0 * q_1) + (1 - q_0^2 * q_1) * HR / (2 + HR) ) / (p_0^2 * p_1)
    print(paste0("pi_t2 = ", pi_t2))
    pi_t3 <- (p_0 * q_1^2 + 2 * q_1 * (q_0 * q_1 - 1) / (1 + HR) - (q_0 * q_1^2 - 1) / (1 + 2*HR)) / (p_0 * p_1^2) 
    print(paste0("pi_t3 = ", pi_t3))
  } else{
    # if there is no mortality risk in at least one group then the terms with 
    # pi_t1, pi_t2, and pi_t3 vanish, i. e. 
    pi_t1 <- 0
    pi_t2 <- 0
    pi_t3 <- 0
  }
  delta <- mu_1 - mu_0
  print(delta)
  
  shift <- delta / sigma
  
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
  print(paste0("pi_x1 = ", pi_x1))
  
  pi_x2 <- integrate(integrand2, -Inf, Inf)$value
  print(paste0("pi_x2 = ", pi_x2))

  pi_x3 <- 2*pi_x1 - 1 + integrate(integrand3, -Inf, Inf)$value
  print(paste0("pi_x3 = ", pi_x3))
  
  pi_u1 <- p_0 * p_1 * pi_t1 + p_0 * q_1 + q_0 * q_1 * pi_x1
  pi_u2 <- p_0^2 * q_1 + p_0^2 * p_1 * pi_t2 + 2 * p_0 * q_0 * q_1 * pi_x1 + q_0^2 * q_1 * pi_x2
  pi_u3 <- p_0 * q_1^2 + p_0 * p_1^2 * pi_t3 + 2 * p_0 * p_1 * q_1 * pi_t1 + q_0 * q_1^2 * pi_x3
  
  mu_u <- pi_u1
  sigma_u <- sqrt((pi_u1 * (1 - pi_u1) + (n_0 - 1)*(pi_u2 - pi_u1^2) + (n_1 - 1)*(pi_u3 - pi_u1^2)) / (n_0 * n_1) )
    
  list(mu_u, sigma_u)
  return(list(mean_u = mu_u, sigma_u = sigma_u))
}

# Moments under the hypothesis that both treatments are identical with zero mortality
mymoments_H_id <- moments(p_0 = 0, p_1 = 0, mu_0 = 0.3, mu_1 = 0.3, sigma = 0.1, tau = 1, 
                          n_0 = 48, n_1 = 96)
print(mymoments_H_id)


z_1_a <- qnorm(1 - 0.0294)


# Moments under the hypothesis that both treatments have zero mortality but 
# differ with respect to quantitative endpoint
mymoments_H_diff1 <- moments(p_0 = 0, p_1 = 0, mu_0 = 0.25, mu_1 = 0.3, sigma = 0.1, tau = 1, 
                             n_0 = 48, n_1 = 96)
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
mymoments_H_diff3 <- moments(p_0 = 0.1, p_1 = 0.137, mu_0 = 0.3, mu_1 = 0.25, sigma = 0.1, tau = 1, 
                             n_0 = 48, n_1 = 96)
print(mymoments_H_diff3)


