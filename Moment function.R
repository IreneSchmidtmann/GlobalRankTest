moments <- function(p_0 = 0.2, p_1 = 0.3, mu_0 = 1, mu_1 = 1, sigma = 1, tau=1){
  if (p_0 > 0 && p_1 > 0) {
    q_0 <- 1 - p_0
    q_1 <- 1 - p_1
    lambda_0 <- -log(q_0) / tau
    lambda_1 <- -log(q_1) / tau
    HR <- lambda_1 / lambda_0
    pi_t1 <- (p_0 * (p_1  - 1) * HR / (1 + HR) + p_1 / (1 + HR)) / (p_0 * p_1)
    print(paste0("pi_t1 = ", pi_t1))
    pi_t2 <- (p_1 + 2/(1 + HR) * (q_0 * q_1 - 1) - (q_0^2 * q_1 - 1) / (2 + HR) ) / (p_0^2 * p_1)
    print(paste0("pi_t2 = ", pi_t2))
    pi_t3 <- (p_0 * q_1^2 + 2 * q_1 * (q_0 * q_1 - 1) / (1 + HR) - (q_0 * q_1^2 - 1) / (2 + HR)) / (p_0 * p_1^2) 
    print(paste0("pi_t3 = ", pi_t3))
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
  
  mu_u <- mu_0 + mu_1
  sigma_u <- 2 * sigma
  list(mu_u, sigma_u)
  return(list(mean_u = mu_u, sigma_u = sigma_u))
}

mymoments <- moments(p_0 = 0.2, p_1 = 0.2, mu_0 = 0.3, mu_1 = 0.25, sigma = 1, tau = 1)
print(mymoments)




