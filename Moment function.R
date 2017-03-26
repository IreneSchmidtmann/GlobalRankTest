moments <- function(p_0 = 0.2, p_1 = 0.3, mu_0 = 1, mu_1 = 1, sigma = 1, tau = 1, 
                    n_0 = 50, n_1 = 50){
  if (p_0 > 0 && p_1 > 0) {
    q_0 <- 1 - p_0
    q_1 <- 1 - p_1
    lambda_0 <- -log(q_0) / tau
    lambda_1 <- -log(q_1) / tau
    HR <- lambda_1 / lambda_0
    pi_t1 <- (p_0 * (p_1  - 1) * HR / (1 + HR) + p_1 / (1 + HR)) / (p_0 * p_1)
    print(paste0("pi_t1 = ", pi_t1))
    pi_t2 <- (p_1 - 2 * HR / (1 + HR) * (1 - q_0 * q_1) + (1 - q_0^2 * q_1) * HR / (2 + HR) ) / (p_0^2 * p_1)
    print(paste0("pi_t2 = ", pi_t2))
    pi_t3 <- (p_0 * q_1^2 + 2 * q_1 * (q_0 * q_1 - 1) / (1 + HR) - (q_0 * q_1^2 - 1) / (1 + 2*HR)) / (p_0 * p_1^2) 
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
  
  pi_u1 <- p_0 * p_1 * pi_t1 + p_0 * q_1 + q_0 * q_1 * pi_x1
  pi_u2 <- p_0^2 * q_1 + p_0^2 * p_1 * pi_t2 + 2 * p_0 * q_0 * q_1 * pi_x1 + q_0^2 * q_1 * pi_x2
  pi_u3 <- p_0 * q_1^2 + p_0 * p_1^2 * pi_t3 + 2 * p_0 * p_1 * q_1 * pi_t1 + q_0 * q_1^2 * pi_x3
  
  mu_u <- pi_u1
  sigma_u <- sqrt((pi_u1 * (1 - pi_u1) + (n_0 - 1)*(pi_u2 - pi_u1^2) + (n_1 - 1)*(pi_u3 - pi_u1^2)) / (n_0 * n_1) )
    
  list(mu_u, sigma_u)
  return(list(mean_u = mu_u, sigma_u = sigma_u))
}

mymoments <- moments(p_0 = 0.0000001, p_1 = 0.0000001, mu_0 = 0.3, mu_1 = 0.3, sigma = 1, tau = 1, 
                     n_0 = 48, n_1 = 96)
print(mymoments)

z_1 <- (0.5 - mymoments$mean_u) / mymoments$sigma_u
print(z_1)

## for evaluating pi_t2 - Wolfram alpha syntax
## integrate (1 - exp(-a_0 v)) a_1 exp(-a_1 v) from v=u to t
