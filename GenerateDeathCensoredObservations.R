library(flexsurv)

# ------------------------------------------------------------------------------
# Function generate_dco generates death censored observations for two groups.
#
# The time to death follows a log-logistic function with scale parameters
# $a_0, a_1$ and shape parameters $b_0, b_1$ respectively. Scale parameters are
# determined based on probability of death until time $\tau$
# The quantitative parameter follows a normal distribution with means
# $\mu_0, mu_1$ respectively, and commmon variance $\sigma^2$
#
# Time to death and the quantitative endpoint are here assumed to be independent
# ------------------------------------------------------------------------------

# it returns a list vectors consisting of
# $x_0$ uncensored values for the quantitative variable in the reference group
# $x_1$ uncensored values for the quantitative variable in the treatment group
# $t_0$ time to death in the reference group
# $t_1$ time to death in the treatment group
# $xx_0 = x_0$ for uncensored observations,
#              $min(x_0, x_1) - 1 + t_0$ for the censored observations in the
#              reference group in the untied case
#              $min(x_0, x_1) - 1$ for the censored observations in the
#              reference group in the tied case
# $xx_1 = x_1$ for uncensored observations,
#              $min(x_0, x_1) - 1 + t_1$ for the censored observations in the
#              reference group in the untied case
#              $min(x_0, x_1) - 1$ for the censored observations in the
#              reference group in the tied case


# function parameters
# $n_0 = $ sample size in reference group \\
# $n_1 = $ sample size in treatment group \\
# $p_0 = $ probability of censoring by death in reference group \\
# $p_1 = $ probability of censoring by death in treatment group\\
# $\mu_0 = $, mean in reference group \\
# $\mu_1 = $, mean in treatment group \\
# $\sigma =  SD(X_i), i=0, 1$ \\
# $\tau$ = time when quantitative endpoint is determined, this is needed to obtain
# the scale parameters of the log-logistic distribution
# $b_0 = $ shape parameter in reference group \\
# $b_1 = $ shape parameter in treatment group \\
# verbose = TRUE / FALSE controls output
# tied = TRUE / FALSE determines whether untied or tied version should be computed
#

generate_dco <- function(n_0 = 50, n_1 = 50,
                         p_0 = 0.2, p_1 = 0.2,
                         mu_0 = 1, mu_1 = 1, sigma = 1, tau = 1,
                         b_0 = 1, b_1 = 1,
                         verbose = FALSE, tied = FALSE){
# $a_0 = $ scale parameter in reference group, derived from p_0
# $a_1 = $ scale parameter in treatment group, derived from p_1
  q_0 <- 1 - p_0
  q_1 <- 1 - p_1
  t_0 <- rep(NA, n_0)
  t_1 <- rep(NA, n_1)
  if (p_0 > 0) {
    # determine scale parameter
    a_0 <- tau * (q_0 / p_0)^(1/b_0)
    # generate time to death variable
    t_0 <- rllogis(n_0, shape = b_0, scale = a_0)
    # generate death indicators
    d_0 <- ifelse(t_0 <= tau, 1, 0)
  }else{
    d_0 <- rep(0, n_0)
  }
  if (p_1 > 0) {
    # determine scale parameter
    a_1 <- tau * (q_1 / p_1)^(1/b_1)
    # generate time to death variable
    t_1 <- rllogis(n_1, shape = b_1, scale = a_1)
    # generate death indicators
    d_1 <- ifelse(t_1 <= tau, 1, 0)
  }else{
    d_1 <- rep(0, n_1)
  }
  # generate quantiative endpoints
  x_0 <- rnorm(n_0, mean = mu_0, sd = sigma)
  x_1 <- rnorm(n_1, mean = mu_1, sd = sigma)
  eta <- min(c(x_0, x_1)) - 1 - tau
  if (!tied) {
    xx_0 <- ifelse(d_0 == 1, eta + t_0, x_0)
    xx_1 <- ifelse(d_1 == 1, eta + t_1, x_1)
  }else {
    xx_0 <- ifelse(d_0 == 1, eta, x_0)
    xx_1 <- ifelse(d_1 == 1, eta, x_1)
  }
  if (verbose) {
    print(paste0("t_0 = ", round(t_0, digits = 3)))
    print(paste0("t_1 = ", round(t_1, digits = 3)))
    print(paste0("d_0 = ", round(d_0, digits = 3)))
    print(paste0("d_1 = ", round(d_1, digits = 3)))
    print(paste0("x_0 = ", round(x_0, digits = 3)))
    print(paste0("x_1 = ", round(x_1, digits = 3)))
    print(paste0("eta = ", round(eta, digits = 3)))
    print(paste0("xx_0 = ", round(xx_0, digits = 3)))
    print(paste0("xx_1 = ", round(xx_1, digits = 3)))
  }
  list(xx_0, xx_1, x_0, x_1, t_0, t_1)
  return(list(xx_0 = xx_0, xx_1 = xx_1,
              x_0 = x_0, x_1 = x_1,
              t_0 = t_0, t_1 = t_1))
}

# set seed for random number generation for reproducibility
set.seed(270318)
no_scenarios <- 3
n_0 <- c(49, 51, 54)
n_1 <- 2 * n_0
p_0 <- c(0.00, 0.01, 0.02)
RR <- c(1, 1, 1)
p_1 <- RR * p_0
tied <- FALSE
mu_0_0 <- 0.3
mu_1_0 <- 0.25
mu_0_1 <- 0.3
mu_1_1 <- 0.3
sigma <- 0.1
tau <- 1
verbose <- FALSE

r <- 10
wx <- matrix(nrow = r, ncol = no_scenarios)
mw <- matrix(nrow = r, ncol = no_scenarios)
for (sc in (1:no_scenarios)) {
  for (i in seq(1:r)) {
    results_gen_dco <- generate_dco(n_0 = n_0, n_1 = n_1,
                                    p_0 = p_0, p_1 = p_1,
                                    mu_0 = mu_0_1, mu_1 = mu_1_1,
                                    sigma = sigma,
                                    tau = 1,
                                    b_0 = 1, b_1 = 1,
                                    verbose = verbose, tied = FALSE)
    wx[i, sc] <- wilcox.test(results_gen_dco$xx_1, results_gen_dco$xx_0)$statistic
  }
}

mw <- wx / (n_0 * n_1)

moments_H_0 <- moments(p_0 = p_0, p_1 = p_1, mu_0 = mu_0_0, mu_1 = mu_1_0,
                       sigma = sigma, tau = tau, n_0 = n_0, n_1 = n_1,
                       verbose = verbose, tied = tied)
mw_standardized <- (mw - moments_H_0$mu_u) / moments_H_0$sigma_u
p_values <- pnorm(mw_standardized, lower.tail = FALSE)
table(p_values <= 0.025)
