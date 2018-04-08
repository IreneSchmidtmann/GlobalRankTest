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
# define scenarios
n_total <- c(147, 153, 162, 186, 237, 390,
             147, 153, 156, 174, 204, 276,
             147, 147, 147, 147, 144, 129,
             147, 141, 135, 120, 96, 60)
n_0 <- n_total / 3
n_1 <- 2 * n_0
# number of scenarios considered
no_scenarios <- length(n_total)
# time at which quantitative endpoint is determined
tau <- 1
# risk of death by time \tau in reference group p_0 = 0, 0.01, 0.02, 0.05, 0.1, 0.2
p_0 <- rep(c(0, 0.01, 0.02, 0.05, 0.1, 0.2), 4)
# set vector of relative risks
RR <-  rep(c(1, 1.2, 1.75, 2.5), rep(6, 4))
# risk of death by time \tau in new treatment group
p_1 <- RR * p_0
# means of quantitative variable under null hypothesis
mu_0_0 <- 0.3
mu_1_0 <- 0.25
# means of quantitative variable under alternative hypothesis
mu_0_1 <- 0.3
mu_1_1 <- 0.3
# variance of quantitative variable
sigma <- 0.1
# consider untied case
tied <- FALSE
verbose <- FALSE

# number of replicates
r <- 100
# define column names for result matrices
my_colnames <- paste0("RR = ", RR, ", p_0 = ", p_0)
# matrix to hold test statistics returned by function wilcox.test
wx <- matrix(nrow = r, ncol = no_scenarios)
# matrix to hold Mann-Whitney statistics (derived from rank sum)
mw <- matrix(nrow = r, ncol = no_scenarios)
# matrix to hold standardized Mann-Whitney statistics
mw_standardized <- matrix(nrow = r, ncol = no_scenarios)
# matrix to hold p-values corresponding to null hypotheses given by parameters
# n_0, n_1, p_0, p_1, mu_0_0, mu_0_1, sigma, tau
p_values <- matrix(nrow = r, ncol = no_scenarios)
# assign column names to result matrices
colnames(wx) <- my_colnames
colnames(mw) <- my_colnames
colnames(p_values) <- my_colnames

# choose shape parameters for log-logistic distribution
beta <- c(0.8, 1.0, 1.2)
# matrix to hold power obtained in each scenario (defined by p_0, RR and n_total)
# for each value of beta
power_beta <- matrix(nrow = no_scenarios, ncol = length(beta))
colnames(power_beta) <- beta

# beta_loop
for (b in seq_along(beta)) {
  # scenario loop
  for (sc in seq(1:no_scenarios)) {
    # replicate loop
    for (i in seq(1:r)) {
      # generate death censored observations
      results_gen_dco <- generate_dco(n_0 = n_0[sc], n_1 = n_1[sc],
                                      p_0 = p_0[sc], p_1 = p_1[sc],
                                      mu_0 = mu_0_1, mu_1 = mu_1_1,
                                      sigma = sigma,
                                      tau = 1,
                                      b_0 = beta, b_1 = beta,
                                      verbose = FALSE, tied = FALSE)
      # compute rank sum statistics based on $\tilde{X_{0k}}$ and $\tilde{X_{1l}}$
      wx[i, sc] <- wilcox.test(results_gen_dco$xx_1, results_gen_dco$xx_0)$statistic
    }
    # obtain Mann-Whitney statistic
    mw[ ,sc] <- wx[ ,sc] / (n_0[sc] * n_1[sc])
    # compute moments under H_0
    moments_H_0 <- moments(n_0 = n_0[sc], n_1 = n_1[sc],
                           p_0 = p_0[sc], p_1 = p_1[sc],
                           mu_0 = mu_0_0, mu_1 = mu_1_0,
                           sigma = sigma, tau = tau,
                           verbose = verbose, tied = tied)
    # obtain standardized Mann-Whitney statistic
    mw_standardized[ ,sc] <- (mw[ ,sc] - moments_H_0$mu_u) / moments_H_0$sigma_u
  }
  # compute p_values for each replicate in each scenario
  p_values <- pnorm(mw_standardized, lower.tail = FALSE)
  # obtain power by counting how many time a p-value is <= 0.025
  power_beta[, b] <- apply(p_values <= 0.025, 2, sum) / r
}

# data frame to hold power determined in each scenario for each beta
power_results <- data.frame(
  RR_column,
  p_0,
  n_total,
  power_beta
)