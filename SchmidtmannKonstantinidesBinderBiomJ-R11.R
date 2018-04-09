## ----init, include=FALSE-------------------------------------------------
# load all libraries and functions
library(knitr)
library(xtable)
library(latex2exp)
library(flexsurv)


# ------------------------------------------------------------------------------
# function moments computes mean and variance of Wilcoxon-Mann-Whitney statistic
# in the presence of death-censored observations
# ------------------------------------------------------------------------------
# Assumptions:
# - normal one-sided test for difference $\epsilon = 0$
# - non-inferiority margin $\epsilon > 0$

# it returns $\mu_u$ and $\sigma_u$

# function parameters
# $p_0 = $ probability of censoring by death in reference group \\
# $p_1 = $ probability of censoring by death in treatment group\\
# $\mu_0 = $, mean in reference group \\
# $\mu_1 = $, mean in treatment group \\
# $\sigma =  SD(X_i), i=0, 1$ \\
# $\tau$ = time when quantitative endpoint is determined, this is needed to calculate
# the hazards \\
# $n_0 = $ sample size in reference group \\
# $n_1 = $ sample size in treatment group \\
# verbose = TRUE / FALSE controls output
# tied = TRUE / FALSE determines whether untied or tied version should be computed
# ------------------------------------------------------------------------------

moments <- function(p_0 = 0.2, p_1 = 0.3, mu_0 = 1, mu_1 = 1, sigma = 1, tau = 1,
                    n_0 = 50, n_1 = 50, verbose = FALSE, tied = FALSE){
  q_0 <- 1 - p_0
  q_1 <- 1 - p_1
  # pi_t1, pi_t2 and pi_t3 are only needed in the untied case
  if (!tied) {
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
  if (!tied) {
    pi_u1 <- p_0 * p_1 * pi_t1 + p_0 * q_1 + q_0 * q_1 * pi_x1
    pi_u2 <- p_0^2 * q_1 + p_0^2 * p_1 * pi_t2 + 2 * p_0 * q_0 * q_1 * pi_x1 + q_0^2 * q_1 * pi_x2
    pi_u3 <- p_0 * q_1^2 + p_0 * p_1^2 * pi_t3 + 2 * p_0 * p_1 * q_1 * pi_t1 + q_0 * q_1^2 * pi_x3
  }else{
    pi_u1 <- p_0 * p_1 / 2 + p_0 * q_1 + q_0 * q_1 * pi_x1
    pi_u2 <- p_0^2 * q_1 + p_0^2 * p_1 / 3 + 2 * p_0 * q_0 * q_1 * pi_x1 + q_0^2 * q_1 * pi_x2
    pi_u3 <- p_0 * q_1^2 + p_0 * p_1^2 / 3 + p_0 * p_1 * q_1  + q_0 * q_1^2 * pi_x3
  }
  if (verbose) {
    print(paste0("pi_u1 = ", pi_u1))
    print(paste0("pi_u2 = ", pi_u2))
    print(paste0("pi_u3 = ", pi_u3))
  }

  mu_u <- pi_u1
  if (!tied) {
    sigma_u <- sqrt((pi_u1 * (1 - pi_u1) + (n_0 - 1)*(pi_u2 - pi_u1^2) + (n_1 - 1)*(pi_u3 - pi_u1^2)) / (n_0 * n_1) )
  }else{
    sigma_u <- sqrt((pi_u1 * (1 - pi_u1) + (n_0 - 1)*(pi_u2 - pi_u1^2 - p_0^2 * p_1 / 12)
                                         + (n_1 - 1)*(pi_u3 - pi_u1^2 - p_0 * p_1^2 / 12) - p_0 * p_1 / 4) / (n_0 * n_1) )
  }
  if (verbose) {
    print(paste0("mu_u = ", mu_u))
    print(paste0("sigma_u = ", sigma_u))
  }
  list(mu_u, sigma_u)
  return(list(mu_u = mu_u, sigma_u = sigma_u))
}

# ------------------------------------------------------------------------------
# Function power_gr determines power for the Wilcoxon-Mann-Whitney test for non-inferiority in
# the presence of death-censored observations
# ------------------------------------------------------------------------------
# Assumptions
# - normal one-sided test for difference $\epsilon = 0$
# - non-inferiority margin $\epsilon > 0$

# epsilon is computed from \mu_0$ under $H_0$
# $H_0: \mu(U) = P(\tilde{X_{0k}} < \tilde{X_{1l}}) \leq \frac{1}{2} - \epsilon$ \\
# $H_1: \mu(U) = P(\tilde{X_{0k}} < \tilde{X_{1l}}) > \frac{1}{2} - \epsilon $ \\

# Function parameters.
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
# verbose = TRUE / FALSE controls output
# tied = TRUE / FALSE determines whether untied or tied version should be computed

power_gr <- function(alpha,
                     p_0_0 = 0.2, p_1_0 = 0.2, p_0_1 = 0.2, p_1_1 = 0.2,
                     mu_0_0 = 1, mu_1_0 = 1, mu_0_1 = 1, mu_1_1 = 1,
                     sigma = 1, tau = 1,
                     n_0 = 50, n_1 = 50, verbose = FALSE, tied = FALSE){
  # determine moments unter $H_0$
  moments_H_0 <- moments(p_0 = p_0_0, p_1 = p_1_0, mu_0 = mu_0_0, mu_1 = mu_1_0,
                         sigma = sigma, tau = tau, n_0 = n_0, n_1 = n_1,
                         verbose = verbose, tied = tied)
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

# set global chunk options
opts_chunk$set(warning = FALSE, message = FALSE)

## ----powerplot_untied, fig.lp="fig:", fig.cap = 'Power and sample size for untied case with sample size allocation for reference to new treatment 1:2. In each panel a  different relative risk (RR) of death in the new treatment group (risk $p_1$) compared to the reference group (risk $p_0$) is assumed under the null hypothesis. The power is computed for the alternative of equivalence of both treatments. Different line types correspond to different risks in the reference group.', echo=FALSE, eval=TRUE, out.width='0.8\\linewidth',out.height='0.9\\linewidth'----
# source("PowerAndMoments.R") #include functions to compute moments and power
# explore power in untied case for specific situation described above, varying the
# probabiliy of death
# sample size in each group
n_0 <- seq(1:150)
n_1 <- 2*n_0

# chosen probabilities of death
chosen_ps <- c(0, 0.01, 0.02, 0.05, 0.1, 0.2)
l_chosen_ps <- length(chosen_ps) # number of chosen probabilities of death

# chosen relative risks (treatment(1)  vs control(0) under H0)
chosen_risks <- c(1, 1.2, 1.75, 2.5)
l_chosen_risks <- length(chosen_risks) # number of chosen risks
figure_elements <- LETTERS[1:l_chosen_risks]

# vector to hold sample sizes
samplesize_80 <- c(rep(NA, l_chosen_ps))
names(samplesize_80) <- as.character(chosen_ps)

# matrix to hold power for all combinations of probabilities of death and sample size
MyPower.Matrix <- matrix(rep(NA, length(n_0)*l_chosen_ps), ncol = l_chosen_ps)
dimnames(MyPower.Matrix)[[2]] <- as.character(chosen_ps)

# vector to hold non-inferiority margins corresponding to chosen situation
epsilon_ps <- c(rep(NA, l_chosen_ps))
names(epsilon_ps) <- as.character(chosen_ps)

# vector of datasets to hold combination of chosen p's, computed epsilon's and sample
# sizes for power 80%
results <- list()

par(mfrow = c(2, 2)) # construct panel of 2x2 graphs, this only works fine if four graphs are produced

for (RR in seq_along(chosen_risks)) {
  for (p in chosen_ps) {
    MyPower <- power_gr(alpha = 0.025,
                        p_0_0 = p, p_1_0 = chosen_risks[RR]*p, p_0_1 = p, p_1_1 = p,
                        mu_0_0 = 0.3, mu_1_0 = 0.25, mu_0_1 = 0.3, mu_1_1 = 0.3,
                        sigma = 0.1, tau = 1,
                        n_0 = n_0, n_1 = n_1,
                        verbose = FALSE, tied = FALSE)
     MyPower.Matrix[, as.character(p)] <- MyPower$power
     epsilon_ps[as.character(p)] <- MyPower$epsilon
     index <- min(which(MyPower$power >= 0.8))
     samplesize_80[as.character(p)] <- n_0[index] + n_1[index]
  }
  results[[RR]] <- data.frame(chosen_ps, epsilon_ps, samplesize_80)
  if (RR == 1) {
    SummaryPower.Matrix <- MyPower.Matrix

  }else{
    SummaryPower.Matrix <- rbind(SummaryPower.Matrix, MyPower.Matrix)
  }

  # plot power as function of sample size and p
  # first establish empty plot with axes
  plot(NULL, xlab = "Total sample size", ylab = "Power", xlim = c(0, max(n_0) + max(n_1)), ylim = c(0, 1))
  my.lt <- c(1:6)     # define line types for graph
  my.col <- c(grey.colors(l_chosen_ps, start = 0, end = 0.6, gamma = 2.2))
                       # define 'colors' for graph
  # plots for p_0 ^= 0
  for (i in seq_along(chosen_ps[2:l_chosen_ps])) {
    lines(x = n_0 + n_1, y = MyPower.Matrix[, as.character(chosen_ps[i + 1])],
          lwd = 2, lty = my.lt[i + 1], col = my.col[i + 1])
  }
  # plots for p_0 = 0 (black line should always be visible)
  lines(x = n_0 + n_1, y = MyPower.Matrix[, "0"], type = "l", lwd = 2, lty = 1, col = "black")

  legend("bottomright", title = TeX("$p_1 = $ RR $ \\cdot p_0$"),
         legend = as.character(chosen_ps), inset = 0.05, bty = "n", lwd = 2,
         lty = my.lt, col = my.col )
  legend("topleft", legend = figure_elements[RR], cex = 1.5, inset = -0.05, bty = "n")
  title(main = paste0("RR = ", as.character(chosen_risks[RR])))
}
# selects every 10th element of the sample power matrix
ii <- seq(10,(4*length(n_0)),10)
TableAppendix <-  SummaryPower.Matrix[ii, ]

ii <- 3 * seq(10,length(n_0),10)
RR_label <- c(chosen_risks[1], rep("", length(ii) - 1),
              chosen_risks[2], rep("", length(ii) - 1),
              chosen_risks[3], rep("", length(ii) - 1),
              chosen_risks[4], rep("", length(ii) - 1))
TableAppendix <- data.frame(RR_label,
                            rep(ii, 4),
                            TableAppendix)


## ----power_tied, echo=FALSE, eval=TRUE-----------------------------------
# source("PowerAndMoments.R") #include functions to compute moments and power
# explore power in tied case for specific situation described above, varying the
# probabiliy of death
# constants as above
# sample size in each group: n_0, n_1
# chosen probabilities of death: chosen_ps and l_chosen_ps
# chosen relative risks: chosen_risks, l_chosen_risks

# vector to hold sample sizes
samplesize_80_tied <- c(rep(NA, l_chosen_ps))
names(samplesize_80_tied) <- as.character(chosen_ps)

# matrix to hold power for all combinations of probabilities of death and sample size
MyPower.Matrix.Tied <- matrix(rep(NA, length(n_0)*l_chosen_ps), ncol = l_chosen_ps)
dimnames(MyPower.Matrix.Tied)[[2]] <- as.character(chosen_ps)

# vector to hold non-inferiority margins corresponding to chosen situation
epsilon_ps_tied <- c(rep(NA, l_chosen_ps))
names(epsilon_ps_tied) <- as.character(chosen_ps)

# vector of datasets to hold combination of chosen p's, computed epsilon's and sample
# sizes for power 80%
results_tied <- list()

for (RR in seq_along(chosen_risks)) {
  for (p in chosen_ps) {
    MyPower.Tied <- power_gr(alpha = 0.025,
                        p_0_0 = p, p_1_0 = chosen_risks[RR]*p, p_0_1 = p, p_1_1 = p,
                        mu_0_0 = 0.3, mu_1_0 = 0.25, mu_0_1 = 0.3, mu_1_1 = 0.3,
                        sigma = 0.1, tau = 1,
                        n_0 = n_0, n_1 = n_1,
                        verbose = FALSE, tied = TRUE)
     MyPower.Matrix.Tied[, as.character(p)] <- MyPower.Tied$power
     epsilon_ps_tied[as.character(p)] <- MyPower.Tied$epsilon
     index <- min(which(MyPower.Tied$power >= 0.8))
     samplesize_80_tied[as.character(p)] <- n_0[index] + n_1[index]
  }
  results_tied[[RR]] <- data.frame(chosen_ps, epsilon_ps_tied, samplesize_80_tied)
  if (RR == 1) {
    SummaryPower.Matrix.Tied <- MyPower.Matrix.Tied

  }else{
    SummaryPower.Matrix.Tied <- rbind(SummaryPower.Matrix.Tied, MyPower.Matrix.Tied)
  }
}

# selects every 10th element of the sample power matrix
ii <- seq(10,(4*length(n_0)),10)
TableAppendix.Tied <-  SummaryPower.Matrix.Tied[ii, ]

ii <- 3 * seq(10,length(n_0),10)
TableAppendix.Tied <- data.frame(RR_label,
                            rep(ii, 4),
                            TableAppendix.Tied)

## ----xtableA, echo=FALSE, results="asis"---------------------------------
results_for_table <- data.frame(
  chosen_ps = double(),
  epsilon_ps = double(),
  samplesize_80 = double(),
  epsilon_ps_tied = double(),
  samplesize_80_tied = double()
)
RR_column <- character()
for (RR in seq_along(chosen_risks)) {
  results_for_table <- rbind(results_for_table,
                             cbind(results[[RR]], results_tied[[RR]]$epsilon_ps_tied, results_tied[[RR]]$samplesize_80_tied))
  RR_column <- c(RR_column, as.character(chosen_risks[RR]), rep(" ", l_chosen_ps - 1))
}
results_for_table <- cbind(RR_column, results_for_table)

colnames(results_for_table) <- c("RR", "$p_0$", "$\\varepsilon^{\\text{(untied)}}$",
                                 "$n^{\\text{(untied)}}_{\\text{total}}$",
                                 "$\\varepsilon^{\\text{(tied)}}$",
                                 "$n^{\\text{(tied)}}_{\\text{total}}$")
print(xtable(results_for_table,caption = 'Effective non-inferiority margin and total samples size as function of relative risk (RR) of death in the new treatment group (risk $p_1$) compared to the reference group (risk $p_0$) for given $\\alpha =0.025$ and power of 80\\% with sample size allocation for reference to new treatment 1:2.',label = "tab:test",digits = c(0, 0, 2, 3, 0, 3, 0)),
      include.rownames = FALSE,
      caption.placement = "top",
      sanitize.colnames.function = function(x){x})

## ----sim, echo=FALSE, eval=TRUE, cache = TRUE----------------------------
# set seed for random number generation for reproducibility
set.seed(270318)
# define scenarios - use parameters from above
n_total <- c(results[[1]]$samplesize_80,
             results[[2]]$samplesize_80,
             results[[3]]$samplesize_80,
             results[[4]]$samplesize_80)
n_0 <- n_total / 3
n_1 <- 2 * n_0
# number of scenarios considered
no_scenarios <- length(n_total)
# time at which quantitative endpoint is determined
tau <- 1
# risk of death by time \tau in reference group p_0 = 0, 0.01, 0.02, 0.05, 0.1, 0.2
p_0 <- rep(chosen_ps, 4)
# set vector of relative risks
RR <-  rep(chosen_risks, rep(6, 4))
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
r <- 1000
# define column names for result matrices
my_colnames <- paste0("RR = ", RR, ", p_0 = ", p_0)
my_colnames_tied <- paste0("RR = ", RR, ", p_0 = ", p_0, " (tied)")
# matrix to hold test statistics returned by function wilcox.test (untied and
# tied case)
wx <- matrix(nrow = r, ncol = no_scenarios)
wx_tied <- matrix(nrow = r, ncol = no_scenarios)
# matrix to hold Mann-Whitney statistics (derived from rank sum) (untied and
# tied case)
mw <- matrix(nrow = r, ncol = no_scenarios)
mw_tied <- matrix(nrow = r, ncol = no_scenarios)
# matrix to hold standardized Mann-Whitney statistics (untied and tied case)
mw_standardized <- matrix(nrow = r, ncol = no_scenarios)
mw_standardized_tied <- matrix(nrow = r, ncol = no_scenarios)
# matrix to hold p-values corresponding to null hypotheses given by parameters
# n_0, n_1, p_0, p_1, mu_0_0, mu_0_1, sigma, tau  (untied and tied case)
p_values <- matrix(nrow = r, ncol = no_scenarios)
p_values_tied <- matrix(nrow = r, ncol = no_scenarios)
# assign column names to result matrices
colnames(wx) <- my_colnames
colnames(wx_tied) <- my_colnames_tied
colnames(mw) <- my_colnames
colnames(mw_tied) <- my_colnames_tied
colnames(mw_standardized) <- my_colnames
colnames(mw_standardized_tied) <- my_colnames_tied
colnames(p_values) <- my_colnames
colnames(p_values_tied) <- my_colnames_tied


# choose shape parameters for log-logistic distribution
beta <- c(0.8, 1.0, 1.2)
# matrix to hold power obtained in each scenario (defined by p_0, RR and n_total)
# for each value of beta (untied and tied case)
power_beta <- matrix(nrow = no_scenarios, ncol = length(beta))
power_beta_tied <- matrix(nrow = no_scenarios, ncol = length(beta))
colnames(power_beta) <- beta
colnames(power_beta_tied) <- beta

# beta_loop
for (b in seq_along(beta)) {
  # scenario loop
  for (sc in seq(1:no_scenarios)) {
    # replicate loop
    for (i in seq(1:r)) {
      # generate death censored observations - untied
      results_gen_dco <- generate_dco(n_0 = n_0[sc], n_1 = n_1[sc],
                                      p_0 = p_0[sc], p_1 = p_1[sc],
                                      mu_0 = mu_0_1, mu_1 = mu_1_1,
                                      sigma = sigma,
                                      tau = 1,
                                      b_0 = beta, b_1 = beta,
                                      verbose = FALSE, tied = FALSE)
      # compute rank sum statistics based on $\tilde{X_{0k}}$ and $\tilde{X_{1l}}$
      wx[i, sc] <- wilcox.test(results_gen_dco$xx_1, results_gen_dco$xx_0)$statistic
      # generate death censored observations - tied
      results_gen_dco_tied <- generate_dco(n_0 = n_0[sc], n_1 = n_1[sc],
                                           p_0 = p_0[sc], p_1 = p_1[sc],
                                           mu_0 = mu_0_1, mu_1 = mu_1_1,
                                           sigma = sigma,
                                           tau = 1,
                                           b_0 = beta, b_1 = beta,
                                           verbose = FALSE, tied = TRUE)
      # compute rank sum statistics based on $\tilde{X_{0k}}$ and $\tilde{X_{1l}}$
      wx_tied[i, sc] <- wilcox.test(results_gen_dco_tied$xx_1, results_gen_dco_tied$xx_0,
                                    exact = FALSE)$statistic
    }
    # obtain Mann-Whitney statistics
    mw[ ,sc] <- wx[ ,sc] / (n_0[sc] * n_1[sc])
    mw_tied[ ,sc] <- wx_tied[ ,sc] / (n_0[sc] * n_1[sc])
    # compute moments under H_0
    moments_H_0 <- moments(n_0 = n_0[sc], n_1 = n_1[sc],
                           p_0 = p_0[sc], p_1 = p_1[sc],
                           mu_0 = mu_0_0, mu_1 = mu_1_0,
                           sigma = sigma, tau = tau,
                           verbose = verbose, tied = tied)
    moments_H_0_tied <- moments(n_0 = n_0[sc], n_1 = n_1[sc],
                                p_0 = p_0[sc], p_1 = p_1[sc],
                                mu_0 = mu_0_0, mu_1 = mu_1_0,
                                sigma = sigma, tau = tau,
                                verbose = verbose, tied = tied)
    # obtain standardized Mann-Whitney statistic
    mw_standardized[ ,sc] <- (mw[ ,sc] - moments_H_0$mu_u) / moments_H_0$sigma_u
    mw_standardized_tied[ ,sc] <- (mw_tied[ ,sc] - moments_H_0_tied$mu_u) / moments_H_0_tied$sigma_u
  }
  # compute p_values for each replicate in each scenario
  p_values <- pnorm(mw_standardized, lower.tail = FALSE)
  p_values_tied <- pnorm(mw_standardized_tied, lower.tail = FALSE)
  # obtain power by counting how many time a p-value is <= 0.025
  power_beta[, b] <- apply(p_values <= 0.025, 2, sum) / r
  power_beta_tied[, b] <- apply(p_values_tied <= 0.025, 2, sum) / r
}

# data frame to hold power determined in each scenario for each beta
power_results <- data.frame(
  RR_column,
  p_0,
  n_total,
  power_beta,
  power_beta_tied
)

## ----xtableB, echo=FALSE, results="asis"---------------------------------
colnames(power_results) <- c("RR", "$p_0$",
                                 "$n$",
                                 "untied: $\\beta = 0.8$",
                                 "$\\beta = 1.0$",
                                 "$\\beta = 1.2$",
                                 "tied: $\\beta = 0.8$",
                                 "$\\beta = 1.0$",
                                 "$\\beta = 1.2$"
                                 )
print(xtable(power_results,caption = 'Observed power as function of relative risk
             (RR) of death in the new treatment group (risk $p_1$) compared to
             the reference group (risk $p_0$), for given $\\alpha =0.025$,
             total sample size $n$, shape parameter $\\beta$ of the log-logistic time
             to event distribution. Sample size allocation for reference to new
             treatment 1:2. ',
             label = "tab:simresult",digits = c(0, 0, 2, 0, 2, 2, 2, 2, 2, 2)),
      include.rownames = FALSE,
      caption.placement = "top",
      sanitize.colnames.function = function(x){x})

## ----xtableC, echo=FALSE, results="asis"---------------------------------
colnames(TableAppendix) <- c("RR", "$n_{\\text{total}}$",
                             "$p_0 = 0$",
                             "$p_0 = 0.01$",
                             "$p_0 = 0.02$",
                             "$p_0 = 0.05$",
                             "$p_0 = 0.1$",
                             "$p_0 = 0.2$"
                             )
# to fit table on page: use abridged version for table, omit all lines
# where power is > 0.95 for all p_0's
TableAppendix <- TableAppendix[(TableAppendix[, 3] <= 0.95 |
                                  TableAppendix[, 4] <= 0.95 |
                                  TableAppendix[, 5] <= 0.95 |
                                  TableAppendix[, 6] <= 0.95 |
                                  TableAppendix[, 7] <= 0.95 |
                                  TableAppendix[, 8] <= 0.95), ]
print(xtable(TableAppendix,caption = 'Computed power as function of relative risk
             (RR) of death in the new treatment group (risk $p_1$) compared to
             the reference group (risk $p_0$), for given $\\alpha =0.025$,
             total sample size $n_{\\text{total}}$, Sample size allocation for
             reference to new treatment 1:2. Untied case.',
             label = "tab:AppendixUntied",digits = c(0, 0, 0, 3, 3, 3, 3, 3, 3)),
      include.rownames = FALSE,
      caption.placement = "top",
      sanitize.colnames.function = function(x){x})


## ----xtableD, echo=FALSE, results="asis"---------------------------------
colnames(TableAppendix.Tied) <- c("RR", "$n_{\\text{total}}$",
                             "$p_0 = 0$",
                             "$p_0 = 0.01$",
                             "$p_0 = 0.02$",
                             "$p_0 = 0.05$",
                             "$p_0 = 0.1$",
                             "$p_0 = 0.2$"
                             )
# to fit table on page: use abridged version for table, omit all lines
# where power is > 0.95 for all p_0's
TableAppendix.Tied <- TableAppendix.Tied[(TableAppendix.Tied[, 3] <= 0.95 |
                                  TableAppendix.Tied[, 4] <= 0.95 |
                                  TableAppendix.Tied[, 5] <= 0.95 |
                                  TableAppendix.Tied[, 6] <= 0.95 |
                                  TableAppendix.Tied[, 7] <= 0.95 |
                                  TableAppendix.Tied[, 8] <= 0.95), ]
print(xtable(TableAppendix.Tied,caption = 'Computed power as function of relative risk
             (RR) of death in the new treatment group (risk $p_1$) compared to
             the reference group (risk $p_0$), for given $\\alpha =0.025$,
             total sample size $n_{\\text{total}}$, Sample size allocation for
             reference to new treatment 1:2. Tied case.',
             label = "tab:AppendixTied",digits = c(0, 0, 0, 3, 3, 3, 3, 3, 3)),
      include.rownames = FALSE,
      caption.placement = "top",
      sanitize.colnames.function = function(x){x})


## ----sessioninfo, echo=FALSE---------------------------------------------
# provide session information
sessionInfo.txt <- sessionInfo()

