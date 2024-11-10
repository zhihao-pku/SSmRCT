#' Power of mRCT using Japan's Method 1 for count endpoints
#'
#' Based on Japan's Method 1, given the global and target region sample sizes, calculate and simulate the marginal probabilities, conditional probabilities, and joint probabilities of global success and efficacy consistency between target region and globally, in clinical trials using superiority, non-inferiority, and equivalence designs with count endpoints.
#'
#' @rdname getPwr_Count_JM1
#'
#' @name getPwr_Count_JM1
#'
#' @param delta_j log(RR) between treatment and control groups for target region.
#' @param delta_nj log(RR) between treatment and control groups for other regions. When \code{delta_nj} is not \code{NA}, \code{delta_a} will be calculated automatically.
#' @param delta_a log(RR) between treatment and control groups globally.
#' @param lambda0_j Baseline hazard of control group for target region.
#' @param lambda0_nj Baseline hazard of control group for other regions. When \code{lambda0_nj} is not \code{NA}, \code{lambda0_a} will be calculated automatically.
#' @param lambda0_a Baseline hazard of control group globally.
#' @param t Average exposure time.
#' @param k The over-dispersion parameter for negative binomial distribution, which is 0 for poisson distribution.
#' @param f Proportion of sample size allocated to target region.
#' @param pi Proportion of global efficacy to retain. The default value is 0.5, which means retaining half of the efficacy.
#' @param cut A positive value for non-inferiority or equivalence margin. For example, if the non-inferiority margin for RR is 0.6, then \code{cut = -log(0.6)}. If the non-inferiority margin for RR is 1.3, then \code{cut = log(1.3)}.
#' @param alpha One-sided type I error rate for global success. The default value is 0.025.
#' @param beta Type II error rate for global success, which is used to calculate global sample size only when \code{N} is \code{NA}.
#' @param N Global sample size. When \code{N} is \code{NA} and \code{beta} is not \code{NA}, \code{N} will be calculated automatically.
#' @param r Ratio of the sample sizes of the treatment group to the control group. The default value is 1.
#' @param direct \code{direct = 1} indicates that a larger RR is preferable, while \code{direct = -1} indicates that a smaller RR is preferable.
#' @param sim Logical value. When set to \code{FALSE}, theoretical calculation is performed. When set to \code{TRUE}, simulation is used, which is more time-consuming.
#' @param nsim Number of simulations.
#' @param seed Random seed for simulation.
#' @param numcore Number of CPU cores to use during simulation.
#'
#' @return A data frame containing input parameters and returned power.
#' \describe{
#'   \item{\code{pwr1 }}{The marginal probability of global success.}
#'   \item{\code{pwr2 }}{The marginal probability that the target region efficacy is consistent with the global efficacy.}
#'   \item{\code{pwr3 }}{The conditional probability that the target region efficacy is consistent with the global efficacy given global success.}
#'   \item{\code{pwr4 }}{The joint probability of global success and the target region efficacy being consistent with the global efficacy.}
#' }
#'
#' @details
#' Taking the larger RR is preferable as an example. The global success criterion and the efficacy consistency criterion between target region and globally
#'
#' in superiority design:
#' \deqn{Z_a = \frac{\hat \delta_a}{\sqrt{Var(\hat \delta_a)}} > \Phi^{-1}(1 - \alpha)}
#' \deqn{\hat \delta_j - \pi\hat \delta_a > 0}
#'
#' in non-inferiority design:
#' \deqn{Z_a = \frac{\hat \delta_a + \Delta}{\sqrt{Var(\hat \delta_a)}} > \Phi^{-1}(1 - \alpha)}
#' \deqn{\hat \delta_j - \hat \delta_a + \pi\Delta > 0}
#'
#' in equivalence design:
#' \deqn{Z_{a_u} = \frac{\hat \delta_a + \Delta}{\sqrt{Var(\hat \delta_a)}} > \Phi^{-1}(1 - \alpha)\text{ and }Z_{a_l} = \frac{\hat \delta_a - \Delta}{\sqrt{Var(\hat \delta_a)}} < \Phi^{-1}(\alpha)}
#' \deqn{\hat \delta_j - \hat \delta_a + \pi\Delta > 0\text{ and }\hat \delta_j - \hat \delta_a - \pi\Delta < 0}
#'
#' Where \eqn{\hat \delta = log(\hat {RR})} between treatment and control groups, and \eqn{\Delta} is the non-inferiority or equivalence margin (\code{cut}).
#'
#' @references
#' 1. Quan H, Li M, Chen J, et al. Assessment of Consistency of Treatment Effects in Multiregional Clinical Trials. Drug Information J. 2010;44(5):617-632. doi:10.1177/009286151004400509
#'
#' 2. Liao JJZ, Yu Z, Li Y. Sample size allocation in multiregional equivalence studies. Pharm Stat. 2018;17(5):570-577. doi:10.1002/pst.1871
#'
#' @export
#'
#' @examples
#' getPwr_Count_Super_JM1(
#'   delta_j = log(1.2),
#'   delta_a = log(1.3),
#'   lambda0_j = 0.1, lambda0_a = 0.1, t = 5, k = 0, f = seq(0.1, 0.9, 0.1),
#'   pi = 0.5, alpha = 0.025, beta = NA, N = 300, r = 1, sim = FALSE
#' )
#'
#' # delta_a will be calculated based on delta_j and delta_nj,
#' # and lambda0_a will be calculated based on lambda0_j and lambda0_nj.
#' # Global sample size will be calculated based on beta.
#' getPwr_Count_Noninf_JM1(
#'   delta_j = log(1.1),
#'   delta_nj = log(1.0),
#'   lambda0_j = 0.1, lambda0_nj = 0.1, t = 5, k = 0, f = seq(0.1, 0.9, 0.1),
#'   pi = 0.5, cut = log(1.3),
#'   alpha = 0.025, beta = 0.2, N = NA, r = 1, direct = -1, sim = FALSE
#' )
getPwr_Count_Super_JM1 <- function(delta_j, delta_nj = NA, delta_a = NA, lambda0_j, lambda0_nj = NA, lambda0_a = NA, t, k = 0, f, pi = 0.5, alpha = 0.025, beta = NA, N = NA, r = 1, sim = FALSE, nsim = 1000, seed = 0, numcore = 4) {
  isNA_delta_nj <- is.na(delta_nj)
  isNA_delta_a <- is.na(delta_a)
  isNA_lambda0_nj <- is.na(lambda0_nj)
  isNA_lambda0_a <- is.na(lambda0_a)
  if (!sim) {
    eg <- as.data.frame(expand.grid(delta_j = delta_j, delta_nj = delta_nj, delta_a = delta_a, lambda0_j = lambda0_j, lambda0_nj = lambda0_nj, lambda0_a = lambda0_a, t = t, k = k, f = f, pi = pi, alpha = alpha, beta = beta, N = N, r = r, stringsAsFactors = FALSE))
    res <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      delta_j <- R$delta_j
      delta_nj <- R$delta_nj
      delta_a <- R$delta_a
      lambda0_j <- R$lambda0_j
      lambda0_nj <- R$lambda0_nj
      lambda0_a <- R$lambda0_a
      t <- R$t
      k <- R$k
      f <- R$f
      pi <- R$pi
      alpha <- R$alpha
      beta <- R$beta
      N <- R$N
      r <- R$r
      if (isNA_delta_nj & (!isNA_delta_a)) {
        delta_nj <- (delta_a - delta_j * f) / (1 - f)
      }
      if (isNA_delta_a & (!isNA_delta_nj)) {
        delta_a <- delta_j * f + delta_nj * (1 - f)
      }
      if (isNA_lambda0_nj & (!isNA_lambda0_a)) {
        lambda0_nj <- (lambda0_a - lambda0_j * f) / (1 - f)
      }
      if (isNA_lambda0_a & (!isNA_lambda0_nj)) {
        lambda0_a <- lambda0_j * f + lambda0_nj * (1 - f)
      }
      if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
        N <- getN_Count_Super(delta = delta_a, lambda0 = lambda0_a, t = t, k = k, alpha = alpha, beta = beta, N = NA, r = r)$N
      }
      Nj <- N * f
      lambda1_j <- exp(delta_j) * lambda0_j
      lambda1_a <- exp(delta_a) * lambda0_a
      sigma1_j <- sqrt(1 / lambda1_j / t + k)
      sigma0_j <- sqrt(1 / lambda0_j / t + k)
      sigma1_a <- sqrt(1 / lambda1_a / t + k)
      sigma0_a <- sqrt(1 / lambda0_a / t + k)
      var_j <- sigma1_j^2 / (r * Nj / (1 + r)) + sigma0_j^2 / (Nj / (1 + r))
      var_a <- sigma1_a^2 / (r * N / (1 + r)) + sigma0_a^2 / (N / (1 + r))
      sej <- sqrt(var_j + pi^2 * var_a - 2 * pi * sqrt(f) * sqrt(var_j * var_a))
      uj <- (delta_j - pi * delta_a) / sej
      se <- sqrt(var_a)
      u <- delta_a / se
      cov <- sqrt(f) * sqrt(var_a * var_j) - pi * var_a
      corr <- cov / (sej * se)
      M <- matrix(c(1, corr, corr, 1), nrow = 2, byrow = T)
      if (delta_a < 0) {
        uj <- (-1) * uj
        u <- (-1) * u
      }
      pwr1 <- pmvnorm(lower = c(qnorm(1 - alpha), -Inf), upper = c(Inf, Inf), mean = c(u, uj), corr = M)
      pwr2 <- pmvnorm(lower = c(-Inf, 0), upper = c(Inf, Inf), mean = c(u, uj), corr = M)
      pwr3 <- pmvnorm(lower = c(qnorm(1 - alpha), 0), upper = c(Inf, Inf), mean = c(u, uj), corr = M)
      pwr4 <- pwr3 / pwr1
      data.frame(delta_a, delta_j, delta_nj, lambda0_a, lambda0_j, lambda0_nj, t, k, f, pi, alpha, beta, N, r, pwr1, pwr2, pwr3, pwr4)
    }, .options = furrr_options(seed = TRUE))
  }
  if (sim) {
    if (isNA_delta_nj & (!isNA_delta_a)) {
      warning("When given global and target region parameters, the results of the data simulation may be incorrect; please refer to the theoretical calculation results instead.")
    }
    eg <- as.data.frame(expand.grid(delta_j = delta_j, delta_nj = delta_nj, delta_a = delta_a, lambda0_j = lambda0_j, lambda0_nj = lambda0_nj, lambda0_a = lambda0_a, t = t, k = k, f = f, pi = pi, alpha = alpha, beta = beta, N = N, r = r, nsim = 1:nsim, stringsAsFactors = FALSE))
    if (numcore >= 2) {
      plan(multisession, workers = numcore)
    }
    ss <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      delta_j <- R$delta_j
      delta_nj <- R$delta_nj
      delta_a <- R$delta_a
      lambda0_j <- R$lambda0_j
      lambda0_nj <- R$lambda0_nj
      lambda0_a <- R$lambda0_a
      t <- R$t
      k <- R$k
      f <- R$f
      pi <- R$pi
      alpha <- R$alpha
      beta <- R$beta
      N <- R$N
      r <- R$r
      if (isNA_delta_nj & (!isNA_delta_a)) {
        delta_nj <- (delta_a - delta_j * f) / (1 - f)
      }
      if (isNA_delta_a & (!isNA_delta_nj)) {
        delta_a <- delta_j * f + delta_nj * (1 - f)
      }
      if (isNA_lambda0_nj & (!isNA_lambda0_a)) {
        lambda0_nj <- (lambda0_a - lambda0_j * f) / (1 - f)
      }
      if (isNA_lambda0_a & (!isNA_lambda0_nj)) {
        lambda0_a <- lambda0_j * f + lambda0_nj * (1 - f)
      }
      if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
        N <- getN_Count_Super(delta = delta_a, lambda0 = lambda0_a, t = t, k = k, alpha = alpha, beta = beta, N = NA, r = r)$N
      }
      Nj <- N * f
      set.seed(i + seed)
      if (k == 0) {
        xt_j <- rpois(n = Nj * r / (1 + r), lambda = lambda0_j * exp(delta_j) * t)
        xc_j <- rpois(n = Nj / (1 + r), lambda = lambda0_j * t)
        xt_nj <- rpois(n = (N - Nj) * r / (1 + r), lambda = lambda0_nj * exp(delta_nj) * t)
        xc_nj <- rpois(n = (N - Nj) / (1 + r), lambda = lambda0_nj * t)
      }
      if (k > 0) {
        xt_j <- rnbinom(n = Nj * r / (1 + r), size = k, mu = lambda0_j * exp(delta_j) * t)
        xc_j <- rnbinom(n = Nj / (1 + r), size = k, mu = lambda0_j * t)
        xt_nj <- rnbinom(n = (N - Nj) * r / (1 + r), size = k, mu = lambda0_nj * exp(delta_nj) * t)
        xc_nj <- rnbinom(n = (N - Nj) / (1 + r), size = k, mu = lambda0_nj * t)
      }
      xt <- c(xt_j, xt_nj)
      xc <- c(xc_j, xc_nj)
      dat_a <- data.frame(x = c(xt, xc), trt = c(rep(1, length(xt)), rep(0, length(xc))))
      dat_j <- data.frame(x = c(xt_j, xc_j), trt = c(rep(1, length(xt_j)), rep(0, length(xc_j))))
      if (k == 0) {
        fit_j <- glm(x ~ trt, dat = dat_j, family = poisson(link = "log"))
        fit_a <- glm(x ~ trt, dat = dat_a, family = poisson(link = "log"))
      }
      if (k > 0) {
        fit_j <- glm.nb(x ~ trt, dat = dat_j)
        fit_a <- glm.nb(x ~ trt, dat = dat_a)
      }
      coef_j <- coef(fit_j)[2]
      coef_a <- coef(fit_a)[2]
      za <- summary(fit_a)$coefficients[2, 3]
      zj <- coef_j - pi * coef_a
      if (delta_a < 0) {
        za <- (-1) * za
        zj <- (-1) * zj
      }
      succ_a <- if_else(za > qnorm(1 - alpha), 1, 0)
      succ_j <- if_else(zj > 0, 1, 0)
      data.frame(delta_a, delta_j, delta_nj, lambda0_a, lambda0_j, lambda0_nj, t, k, f, pi, alpha, beta, N, r, succ_a, succ_j)
    }, .progress = TRUE, .options = furrr_options(seed = TRUE))
    if (numcore >= 2) {
      plan(sequential)
    }
    res <- suppressMessages({
      ss %>%
        group_by(delta_a, delta_j, delta_nj, lambda0_a, lambda0_j, lambda0_nj, t, k, f, pi, alpha, beta, N, r) %>%
        summarise(pwr1 = mean(succ_a, na.rm = TRUE), pwr2 = mean(succ_j, na.rm = TRUE), pwr3 = mean(succ_a & succ_j, na.rm = TRUE), pwr4 = mean(succ_j[succ_a == 1], na.rm = TRUE)) %>%
        arrange(f) %>%
        as.data.frame()
    })
  }
  return(res)
}


#' @rdname getPwr_Count_JM1
#' @export
getPwr_Count_Noninf_JM1 <- function(delta_j, delta_nj = NA, delta_a = NA, lambda0_j, lambda0_nj = NA, lambda0_a = NA, t, k = 0, f, pi = 0.5, cut, alpha = 0.025, beta = NA, N = NA, r = 1, direct = 1, sim = FALSE, nsim = 1000, seed = 0, numcore = 4) {
  isNA_delta_nj <- is.na(delta_nj)
  isNA_delta_a <- is.na(delta_a)
  isNA_lambda0_nj <- is.na(lambda0_nj)
  isNA_lambda0_a <- is.na(lambda0_a)
  if (!sim) {
    eg <- as.data.frame(expand.grid(delta_j = delta_j, delta_nj = delta_nj, delta_a = delta_a, lambda0_j = lambda0_j, lambda0_nj = lambda0_nj, lambda0_a = lambda0_a, t = t, k = k, f = f, pi = pi, cut = cut, alpha = alpha, beta = beta, N = N, r = r, stringsAsFactors = FALSE))
    res <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      delta_j <- R$delta_j
      delta_nj <- R$delta_nj
      delta_a <- R$delta_a
      lambda0_j <- R$lambda0_j
      lambda0_nj <- R$lambda0_nj
      lambda0_a <- R$lambda0_a
      t <- R$t
      k <- R$k
      f <- R$f
      pi <- R$pi
      cut <- R$cut
      alpha <- R$alpha
      beta <- R$beta
      N <- R$N
      r <- R$r
      if (isNA_delta_nj & (!isNA_delta_a)) {
        delta_nj <- (delta_a - delta_j * f) / (1 - f)
      }
      if (isNA_delta_a & (!isNA_delta_nj)) {
        delta_a <- delta_j * f + delta_nj * (1 - f)
      }
      if (isNA_lambda0_nj & (!isNA_lambda0_a)) {
        lambda0_nj <- (lambda0_a - lambda0_j * f) / (1 - f)
      }
      if (isNA_lambda0_a & (!isNA_lambda0_nj)) {
        lambda0_a <- lambda0_j * f + lambda0_nj * (1 - f)
      }
      if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
        N <- getN_Count_Noninf(delta = delta_a, lambda0 = lambda0_a, t = t, k = k, cut = cut, alpha = alpha, beta = beta, N = NA, r = r, direct = direct)$N
      }
      Nj <- N * f
      lambda1_j <- exp(delta_j) * lambda0_j
      lambda1_a <- exp(delta_a) * lambda0_a
      sigma1_j <- sqrt(1 / lambda1_j / t + k)
      sigma0_j <- sqrt(1 / lambda0_j / t + k)
      sigma1_a <- sqrt(1 / lambda1_a / t + k)
      sigma0_a <- sqrt(1 / lambda0_a / t + k)
      var_j <- sigma1_j^2 / (r * Nj / (1 + r)) + sigma0_j^2 / (Nj / (1 + r))
      var_a <- sigma1_a^2 / (r * N / (1 + r)) + sigma0_a^2 / (N / (1 + r))
      sej <- sqrt(var_j + var_a - 2 * sqrt(f) * sqrt(var_j * var_a))
      uj <- if_else(direct == 1, (delta_j - delta_a + pi * cut) / sej, (delta_j - delta_a - pi * cut) / sej)
      se <- sqrt(var_a)
      u <- if_else(direct == 1, (delta_a + cut) / se, (delta_a - cut) / se)
      cov <- sqrt(f) * sqrt(var_a * var_j) - var_a
      corr <- cov / (sej * se)
      M <- matrix(c(1, corr, corr, 1), nrow = 2, byrow = T)
      if (direct == -1) {
        uj <- (-1) * uj
        u <- (-1) * u
      }
      pwr1 <- pmvnorm(lower = c(qnorm(1 - alpha), -Inf), upper = c(Inf, Inf), mean = c(u, uj), corr = M)
      pwr2 <- pmvnorm(lower = c(-Inf, 0), upper = c(Inf, Inf), mean = c(u, uj), corr = M)
      pwr3 <- pmvnorm(lower = c(qnorm(1 - alpha), 0), upper = c(Inf, Inf), mean = c(u, uj), corr = M)
      pwr4 <- pwr3 / pwr1
      data.frame(delta_a, delta_j, delta_nj, lambda0_a, lambda0_j, lambda0_nj, t, k, f, pi, cut, alpha, beta, N, r, direct, pwr1, pwr2, pwr3, pwr4)
    }, .options = furrr_options(seed = TRUE))
  }
  if (sim) {
    if (isNA_delta_nj & (!isNA_delta_a)) {
      warning("When given global and target region parameters, the data simulation results may be incorrect; please refer to the theoretical calculation results instead.")
    }
    eg <- as.data.frame(expand.grid(delta_j = delta_j, delta_nj = delta_nj, delta_a = delta_a, lambda0_j = lambda0_j, lambda0_nj = lambda0_nj, lambda0_a = lambda0_a, t = t, k = k, f = f, pi = pi, cut = cut, alpha = alpha, beta = beta, N = N, r = r, nsim = 1:nsim, stringsAsFactors = FALSE))
    if (numcore >= 2) {
      plan(multisession, workers = numcore)
    }
    ss <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      delta_j <- R$delta_j
      delta_nj <- R$delta_nj
      delta_a <- R$delta_a
      lambda0_j <- R$lambda0_j
      lambda0_nj <- R$lambda0_nj
      lambda0_a <- R$lambda0_a
      t <- R$t
      k <- R$k
      f <- R$f
      pi <- R$pi
      cut <- R$cut
      alpha <- R$alpha
      beta <- R$beta
      N <- R$N
      r <- R$r
      if (isNA_delta_nj & (!isNA_delta_a)) {
        delta_nj <- (delta_a - delta_j * f) / (1 - f)
      }
      if (isNA_delta_a & (!isNA_delta_nj)) {
        delta_a <- delta_j * f + delta_nj * (1 - f)
      }
      if (isNA_lambda0_nj & (!isNA_lambda0_a)) {
        lambda0_nj <- (lambda0_a - lambda0_j * f) / (1 - f)
      }
      if (isNA_lambda0_a & (!isNA_lambda0_nj)) {
        lambda0_a <- lambda0_j * f + lambda0_nj * (1 - f)
      }
      if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
        N <- getN_Count_Noninf(delta = delta_a, lambda0 = lambda0_a, t = t, k = k, cut = cut, alpha = alpha, beta = beta, N = NA, r = r, direct = direct)$N
      }
      Nj <- N * f
      set.seed(i + seed)
      if (k == 0) {
        xt_j <- rpois(n = Nj * r / (1 + r), lambda = lambda0_j * exp(delta_j) * t)
        xc_j <- rpois(n = Nj / (1 + r), lambda = lambda0_j * t)
        xt_nj <- rpois(n = (N - Nj) * r / (1 + r), lambda = lambda0_nj * exp(delta_nj) * t)
        xc_nj <- rpois(n = (N - Nj) / (1 + r), lambda = lambda0_nj * t)
      }
      if (k > 0) {
        xt_j <- rnbinom(n = Nj * r / (1 + r), size = k, mu = lambda0_j * exp(delta_j) * t)
        xc_j <- rnbinom(n = Nj / (1 + r), size = k, mu = lambda0_j * t)
        xt_nj <- rnbinom(n = (N - Nj) * r / (1 + r), size = k, mu = lambda0_nj * exp(delta_nj) * t)
        xc_nj <- rnbinom(n = (N - Nj) / (1 + r), size = k, mu = lambda0_nj * t)
      }
      xt <- c(xt_j, xt_nj)
      xc <- c(xc_j, xc_nj)
      dat_a <- data.frame(x = c(xt, xc), trt = c(rep(1, length(xt)), rep(0, length(xc))))
      dat_j <- data.frame(x = c(xt_j, xc_j), trt = c(rep(1, length(xt_j)), rep(0, length(xc_j))))
      if (k == 0) {
        fit_j <- glm(x ~ trt, dat = dat_j, family = poisson(link = "log"))
        fit_a <- glm(x ~ trt, dat = dat_a, family = poisson(link = "log"))
      }
      if (k > 0) {
        fit_j <- glm.nb(x ~ trt, dat = dat_j)
        fit_a <- glm.nb(x ~ trt, dat = dat_a)
      }
      coef_j <- coef(fit_j)[2]
      coef_a <- coef(fit_a)[2]
      za <- if_else(direct == 1, summary(fit_a)$coefficients[2, 3] + cut / summary(fit_a)$coefficients[2, 2], summary(fit_a)$coefficients[2, 3] - cut / summary(fit_a)$coefficients[2, 2])
      zj <- if_else(direct == 1, coef_j - coef_a + pi * cut, coef_j - coef_a - pi * cut)
      if (direct == -1) {
        zj <- (-1) * zj
        za <- (-1) * za
      }
      succ_a <- if_else(za > qnorm(1 - alpha), 1, 0)
      succ_j <- if_else(zj > 0, 1, 0)
      data.frame(delta_a, delta_j, delta_nj, lambda0_a, lambda0_j, lambda0_nj, t, k, f, pi, cut, alpha, beta, N, r, direct, succ_a, succ_j)
    }, .progress = TRUE, .options = furrr_options(seed = TRUE))
    if (numcore >= 2) {
      plan(sequential)
    }
    res <- suppressMessages({
      ss %>%
        group_by(delta_a, delta_j, delta_nj, lambda0_a, lambda0_j, lambda0_nj, t, k, f, pi, alpha, beta, N, r, direct) %>%
        summarise(pwr1 = mean(succ_a, na.rm = TRUE), pwr2 = mean(succ_j, na.rm = TRUE), pwr3 = mean(succ_a & succ_j, na.rm = TRUE), pwr4 = mean(succ_j[succ_a == 1], na.rm = TRUE)) %>%
        arrange(f) %>%
        as.data.frame()
    })
  }
  return(res)
}

#' @rdname getPwr_Count_JM1
#' @export
getPwr_Count_Equi_JM1 <- function(delta_j, delta_nj = NA, delta_a = NA, lambda0_j, lambda0_nj = NA, lambda0_a = NA, t, k = 0, f, pi = 0.5, cut, alpha = 0.025, beta = NA, N = NA, r = 1, sim = FALSE, nsim = 1000, seed = 0, numcore = 4) {
  isNA_delta_nj <- is.na(delta_nj)
  isNA_delta_a <- is.na(delta_a)
  isNA_lambda0_nj <- is.na(lambda0_nj)
  isNA_lambda0_a <- is.na(lambda0_a)
  if (!sim) {
    eg <- as.data.frame(expand.grid(delta_j = delta_j, delta_nj = delta_nj, delta_a = delta_a, lambda0_j = lambda0_j, lambda0_nj = lambda0_nj, lambda0_a = lambda0_a, t = t, k = k, f = f, pi = pi, cut = cut, alpha = alpha, beta = beta, N = N, r = r, stringsAsFactors = FALSE))
    res <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      delta_j <- R$delta_j
      delta_nj <- R$delta_nj
      delta_a <- R$delta_a
      lambda0_j <- R$lambda0_j
      lambda0_nj <- R$lambda0_nj
      lambda0_a <- R$lambda0_a
      t <- R$t
      k <- R$k
      f <- R$f
      pi <- R$pi
      cut <- R$cut
      alpha <- R$alpha
      beta <- R$beta
      N <- R$N
      r <- R$r
      if (isNA_delta_nj & (!isNA_delta_a)) {
        delta_nj <- (delta_a - delta_j * f) / (1 - f)
      }
      if (isNA_delta_a & (!isNA_delta_nj)) {
        delta_a <- delta_j * f + delta_nj * (1 - f)
      }
      if (isNA_lambda0_nj & (!isNA_lambda0_a)) {
        lambda0_nj <- (lambda0_a - lambda0_j * f) / (1 - f)
      }
      if (isNA_lambda0_a & (!isNA_lambda0_nj)) {
        lambda0_a <- lambda0_j * f + lambda0_nj * (1 - f)
      }
      if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
        N <- getN_Count_Equi(delta = delta_a, lambda0 = lambda0_a, t = t, k = k, cut = cut, alpha = alpha, beta = beta, N = NA, r = r)$N
      }
      Nj <- N * f
      lambda1_j <- exp(delta_j) * lambda0_j
      lambda1_a <- exp(delta_a) * lambda0_a
      sigma1_j <- sqrt(1 / lambda1_j / t + k)
      sigma0_j <- sqrt(1 / lambda0_j / t + k)
      sigma1_a <- sqrt(1 / lambda1_a / t + k)
      sigma0_a <- sqrt(1 / lambda0_a / t + k)
      var_j <- sigma1_j^2 / (r * Nj / (1 + r)) + sigma0_j^2 / (Nj / (1 + r))
      var_a <- sigma1_a^2 / (r * N / (1 + r)) + sigma0_a^2 / (N / (1 + r))
      sej <- sqrt(var_j + var_a - 2 * sqrt(f) * sqrt(var_j * var_a))
      uj1 <- (delta_j - delta_a + pi * cut) / sej
      uj2 <- (delta_j - delta_a - pi * cut) / sej
      se <- sqrt(var_a)
      u1 <- (delta_a + cut) / se
      u2 <- (delta_a - cut) / se
      cov <- sqrt(f) * sqrt(var_a * var_j) - var_a
      corr <- cov / (sej * se)
      M <- matrix(c(1, 1, corr, corr, 1, 1, corr, corr, corr, corr, 1, 1, corr, corr, 1, 1), nrow = 4, byrow = T)
      pwr1 <- pmvnorm(lower = c(qnorm(1 - alpha), -Inf, -Inf, -Inf), upper = c(Inf, -qnorm(1 - alpha), Inf, Inf), mean = c(u1, u2, uj1, uj2), corr = M)
      pwr2 <- pmvnorm(lower = c(-Inf, -Inf, 0, -Inf), upper = c(Inf, Inf, Inf, 0), mean = c(u1, u2, uj1, uj2), corr = M)
      pwr3 <- pmvnorm(lower = c(qnorm(1 - alpha), -Inf, 0, -Inf), upper = c(Inf, -qnorm(1 - alpha), Inf, 0), mean = c(u1, u2, uj1, uj2), corr = M)
      pwr4 <- pwr3 / pwr1
      data.frame(delta_a, delta_j, delta_nj, lambda0_a, lambda0_j, lambda0_nj, t, k, f, cut, pi, alpha, beta, N, r, pwr1, pwr2, pwr3, pwr4)
    }, .options = furrr_options(seed = TRUE))
  }
  if (sim) {
    if (isNA_delta_nj & (!isNA_delta_a)) {
      warning("When given global and target region parameters, the data simulation results may be incorrect; please refer to the theoretical calculation results instead.")
    }
    eg <- as.data.frame(expand.grid(delta_j = delta_j, delta_nj = delta_nj, delta_a = delta_a, lambda0_j = lambda0_j, lambda0_nj = lambda0_nj, lambda0_a = lambda0_a, t = t, k = k, f = f, pi = pi, cut = cut, alpha = alpha, beta = beta, N = N, r = r, nsim = 1:nsim, stringsAsFactors = FALSE))
    if (numcore >= 2) {
      plan(multisession, workers = numcore)
    }
    ss <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      delta_j <- R$delta_j
      delta_nj <- R$delta_nj
      delta_a <- R$delta_a
      lambda0_j <- R$lambda0_j
      lambda0_nj <- R$lambda0_nj
      lambda0_a <- R$lambda0_a
      t <- R$t
      k <- R$k
      f <- R$f
      pi <- R$pi
      cut <- R$cut
      alpha <- R$alpha
      beta <- R$beta
      N <- R$N
      r <- R$r
      if (isNA_delta_nj & (!isNA_delta_a)) {
        delta_nj <- (delta_a - delta_j * f) / (1 - f)
      }
      if (isNA_delta_a & (!isNA_delta_nj)) {
        delta_a <- delta_j * f + delta_nj * (1 - f)
      }
      if (isNA_lambda0_nj & (!isNA_lambda0_a)) {
        lambda0_nj <- (lambda0_a - lambda0_j * f) / (1 - f)
      }
      if (isNA_lambda0_a & (!isNA_lambda0_nj)) {
        lambda0_a <- lambda0_j * f + lambda0_nj * (1 - f)
      }
      if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
        N <- getN_Count_Equi(delta = delta_a, lambda0 = lambda0_a, t = t, k = k, cut = cut, alpha = alpha, beta = beta, N = NA, r = r)$N
      }
      Nj <- N * f
      set.seed(i + seed)
      if (k == 0) {
        xt_j <- rpois(n = Nj * r / (1 + r), lambda = lambda0_j * exp(delta_j) * t)
        xc_j <- rpois(n = Nj / (1 + r), lambda = lambda0_j * t)
        xt_nj <- rpois(n = (N - Nj) * r / (1 + r), lambda = lambda0_nj * exp(delta_nj) * t)
        xc_nj <- rpois(n = (N - Nj) / (1 + r), lambda = lambda0_nj * t)
      }
      if (k > 0) {
        xt_j <- rnbinom(n = Nj * r / (1 + r), size = k, mu = lambda0_j * exp(delta_j) * t)
        xc_j <- rnbinom(n = Nj / (1 + r), size = k, mu = lambda0_j * t)
        xt_nj <- rnbinom(n = (N - Nj) * r / (1 + r), size = k, mu = lambda0_nj * exp(delta_nj) * t)
        xc_nj <- rnbinom(n = (N - Nj) / (1 + r), size = k, mu = lambda0_nj * t)
      }
      xt <- c(xt_j, xt_nj)
      xc <- c(xc_j, xc_nj)
      dat_a <- data.frame(x = c(xt, xc), trt = c(rep(1, length(xt)), rep(0, length(xc))))
      dat_j <- data.frame(x = c(xt_j, xc_j), trt = c(rep(1, length(xt_j)), rep(0, length(xc_j))))
      if (k == 0) {
        fit_j <- glm(x ~ trt, dat = dat_j, family = poisson(link = "log"))
        fit_a <- glm(x ~ trt, dat = dat_a, family = poisson(link = "log"))
      }
      if (k > 0) {
        fit_j <- glm.nb(x ~ trt, dat = dat_j)
        fit_a <- glm.nb(x ~ trt, dat = dat_a)
      }
      coef_j <- coef(fit_j)[2]
      coef_a <- coef(fit_a)[2]
      za1 <- summary(fit_a)$coefficients[2, 3] + cut / summary(fit_a)$coefficients[2, 2]
      za2 <- summary(fit_a)$coefficients[2, 3] - cut / summary(fit_a)$coefficients[2, 2]
      zj1 <- coef_j - coef_a + pi * cut
      zj2 <- coef_j - coef_a - pi * cut
      succ_a <- if_else(za1 > qnorm(1 - alpha) & za2 < -qnorm(1 - alpha), 1, 0)
      succ_j <- if_else(zj1 > 0 & zj2 < 0, 1, 0)
      data.frame(delta_a, delta_j, delta_nj, lambda0_a, lambda0_j, lambda0_nj, t, k, f, pi, cut, alpha, beta, N, r, succ_a, succ_j)
    }, .progress = TRUE, .options = furrr_options(seed = TRUE))
    if (numcore >= 2) {
      plan(sequential)
    }
    res <- suppressMessages({
      ss %>%
        group_by(delta_a, delta_j, delta_nj, lambda0_a, lambda0_j, lambda0_nj, t, k, f, pi, cut, alpha, beta, N, r) %>%
        summarise(pwr1 = mean(succ_a, na.rm = TRUE), pwr2 = mean(succ_j, na.rm = TRUE), pwr3 = mean(succ_a & succ_j, na.rm = TRUE), pwr4 = mean(succ_j[succ_a == 1], na.rm = TRUE)) %>%
        arrange(f) %>%
        as.data.frame()
    })
  }
  return(res)
}
