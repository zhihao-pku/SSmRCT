#' Power of mRCT using Japan's Method 1 for continuous endpoints
#'
#' Based on Japan's Method 1, given the global and target region sample sizes, calculate and simulate the marginal probabilities, conditional probabilities, and joint probabilities of global success and efficacy consistency between target region and globally, in clinical trials using superiority, non-inferiority, and  equivalence designs with continuous endpoints.
#'
#' @rdname getPwr_Con_JM1
#'
#' @name getPwr_Con_JM1
#'
#' @param delta_j Mean difference between treatment and control groups in target region.
#' @param delta_nj Mean difference between treatment and control groups in other regions. When \code{delta_nj} is not \code{NA}, \code{delta_a} will be calculated automatically.
#' @param delta_a Mean difference between treatment and control groups globally.
#' @param sigma Common standard deviation.
#' @param f Proportion of sample size allocated to target region.
#' @param pi Proportion of global efficacy to retain. The default value is 0.5, which means retaining half of the efficacy.
#' @param cut A positive value for non-inferiority or equivalence margin.
#' @param alpha One-sided type I error rate for global success. The default value is 0.025.
#' @param beta Type II error rate for global success, which is used to calculate global sample size only when \code{N} is \code{NA}.
#' @param N Global sample size. When \code{N} is \code{NA} and \code{beta} is not \code{NA}, \code{N} will be calculated automatically.
#' @param r Ratio of the sample sizes of the treatment group to the control group. The default value is 1.
#' @param direct \code{direct = 1} indicates that a larger mean difference is preferable, while \code{direct = -1} indicates that a smaller mean difference is preferable.
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
#' Taking the larger mean difference is preferable as an example. The global success criterion and the efficacy consistency criterion between target region and globally
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
#' Where \eqn{\hat \delta} is the mean difference between treatment and control groups, and \eqn{\Delta} is the non-inferiority or equivalence margin.
#'
#' @references
#' 1. Quan H, Li M, Chen J, et al. Assessment of Consistency of Treatment Effects in Multiregional Clinical Trials. Drug Information J. 2010;44(5):617-632. doi:10.1177/009286151004400509
#'
#' 2. Liao JJZ, Yu Z, Li Y. Sample size allocation in multiregional equivalence studies. Pharm Stat. 2018;17(5):570-577. doi:10.1002/pst.1871
#'
#' @export
#'
#' @examples
#' getPwr_Con_Super_JM1(
#'   delta_j = 0.5, delta_a = 0.7, sigma = 1, f = seq(0.1, 0.9, 0.1),
#'   pi = 0.5, alpha = 0.025, beta = NA, N = 100, r = 1, sim = FALSE
#' )
#'
#' # delta_a will be calculated based on delta_j and delta_nj.
#' # Global sample size will be calculated based on beta.
#' getPwr_Con_Noninf_JM1(
#'   delta_j = 0.2, delta_nj = 0.1, sigma = 1, f = seq(0.1, 0.9, 0.1),
#'   pi = 0.5, cut = 0.4, alpha = 0.025, beta = 0.2, N = NA, r = 1, direct = -1,
#'   sim = FALSE
#' )
getPwr_Con_Super_JM1 <- function(delta_j, delta_nj = NA, delta_a = NA, sigma, f, pi = 0.5, alpha = 0.025, beta = NA, N = NA, r = 1, sim = FALSE, nsim = 1000, seed = 0, numcore = 4) {
  if (!sim) {
    eg <- as.data.frame(expand.grid(delta_j = delta_j, delta_nj = delta_nj, delta_a = delta_a, sigma = sigma, f = f, pi = pi, alpha = alpha, beta = beta, N = N, r = r, stringsAsFactors = FALSE))
    res <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      delta_j <- R$delta_j
      delta_nj <- R$delta_nj
      delta_a <- R$delta_a
      sigma <- R$sigma
      f <- R$f
      pi <- R$pi
      alpha <- R$alpha
      beta <- R$beta
      N <- R$N
      r <- R$r
      gr <- 2 + r + 1 / r
      if (is.na(delta_nj) & (!is.na(delta_a))) {
        delta_nj <- (delta_a - delta_j * f) / (1 - f)
      }
      if (is.na(delta_a) & (!is.na(delta_nj))) {
        delta_a <- delta_j * f + delta_nj * (1 - f)
      }
      if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
        N <- getN_Con_Super(delta = delta_a, sigma = sigma, alpha = alpha, beta = beta, N = NA, r = r)$N
      }
      Nj <- N * f
      sej <- sqrt(gr * sigma^2 / Nj + pi^2 * gr * sigma^2 / N - 2 * pi * sqrt(f) * sqrt(gr * sigma^2 / Nj * gr * sigma^2 / N))
      uj <- (delta_j - pi * delta_a) / sej
      se <- sqrt(gr * sigma^2 / N)
      u <- delta_a / se
      cov <- sqrt(f) * sqrt(gr * sigma^2 / N * gr * sigma^2 / Nj) - pi * gr * sigma^2 / N
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
      data.frame(delta_a, delta_j, delta_nj, sigma, f, pi, alpha, beta, N, r, pwr1, pwr2, pwr3, pwr4)
    }, .options = furrr_options(seed = TRUE))
  }
  if (sim) {
    eg <- as.data.frame(expand.grid(delta_j = delta_j, delta_nj = delta_nj, delta_a = delta_a, sigma = sigma, f = f, pi = pi, alpha = alpha, beta = beta, N = N, r = r, nsim = 1:nsim, stringsAsFactors = FALSE))
    if (numcore >= 2) {
      plan(multisession, workers = numcore)
    }
    ss <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      delta_j <- R$delta_j
      delta_nj <- R$delta_nj
      delta_a <- R$delta_a
      sigma <- R$sigma
      f <- R$f
      pi <- R$pi
      alpha <- R$alpha
      beta <- R$beta
      N <- R$N
      r <- R$r
      if (is.na(delta_nj) & (!is.na(delta_a))) {
        delta_nj <- (delta_a - delta_j * f) / (1 - f)
      }
      if (is.na(delta_a) & (!is.na(delta_nj))) {
        delta_a <- delta_j * f + delta_nj * (1 - f)
      }
      if (is.na(N) & (!is.na(beta))) {
        N <- getN_Con_Super(delta = delta_a, sigma = sigma, alpha = alpha, beta = beta, N = NA, r = r)$N
      }
      Nj <- N * f
      set.seed(i + seed)
      xt_j <- rnorm(n = Nj * r / (1 + r), mean = delta_j, sd = sigma)
      xc_j <- rnorm(n = Nj / (1 + r), mean = 0, sd = sigma)
      xt_nj <- rnorm(n = (N - Nj) * r / (1 + r), mean = delta_nj, sd = sigma)
      xc_nj <- rnorm(n = (N - Nj) / (1 + r), mean = 0, sd = sigma)
      xt <- c(xt_j, xt_nj)
      xc <- c(xc_j, xc_nj)
      za <- (mean(xt) - mean(xc)) / (sigma * sqrt(1 / length(xt) + 1 / length(xc)))
      zj <- mean(xt_j) - mean(xc_j) - pi * (mean(xt) - mean(xc))
      if (delta_a < 0) {
        za <- (-1) * za
        zj <- (-1) * zj
      }
      succ_a <- if_else(za > qnorm(1 - alpha), 1, 0)
      succ_j <- if_else(zj > 0, 1, 0)
      data.frame(delta_a, delta_j, delta_nj, sigma, f, pi, alpha, beta, N, r, succ_a, succ_j)
    }, .progress = TRUE, .options = furrr_options(seed = TRUE))
    if (numcore >= 2) {
      plan(sequential)
    }
    res <- suppressMessages({
      ss %>%
        group_by(delta_a, delta_j, delta_nj, sigma, f, pi, alpha, beta, N, r) %>%
        summarise(pwr1 = mean(succ_a, na.rm = TRUE), pwr2 = mean(succ_j, na.rm = TRUE), pwr3 = mean(succ_a & succ_j, na.rm = TRUE), pwr4 = mean(succ_j[succ_a == 1], na.rm = TRUE)) %>%
        arrange(f) %>%
        as.data.frame()
    })
  }
  return(res)
}

#' @rdname getPwr_Con_JM1
#' @export
getPwr_Con_Noninf_JM1 <- function(delta_j, delta_nj = NA, delta_a = NA, sigma, f, pi = 0.5, cut, alpha = 0.025, beta = NA, N = NA, r = 1, direct = 1, sim = FALSE, nsim = 1000, seed = 0, numcore = 4) {
  if (!sim) {
    eg <- as.data.frame(expand.grid(delta_j = delta_j, delta_nj = delta_nj, delta_a = delta_a, sigma = sigma, f = f, pi = pi, cut = cut, alpha = alpha, beta = beta, N = N, r = r, stringsAsFactors = FALSE))
    res <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      delta_j <- R$delta_j
      delta_nj <- R$delta_nj
      delta_a <- R$delta_a
      sigma <- R$sigma
      f <- R$f
      pi <- R$pi
      cut <- R$cut
      alpha <- R$alpha
      beta <- R$beta
      N <- R$N
      r <- R$r
      gr <- 2 + r + 1 / r
      if (is.na(delta_nj) & (!is.na(delta_a))) {
        delta_nj <- (delta_a - delta_j * f) / (1 - f)
      }
      if (is.na(delta_a) & (!is.na(delta_nj))) {
        delta_a <- delta_j * f + delta_nj * (1 - f)
      }
      if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
        N <- getN_Con_Noninf(delta = delta_a, sigma = sigma, cut = cut, alpha = alpha, beta = beta, N = NA, r = r)$N
      }
      Nj <- N * f
      se <- sqrt(gr * sigma^2 / N)
      u <- if_else(direct == 1, (delta_a + cut) / se, (delta_a - cut) / se)
      sej <- sqrt(gr * sigma^2 / Nj + gr * sigma^2 / N - 2 * sqrt(f) * sqrt(gr * sigma^2 / Nj * gr * sigma^2 / N))
      uj <- if_else(direct == 1, (delta_j - delta_a + pi * cut) / sej, (delta_j - delta_a - pi * cut) / sej)
      cov <- sqrt(f) * sqrt(gr * sigma^2 / N * gr * sigma^2 / Nj) - gr * sigma^2 / N
      corr <- cov / (sej * se)
      M <- matrix(c(1, corr, corr, 1), nrow = 2, byrow = T)
      if (direct == -1) {
        u <- (-1) * u
        uj <- (-1) * uj
      }
      pwr1 <- pmvnorm(lower = c(qnorm(1 - alpha), -Inf), upper = c(Inf, Inf), mean = c(u, uj), corr = M)
      pwr2 <- pmvnorm(lower = c(-Inf, 0), upper = c(Inf, Inf), mean = c(u, uj), corr = M)
      pwr3 <- pmvnorm(lower = c(qnorm(1 - alpha), 0), upper = c(Inf, Inf), mean = c(u, uj), corr = M)
      pwr4 <- pwr3 / pwr1
      data.frame(delta_a, delta_j, delta_nj, sigma, f, pi, cut, alpha, beta, N, r, direct, pwr1, pwr2, pwr3, pwr4)
    }, .options = furrr_options(seed = TRUE))
  }
  if (sim) {
    eg <- as.data.frame(expand.grid(delta_j = delta_j, delta_nj = delta_nj, delta_a = delta_a, sigma = sigma, f = f, pi = pi, cut = cut, alpha = alpha, beta = beta, N = N, r = r, nsim = 1:nsim, stringsAsFactors = FALSE))
    if (numcore >= 2) {
      plan(multisession, workers = numcore)
    }
    ss <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      delta_j <- R$delta_j
      delta_nj <- R$delta_nj
      delta_a <- R$delta_a
      sigma <- R$sigma
      f <- R$f
      pi <- R$pi
      cut <- R$cut
      alpha <- R$alpha
      beta <- R$beta
      N <- R$N
      r <- R$r
      if (is.na(delta_nj) & (!is.na(delta_a))) {
        delta_nj <- (delta_a - delta_j * f) / (1 - f)
      }
      if (is.na(delta_a) & (!is.na(delta_nj))) {
        delta_a <- delta_j * f + delta_nj * (1 - f)
      }
      if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
        N <- getN_Con_Noninf(delta = delta_a, sigma = sigma, cut = cut, alpha = alpha, beta = beta, N = NA, r = r)$N
      }
      Nj <- N * f
      set.seed(i + seed)
      xt_j <- rnorm(n = Nj * r / (1 + r), mean = delta_j, sd = sigma)
      xc_j <- rnorm(n = Nj / (1 + r), mean = 0, sd = sigma)
      xt_nj <- rnorm(n = (N - Nj) * r / (1 + r), mean = delta_nj, sd = sigma)
      xc_nj <- rnorm(n = (N - Nj) / (1 + r), mean = 0, sd = sigma)
      xt <- c(xt_j, xt_nj)
      xc <- c(xc_j, xc_nj)
      za <- if_else(direct == 1, (mean(xt) - mean(xc) + cut) / (sigma * sqrt(1 / length(xt) + 1 / length(xc))), (mean(xt) - mean(xc) - cut) / (sigma * sqrt(1 / length(xt) + 1 / length(xc))))
      zj <- if_else(direct == 1, mean(xt_j) - mean(xc_j) - (mean(xt) - mean(xc)) + pi * cut, mean(xt_j) - mean(xc_j) - (mean(xt) - mean(xc)) - pi * cut)
      if (direct == -1) {
        za <- (-1) * za
        zj <- (-1) * zj
      }
      succ_a <- if_else(za > qnorm(1 - alpha), 1, 0)
      succ_j <- if_else(zj > 0, 1, 0)
      data.frame(delta_a, delta_j, delta_nj, sigma, f, pi, cut, alpha, beta, N, r, direct, succ_a, succ_j)
    }, .progress = TRUE, .options = furrr_options(seed = TRUE))
    if (numcore >= 2) {
      plan(sequential)
    }
    res <- suppressMessages({
      ss %>%
        group_by(delta_a, delta_j, delta_nj, sigma, f, pi, cut, alpha, beta, N, r, direct, ) %>%
        summarise(pwr1 = mean(succ_a, na.rm = TRUE), pwr2 = mean(succ_j, na.rm = TRUE), pwr3 = mean(succ_a & succ_j, na.rm = TRUE), pwr4 = mean(succ_j[succ_a == 1], na.rm = TRUE)) %>%
        arrange(f) %>%
        as.data.frame()
    })
  }
  return(res)
}

#' @rdname getPwr_Con_JM1
#' @export
getPwr_Con_Equi_JM1 <- function(delta_j, delta_nj = NA, delta_a = NA, sigma, f, pi = 0.5, cut, alpha = 0.025, beta = NA, N = NA, r = 1, sim = FALSE, nsim = 1000, seed = 0, numcore = 4) {
  if (!sim) {
    eg <- as.data.frame(expand.grid(delta_j = delta_j, delta_nj = delta_nj, delta_a = delta_a, sigma = sigma, f = f, pi = pi, cut = cut, alpha = alpha, beta = beta, N = N, r = r, stringsAsFactors = FALSE))
    res <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      delta_j <- R$delta_j
      delta_nj <- R$delta_nj
      delta_a <- R$delta_a
      sigma <- R$sigma
      f <- R$f
      pi <- R$pi
      cut <- R$cut
      alpha <- R$alpha
      beta <- R$beta
      N <- R$N
      r <- R$r
      gr <- 2 + r + 1 / r
      if (is.na(delta_nj) & (!is.na(delta_a))) {
        delta_nj <- (delta_a - delta_j * f) / (1 - f)
      }
      if (is.na(delta_a) & (!is.na(delta_nj))) {
        delta_a <- delta_j * f + delta_nj * (1 - f)
      }
      if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
        N <- getN_Con_Equi(delta = delta_a, sigma = sigma, cut = cut, alpha = alpha, beta = beta, N = NA, r = r)$N
      }
      Nj <- N * f
      se <- sqrt(gr * sigma^2 / N)
      u1 <- (delta_a + cut) / se
      u2 <- (delta_a - cut) / se
      sej <- sqrt(gr * sigma^2 / Nj + gr * sigma^2 / N - 2 * sqrt(f) * sqrt(gr * sigma^2 / Nj * gr * sigma^2 / N))
      uj1 <- (delta_j - delta_a + pi * cut) / sej
      uj2 <- (delta_j - delta_a - pi * cut) / sej
      cov <- sqrt(f) * sqrt(gr * sigma^2 / N * gr * sigma^2 / Nj) - gr * sigma^2 / N
      corr <- cov / (sej * se)
      M <- matrix(c(1, 1, corr, corr, 1, 1, corr, corr, corr, corr, 1, 1, corr, corr, 1, 1), nrow = 4, byrow = T)
      pwr1 <- pmvnorm(lower = c(qnorm(1 - alpha), -Inf, -Inf, -Inf), upper = c(Inf, -qnorm(1 - alpha), Inf, Inf), mean = c(u1, u2, uj1, uj2), corr = M)
      pwr2 <- pmvnorm(lower = c(-Inf, -Inf, 0, -Inf), upper = c(Inf, Inf, Inf, 0), mean = c(u1, u2, uj1, uj2), corr = M)
      pwr3 <- pmvnorm(lower = c(qnorm(1 - alpha), -Inf, 0, -Inf), upper = c(Inf, -qnorm(1 - alpha), Inf, 0), mean = c(u1, u2, uj1, uj2), corr = M)
      pwr4 <- pwr3 / pwr1
      data.frame(delta_a, delta_j, delta_nj, sigma, f, pi, cut, alpha, beta, N, r, pwr1, pwr2, pwr3, pwr4)
    }, .options = furrr_options(seed = TRUE))
  }
  if (sim) {
    eg <- as.data.frame(expand.grid(delta_j = delta_j, delta_nj = delta_nj, delta_a = delta_a, sigma = sigma, f = f, pi = pi, cut = cut, alpha = alpha, beta = beta, N = N, r = r, nsim = 1:nsim, stringsAsFactors = FALSE))
    if (numcore >= 2) {
      plan(multisession, workers = numcore)
    }
    ss <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      delta_j <- R$delta_j
      delta_nj <- R$delta_nj
      delta_a <- R$delta_a
      sigma <- R$sigma
      f <- R$f
      pi <- R$pi
      cut <- R$cut
      alpha <- R$alpha
      beta <- R$beta
      N <- R$N
      r <- R$r
      if (is.na(delta_nj) & (!is.na(delta_a))) {
        delta_nj <- (delta_a - delta_j * f) / (1 - f)
      }
      if (is.na(delta_a) & (!is.na(delta_nj))) {
        delta_a <- delta_j * f + delta_nj * (1 - f)
      }
      if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
        N <- getN_Con_Equi(delta = delta_a, sigma = sigma, cut = cut, alpha = alpha, beta = beta, N = NA, r = r)$N
      }
      Nj <- N * f
      set.seed(i + seed)
      xt_j <- rnorm(n = Nj * r / (1 + r), mean = delta_j, sd = sigma)
      xc_j <- rnorm(n = Nj / (1 + r), mean = 0, sd = sigma)
      xt_nj <- rnorm(n = (N - Nj) * r / (1 + r), mean = delta_nj, sd = sigma)
      xc_nj <- rnorm(n = (N - Nj) / (1 + r), mean = 0, sd = sigma)
      xt <- c(xt_j, xt_nj)
      xc <- c(xc_j, xc_nj)
      za1 <- (mean(xt) - mean(xc) + cut) / (sigma * sqrt(1 / length(xt) + 1 / length(xc)))
      za2 <- (mean(xt) - mean(xc) - cut) / (sigma * sqrt(1 / length(xt) + 1 / length(xc)))
      zj1 <- mean(xt_j) - mean(xc_j) - (mean(xt) - mean(xc)) + pi * cut
      zj2 <- mean(xt_j) - mean(xc_j) - (mean(xt) - mean(xc)) - pi * cut
      succ_a <- if_else(za1 > qnorm(1 - alpha) & za2 < -qnorm(1 - alpha), 1, 0)
      succ_j <- if_else(zj1 > 0 & zj2 < 0, 1, 0)
      data.frame(delta_a, delta_j, delta_nj, sigma, f, pi, cut, alpha, beta, N, r, succ_a, succ_j)
    }, .progress = TRUE, .options = furrr_options(seed = TRUE))
    if (numcore >= 2) {
      plan(sequential)
    }
    res <- suppressMessages({
      ss %>%
        group_by(delta_a, delta_j, delta_nj, sigma, f, pi, cut, alpha, beta, N, r, ) %>%
        summarise(pwr1 = mean(succ_a, na.rm = TRUE), pwr2 = mean(succ_j, na.rm = TRUE), pwr3 = mean(succ_a & succ_j, na.rm = TRUE), pwr4 = mean(succ_j[succ_a == 1], na.rm = TRUE)) %>%
        arrange(f) %>%
        as.data.frame()
    })
  }
  return(res)
}
