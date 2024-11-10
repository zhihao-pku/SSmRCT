#' Power of mRCT using Japan's Method 1 for binary endpoints
#'
#' Based on Japan's Method 1, given the global and target region sample sizes, calculate and simulate the marginal probabilities, conditional probabilities, and joint probabilities of global success and efficacy consistency between target region and globally, in clinical trials using superiority, non-inferiority, and equivalence designs with binary endpoints.
#'
#' @rdname getPwr_Bin_JM1
#'
#' @name getPwr_Bin_JM1
#'
#' @param p1_j Rate of treatment group in target region.
#' @param p0_j Rate of control group in target region.
#' @param p1_nj Rate of treatment group in other regions. When \code{p1_nj} is not \code{NA}, \code{p1_a} will be calculated automatically.
#' @param p0_nj Rate of control group in other regions. When \code{p0_nj} is not \code{NA}, \code{p0_a} will be calculated automatically.
#' @param p1_a Rate of treatment group globally.
#' @param p0_a Rate of the control group globally.
#' @param f Proportion of sample size allocated to target region.
#' @param pi Proportion of global efficacy to retain. The default value is 0.5, which means retaining half of the efficacy.
#' @param cut A positive value for non-inferiority or equivalence margin. For RD, the margin is on the original scale. For RR and OR, the margin is on the log scale. For example, if the non-inferiority margins for RD, RR, OR are 0.2, 0.6, and 1.3, then the \code{cut = 0.2}, \code{cut = -log(0.6)}, and \code{cut = log(1.3)}, respectively.
#' @param alpha One-sided type I error rate for global success. The default value is 0.025.
#' @param beta Type II error rate for global success, which is used to calculate global sample size only when \code{N} is \code{NA}.
#' @param N Global sample size. When N is \code{NA} and \code{alpha} and \code{beta} are not \code{NA}, \code{N} will be calculated automatically.
#' @param r Ratio of the sample sizes of the treatment group to the control group. The default value is 1.
#' @param scale Optional values are "RD" for rate difference, "RR" for relative risk, and "OR" for odds ratio. Default value is "RD".
#' @param direct If \code{direct = 1}, larger values of RD/RR/OR are preferable. If \code{direct = -1}, smaller values of RD/RR/OR are preferable.
#' @param sim Logical value. When set to \code{FALSE}, theoretical calculation is performed. When set to \code{TRUE}, simulation is used, which is more time-consuming.
#' @param nsim Number of simulations.
#' @param seed Random seed for simulation.
#' @param numcore Number of CPU cores to use during simulation.
#'
#' @return A data frame containing input parameters and returned power.
#' \describe{
#'   \item{\code{pwr1}}{The marginal probability of global success.}
#'   \item{\code{pwr2}}{The marginal probability that the target region efficacy is consistent with the global efficacy.}
#'   \item{\code{pwr3}}{The conditional probability that the target region efficacy is consistent with the global efficacy given global success.}
#'   \item{\code{pwr4}}{The joint probability of global success and the target region efficacy being consistent with the global efficacy.}
#' }
#'
#' @details
#' Taking the larger RD/RR/OR is preferable as an example. The global success criterion and the efficacy consistency criterion between target region and globally
#'
#' in superiority design:
#' \deqn{Z_a = \frac{\hat \delta_a}{\sqrt{Var(\hat \delta_a)}}> \Phi^{-1}(1 - \alpha)}
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
#' Where \eqn{\hat \delta = \hat p_1 - \hat p_0} for RD, \eqn{\hat \delta = log(\frac{\hat p_1}{\hat p_0})} for RR, and  \eqn{\hat \delta = log(\frac{\hat p_1 / (1 - \hat p_1)}{\hat p_0 / (1 - \hat p_0)})}  for OR.  \eqn{\Delta} is the non-inferiority or equivalence margin (\code{cut}).
#'
#' @references
#' 1. Quan H, Li M, Chen J, et al. Assessment of Consistency of Treatment Effects in Multiregional Clinical Trials. Drug Information J. 2010;44(5):617-632. doi:10.1177/009286151004400509
#'
#' 2. Liao JJZ, Yu Z, Li Y. Sample size allocation in multiregional equivalence studies. Pharm Stat. 2018;17(5):570-577. doi:10.1002/pst.1871
#'
#' @export
#'
#' @examples
#' getPwr_Bin_Super_JM1(
#'   p1_j = 0.7, p0_j = 0.5, p1_a = 0.75, p0_a = 0.5, f = seq(0.1, 0.9, 0.1),
#'   pi = 0.5, alpha = 0.025, beta = NA, N = 200, r = 1, scale = "RD", sim = FALSE
#' )
#'
#' # p1_a and p0_a will be calculated based on p1_j and p1_nj, p0_j and p0_nj, respectively.
#' # Global sample size will be calculated based on beta.
#' getPwr_Bin_Noninf_JM1(
#'   p1_j = 0.6, p0_j = 0.5, p1_nj = 0.5, p0_nj = 0.5, f = seq(0.1, 0.9, 0.1),
#'   pi = 0.5, cut = log(1.4),
#'   alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RR", direct = -1,
#'   sim = FALSE
#' )
getPwr_Bin_Super_JM1 <- function(p1_j, p0_j, p1_nj = NA, p0_nj = NA, p1_a = NA, p0_a = NA, f, pi = 0.5, alpha = 0.025, beta = NA, N = NA, r = 1, scale = "RD", sim = FALSE, nsim = 1000, seed = 0, numcore = 4) {
  if (!sim) {
    eg <- as.data.frame(expand.grid(p1_j = p1_j, p0_j = p0_j, p1_nj = p1_nj, p0_nj = p0_nj, p1_a = p1_a, p0_a = p0_a, f = f, pi = pi, alpha = alpha, beta = beta, N = N, r = r, scale = scale, stringsAsFactors = FALSE))
    res <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      p1_j <- R$p1_j
      p0_j <- R$p0_j
      p1_nj <- R$p1_nj
      p0_nj <- R$p0_nj
      p1_a <- R$p1_a
      p0_a <- R$p0_a
      f <- R$f
      pi <- R$pi
      alpha <- R$alpha
      beta <- R$beta
      N <- R$N
      r <- R$r
      scale <- R$scale
      if (is.na(p1_nj) & (!is.na(p1_a))) {
        p1_nj <- (p1_a - p1_j * f) / (1 - f)
        p0_nj <- (p0_a - p0_j * f) / (1 - f)
      }
      if (is.na(p1_a) & (!is.na(p1_nj))) {
        p1_a <- p1_j * f + p1_nj * (1 - f)
        p0_a <- p0_j * f + p0_nj * (1 - f)
      }
      if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
        N <- getN_Bin_Super(p1 = p1_a, p0 = p0_a, alpha = alpha, beta = beta, N = NA, r = r, scale = scale)$N
      }
      Nj <- N * f
      if (scale == "RD") {
        delta_a <- p1_a - p0_a
        delta_j <- p1_j - p0_j
        var_j <- p1_j * (1 - p1_j) / (r * Nj / (1 + r)) + p0_j * (1 - p0_j) / (Nj / (1 + r))
        var_a <- p1_a * (1 - p1_a) / (r * N / (1 + r)) + p0_a * (1 - p0_a) / (N / (1 + r))
      }
      if (scale == "RR") {
        delta_a <- log(p1_a / p0_a)
        delta_j <- log(p1_j / p0_j)
        var_j <- (1 - p1_j) / p1_j / (r * Nj / (1 + r)) + (1 - p0_j) / p0_j / (Nj / (1 + r))
        var_a <- (1 - p1_a) / p1_a / (r * N / (1 + r)) + (1 - p0_a) / p0_a / (N / (1 + r))
      }
      if (scale == "OR") {
        delta_a <- log(p1_a / (1 - p1_a) / p0_a * (1 - p0_a))
        delta_j <- log(p1_j / (1 - p1_j) / p0_j * (1 - p0_j))
        var_j <- 1 / (1 - p1_j) / p1_j / (r * Nj / (1 + r)) + 1 / (1 - p0_j) / p0_j / (Nj / (1 + r))
        var_a <- 1 / (1 - p1_a) / p1_a / (r * N / (1 + r)) + 1 / (1 - p0_a) / p0_a / (N / (1 + r))
      }
      sej <- sqrt(var_j + pi^2 * var_a - 2 * pi * sqrt(f) * sqrt(var_j * var_a))
      uj <- (delta_j - pi * delta_a) / sej
      se <- sqrt(var_a)
      u <- delta_a / se
      cov <- sqrt(f) * sqrt(var_a * var_j) - pi * var_a
      corr <- cov / (sej * se)
      M <- matrix(c(1, corr, corr, 1), nrow = 2, byrow = T)
      if (delta_a < 0) {
        u <- (-1) * u
        uj <- (-1) * uj
      }
      pwr1 <- pmvnorm(lower = c(qnorm(1 - alpha), -Inf), upper = c(Inf, Inf), mean = c(u, uj), corr = M)
      pwr2 <- pmvnorm(lower = c(-Inf, 0), upper = c(Inf, Inf), mean = c(u, uj), corr = M)
      pwr3 <- pmvnorm(lower = c(qnorm(1 - alpha), 0), upper = c(Inf, Inf), mean = c(u, uj), corr = M)
      pwr4 <- pwr3 / pwr1
      data.frame(p1_a, p0_a, p1_j, p0_j, p1_nj, p0_nj, delta_a, delta_j, scale, f, pi, alpha, beta, N, r, pwr1, pwr2, pwr3, pwr4)
    }, .options = furrr_options(seed = TRUE))
  }
  if (sim) {
    eg <- as.data.frame(expand.grid(p1_j = p1_j, p0_j = p0_j, p1_nj = p1_nj, p0_nj = p0_nj, p1_a = p1_a, p0_a = p0_a, f = f, pi = pi, alpha = alpha, beta = beta, N = N, r = r, scale = scale, nsim = 1:nsim, stringsAsFactors = FALSE))
    if (numcore >= 2) {
      plan(multisession, workers = numcore)
    }
    ss <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      p1_j <- R$p1_j
      p0_j <- R$p0_j
      p1_nj <- R$p1_nj
      p0_nj <- R$p0_nj
      p1_a <- R$p1_a
      p0_a <- R$p0_a
      f <- R$f
      pi <- R$pi
      alpha <- R$alpha
      beta <- R$beta
      N <- R$N
      r <- R$r
      scale <- R$scale
      if (is.na(p1_nj) & (!is.na(p1_a))) {
        p1_nj <- (p1_a - p1_j * f) / (1 - f)
        p0_nj <- (p0_a - p0_j * f) / (1 - f)
      }
      if (is.na(p1_a) & (!is.na(p1_nj))) {
        p1_a <- p1_j * f + p1_nj * (1 - f)
        p0_a <- p0_j * f + p0_nj * (1 - f)
      }
      if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
        N <- getN_Bin_Super(p1 = p1_a, p0 = p0_a, alpha = alpha, beta = beta, N = NA, r = r, scale = scale)$N
      }
      Nj <- N * f
      set.seed(i + seed)
      xt_j <- rbinom(n = Nj * r / (1 + r), size = 1, prob = p1_j)
      xc_j <- rbinom(n = Nj / (1 + r), size = 1, prob = p0_j)
      xt_nj <- rbinom(n = (N - Nj) * r / (1 + r), size = 1, prob = p1_nj)
      xc_nj <- rbinom(n = (N - Nj) / (1 + r), size = 1, prob = p0_nj)
      xt <- c(xt_j, xt_nj)
      xc <- c(xc_j, xc_nj)
      p1_a_est <- mean(xt)
      p0_a_est <- mean(xc)
      p1_j_est <- mean(xt_j)
      p0_j_est <- mean(xc_j)
      if (scale == "RD") {
        delta_a <- p1_a - p0_a
        delta_j <- p1_j - p0_j
        delta_a_est <- p1_a_est - p0_a_est
        delta_j_est <- p1_j_est - p0_j_est
        var_a <- p1_a_est * (1 - p1_a_est) / (r * N / (1 + r)) + p0_a_est * (1 - p0_a_est) / (N / (1 + r))
      }
      if (scale == "RR") {
        delta_a <- log(p1_a / p0_a)
        delta_j <- log(p1_j / p0_j)
        delta_a_est <- log(p1_a_est / p0_a_est)
        delta_j_est <- log(p1_j_est / p0_j_est)
        var_a <- (1 - p1_a_est) / p1_a_est / (r * N / (1 + r)) + (1 - p0_a_est) / p0_a_est / (N / (1 + r))
      }
      if (scale == "OR") {
        delta_a <- log(p1_a / (1 - p1_a) / p0_a * (1 - p0_a))
        delta_j <- log(p1_j / (1 - p1_j) / p0_j * (1 - p0_j))
        delta_a_est <- log(p1_a_est / (1 - p1_a_est) / p0_a_est * (1 - p0_a_est))
        delta_j_est <- log(p1_j_est / (1 - p1_j_est) / p0_j_est * (1 - p0_j_est))
        var_a <- 1 / (1 - p1_a_est) / p1_a_est / (r * N / (1 + r)) + 1 / (1 - p0_a_est) / p0_a_est / (N / (1 + r))
      }
      za <- delta_a_est / sqrt(var_a)
      zj <- delta_j_est - pi * delta_a_est
      if (delta_a < 0) {
        za <- (-1) * za
        zj <- (-1) * zj
      }
      succ_a <- if_else(za > qnorm(1 - alpha), 1, 0)
      succ_j <- if_else(zj > 0, 1, 0)
      data.frame(p1_a, p0_a, p1_j, p0_j, p1_nj, p0_nj, delta_a, delta_j, scale, f, pi, alpha, beta, N, r, succ_a, succ_j)
    }, .progress = TRUE, .options = furrr_options(seed = TRUE))
    if (numcore >= 2) {
      plan(sequential)
    }
    res <- suppressMessages({
      ss %>%
        group_by(p1_a, p0_a, p1_j, p0_j, p1_nj, p0_nj, delta_a, delta_j, scale, f, pi, alpha, beta, N, r) %>%
        summarise(pwr1 = mean(succ_a, na.rm = TRUE), pwr2 = mean(succ_j, na.rm = TRUE), pwr3 = mean(succ_a & succ_j, na.rm = TRUE), pwr4 = mean(succ_j[succ_a == 1], na.rm = TRUE)) %>%
        arrange(f) %>%
        as.data.frame()
    })
  }
  return(res)
}

#' @rdname getPwr_Bin_JM1
#' @export
getPwr_Bin_Noninf_JM1 <- function(p1_j, p0_j, p1_nj = NA, p0_nj = NA, p1_a = NA, p0_a = NA, f, pi = 0.5, cut, alpha = 0.025, beta = NA, N = NA, r = 1, scale = "RD", direct = 1, sim = FALSE, nsim = 1000, seed = 0, numcore = 4) {
  if (!sim) {
    eg <- as.data.frame(expand.grid(p1_j = p1_j, p0_j = p0_j, p1_nj = p1_nj, p0_nj = p0_nj, p1_a = p1_a, p0_a = p0_a, f = f, pi = pi, cut = cut, alpha = alpha, beta = beta, N = N, r = r, scale = scale, stringsAsFactors = FALSE))
    res <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      p1_j <- R$p1_j
      p0_j <- R$p0_j
      p1_nj <- R$p1_nj
      p0_nj <- R$p0_nj
      p1_a <- R$p1_a
      p0_a <- R$p0_a
      f <- R$f
      pi <- R$pi
      cut <- R$cut
      alpha <- R$alpha
      beta <- R$beta
      N <- R$N
      r <- R$r
      scale <- R$scale
      if (is.na(p1_nj) & (!is.na(p1_a))) {
        p1_nj <- (p1_a - p1_j * f) / (1 - f)
        p0_nj <- (p0_a - p0_j * f) / (1 - f)
      }
      if (is.na(p1_a) & (!is.na(p1_nj))) {
        p1_a <- p1_j * f + p1_nj * (1 - f)
        p0_a <- p0_j * f + p0_nj * (1 - f)
      }
      if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
        N <- getN_Bin_Noninf(p1 = p1_a, p0 = p0_a, cut = cut, alpha = alpha, beta = beta, N = NA, r = r, scale = scale, direct = direct)$N
      }
      Nj <- N * f
      if (scale == "RD") {
        delta_a <- p1_a - p0_a
        delta_j <- p1_j - p0_j
        var_j <- p1_j * (1 - p1_j) / (r * Nj / (1 + r)) + p0_j * (1 - p0_j) / (Nj / (1 + r))
        var_a <- p1_a * (1 - p1_a) / (r * N / (1 + r)) + p0_a * (1 - p0_a) / (N / (1 + r))
      }
      if (scale == "RR") {
        delta_a <- log(p1_a / p0_a)
        delta_j <- log(p1_j / p0_j)
        var_j <- (1 - p1_j) / p1_j / (r * Nj / (1 + r)) + (1 - p0_j) / p0_j / (Nj / (1 + r))
        var_a <- (1 - p1_a) / p1_a / (r * N / (1 + r)) + (1 - p0_a) / p0_a / (N / (1 + r))
      }
      if (scale == "OR") {
        delta_a <- log(p1_a / (1 - p1_a) / p0_a * (1 - p0_a))
        delta_j <- log(p1_j / (1 - p1_j) / p0_j * (1 - p0_j))
        var_j <- 1 / (1 - p1_j) / p1_j / (r * Nj / (1 + r)) + 1 / (1 - p0_j) / p0_j / (Nj / (1 + r))
        var_a <- 1 / (1 - p1_a) / p1_a / (r * N / (1 + r)) + 1 / (1 - p0_a) / p0_a / (N / (1 + r))
      }
      sej <- sqrt(var_j + var_a - 2 * sqrt(f) * sqrt(var_j * var_a))
      uj <- if_else(direct == 1, (delta_j - delta_a + pi * cut) / sej, (delta_j - delta_a - pi * cut) / sej)
      se <- sqrt(var_a)
      u <- if_else(direct == 1, (delta_a + cut) / se, (delta_a - cut) / se)
      cov <- sqrt(f) * sqrt(var_j * var_a) - var_a
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
      data.frame(p1_a, p0_a, p1_j, p0_j, p1_nj, p0_nj, delta_a, delta_j, scale, f, pi, cut, alpha, beta, N, r, direct, pwr1, pwr2, pwr3, pwr4)
    }, .options = furrr_options(seed = TRUE))
  }
  if (sim) {
    eg <- as.data.frame(expand.grid(p1_j = p1_j, p0_j = p0_j, p1_nj = p1_nj, p0_nj = p0_nj, p1_a = p1_a, p0_a = p0_a, f = f, pi = pi, cut = cut, alpha = alpha, beta = beta, N = N, r = r, scale = scale, nsim = 1:nsim, stringsAsFactors = FALSE))
    if (numcore >= 2) {
      plan(multisession, workers = numcore)
    }
    ss <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      p1_j <- R$p1_j
      p0_j <- R$p0_j
      p1_nj <- R$p1_nj
      p0_nj <- R$p0_nj
      p1_a <- R$p1_a
      p0_a <- R$p0_a
      f <- R$f
      pi <- R$pi
      cut <- R$cut
      alpha <- R$alpha
      beta <- R$beta
      N <- R$N
      r <- R$r
      scale <- R$scale
      if (is.na(p1_nj) & (!is.na(p1_a))) {
        p1_nj <- (p1_a - p1_j * f) / (1 - f)
        p0_nj <- (p0_a - p0_j * f) / (1 - f)
      }
      if (is.na(p1_a) & (!is.na(p1_nj))) {
        p1_a <- p1_j * f + p1_nj * (1 - f)
        p0_a <- p0_j * f + p0_nj * (1 - f)
      }
      if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
        N <- getN_Bin_Noninf(p1 = p1_a, p0 = p0_a, cut = cut, alpha = alpha, beta = beta, N = NA, r = r, scale = scale, direct = direct)$N
      }
      Nj <- N * f
      set.seed(i + seed)
      xt_j <- rbinom(n = Nj * r / (1 + r), size = 1, prob = p1_j)
      xc_j <- rbinom(n = Nj / (1 + r), size = 1, prob = p0_j)
      xt_nj <- rbinom(n = (N - Nj) * r / (1 + r), size = 1, prob = p1_nj)
      xc_nj <- rbinom(n = (N - Nj) / (1 + r), size = 1, prob = p0_nj)
      xt <- c(xt_j, xt_nj)
      xc <- c(xc_j, xc_nj)
      p1_a_est <- mean(xt)
      p0_a_est <- mean(xc)
      p1_j_est <- mean(xt_j)
      p0_j_est <- mean(xc_j)
      if (scale == "RD") {
        delta_a <- p1_a - p0_a
        delta_j <- p1_j - p0_j
        delta_a_est <- p1_a_est - p0_a_est
        delta_j_est <- p1_j_est - p0_j_est
        var_a <- p1_a_est * (1 - p1_a_est) / (r * N / (1 + r)) + p0_a_est * (1 - p0_a_est) / (N / (1 + r))
      }
      if (scale == "RR") {
        delta_a <- log(p1_a / p0_a)
        delta_j <- log(p1_j / p0_j)
        delta_a_est <- log(p1_a_est / p0_a_est)
        delta_j_est <- log(p1_j_est / p0_j_est)
        var_a <- (1 - p1_a_est) / p1_a_est / (r * N / (1 + r)) + (1 - p0_a_est) / p0_a_est / (N / (1 + r))
      }
      if (scale == "OR") {
        delta_a <- log(p1_a / (1 - p1_a) / p0_a * (1 - p0_a))
        delta_j <- log(p1_j / (1 - p1_j) / p0_j * (1 - p0_j))
        delta_a_est <- log(p1_a_est / (1 - p1_a_est) / p0_a_est * (1 - p0_a_est))
        delta_j_est <- log(p1_j_est / (1 - p1_j_est) / p0_j_est * (1 - p0_j_est))
        var_a <- 1 / (1 - p1_a_est) / p1_a_est / (r * N / (1 + r)) + 1 / (1 - p0_a_est) / p0_a_est / (N / (1 + r))
      }
      za <- if_else(direct == 1, (delta_a_est + cut) / sqrt(var_a), (delta_a_est - cut) / sqrt(var_a))
      zj <- if_else(direct == 1, delta_j_est - delta_a_est + pi * cut, delta_j_est - delta_a_est - pi * cut)
      if (direct == -1) {
        za <- (-1) * za
        zj <- (-1) * zj
      }
      succ_a <- if_else(za > qnorm(1 - alpha), 1, 0)
      succ_j <- if_else(zj > 0, 1, 0)
      data.frame(p1_a, p0_a, p1_j, p0_j, p1_nj, p0_nj, delta_a, delta_j, scale, f, pi, cut, alpha, beta, N, r, direct, succ_a, succ_j)
    }, .progress = TRUE, .options = furrr_options(seed = TRUE))
    if (numcore >= 2) {
      plan(sequential)
    }
    res <- suppressMessages({
      ss %>%
        group_by(p1_a, p0_a, p1_j, p0_j, p1_nj, p0_nj, delta_a, delta_j, scale, f, pi, cut, alpha, beta, N, r, direct, ) %>%
        summarise(pwr1 = mean(succ_a, na.rm = TRUE), pwr2 = mean(succ_j, na.rm = TRUE), pwr3 = mean(succ_a & succ_j, na.rm = TRUE), pwr4 = mean(succ_j[succ_a == 1], na.rm = TRUE)) %>%
        arrange(f) %>%
        as.data.frame()
    })
  }
  return(res)
}

#' @rdname getPwr_Bin_JM1
#' @export
getPwr_Bin_Equi_JM1 <- function(p1_j, p0_j, p1_nj = NA, p0_nj = NA, p1_a = NA, p0_a = NA, f, pi = 0.5, cut, alpha, beta = NA, N = NA, r = 1, scale = "RD", sim = FALSE, nsim = 1000, seed = 0, numcore = 4) {
  if (!sim) {
    eg <- as.data.frame(expand.grid(p1_j = p1_j, p0_j = p0_j, p1_nj = p1_nj, p0_nj = p0_nj, p1_a = p1_a, p0_a = p0_a, f = f, pi = pi, cut = cut, alpha = alpha, beta = beta, N = N, r = r, scale = scale, stringsAsFactors = FALSE))
    res <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      p1_j <- R$p1_j
      p0_j <- R$p0_j
      p1_nj <- R$p1_nj
      p0_nj <- R$p0_nj
      p1_a <- R$p1_a
      p0_a <- R$p0_a
      f <- R$f
      pi <- R$pi
      cut <- R$cut
      alpha <- R$alpha
      beta <- R$beta
      N <- R$N
      r <- R$r
      scale <- R$scale
      if (is.na(p1_nj) & (!is.na(p1_a))) {
        p1_nj <- (p1_a - p1_j * f) / (1 - f)
        p0_nj <- (p0_a - p0_j * f) / (1 - f)
      }
      if (is.na(p1_a) & (!is.na(p1_nj))) {
        p1_a <- p1_j * f + p1_nj * (1 - f)
        p0_a <- p0_j * f + p0_nj * (1 - f)
      }
      if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
        N <- getN_Bin_Equi(p1 = p1_a, p0 = p0_a, cut = cut, alpha = alpha, beta = beta, N = NA, r = r, scale = scale)$N
      }
      Nj <- N * f
      if (scale == "RD") {
        delta_a <- p1_a - p0_a
        delta_j <- p1_j - p0_j
        var_j <- p1_j * (1 - p1_j) / (r * Nj / (1 + r)) + p0_j * (1 - p0_j) / (Nj / (1 + r))
        var_a <- p1_a * (1 - p1_a) / (r * N / (1 + r)) + p0_a * (1 - p0_a) / (N / (1 + r))
      }
      if (scale == "RR") {
        delta_a <- log(p1_a / p0_a)
        delta_j <- log(p1_j / p0_j)
        var_j <- (1 - p1_j) / p1_j / (r * Nj / (1 + r)) + (1 - p0_j) / p0_j / (Nj / (1 + r))
        var_a <- (1 - p1_a) / p1_a / (r * N / (1 + r)) + (1 - p0_a) / p0_a / (N / (1 + r))
      }
      if (scale == "OR") {
        delta_a <- log(p1_a / (1 - p1_a) / p0_a * (1 - p0_a))
        delta_j <- log(p1_j / (1 - p1_j) / p0_j * (1 - p0_j))
        var_j <- 1 / (1 - p1_j) / p1_j / (r * Nj / (1 + r)) + 1 / (1 - p0_j) / p0_j / (Nj / (1 + r))
        var_a <- 1 / (1 - p1_a) / p1_a / (r * N / (1 + r)) + 1 / (1 - p0_a) / p0_a / (N / (1 + r))
      }
      sej <- sqrt(var_j + var_a - 2 * sqrt(f) * sqrt(var_j * var_a))
      uj1 <- (delta_j - delta_a + pi * cut) / sej
      uj2 <- (delta_j - delta_a - pi * cut) / sej
      se <- sqrt(var_a)
      u1 <- (delta_a + cut) / se
      u2 <- (delta_a - cut) / se
      cov <- sqrt(f) * sqrt(var_j * var_a) - var_a
      corr <- cov / (sej * se)
      M <- matrix(c(1, 1, corr, corr, 1, 1, corr, corr, corr, corr, 1, 1, corr, corr, 1, 1), nrow = 4, byrow = T)
      pwr1 <- pmvnorm(lower = c(qnorm(1 - alpha), -Inf, -Inf, -Inf), upper = c(Inf, -qnorm(1 - alpha), Inf, Inf), mean = c(u1, u2, uj1, uj2), corr = M)
      pwr2 <- pmvnorm(lower = c(-Inf, -Inf, 0, -Inf), upper = c(Inf, Inf, Inf, 0), mean = c(u1, u2, uj1, uj2), corr = M)
      pwr3 <- pmvnorm(lower = c(qnorm(1 - alpha), -Inf, 0, -Inf), upper = c(Inf, -qnorm(1 - alpha), Inf, 0), mean = c(u1, u2, uj1, uj2), corr = M)
      pwr4 <- pwr3 / pwr1
      data.frame(p1_a, p0_a, p1_j, p0_j, p1_nj, p0_nj, delta_a, delta_j, scale, f, pi, cut, alpha, beta, N, r, pwr1, pwr2, pwr3, pwr4)
    }, .options = furrr_options(seed = TRUE))
  }
  if (sim) {
    eg <- as.data.frame(expand.grid(p1_j = p1_j, p0_j = p0_j, p1_nj = p1_nj, p0_nj = p0_nj, p1_a = p1_a, p0_a = p0_a, f = f, pi = pi, cut = cut, alpha = alpha, beta = beta, N = N, r = r, scale = scale, nsim = 1:nsim, stringsAsFactors = FALSE))
    if (numcore >= 2) {
      plan(multisession, workers = numcore)
    }
    ss <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      p1_j <- R$p1_j
      p0_j <- R$p0_j
      p1_nj <- R$p1_nj
      p0_nj <- R$p0_nj
      p1_a <- R$p1_a
      p0_a <- R$p0_a
      f <- R$f
      pi <- R$pi
      cut <- R$cut
      alpha <- R$alpha
      beta <- R$beta
      N <- R$N
      r <- R$r
      scale <- R$scale
      if (is.na(p1_nj) & (!is.na(p1_a))) {
        p1_nj <- (p1_a - p1_j * f) / (1 - f)
        p0_nj <- (p0_a - p0_j * f) / (1 - f)
      }
      if (is.na(p1_a) & (!is.na(p1_nj))) {
        p1_a <- p1_j * f + p1_nj * (1 - f)
        p0_a <- p0_j * f + p0_nj * (1 - f)
      }
      if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
        N <- getN_Bin_Equi(p1 = p1_a, p0 = p0_a, cut = cut, alpha = alpha, beta = beta, N = NA, r = r, scale = scale)$N
      }
      Nj <- N * f
      set.seed(i + seed)
      xt_j <- rbinom(n = Nj * r / (1 + r), size = 1, prob = p1_j)
      xc_j <- rbinom(n = Nj / (1 + r), size = 1, prob = p0_j)
      xt_nj <- rbinom(n = (N - Nj) * r / (1 + r), size = 1, prob = p1_nj)
      xc_nj <- rbinom(n = (N - Nj) / (1 + r), size = 1, prob = p0_nj)
      xt <- c(xt_j, xt_nj)
      xc <- c(xc_j, xc_nj)
      p1_a_est <- mean(xt)
      p0_a_est <- mean(xc)
      p1_j_est <- mean(xt_j)
      p0_j_est <- mean(xc_j)
      if (scale == "RD") {
        delta_a <- p1_a - p0_a
        delta_j <- p1_j - p0_j
        delta_a_est <- p1_a_est - p0_a_est
        delta_j_est <- p1_j_est - p0_j_est
        var_a <- p1_a_est * (1 - p1_a_est) / (r * N / (1 + r)) + p0_a_est * (1 - p0_a_est) / (N / (1 + r))
      }
      if (scale == "RR") {
        delta_a <- log(p1_a / p0_a)
        delta_j <- log(p1_j / p0_j)
        delta_a_est <- log(p1_a_est / p0_a_est)
        delta_j_est <- log(p1_j_est / p0_j_est)
        var_a <- (1 - p1_a_est) / p1_a_est / (r * N / (1 + r)) + (1 - p0_a_est) / p0_a_est / (N / (1 + r))
      }
      if (scale == "OR") {
        delta_a <- log(p1_a / (1 - p1_a) / p0_a * (1 - p0_a))
        delta_j <- log(p1_j / (1 - p1_j) / p0_j * (1 - p0_j))
        delta_a_est <- log(p1_a_est / (1 - p1_a_est) / p0_a_est * (1 - p0_a_est))
        delta_j_est <- log(p1_j_est / (1 - p1_j_est) / p0_j_est * (1 - p0_j_est))
        var_a <- 1 / (1 - p1_a_est) / p1_a_est / (r * N / (1 + r)) + 1 / (1 - p0_a_est) / p0_a_est / (N / (1 + r))
      }
      za1 <- (delta_a_est + cut) / sqrt(var_a)
      za2 <- (delta_a_est - cut) / sqrt(var_a)
      zj1 <- delta_j_est - delta_a_est + pi * cut
      zj2 <- delta_j_est - delta_a_est - pi * cut
      succ_a <- if_else(za1 > qnorm(1 - alpha) & za2 < -qnorm(1 - alpha), 1, 0)
      succ_j <- if_else(zj1 > 0 & zj2 < 0, 1, 0)
      data.frame(p1_a, p0_a, p1_j, p0_j, p1_nj, p0_nj, delta_a, delta_j, scale, f, pi, cut, alpha, beta, N, r, succ_a, succ_j)
    }, .progress = TRUE, .options = furrr_options(seed = TRUE))
    if (numcore >= 2) {
      plan(sequential)
    }
    res <- suppressMessages({
      ss %>%
        group_by(p1_a, p0_a, p1_j, p0_j, p1_nj, p0_nj, delta_a, delta_j, scale, f, pi, cut, alpha, beta, N, r) %>%
        summarise(pwr1 = mean(succ_a, na.rm = TRUE), pwr2 = mean(succ_j, na.rm = TRUE), pwr3 = mean(succ_a & succ_j, na.rm = TRUE), pwr4 = mean(succ_j[succ_a == 1], na.rm = TRUE)) %>%
        arrange(f) %>%
        as.data.frame()
    })
  }
  return(res)
}
