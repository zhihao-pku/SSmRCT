#' Power of mRCT using Japan's Method 1 for survival endpoints
#'
#' Based on Japan's Method 1, given the global and target region sample sizes, calculate and simulate the marginal probabilities, conditional probabilities, and joint probabilities of global success and efficacy consistency between target region and global, in clinical trials using superiority, non-inferiority, and equivalence designs with survival endpoints.
#'
#' @rdname getPwr_Surv_JM1
#'
#' @name getPwr_Surv_JM1
#'
#' @param delta_j log(HR) between treatment and control groups in target region.
#' @param delta_nj log(HR) between treatment and control groups in other regions. When \code{delta_nj} is not \code{NA}, \code{delta_a} will be calculated automatically.
#' @param delta_a log(HR) between treatment and control groups globally.
#' @param f Proportion of sample size allocated to target region.
#' @param pi Proportion of global efficacy to retain. The default value is 0.5, which means retaining half of the efficacy.
#' @param cut A positive value for non-inferiority or equivalence margin. For example, if the non-inferiority margin for HR is 0.6, then \code{cut = -log(0.6)}. If the non-inferiority margin for HR is 1.3, then \code{cut = log(1.3)}.
#' @param alpha One-sided type I error rate for global success. The default value is 0.025.
#' @param beta Type II error rate for global success, which is used to calculate global sample size only when \code{N} is \code{NA}.
#' @param N Global sample size. When \code{N} is \code{NA} and \code{beta} is not \code{NA}, \code{N} will be calculated automatically.
#' @param r Ratio of the sample sizes of the treatment group to the control group. The default value is 1.
#' @param criterion If \code{criterion = 1}, the consistency criterion defined on the log(HR) scale will be used. If \code{criterion = 2}, the consistency criterion defined on the HR scale will be used. See \code{details} for more information.
#' @param direct \code{direct = 1} indicates that a larger HR is preferable, while \code{direct = -1} indicates that a smaller HR is preferable.
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
#' Taking the larger HR is preferable as an example. The global success criterion and the efficacy consistency criterion between target region and globally
#'
#' \code{criterion = 1} in superiority design:
#' \deqn{Z_a = \frac{\hat \delta_a}{\sqrt{Var(\hat \delta_a)}} > \Phi^{-1}(1 - \alpha)}
#' \deqn{\hat \delta_j - \pi\hat \delta_a > 0}
#'
#' \code{criterion = 2} in superiority design:
#' \deqn{\text{uper limit of 95\% CI for }(1 - e^{\hat \delta_a}) < 0}
#' \deqn{1 - e^{\hat \delta_j} - \pi(1 - e^{\hat \delta_a}) < 0}
#'
#' \code{criterion = 1} in non-inferiority design:
#' \deqn{Z_a = \frac{\hat \delta_a + \Delta}{\sqrt{Var(\hat \delta_a)}} > \Phi^{-1}(1 - \alpha)}
#' \deqn{\hat \delta_j - \hat \delta_a + \pi\Delta > 0}
#'
#' \code{criterion = 2} in non-inferiority design:
#' \deqn{\text{uper limit of 95\% CI for }(\frac{1}{e^{\Delta}} - e^{\hat \delta_a})<0}
#' \deqn{\frac{1}{e^{\Delta}} - e^{\hat \delta_j} - \pi(\frac{1}{e^{\Delta}} - e^{\hat \delta_a}) < 0}
#'
#' \code{criterion = 1} in equivalence design:
#' \deqn{Z_{a_u} = \frac{\hat \delta_a + \Delta}{\sqrt{Var(\hat \delta_a)}} > \Phi^{-1}(1 - \alpha)\text{ and }Z_{a_l} = \frac{\hat \delta_a - \Delta}{\sqrt{Var(\hat \delta_a)}} < \Phi^{-1}(\alpha)}
#' \deqn{\hat \delta_j - \hat \delta_a + \pi\Delta > 0\text{ and }\hat \delta_j - \hat \delta_a - \pi\Delta < 0}
#'
#' \code{criterion = 2} in equivalence design:
#' \deqn{\text{uper limit of 95\% CI for }(\frac{1}{e^{\Delta}} - e^{\hat \delta_a}) < 0\text{ and lower limit of 95\% CI for }e^{\Delta} - e^{\hat \delta_a} > 0}
#' \deqn{\Delta - e^{\hat \delta_j} - \pi(\Delta - e^{\hat \delta_a}) < 0 \text{ and } \frac{1}{\Delta} - e^{\hat \delta_j} - \pi(\frac{1}{\Delta} - e^{\hat \delta_a}) > 0}
#'
#' Where \eqn{\hat \delta = log(\hat {HR})} between treatment and control groups, and \eqn{\Delta} is the non-inferiority or equivalence margin (\code{cut}).
#'
#' For \code{criterion = 2}, by delta method, \eqn{e^{\hat \delta}} approximately follows a distribution of \eqn{N(e^{\delta}, (e^{\delta}\sqrt{\frac{1}{N / (r + 1)} + \frac{1}{Nr / (r + 1)}})^2)}, used in theoretical calculation (\code{sim=FALSE}).
#'
#' @references
#' 1. Quan H, Li M, Chen J, et al. Assessment of Consistency of Treatment Effects in Multiregional Clinical Trials. Drug Information J. 2010;44(5):617-632. doi:10.1177/009286151004400509
#'
#' 2. Liao JJZ, Yu Z, Li Y. Sample size allocation in multiregional equivalence studies. Pharm Stat. 2018;17(5):570-577. doi:10.1002/pst.1871
#'
#' @export
#'
#' @examples
#' getPwr_Surv_Super_JM1(
#'   delta_j = log(1.3),
#'   delta_a = log(1.4),
#'   f = seq(0.1, 0.9, 0.1),
#'   pi = 0.5, alpha = 0.025, beta = NA, N = 200, r = 1, criterion = 1,
#'   sim = FALSE
#' )
#'
#' # delta_a will be calculated based on delta_j and delta_nj.
#' # Global sample size will be calculated based on beta.
#' getPwr_Surv_Noninf_JM1(
#'   delta_j = log(1.1),
#'   delta_nj = log(1.0),
#'   f = seq(0.1, 0.9, 0.1),
#'   cut = log(1.3),
#'   pi = 0.5, alpha = 0.025, beta = 0.2, N = NA, r = 1, criterion = 2,
#'   direct = -1, sim = FALSE
#' )
getPwr_Surv_Super_JM1 <- function(delta_j, delta_nj = NA, delta_a = NA, f, pi = 0.5, alpha = 0.025, beta = NA, N = NA, r = 1, criterion = 1, sim = FALSE, nsim = 1000, seed = 0, numcore = 4) {
  isNA_delta_nj <- is.na(delta_nj)
  isNA_delta_a <- is.na(delta_a)
  if (!sim) {
    eg <- as.data.frame(expand.grid(delta_j = delta_j, delta_nj = delta_nj, delta_a = delta_a, f = f, pi = pi, alpha = alpha, beta = beta, N = N, r = r, criterion = criterion, stringsAsFactors = FALSE))
    res <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      delta_j <- R$delta_j
      delta_nj <- R$delta_nj
      delta_a <- R$delta_a
      f <- R$f
      pi <- R$pi
      alpha <- R$alpha
      beta <- R$beta
      N <- R$N
      r <- R$r
      criterion <- R$criterion
      gr <- 2 + r + 1 / r
      if (isNA_delta_nj & (!isNA_delta_a)) {
        delta_nj <- (delta_a - delta_j * f) / (1 - f)
      }
      if (isNA_delta_a & (!isNA_delta_nj)) {
        delta_a <- delta_j * f + delta_nj * (1 - f)
      }
      if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
        N <- getN_Surv_Super(delta = delta_a, alpha = alpha, beta = beta, N = NA, r = r, criterion = criterion)$N
      }
      Nj <- N * f
      if (criterion == 1) {
        se <- sqrt(gr / N)
        u <- delta_a / se
        sej <- sqrt(gr / Nj + pi^2 * gr / N - 2 * pi * sqrt(f) * sqrt(gr / Nj * gr / N))
        uj <- (delta_j - pi * delta_a) / sej
        cov <- sqrt(f) * sqrt(gr / N * gr / Nj) - pi * gr / N
      }
      if (criterion == 2) {
        se <- sqrt(gr / N * exp(delta_a)^2)
        u <- -(1 - exp(delta_a)) / se
        sej <- sqrt(gr / Nj * exp(delta_j)^2 + pi^2 * gr / N * exp(delta_a)^2 - 2 * pi * gr / N * exp(delta_j) * exp(delta_a))
        uj <- -(1 - exp(delta_j) - pi * (1 - exp(delta_a))) / sej
        cov <- gr / N * exp(delta_j) * exp(delta_a) - pi * gr / N * exp(delta_a)^2
      }
      corr <- cov / (sej * se)
      if (delta_a < 0) {
        u <- (-1) * u
        uj <- (-1) * uj
      }
      M <- matrix(c(1, corr, corr, 1), nrow = 2, byrow = T)
      pwr1 <- pmvnorm(lower = c(qnorm(1 - alpha), -Inf), upper = c(Inf, Inf), mean = c(u, uj), corr = M)
      pwr2 <- pmvnorm(lower = c(-Inf, 0), upper = c(Inf, Inf), mean = c(u, uj), corr = M)
      pwr3 <- pmvnorm(lower = c(qnorm(1 - alpha), 0), upper = c(Inf, Inf), mean = c(u, uj), corr = M)
      pwr4 <- pwr3 / pwr1
      data.frame(delta_a, delta_j, delta_nj, f, pi, alpha, beta, N, r, criterion, pwr1, pwr2, pwr3, pwr4)
    }, .options = furrr_options(seed = TRUE))
  }
  if (sim) {
    eg <- as.data.frame(expand.grid(delta_j = delta_j, delta_nj = delta_nj, delta_a = delta_a, f = f, pi = pi, alpha = alpha, beta = beta, N = N, r = r, criterion = criterion, nsim = 1:nsim, stringsAsFactors = FALSE))
    if (numcore >= 2) {
      plan(multisession, workers = numcore)
    }
    ss <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      delta_j <- R$delta_j
      delta_nj <- R$delta_nj
      delta_a <- R$delta_a
      f <- R$f
      pi <- R$pi
      alpha <- R$alpha
      beta <- R$beta
      N <- R$N
      r <- R$r
      criterion <- R$criterion
      if (isNA_delta_nj & (!isNA_delta_a)) {
        delta_nj <- (delta_a - delta_j * f) / (1 - f)
      }
      if (isNA_delta_a & (!isNA_delta_nj)) {
        delta_a <- delta_j * f + delta_nj * (1 - f)
      }
      if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
        N <- getN_Surv_Super(delta = delta_a, alpha = alpha, beta = beta, N = NA, r = r, criterion = criterion)$N
      }
      Nj <- N * f
      set.seed(i + seed)
      xt_j <- rexp(n = Nj * r / (1 + r), rate = exp(delta_j))
      xc_j <- rexp(n = Nj / (1 + r), rate = 1)
      dat_j <- data.frame(time = c(xt_j, xc_j), status = 1, trt = c(rep(1, length(xt_j)), rep(0, length(xc_j))))
      if (isNA_delta_nj & (!isNA_delta_a)) {
        ut_nj <- runif(n = (N - Nj) * r / (1 + r), min = 0, max = 1)
        uc_nj <- runif(n = (N - Nj) / (1 + r), min = 0, max = 1)
        xc_nj <- -log(1 - uc_nj)
        xt_nj <- c()
        g <- 1
        while (g < (length(ut_nj) + 1)) {
          u <- ut_nj[g]
          solve_t <- function(t) {
            exp(-exp(delta_a) * t) - exp(-exp(delta_j) * t) * f - (1 - u) * (1 - f)
          }
          t <- uniroot(f = solve_t, interval = c(0, 1e+06))$root
          xt_nj <- c(xt_nj, t)
          g <- g + 1
        }
        xt <- c(xt_j, xt_nj)
        xc <- c(xc_j, xc_nj)
      }
      if (isNA_delta_a & (!isNA_delta_nj)) {
        xt_nj <- rexp(n = (N - Nj) * r / (1 + r), rate = exp(delta_nj))
        xc_nj <- rexp(n = (N - Nj) / (1 + r), rate = 1)
        xt <- c(xt_j, xt_nj)
        xc <- c(xc_j, xc_nj)
      }
      dat_a <- data.frame(time = c(xt, xc), status = 1, trt = c(rep(1, length(xt)), rep(0, length(xc))))
      fit_j <- coxph(Surv(time, status) ~ trt, dat = dat_j)
      coef_j <- coef(fit_j)
      fit_a <- coxph(Surv(time, status) ~ trt, dat = dat_a)
      coef_a <- coef(fit_a)
      za <- summary(fit_a)$coefficients[4]
      if (criterion == 1) {
        zj <- coef_j - pi * coef_a
      }
      if (criterion == 2) {
        zj <- -(1 - exp(coef_j) - pi * (1 - exp(coef_a)))
      }
      if (delta_a < 0) {
        za <- (-1) * za
        zj <- (-1) * zj
      }
      succ_a <- if_else(za > qnorm(1 - alpha), 1, 0)
      succ_j <- if_else(zj > 0, 1, 0)
      data.frame(delta_a, delta_j, delta_nj, f, pi, alpha, beta, N, r, criterion, succ_a, succ_j)
    }, .progress = TRUE, .options = furrr_options(seed = TRUE))
    if (numcore >= 2) {
      plan(sequential)
    }
    res <- suppressMessages({
      if (isNA_delta_nj & (!isNA_delta_a)) {
        df <- ss %>%
          group_by(delta_a, delta_j, delta_nj, f, pi, alpha, beta, N, r, criterion) %>%
          summarise(pwr1 = mean(succ_a, na.rm = TRUE), pwr2 = mean(succ_j, na.rm = TRUE), pwr3 = mean(succ_a & succ_j, na.rm = TRUE), pwr4 = mean(succ_j[succ_a == 1], na.rm = TRUE)) %>%
          arrange(f) %>%
          as.data.frame()
      }
      if (isNA_delta_a & (!isNA_delta_nj)) {
        df <- ss %>%
          group_by(delta_a, delta_j, delta_nj, f, pi, alpha, beta, N, r, criterion) %>%
          summarise(pwr1 = mean(succ_a, na.rm = TRUE), pwr2 = mean(succ_j, na.rm = TRUE), pwr3 = mean(succ_a & succ_j, na.rm = TRUE), pwr4 = mean(succ_j[succ_a == 1], na.rm = TRUE)) %>%
          arrange(f) %>%
          as.data.frame()
      }
      df
    })
  }
  return(res)
}

#' @rdname getPwr_Surv_JM1
#' @export
getPwr_Surv_Noninf_JM1 <- function(delta_j, delta_nj = NA, delta_a = NA, f, pi = 0.5, cut, alpha = 0.025, beta = NA, N = NA, r = 1, criterion = 1, direct = 1, sim = FALSE, nsim = 1000, seed = 0, numcore = 4) {
  isNA_delta_nj <- is.na(delta_nj)
  isNA_delta_a <- is.na(delta_a)
  if (!sim) {
    eg <- as.data.frame(expand.grid(delta_j = delta_j, delta_nj = delta_nj, delta_a = delta_a, f = f, pi = pi, cut = cut, alpha = alpha, beta = beta, N = N, r = r, criterion = criterion, stringsAsFactors = FALSE))
    res <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      delta_j <- R$delta_j
      delta_nj <- R$delta_nj
      delta_a <- R$delta_a
      f <- R$f
      pi <- R$pi
      cut <- R$cut
      alpha <- R$alpha
      beta <- R$beta
      N <- R$N
      r <- R$r
      criterion <- R$criterion
      gr <- 2 + r + 1 / r
      if (isNA_delta_nj & (!isNA_delta_a)) {
        delta_nj <- (delta_a - delta_j * f) / (1 - f)
      }
      if (isNA_delta_a & (!isNA_delta_nj)) {
        delta_a <- delta_j * f + delta_nj * (1 - f)
      }
      if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
        N <- getN_Surv_Noninf(delta = delta_a, cut = cut, alpha = alpha, beta = beta, N = NA, r = r, criterion = criterion, direct = direct)$N
      }
      Nj <- N * f
      if (criterion == 1) {
        sej <- sqrt(gr / Nj + gr / N - 2 * sqrt(f) * sqrt(gr / Nj * gr / N))
        uj <- if_else(direct == 1, (delta_j - delta_a + pi * cut) / sej, (delta_j - delta_a - pi * cut) / sej)
        se <- sqrt(gr / N)
        u <- if_else(direct == 1, (delta_a + cut) / se, (delta_a - cut) / se)
        cov <- sqrt(f) * sqrt(gr / N * gr / Nj) - gr / N
      }
      if (criterion == 2) {
        sej <- sqrt(gr / Nj * exp(delta_j)^2 + pi^2 * gr / N * exp(delta_a)^2 - 2 * pi * gr / N * exp(delta_j) * exp(delta_a))
        uj <- if_else(direct == 1, -(1 / exp(cut) - exp(delta_j) - pi * (1 / exp(cut) - exp(delta_a))) / sej, -(exp(cut) - exp(delta_j) - pi * (exp(cut) - exp(delta_a))) / sej)
        se <- sqrt(gr / N * exp(delta_a)^2)
        u <- if_else(direct == 1, -(1 / exp(cut) - exp(delta_a)) / se, -(exp(cut) - exp(delta_a)) / se)
        cov <- gr / N * exp(delta_j) * exp(delta_a) - pi * gr / N * exp(delta_a)^2
      }
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
      data.frame(delta_a, delta_j, delta_nj, f, pi, cut, alpha, beta, N, r, direct, pwr1, pwr2, pwr3, pwr4)
    }, .options = furrr_options(seed = TRUE))
  }
  if (sim) {
    eg <- as.data.frame(expand.grid(delta_j = delta_j, delta_nj = delta_nj, delta_a = delta_a, f = f, pi = pi, cut = cut, alpha = alpha, beta = beta, N = N, r = r, criterion = criterion, nsim = 1:nsim, stringsAsFactors = FALSE))
    if (numcore >= 2) {
      plan(multisession, workers = numcore)
    }
    ss <- future_map_dfr(.x = 1:nsim, .f = function(i) {
      R <- eg[i, ]
      delta_j <- R$delta_j
      delta_nj <- R$delta_nj
      delta_a <- R$delta_a
      f <- R$f
      pi <- R$pi
      cut <- R$cut
      alpha <- R$alpha
      beta <- R$beta
      N <- R$N
      r <- R$r
      criterion <- R$criterion
      if (isNA_delta_nj & (!isNA_delta_a)) {
        delta_nj <- (delta_a - delta_j * f) / (1 - f)
      }
      if (isNA_delta_a & (!isNA_delta_nj)) {
        delta_a <- delta_j * f + delta_nj * (1 - f)
      }
      if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
        N <- getN_Surv_Noninf(delta = delta_a, cut = cut, alpha = alpha, beta = beta, N = NA, r = r, criterion = criterion, direct = direct)$N
      }
      Nj <- N * f
      set.seed(i + seed)
      xt_j <- rexp(n = Nj * r / (1 + r), rate = exp(delta_j))
      xc_j <- rexp(n = Nj / (1 + r), rate = 1)
      dat_j <- data.frame(time = c(xt_j, xc_j), status = 1, trt = c(rep(1, length(xt_j)), rep(0, length(xc_j))))
      if (isNA_delta_nj & (!isNA_delta_a)) {
        ut_nj <- runif(n = (N - Nj) * r / (1 + r), min = 0, max = 1)
        uc_nj <- runif(n = (N - Nj) / (1 + r), min = 0, max = 1)
        xc_nj <- -log(1 - uc_nj)
        xt_nj <- c()
        g <- 1
        while (g < (length(ut_nj) + 1)) {
          u <- ut_nj[g]
          solve_t <- function(t) {
            exp(-exp(delta_a) * t) - exp(-exp(delta_j) * t) * f - (1 - u) * (1 - f)
          }
          t <- uniroot(f = solve_t, interval = c(0, 1e+06))$root
          xt_nj <- c(xt_nj, t)
          g <- g + 1
        }
        xt <- c(xt_j, xt_nj)
        xc <- c(xc_j, xc_nj)
      }
      if (isNA_delta_a & (!isNA_delta_nj)) {
        xt_nj <- rexp(n = (N - Nj) * r / (1 + r), rate = exp(delta_nj))
        xc_nj <- rexp(n = (N - Nj) / (1 + r), rate = 1)
        xt <- c(xt_j, xt_nj)
        xc <- c(xc_j, xc_nj)
      }
      dat_a <- data.frame(time = c(xt, xc), status = 1, trt = c(rep(1, length(xt)), rep(0, length(xc))))
      fit_j <- coxph(Surv(time, status) ~ trt, dat = dat_j)
      coef_j <- coef(fit_j)
      fit_a <- coxph(Surv(time, status) ~ trt, dat = dat_a)
      coef_a <- coef(fit_a)
      za <- if_else(direct == 1, summary(fit_a)$coefficients[4] + cut / summary(fit_a)$coefficients[3], summary(fit_a)$coefficients[4] - cut / summary(fit_a)$coefficients[3])
      if (criterion == 1) {
        zj <- if_else(direct == 1, coef_j - coef_a + pi * cut, coef_j - coef_a - pi * cut)
      }
      if (criterion == 2) {
        zj <- if_else(direct == 1, -(1 / exp(cut) - exp(coef_j) - pi * (1 / exp(cut) - exp(coef_a))), -(exp(cut) - exp(coef_j) - pi * (exp(cut) - exp(coef_a))))
      }
      if (direct == -1) {
        zj <- (-1) * zj
        za <- (-1) * za
      }
      succ_a <- if_else(za > qnorm(1 - alpha), 1, 0)
      succ_j <- if_else(zj > 0, 1, 0)
      data.frame(delta_a, delta_j, delta_nj, f, pi, cut, alpha, beta, N, r, criterion, direct, succ_a, succ_j)
    }, .progress = TRUE, .options = furrr_options(seed = TRUE))
    if (numcore >= 2) {
      plan(sequential)
    }
    res <- suppressMessages({
      if (isNA_delta_nj & (!isNA_delta_a)) {
        df <- ss %>%
          group_by(delta_a, delta_j, delta_nj, f, pi, cut, alpha, beta, N, r, criterion, direct) %>%
          summarise(pwr1 = mean(succ_a, na.rm = TRUE), pwr2 = mean(succ_j, na.rm = TRUE), pwr3 = mean(succ_a & succ_j, na.rm = TRUE), pwr4 = mean(succ_j[succ_a == 1], na.rm = TRUE)) %>%
          arrange(f) %>%
          as.data.frame()
      }
      if (isNA_delta_a & (!isNA_delta_nj)) {
        df <- ss %>%
          group_by(delta_a, delta_j, delta_nj, f, pi, cut, alpha, beta, N, r, criterion, direct) %>%
          summarise(pwr1 = mean(succ_a, na.rm = TRUE), pwr2 = mean(succ_j, na.rm = TRUE), pwr3 = mean(succ_a & succ_j, na.rm = TRUE), pwr4 = mean(succ_j[succ_a == 1], na.rm = TRUE)) %>%
          arrange(f) %>%
          as.data.frame()
      }
      df
    })
  }
  return(res)
}

#' @rdname getPwr_Surv_JM1
#' @export
getPwr_Surv_Equi_JM1 <- function(delta_j, delta_nj = NA, delta_a = NA, f, pi = 0.5, cut, alpha = 0.025, beta = NA, N = NA, r = 1, criterion = 1, sim = FALSE, nsim = 1000, seed = 0, numcore = 4) {
  isNA_delta_nj <- is.na(delta_nj)
  isNA_delta_a <- is.na(delta_a)
  if (!sim) {
    eg <- as.data.frame(expand.grid(delta_j = delta_j, delta_nj = delta_nj, delta_a = delta_a, f = f, pi = pi, cut = cut, alpha = alpha, beta = beta, N = N, r = r, criterion = criterion, stringsAsFactors = FALSE))
    res <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
      R <- eg[i, ]
      delta_j <- R$delta_j
      delta_nj <- R$delta_nj
      delta_a <- R$delta_a
      f <- R$f
      pi <- R$pi
      cut <- R$cut
      alpha <- R$alpha
      beta <- R$beta
      N <- R$N
      r <- R$r
      criterion <- R$criterion
      gr <- 2 + r + 1 / r
      if (isNA_delta_nj & (!isNA_delta_a)) {
        delta_nj <- (delta_a - delta_j * f) / (1 - f)
      }
      if (isNA_delta_a & (!isNA_delta_nj)) {
        delta_a <- delta_j * f + delta_nj * (1 - f)
      }
      if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
        N <- getN_Surv_Equi(delta = delta_a, cut = cut, alpha = alpha, beta = beta, N = NA, r = r, criterion = criterion)$N
      }
      Nj <- N * f
      if (criterion == 1) {
        sej <- sqrt(gr / Nj + gr / N - 2 * sqrt(f) * sqrt(gr / Nj * gr / N))
        uj1 <- (delta_j - delta_a + pi * cut) / sej
        uj2 <- (delta_j - delta_a - pi * cut) / sej
        se <- sqrt(gr / N)
        u1 <- (delta_a + cut) / se
        u2 <- (delta_a - cut) / se
        cov <- sqrt(f) * sqrt(gr / N * gr / Nj) - gr / N
      }
      if (criterion == 2) {
        sej <- sqrt(gr / Nj * exp(delta_j)^2 + pi^2 * gr / N * exp(delta_a)^2 - 2 * pi * gr / N * exp(delta_j) * exp(delta_a))
        uj1 <- -(1 / exp(cut) - exp(delta_j) - pi * (1 / exp(cut) - exp(delta_a))) / sej
        uj2 <- -(exp(cut) - exp(delta_j) - pi * (exp(cut) - exp(delta_a))) / sej
        se <- sqrt(gr / N * exp(delta_a)^2)
        u1 <- -(1 / exp(cut) - exp(delta_a)) / se
        u2 <- -(exp(cut) - exp(delta_a)) / se
        cov <- gr / N * exp(delta_j) * exp(delta_a) - pi * gr / N * exp(delta_a)^2
      }
      corr <- cov / (sej * se)
      M <- matrix(c(1, 1, corr, corr, 1, 1, corr, corr, corr, corr, 1, 1, corr, corr, 1, 1), nrow = 4, byrow = T)
      pwr1 <- pmvnorm(lower = c(qnorm(1 - alpha), -Inf, -Inf, -Inf), upper = c(Inf, -qnorm(1 - alpha), Inf, Inf), mean = c(u1, u2, uj1, uj2), corr = M)
      pwr2 <- pmvnorm(lower = c(-Inf, -Inf, 0, -Inf), upper = c(Inf, Inf, Inf, 0), mean = c(u1, u2, uj1, uj2), corr = M)
      pwr3 <- pmvnorm(lower = c(qnorm(1 - alpha), -Inf, 0, -Inf), upper = c(Inf, -qnorm(1 - alpha), Inf, 0), mean = c(u1, u2, uj1, uj2), corr = M)
      pwr4 <- pwr3 / pwr1
      data.frame(delta_a, delta_j, delta_nj, f, pi, cut, alpha, beta, N, r, pwr1, pwr2, pwr3, pwr4)
    }, .options = furrr_options(seed = TRUE))
  }
  if (sim) {
    eg <- as.data.frame(expand.grid(delta_j = delta_j, delta_nj = delta_nj, delta_a = delta_a, f = f, pi = pi, cut = cut, alpha = alpha, beta = beta, N = N, r = r, criterion = criterion, nsim = 1:nsim, stringsAsFactors = FALSE))
    if (numcore >= 2) {
      plan(multisession, workers = numcore)
    }
    ss <- future_map_dfr(.x = 1:nsim, .f = function(i) {
      R <- eg[i, ]
      delta_j <- R$delta_j
      delta_nj <- R$delta_nj
      delta_a <- R$delta_a
      f <- R$f
      pi <- R$pi
      cut <- R$cut
      alpha <- R$alpha
      beta <- R$beta
      N <- R$N
      r <- R$r
      criterion <- R$criterion
      if (isNA_delta_nj & (!isNA_delta_a)) {
        delta_nj <- (delta_a - delta_j * f) / (1 - f)
      }
      if (isNA_delta_a & (!isNA_delta_nj)) {
        delta_a <- delta_j * f + delta_nj * (1 - f)
      }
      if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
        N <- getN_Surv_Equi(delta = delta_a, cut = cut, alpha = alpha, beta = beta, N = NA, r = r, criterion = criterion)$N
      }
      Nj <- N * f
      set.seed(i + seed)
      xt_j <- rexp(n = Nj * r / (1 + r), rate = exp(delta_j))
      xc_j <- rexp(n = Nj / (1 + r), rate = 1)
      dat_j <- data.frame(time = c(xt_j, xc_j), status = 1, trt = c(rep(1, length(xt_j)), rep(0, length(xc_j))))
      if (isNA_delta_nj & (!isNA_delta_a)) {
        ut_nj <- runif(n = (N - Nj) * r / (1 + r), min = 0, max = 1)
        uc_nj <- runif(n = (N - Nj) / (1 + r), min = 0, max = 1)
        xc_nj <- -log(1 - uc_nj)
        xt_nj <- c()
        g <- 1
        while (g < (length(ut_nj) + 1)) {
          u <- ut_nj[g]
          solve_t <- function(t) {
            exp(-exp(delta_a) * t) - exp(-exp(delta_j) * t) * f - (1 - u) * (1 - f)
          }
          t <- uniroot(f = solve_t, interval = c(0, 1e+06))$root
          xt_nj <- c(xt_nj, t)
          g <- g + 1
        }
        xt <- c(xt_j, xt_nj)
        xc <- c(xc_j, xc_nj)
      }
      if (isNA_delta_a & (!isNA_delta_nj)) {
        xt_nj <- rexp(n = (N - Nj) * r / (1 + r), rate = exp(delta_nj))
        xc_nj <- rexp(n = (N - Nj) / (1 + r), rate = 1)
        xt <- c(xt_j, xt_nj)
        xc <- c(xc_j, xc_nj)
      }
      dat_a <- data.frame(time = c(xt, xc), status = 1, trt = c(rep(1, length(xt)), rep(0, length(xc))))
      fit_j <- coxph(Surv(time, status) ~ trt, dat = dat_j)
      coef_j <- coef(fit_j)
      fit_a <- coxph(Surv(time, status) ~ trt, dat = dat_a)
      coef_a <- coef(fit_a)
      za1 <- summary(fit_a)$coefficients[4] + cut / summary(fit_a)$coefficients[3]
      za2 <- summary(fit_a)$coefficients[4] - cut / summary(fit_a)$coefficients[3]
      if (criterion == 1) {
        zj1 <- coef_j - coef_a + pi * cut
        zj2 <- coef_j - coef_a - pi * cut
      }
      if (criterion == 2) {
        zj1 <- -(1 / exp(cut) - exp(coef_j) - pi * (1 / exp(cut) - exp(coef_a)))
        zj2 <- -(exp(cut) - exp(coef_j) - pi * (exp(cut) - exp(coef_a)))
      }
      succ_a <- if_else(za1 > qnorm(1 - alpha) & za2 < -qnorm(1 - alpha), 1, 0)
      succ_j <- if_else(zj1 > 0 & zj2 < 0, 1, 0)
      data.frame(delta_a, delta_j, delta_nj, f, pi, cut, alpha, beta, N, r, criterion, succ_a, succ_j)
    }, .progress = TRUE, .options = furrr_options(seed = TRUE))
    if (numcore >= 2) {
      plan(sequential)
    }
    res <- suppressMessages({
      if (isNA_delta_nj & (!isNA_delta_a)) {
        df <- ss %>%
          group_by(delta_a, delta_j, delta_nj, f, pi, cut, alpha, beta, N, r, criterion) %>%
          summarise(pwr1 = mean(succ_a, na.rm = TRUE), pwr2 = mean(succ_j, na.rm = TRUE), pwr3 = mean(succ_a & succ_j, na.rm = TRUE), pwr4 = mean(succ_j[succ_a == 1], na.rm = TRUE)) %>%
          arrange(f) %>%
          as.data.frame()
      }
      if (isNA_delta_a & (!isNA_delta_nj)) {
        df <- ss %>%
          group_by(delta_a, delta_j, delta_nj, f, pi, cut, alpha, beta, N, r, criterion) %>%
          summarise(pwr1 = mean(succ_a, na.rm = TRUE), pwr2 = mean(succ_j, na.rm = TRUE), pwr3 = mean(succ_a & succ_j, na.rm = TRUE), pwr4 = mean(succ_j[succ_a == 1], na.rm = TRUE)) %>%
          arrange(f) %>%
          as.data.frame()
      }
      df
    })
  }
  return(res)
}
