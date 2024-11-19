#' Power of mRCT using Japan's Method 1 for survival endpoints
#'
#' Based on Japan's Method 1, given the global and target region number of events, calculate and simulate the marginal probabilities, conditional probabilities, and joint probabilities of global success and efficacy consistency between target region and global, in clinical trials using superiority, non-inferiority, and equivalence designs with survival endpoints.
#'
#' @rdname getPwr_Surv_JM1
#'
#' @name getPwr_Surv_JM1
#'
#' @param delta_j A vector. log(HR) between treatment and control groups in target region.
#' @param delta_nj A vector. log(HR) between treatment and control groups in other regions. When \code{delta_nj} is not \code{NA}, \code{delta_a} will be calculated automatically.
#' @param delta_a A vector. log(HR) between treatment and control groups globally.
#' @param f A vector. Proportion of number of events allocated to target region.
#' @param pi A vector. Proportion of global efficacy to retain. Default value is 0.5, which means retaining half of the efficacy.
#' @param cut A vector. Positive value for non-inferiority or equivalence margin. For example, if the non-inferiority margin for HR is 0.6, then \code{cut = -log(0.6)}. If the non-inferiority margin for HR is 1.3, then \code{cut = log(1.3)}.
#' @param alpha A vector. One-sided type I error rate for global success. Default value is 0.025.
#' @param beta A vector. Type II error rate for global success, which is used to calculate global number of events only when \code{Ne} is \code{NA}.
#' @param Ne A vector. Global number of events. When \code{Ne} is \code{NA} and \code{beta} is not \code{NA}, \code{Ne} will be calculated automatically.
#' @param r A vector. Ratio of the number of events of the treatment group to the control group. Default value is 1.
#' @param criterion A vector. If \code{criterion = 1}, the consistency criterion defined on the log(HR) scale will be used. If \code{criterion = 2}, the consistency criterion defined on the HR scale will be used. See \code{details} for more information.
#' @param direct \code{direct = 1} indicates that a larger HR is preferable, while \code{direct = -1} indicates that a smaller HR is preferable.
#' @param sim Logical value. When set to \code{FALSE}, theoretical calculation is performed. When set to \code{TRUE}, simulation is used, which is more time-consuming.
#' @param nsim Number of simulations.
#' @param seed Random seed for simulation.
#' @param numcore Number of CPU cores to use during simulation. Default value is 2.
#' @param maxNe Maximum possible global number of events (\code{Ne}) in equivalence design. Default value is 1e+06.
#'
#' @return A data frame containing input parameters and returned power.
#' \describe{
#'   \item{\code{pwr1 }}{The marginal probability of global success.}
#'   \item{\code{pwr2 }}{The marginal probability that the target region efficacy is consistent with the global efficacy.}
#'   \item{\code{pwr3 }}{The joint probability of global success and the target region efficacy being consistent with the global efficacy.}
#'   \item{\code{pwr4 }}{The conditional probability that the target region efficacy is consistent with the global efficacy given global success.}
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
#' For \code{criterion = 2}, by delta method, \eqn{1 - e^{\hat \delta}} approximately follows a distribution of \eqn{N(1 - e^{\delta}, (e^{\delta}\sqrt{\frac{1}{N_e / (r + 1)} + \frac{1}{N_e r / (r + 1)}})^2)}, used in theoretical calculation (\code{sim=FALSE}).
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
#'   pi = 0.5, alpha = 0.025, beta = NA, Ne = 200, r = 1, criterion = 1,
#'   sim = FALSE
#' )
#'
#' # Delta_a will be calculated based on delta_j and delta_nj.
#' # Global number of events will be calculated based on alpha and beta.
#' getPwr_Surv_Noninf_JM1(
#'   delta_j = log(1.1),
#'   delta_nj = log(1.0),
#'   f = seq(0.1, 0.9, 0.1),
#'   cut = log(1.3),
#'   pi = 0.5, alpha = 0.025, beta = 0.2, Ne = NA, r = 1, criterion = 2,
#'   direct = -1, sim = FALSE
#' )
getPwr_Surv_Super_JM1 <- function(delta_j, delta_nj = NA, delta_a = NA, f, pi = 0.5, alpha = 0.025, beta = NA, Ne = NA, r = 1, criterion = 1, sim = FALSE, nsim = 1000, seed = 0, numcore = 2) {
  isNA_delta_nj <- is.na(delta_nj)
  isNA_delta_a <- is.na(delta_a)
  eg <- as.data.frame(expand.grid(delta_j = delta_j, delta_nj = delta_nj, delta_a = delta_a, f = f, pi = pi, alpha = alpha, beta = beta, Ne = Ne, r = r, criterion = criterion, stringsAsFactors = FALSE))
  set.seed(seed)
  seed1 <- sample(x = 1:1e8, size = nrow(eg) * nsim, replace = FALSE)
  if (sim & numcore >= 2) {
    future::plan(future::multisession, workers = numcore)
  }
  res <- furrr::future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta_j <- R$delta_j
    delta_nj <- R$delta_nj
    delta_a <- R$delta_a
    f <- R$f
    pi <- R$pi
    alpha <- R$alpha
    beta <- R$beta
    Ne <- R$Ne
    r <- R$r
    criterion <- R$criterion
    if (isNA_delta_nj & isNA_delta_a) {
      stop("Delta_nj and delta_a cannot both be NA.")
    }
    if (!isNA_delta_nj & !isNA_delta_a) {
      warning("When delta_nj is not NA, delta_a will be calculated based on delta_j and delta_nj.")
    }
    if (f < 0 | f > 1) {
      stop("Parameter f should be between 0 and 1.")
    }
    if (pi < 0 | pi > 1) {
      warning("Parameter pi generally is between 0 and 1.")
    }
    if (is.na(beta) & is.na(Ne)) {
      stop("Beta and Ne cannot both be NA.")
    }
    if (!is.na(beta) & (!is.na(Ne))) {
      warning("When beta is not NA, Ne will be automatically calculated.")
    }
    if (!criterion %in% c(1, 2)) {
      stop("Parameter criterion can only be `1` or `2`.")
    }
    if (!is.logical(sim)) {
      stop("Parameter sim should be one of `TRUE` or `FALSE`.")
    }
    if (isNA_delta_nj & (!isNA_delta_a)) {
      delta_nj <- (delta_a - delta_j * f) / (1 - f)
    }
    if (isNA_delta_a & (!isNA_delta_nj)) {
      delta_a <- delta_j * f + delta_nj * (1 - f)
    }
    if (!is.na(beta)) {
      Ne <- getNe_Surv_Super(delta = delta_a, alpha = alpha, beta = beta, Ne = NA, r = r)$Ne
    }
    Nej <- Ne * f
    if (!sim) {
      gr <- 2 + r + 1 / r
      se <- sqrt(gr / Ne)
      u <- delta_a / se
      if (criterion == 1) {
        sej <- sqrt(gr / Nej + pi^2 * gr / Ne - 2 * pi * sqrt(f) * sqrt(gr / Nej * gr / Ne))
        uj <- (delta_j - pi * delta_a) / sej
        cov <- sqrt(f) * sqrt(gr / Ne * gr / Nej) - pi * gr / Ne
        corr <- cov / (sej * se)
      }
      if (criterion == 2) {
        sej <- sqrt(gr / Nej * exp(delta_j)^2 + pi^2 * gr / Ne * exp(delta_a)^2 - 2 * pi * gr / Ne * exp(delta_j) * exp(delta_a))
        uj <- -(1 - exp(delta_j) - pi * (1 - exp(delta_a))) / sej
        cov <- gr / Ne * exp(delta_j) * exp(delta_a) - pi * gr / Ne * exp(delta_a)^2
        se_ <- sqrt(gr / Ne * exp(delta_a)^2)
        corr <- cov / (sej * se_)
      }
      if (delta_a < 0) {
        u <- (-1) * u
        uj <- (-1) * uj
      }
      M <- matrix(c(1, corr, corr, 1), nrow = 2, byrow = T)
      pwr1 <- mvtnorm::pmvnorm(lower = c(qnorm(1 - alpha), -Inf), upper = c(Inf, Inf), mean = c(u, uj), corr = M)
      pwr2 <- mvtnorm::pmvnorm(lower = c(-Inf, 0), upper = c(Inf, Inf), mean = c(u, uj), corr = M)
      pwr3 <- mvtnorm::pmvnorm(lower = c(qnorm(1 - alpha), 0), upper = c(Inf, Inf), mean = c(u, uj), corr = M)
      pwr4 <- pwr3 / pwr1
      df <- data.frame(delta_a, delta_j, delta_nj, f, pi, alpha, beta, Ne, r, criterion, pwr1, pwr2, pwr3, pwr4)
    }
    if (sim) {
      simda <- NULL
      for (ii in 1:nsim) {
        seed2 <- seed1[((i - 1) * nsim + 1):(i * nsim)]
        set.seed(seed2[ii])
        xt_j <- rexp(n = Nej * r / (1 + r), rate = exp(delta_j))
        xc_j <- rexp(n = Nej / (1 + r), rate = 1)
        dat_j <- data.frame(time = c(xt_j, xc_j), status = 1, trt = c(rep(1, length(xt_j)), rep(0, length(xc_j))))
        if (isNA_delta_nj & (!isNA_delta_a)) {
          ut_nj <- runif(n = (Ne - Nej) * r / (1 + r), min = 0, max = 1)
          uc_nj <- runif(n = (Ne - Nej) / (1 + r), min = 0, max = 1)
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
          xt_nj <- rexp(n = (Ne - Nej) * r / (1 + r), rate = exp(delta_nj))
          xc_nj <- rexp(n = (Ne - Nej) / (1 + r), rate = 1)
          xt <- c(xt_j, xt_nj)
          xc <- c(xc_j, xc_nj)
        }
        dat_a <- data.frame(time = c(xt, xc), status = 1, trt = c(rep(1, length(xt)), rep(0, length(xc))))
        fit_j <- survival::coxph(survival::Surv(time, status) ~ trt, dat = dat_j)
        coef_j <- coef(fit_j)
        fit_a <- survival::coxph(survival::Surv(time, status) ~ trt, dat = dat_a)
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
        succ_a <- dplyr::if_else(za > qnorm(1 - alpha), 1, 0)
        succ_j <- dplyr::if_else(zj > 0, 1, 0)
        da <- data.frame(delta_a, delta_j, delta_nj, f, pi, alpha, beta, Ne, r, criterion, succ_a, succ_j)
        simda <- dplyr::bind_rows(simda, da)
      }
      df <- simda %>%
        dplyr::group_by(delta_a, delta_j, delta_nj, f, pi, alpha, beta, Ne, r, criterion) %>%
        dplyr::summarise(pwr1 = mean(succ_a, na.rm = TRUE), pwr2 = mean(succ_j, na.rm = TRUE), pwr3 = mean(succ_a & succ_j, na.rm = TRUE), pwr4 = mean(succ_j[succ_a == 1], na.rm = TRUE), .groups = "keep") %>%
        dplyr::arrange(f) %>%
        as.data.frame()
    }
    df
  }, .progress = TRUE, .options = furrr::furrr_options(seed = TRUE))
  if (sim & numcore >= 2) {
    future::plan(future::sequential)
  }
  return(res)
}


#' @rdname getPwr_Surv_JM1
#' @export
getPwr_Surv_Noninf_JM1 <- function(delta_j, delta_nj = NA, delta_a = NA, f, pi = 0.5, cut, alpha = 0.025, beta = NA, Ne = NA, r = 1, criterion = 1, direct = 1, sim = FALSE, nsim = 1000, seed = 0, numcore = 2) {
  isNA_delta_nj <- is.na(delta_nj)
  isNA_delta_a <- is.na(delta_a)
  eg <- as.data.frame(expand.grid(delta_j = delta_j, delta_nj = delta_nj, delta_a = delta_a, f = f, pi = pi, cut = cut, alpha = alpha, beta = beta, Ne = Ne, r = r, criterion = criterion, stringsAsFactors = FALSE))
  set.seed(seed)
  seed1 <- sample(x = 1:1e8, size = nrow(eg) * nsim, replace = FALSE)
  if (sim & numcore >= 2) {
    future::plan(future::multisession, workers = numcore)
  }
  res <- furrr::future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta_j <- R$delta_j
    delta_nj <- R$delta_nj
    delta_a <- R$delta_a
    f <- R$f
    pi <- R$pi
    cut <- R$cut
    alpha <- R$alpha
    beta <- R$beta
    Ne <- R$Ne
    r <- R$r
    criterion <- R$criterion
    if (isNA_delta_nj & isNA_delta_a) {
      stop("Delta_nj and delta_a cannot both be NA.")
    }
    if (!isNA_delta_nj & !isNA_delta_a) {
      warning("When delta_nj is not NA, delta_a will be calculated based on delta_j and delta_nj.")
    }
    if (f < 0 | f > 1) {
      stop("Parameter f should be between 0 and 1.")
    }
    if (pi < 0 | pi > 1) {
      warning("Parameter pi generally is between 0 and 1.")
    }
    if (is.na(beta) & cut < 0) {
      warning("Parameter cut should be a positive value.")
    }
    if (is.na(beta) & is.na(Ne)) {
      stop("Beta and Ne cannot both be NA.")
    }
    if (!is.na(beta) & (!is.na(Ne))) {
      warning("When beta is not NA, Ne will be automatically calculated.")
    }
    if (!criterion %in% c(1, 2)) {
      stop("Parameter criterion can only be `1` or `2`.")
    }
    if (!direct %in% c(-1, 1)) {
      stop("Parameter direct should be one of `1` or `-1`.")
    }
    if (!is.logical(sim)) {
      stop("Parameter sim should be one of `TRUE` or `FALSE`.")
    }
    if (isNA_delta_nj & (!isNA_delta_a)) {
      delta_nj <- (delta_a - delta_j * f) / (1 - f)
    }
    if (isNA_delta_a & (!isNA_delta_nj)) {
      delta_a <- delta_j * f + delta_nj * (1 - f)
    }
    if (!is.na(beta)) {
      Ne <- getNe_Surv_Noninf(delta = delta_a, cut = cut, alpha = alpha, beta = beta, Ne = NA, r = r, direct = direct)$Ne
    }
    Nej <- Ne * f
    if (!sim) {
      gr <- 2 + r + 1 / r
      se <- sqrt(gr / Ne)
      u <- dplyr::if_else(direct == 1, (delta_a + cut) / se, (delta_a - cut) / se)
      if (criterion == 1) {
        sej <- sqrt(gr / Nej + gr / Ne - 2 * sqrt(f) * sqrt(gr / Nej * gr / Ne))
        uj <- dplyr::if_else(direct == 1, (delta_j - delta_a + pi * cut) / sej, (delta_j - delta_a - pi * cut) / sej)
        cov <- sqrt(f) * sqrt(gr / Ne * gr / Nej) - gr / Ne
        corr <- cov / (sej * se)
      }
      if (criterion == 2) {
        sej <- sqrt(gr / Nej * exp(delta_j)^2 + pi^2 * gr / Ne * exp(delta_a)^2 - 2 * pi * gr / Ne * exp(delta_j) * exp(delta_a))
        uj <- dplyr::if_else(direct == 1, -(1 / exp(cut) - exp(delta_j) - pi * (1 / exp(cut) - exp(delta_a))) / sej, -(exp(cut) - exp(delta_j) - pi * (exp(cut) - exp(delta_a))) / sej)
        cov <- gr / Ne * exp(delta_j) * exp(delta_a) - pi * gr / Ne * exp(delta_a)^2
        se_ <- sqrt(gr / Ne * exp(delta_a)^2)
        corr <- cov / (sej * se_)
      }
      M <- matrix(c(1, corr, corr, 1), nrow = 2, byrow = T)
      if (direct == -1) {
        uj <- (-1) * uj
        u <- (-1) * u
      }
      pwr1 <- mvtnorm::pmvnorm(lower = c(qnorm(1 - alpha), -Inf), upper = c(Inf, Inf), mean = c(u, uj), corr = M)
      pwr2 <- mvtnorm::pmvnorm(lower = c(-Inf, 0), upper = c(Inf, Inf), mean = c(u, uj), corr = M)
      pwr3 <- mvtnorm::pmvnorm(lower = c(qnorm(1 - alpha), 0), upper = c(Inf, Inf), mean = c(u, uj), corr = M)
      pwr4 <- pwr3 / pwr1
      df <- data.frame(delta_a, delta_j, delta_nj, f, pi, cut, alpha, beta, Ne, r, direct, pwr1, pwr2, pwr3, pwr4)
    }
    if (sim) {
      simda <- NULL
      for (ii in 1:nsim) {
        seed2 <- seed1[((i - 1) * nsim + 1):(i * nsim)]
        set.seed(seed2[ii])
        xt_j <- rexp(n = Nej * r / (1 + r), rate = exp(delta_j))
        xc_j <- rexp(n = Nej / (1 + r), rate = 1)
        dat_j <- data.frame(time = c(xt_j, xc_j), status = 1, trt = c(rep(1, length(xt_j)), rep(0, length(xc_j))))
        if (isNA_delta_nj & (!isNA_delta_a)) {
          ut_nj <- runif(n = (Ne - Nej) * r / (1 + r), min = 0, max = 1)
          uc_nj <- runif(n = (Ne - Nej) / (1 + r), min = 0, max = 1)
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
          xt_nj <- rexp(n = (Ne - Nej) * r / (1 + r), rate = exp(delta_nj))
          xc_nj <- rexp(n = (Ne - Nej) / (1 + r), rate = 1)
          xt <- c(xt_j, xt_nj)
          xc <- c(xc_j, xc_nj)
        }
        dat_a <- data.frame(time = c(xt, xc), status = 1, trt = c(rep(1, length(xt)), rep(0, length(xc))))
        fit_j <- survival::coxph(survival::Surv(time, status) ~ trt, dat = dat_j)
        coef_j <- coef(fit_j)
        fit_a <- survival::coxph(survival::Surv(time, status) ~ trt, dat = dat_a)
        coef_a <- coef(fit_a)
        za <- dplyr::if_else(direct == 1, summary(fit_a)$coefficients[4] + cut / summary(fit_a)$coefficients[3], summary(fit_a)$coefficients[4] - cut / summary(fit_a)$coefficients[3])
        if (criterion == 1) {
          zj <- dplyr::if_else(direct == 1, coef_j - coef_a + pi * cut, coef_j - coef_a - pi * cut)
        }
        if (criterion == 2) {
          zj <- dplyr::if_else(direct == 1, -(1 / exp(cut) - exp(coef_j) - pi * (1 / exp(cut) - exp(coef_a))), -(exp(cut) - exp(coef_j) - pi * (exp(cut) - exp(coef_a))))
        }
        if (direct == -1) {
          zj <- (-1) * zj
          za <- (-1) * za
        }
        succ_a <- dplyr::if_else(za > qnorm(1 - alpha), 1, 0)
        succ_j <- dplyr::if_else(zj > 0, 1, 0)
        da <- data.frame(delta_a, delta_j, delta_nj, f, pi, cut, alpha, beta, Ne, r, criterion, direct, succ_a, succ_j)
        simda <- dplyr::bind_rows(simda, da)
      }
      df <- simda %>%
        dplyr::group_by(delta_a, delta_j, delta_nj, f, pi, cut, alpha, beta, Ne, r, criterion, direct) %>%
        dplyr::summarise(pwr1 = mean(succ_a, na.rm = TRUE), pwr2 = mean(succ_j, na.rm = TRUE), pwr3 = mean(succ_a & succ_j, na.rm = TRUE), pwr4 = mean(succ_j[succ_a == 1], na.rm = TRUE), .groups = "keep") %>%
        dplyr::arrange(f) %>%
        as.data.frame()
    }
    df
  }, .progress = TRUE, .options = furrr::furrr_options(seed = TRUE))
  if (sim & numcore >= 2) {
    future::plan(future::sequential)
  }
  return(res)
}


#' @rdname getPwr_Surv_JM1
#' @export
getPwr_Surv_Equi_JM1 <- function(delta_j, delta_nj = NA, delta_a = NA, f, pi = 0.5, cut, alpha = 0.025, beta = NA, Ne = NA, r = 1, criterion = 1, sim = FALSE, nsim = 1000, seed = 0, numcore = 2, maxNe = 1e+06) {
  isNA_delta_nj <- is.na(delta_nj)
  isNA_delta_a <- is.na(delta_a)
  eg <- as.data.frame(expand.grid(delta_j = delta_j, delta_nj = delta_nj, delta_a = delta_a, f = f, pi = pi, cut = cut, alpha = alpha, beta = beta, Ne = Ne, r = r, criterion = criterion, stringsAsFactors = FALSE))
  set.seed(seed)
  seed1 <- sample(x = 1:1e8, size = nrow(eg) * nsim, replace = FALSE)
  if (sim & numcore >= 2) {
    future::plan(future::multisession, workers = numcore)
  }
  res <- furrr::future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta_j <- R$delta_j
    delta_nj <- R$delta_nj
    delta_a <- R$delta_a
    f <- R$f
    pi <- R$pi
    cut <- R$cut
    alpha <- R$alpha
    beta <- R$beta
    Ne <- R$Ne
    r <- R$r
    criterion <- R$criterion
    if (isNA_delta_nj & isNA_delta_a) {
      stop("Delta_nj and delta_a cannot both be NA.")
    }
    if (!isNA_delta_nj & !isNA_delta_a) {
      warning("When delta_nj is not NA, delta_a will be calculated based on delta_j and delta_nj.")
    }
    if (f < 0 | f > 1) {
      stop("Parameter f should be between 0 and 1.")
    }
    if (pi < 0 | pi > 1) {
      warning("Parameter pi generally is between 0 and 1.")
    }
    if (is.na(beta) & cut < 0) {
      warning("Parameter cut should be a positive value.")
    }
    if (is.na(beta) & is.na(Ne)) {
      stop("Beta and Ne cannot both be NA.")
    }
    if (!is.na(beta) & (!is.na(Ne))) {
      warning("When beta is not NA, Ne will be automatically calculated.")
    }
    if (!criterion %in% c(1, 2)) {
      stop("Parameter criterion can only be `1` or `2`.")
    }
    if (!is.logical(sim)) {
      stop("Parameter sim should be one of `TRUE` or `FALSE`.")
    }
    if (isNA_delta_nj & (!isNA_delta_a)) {
      delta_nj <- (delta_a - delta_j * f) / (1 - f)
    }
    if (isNA_delta_a & (!isNA_delta_nj)) {
      delta_a <- delta_j * f + delta_nj * (1 - f)
    }
    if (!is.na(beta)) {
      Ne <- getNe_Surv_Equi(delta = delta_a, cut = cut, alpha = alpha, beta = beta, Ne = NA, r = r, maxNe = maxNe)$Ne
    }
    Nej <- Ne * f
    if (!sim) {
      gr <- 2 + r + 1 / r
      se <- sqrt(gr / Ne)
      u1 <- (delta_a + cut) / se
      u2 <- (delta_a - cut) / se
      if (criterion == 1) {
        sej <- sqrt(gr / Nej + gr / Ne - 2 * sqrt(f) * sqrt(gr / Nej * gr / Ne))
        uj1 <- (delta_j - delta_a + pi * cut) / sej
        uj2 <- (delta_j - delta_a - pi * cut) / sej
        cov <- sqrt(f) * sqrt(gr / Ne * gr / Nej) - gr / Ne
        corr <- cov / (sej * se)
      }
      if (criterion == 2) {
        sej <- sqrt(gr / Nej * exp(delta_j)^2 + pi^2 * gr / Ne * exp(delta_a)^2 - 2 * pi * gr / Ne * exp(delta_j) * exp(delta_a))
        uj1 <- -(1 / exp(cut) - exp(delta_j) - pi * (1 / exp(cut) - exp(delta_a))) / sej
        uj2 <- -(exp(cut) - exp(delta_j) - pi * (exp(cut) - exp(delta_a))) / sej
        cov <- gr / Ne * exp(delta_j) * exp(delta_a) - pi * gr / Ne * exp(delta_a)^2
        se_ <- sqrt(gr / Ne * exp(delta_a)^2)
        corr <- cov / (sej * se_)
      }
      M <- matrix(c(1, 1, corr, corr, 1, 1, corr, corr, corr, corr, 1, 1, corr, corr, 1, 1), nrow = 4, byrow = T)
      pwr1 <- mvtnorm::pmvnorm(lower = c(qnorm(1 - alpha), -Inf, -Inf, -Inf), upper = c(Inf, -qnorm(1 - alpha), Inf, Inf), mean = c(u1, u2, uj1, uj2), corr = M)
      pwr2 <- mvtnorm::pmvnorm(lower = c(-Inf, -Inf, 0, -Inf), upper = c(Inf, Inf, Inf, 0), mean = c(u1, u2, uj1, uj2), corr = M)
      pwr3 <- mvtnorm::pmvnorm(lower = c(qnorm(1 - alpha), -Inf, 0, -Inf), upper = c(Inf, -qnorm(1 - alpha), Inf, 0), mean = c(u1, u2, uj1, uj2), corr = M)
      pwr4 <- pwr3 / pwr1
      df <- data.frame(delta_a, delta_j, delta_nj, f, pi, cut, alpha, beta, Ne, r, pwr1, pwr2, pwr3, pwr4)
    }
    if (sim) {
      simda <- NULL
      for (ii in 1:nsim) {
        seed2 <- seed1[((i - 1) * nsim + 1):(i * nsim)]
        set.seed(seed2[ii])
        xt_j <- rexp(n = Nej * r / (1 + r), rate = exp(delta_j))
        xc_j <- rexp(n = Nej / (1 + r), rate = 1)
        dat_j <- data.frame(time = c(xt_j, xc_j), status = 1, trt = c(rep(1, length(xt_j)), rep(0, length(xc_j))))
        if (isNA_delta_nj & (!isNA_delta_a)) {
          ut_nj <- runif(n = (Ne - Nej) * r / (1 + r), min = 0, max = 1)
          uc_nj <- runif(n = (Ne - Nej) / (1 + r), min = 0, max = 1)
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
          xt_nj <- rexp(n = (Ne - Nej) * r / (1 + r), rate = exp(delta_nj))
          xc_nj <- rexp(n = (Ne - Nej) / (1 + r), rate = 1)
          xt <- c(xt_j, xt_nj)
          xc <- c(xc_j, xc_nj)
        }
        dat_a <- data.frame(time = c(xt, xc), status = 1, trt = c(rep(1, length(xt)), rep(0, length(xc))))
        fit_j <- survival::coxph(survival::Surv(time, status) ~ trt, dat = dat_j)
        coef_j <- coef(fit_j)
        fit_a <- survival::coxph(survival::Surv(time, status) ~ trt, dat = dat_a)
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
        succ_a <- dplyr::if_else(za1 > qnorm(1 - alpha) & za2 < -qnorm(1 - alpha), 1, 0)
        succ_j <- dplyr::if_else(zj1 > 0 & zj2 < 0, 1, 0)
        da <- data.frame(delta_a, delta_j, delta_nj, f, pi, cut, alpha, beta, Ne, r, criterion, succ_a, succ_j)
        simda <- dplyr::bind_rows(simda, da)
      }
      df <- simda %>%
        dplyr::group_by(delta_a, delta_j, delta_nj, f, pi, cut, alpha, beta, Ne, r, criterion) %>%
        dplyr::summarise(pwr1 = mean(succ_a, na.rm = TRUE), pwr2 = mean(succ_j, na.rm = TRUE), pwr3 = mean(succ_a & succ_j, na.rm = TRUE), pwr4 = mean(succ_j[succ_a == 1], na.rm = TRUE), .groups = "keep") %>%
        dplyr::arrange(f) %>%
        as.data.frame()
    }
    df
  }, .progress = TRUE, .options = furrr::furrr_options(seed = TRUE))
  if (sim & numcore >= 2) {
    future::plan(future::sequential)
  }
  return(res)
}
