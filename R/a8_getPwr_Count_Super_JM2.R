#' Power of mRCT using Japan's Method 2 for count endpoints
#'
#' Based on Japan's Method 2, given the global and target region sample sizes, calculate and simulate the marginal probabilities, conditional probabilities, and joint probabilities of global success and efficacy consistency between target region and globally, in clinical trials using superiority, non-inferiority, and equivalence designs with count endpoints.
#'
#' @rdname getPwr_Count_JM2
#'
#' @name getPwr_Count_JM2
#'
#' @param delta_i A vector. log(RR) between treatment and control groups for each region.
#' @param lambda0_i A vector. Baseline hazard of control group for each region.
#' @param t Average exposure time.
#' @param k The over-dispersion parameter for negative binomial distribution, which is 0 for poisson distribution.
#' @param fi A vector. Proportion of sample size allocated to each region.
#' @param cut A positive value for non-inferiority or equivalence margin. For example, if the non-inferiority margin for RR is 0.6, then \code{cut = -log(0.6)}. If the non-inferiority margin for RR is 1.3, then \code{cut = log(1.3)}.
#' @param alpha One-sided type I error rate for global success. The default value is 0.025.
#' @param beta Type II error rate for global success, which is used to calculate global sample size only when \code{N} is \code{NA}.
#' @param N Global sample size. Global sample size. When \code{N} is \code{NA} and \code{beta} is not \code{NA}, \code{N} will be calculated automatically.
#' @param r Ratio of the sample sizes of the treatment group to the control group. The default value is 1.
#' @param direct \code{direct = 1} indicates that a larger RR is preferable, while \code{direct = -1} indicates that a smaller RR is preferable.
#' @param sim Logical value. When set to \code{FALSE}, theoretical calculation is performed. When set to \code{TRUE}, simulation is used, which is more time-consuming.
#' @param nsim Number of simulations.
#' @param seed Random seed for simulation.
#'
#' @return A list where:
#' \itemize{
#'   \item{overall}{
#'     \itemize{
#'       \item{\code{pwr1 }}{The marginal probability of global success.}
#'       \item{\code{pwr2 }}{The marginal probability that all region's efficacy is consistent with the global efficacy.}
#'       \item{\code{pwr3 }}{The conditional probability that all region's efficacy is consistent with the global efficacy given global success.}
#'       \item{\code{pwr4 }}{The joint probability of global success and all region's efficacy being consistent with the global efficacy.}
#'     }
#'   }
#'   \item{\code{pwr_margin }}{The marginal probability that the ith region efficacy is consistent with the global efficacy.}
#'   \item{\code{pwr_condition }}{The conditional probability that the ith region efficacy is consistent with the global efficacy given global success.}
#'   \item{\code{pwr_joint }}{The joint probability of global success and the ith region efficacy being consistent with the global efficacy.}
#' }
#'
#' @details
#' Taking the larger RR is preferable as an example. The global success criterion and the efficacy consistency criterion between target region and globally
#'
#' in superiority design:
#' \deqn{Z_a = \frac{\hat \delta_a}{\sqrt{Var(\hat \delta_a)}} > \Phi^{-1}(1 - \alpha)}
#' \deqn{\hat \delta_i > 0 \text{ for i = 1, 2, .., m}}
#'
#' in non-inferiority design:
#' \deqn{Z_a = \frac{\hat \delta_a + \Delta}{\sqrt{Var(\hat \delta_a)}} > \Phi^{-1}(1 - \alpha)}
#' \deqn{\hat \delta_i + \Delta> 0 \text{ for i = 1, 2, .., m}}
#'
#' in equivalence design:
#' \deqn{Z_{a_u} = \frac{\hat \delta_a + \Delta}{\sqrt{Var(\hat \delta_a)}} > \Phi^{-1}(1 - \alpha)\text{ and }Z_{a_l} = \frac{\hat \delta_a - \Delta}{\sqrt{Var(\hat \delta_a)}} < \Phi^{-1}(\alpha)}
#' \deqn{\hat \delta_i + \Delta> 0 \text{ and } \hat \delta_i - \Delta < 0 \text{ for i = 1, 2, .., m}}
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
#' f_set <- seq(0.1, 0.9, 0.1)
#' map_dfr(
#'   .x = 1:length(f_set),
#'   .f = function(i) {
#'     f <- f_set[i]
#'     res <- getPwr_Count_Super_JM2(
#'       delta_i = c(
#'         log(1.2),
#'         log(1.4)
#'       ),
#'       lambda0_i = c(0.1, 0.1),
#'       t = 5, k = 0, fi = c(f, 1 - f),
#'       alpha = 0.025, beta = NA, N = 300, r = 1, sim = FALSE
#'     )$overall
#'     res$f <- f
#'     res
#'   }
#' )
#'
#' # Global log(RR) will be calculated based on delta_i and fi,
#' # and global lambda0 will be calculated based on lambda0_i and fi
#' # Global sample size will be calculated based on beta.
#' f_set <- seq(0.1, 0.9, 0.1)
#' map_dfr(
#'   .x = 1:length(f_set),
#'   .f = function(i) {
#'     f <- f_set[i]
#'     res <- getPwr_Count_Noninf_JM2(
#'       delta_i = c(
#'         log(1.1),
#'         log(1.0)
#'       ),
#'       lambda0_i = c(0.1, 0.1),
#'       t = 5, k = 0, fi = c(f, 1 - f),
#'       cut = log(1.3),
#'       alpha = 0.025, beta = 0.2, N = NA, r = 1, direct = -1, sim = FALSE
#'     )$overall
#'     res$f <- f
#'     res
#'   }
#' )
getPwr_Count_Super_JM2 <- function(delta_i, lambda0_i, t, k = 0, fi, alpha = 0.025, beta = NA, N = NA, r = 1, sim = FALSE, nsim = 1000, seed = 0) {
  if (!sim) {
    num <- length(delta_i)
    delta_a <- sum(delta_i * fi)
    lambda0_a <- sum(lambda0_i * fi)
    if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
      N <- getN_Count_Super(delta = delta_a, lambda0 = lambda0_a, t = t, k = k, alpha = alpha, beta = beta, N = NA, r = r)$N
    }
    Ni <- N * fi
    lambda1_i <- exp(delta_i) * lambda0_i
    lambda1_a <- exp(delta_a) * lambda0_a
    sigma1_i <- sqrt(1 / lambda1_i / t + k)
    sigma0_i <- sqrt(1 / lambda0_i / t + k)
    sigma1_a <- sqrt(1 / lambda1_a / t + k)
    sigma0_a <- sqrt(1 / lambda0_a / t + k)
    var_i <- sigma1_i^2 / (r * Ni / (1 + r)) + sigma0_i^2 / (Ni / (1 + r))
    var_a <- sigma1_a^2 / (r * N / (1 + r)) + sigma0_a^2 / (N / (1 + r))
    sei <- sqrt(var_i)
    ui <- delta_i / sei
    se <- sqrt(var_a)
    u <- delta_a / se
    M <- diag(num + 1)
    M[1, ] <- c(1, sqrt(fi))
    M[, 1] <- c(1, sqrt(fi))
    if (delta_a < 0) {
      ui <- (-1) * ui
      u <- (-1) * u
    }
    pwr1 <- 1 - pnorm(q = qnorm(1 - alpha), mean = u, sd = 1)
    pwr2 <- prod(1 - pnorm(q = 0, mean = ui, sd = 1))
    pwr3 <- pmvnorm(lower = c(qnorm(1 - alpha), rep(0, num)), upper = rep(Inf, num + 1), mean = c(u, ui), sigma = M)
    pwr4 <- pwr3 / pwr1
    res <- data.frame(delta_a = delta_a, N = N, pwr1 = pwr1, pwr2 = pwr2, pwr3 = pwr3, pwr4 = pwr4)
    pwr_margin <- 1 - pnorm(q = 0, mean = ui, sd = 1)
    pwr_joint <- c()
    for (g in 1:num) {
      pwr_joint_ <- pmvnorm(lower = c(qnorm(1 - alpha), 0), upper = c(Inf, Inf), mean = c(u, ui[g]), sigma = M[c(1, g + 1), c(1, g + 1)])
      pwr_joint <- c(pwr_joint, pwr_joint_)
    }
    pwr_condition <- pwr_joint / pwr1
    L <- list(overall = res, pwr_margin = pwr_margin, pwr_condition = pwr_condition, pwr_joint = pwr_joint)
  }
  if (sim) {
    num <- length(delta_i)
    delta_a <- sum(delta_i * fi)
    lambda0_a <- sum(lambda0_i * fi)
    if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
      N <- getN_Count_Super(delta = delta_a, lambda0 = lambda0_a, t = t, k = k, alpha = alpha, beta = beta, N = NA, r = r)$N
    }
    Ni <- N * fi
    da <- data.frame()
    di <- NULL
    for (j in 1:nsim) {
      set.seed(j + seed)
      Xt <- matrix(NA, nrow = N * r / (1 + r), ncol = num)
      Xc <- matrix(NA, nrow = N / (1 + r), ncol = num)
      coef_i <- c()
      for (g in 1:num) {
        nt_g <- round(Ni[g] * r / (1 + r))
        nc_g <- round(Ni[g] / (1 + r))
        if (k == 0) {
          Xt[1:nt_g, g] <- rpois(n = nt_g, lambda = lambda0_i[g] * exp(delta_i[g]) * t)
          Xc[1:nc_g, g] <- rpois(n = nc_g, lambda = lambda0_i[g] * t)
        }
        if (k > 0) {
          Xt[1:nt_g, g] <- rnbinom(n = nt_g, size = k, mu = lambda0_i[g] * exp(delta_i[g]) * t)
          Xc[1:nc_g, g] <- rnbinom(n = nc_g, size = k, mu = lambda0_i[g] * t)
        }
        dat_i <- data.frame(x = c(Xt[1:nt_g, g], Xc[1:nc_g, g]), trt = c(rep(1, nt_g), rep(0, nc_g)))
        if (k == 0) {
          fit_i <- glm(x ~ trt, dat = dat_i, family = poisson(link = "log"))
        }
        if (k > 0) {
          fit_i <- glm.nb(x ~ trt, dat = dat_i)
        }
        coef_i_ <- coef(fit_i)[2]
        coef_i <- c(coef_i, coef_i_)
      }
      Xt <- as.numeric(Xt)
      Xt <- Xt[!is.na(Xt)]
      Xc <- as.numeric(Xc)
      Xc <- Xc[!is.na(Xc)]
      dat_a <- data.frame(x = c(Xt, Xc), trt = c(rep(1, length(Xt)), rep(0, length(Xc))))
      if (k == 0) {
        fit_a <- glm(x ~ trt, dat = dat_a, family = poisson(link = "log"))
      }
      if (k > 0) {
        fit_a <- glm.nb(x ~ trt, dat = dat_a)
      }
      coef_a <- coef(fit_a)[2]
      za <- summary(fit_a)$coefficients[2, 3]
      zi <- coef_i
      if (delta_a < 0) {
        za <- (-1) * za
        zi <- (-1) * zi
      }
      succ_a <- if_else(za > qnorm(1 - alpha), 1, 0)
      succ_i <- if_else(all(zi > 0), 1, 0)
      da <- bind_rows(da, c(succ_a = succ_a, succ_i = succ_i))
      di <- rbind(di, as.numeric(zi > 0))
    }
    res <- suppressMessages({
      da %>%
        summarise(pwr1 = mean(succ_a, na.rm = TRUE), pwr2 = mean(succ_i, na.rm = TRUE), pwr3 = mean(succ_a & succ_i, na.rm = TRUE), pwr4 = mean(succ_i[succ_a == 1], na.rm = TRUE)) %>%
        mutate(delta_a = delta_a, N = N) %>%
        dplyr::select(delta_a, N, pwr1, pwr2, pwr3, pwr4)
    })
    pwr_margin <- colMeans(di)
    pwr_joint <- c()
    for (g in 1:num) {
      pwr_joint_ <- mean(da$succ_a & di[, g])
      pwr_joint <- c(pwr_joint, pwr_joint_)
    }
    pwr_condition <- colMeans(di[da$succ_a == 1, ])
    L <- list(overall = res, pwr_margin = pwr_margin, pwr_condition = pwr_condition, pwr_joint = pwr_joint)
  }
  return(L)
}

#' @rdname getPwr_Count_JM2
#' @export
getPwr_Count_Noninf_JM2 <- function(delta_i, lambda0_i, t, k = 0, fi, cut, alpha = 0.025, beta = NA, N = NA, r = 1, direct = 1, sim = FALSE, nsim = 1000, seed = 0) {
  if (!sim) {
    num <- length(delta_i)
    delta_a <- sum(delta_i * fi)
    lambda0_a <- sum(lambda0_i * fi)
    if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
      N <- getN_Count_Noninf(delta = delta_a, lambda0 = lambda0_a, t = t, k = k, cut = cut, alpha = alpha, beta = beta, N = NA, r = r, direct = direct)$N
    }
    Ni <- N * fi
    lambda1_i <- exp(delta_i) * lambda0_i
    lambda1_a <- exp(delta_a) * lambda0_a
    sigma1_i <- sqrt(1 / lambda1_i / t + k)
    sigma0_i <- sqrt(1 / lambda0_i / t + k)
    sigma1_a <- sqrt(1 / lambda1_a / t + k)
    sigma0_a <- sqrt(1 / lambda0_a / t + k)
    var_i <- sigma1_i^2 / (r * Ni / (1 + r)) + sigma0_i^2 / (Ni / (1 + r))
    var_a <- sigma1_a^2 / (r * N / (1 + r)) + sigma0_a^2 / (N / (1 + r))
    sei <- sqrt(var_i)
    ui <- if_else(direct == rep(1, num), (delta_i + cut) / sei, (delta_i - cut) / sei)
    se <- sqrt(var_a)
    u <- if_else(direct == 1, (delta_a + cut) / se, (delta_a - cut) / se)
    M <- diag(num + 1)
    M[1, ] <- c(1, sqrt(fi))
    M[, 1] <- c(1, sqrt(fi))
    if (direct == -1) {
      ui <- (-1) * ui
      u <- (-1) * u
    }
    pwr1 <- 1 - pnorm(q = qnorm(1 - alpha), mean = u, sd = 1)
    pwr2 <- prod(1 - pnorm(q = 0, mean = ui, sd = 1))
    pwr3 <- pmvnorm(lower = c(qnorm(1 - alpha), rep(0, num)), upper = rep(Inf, num + 1), mean = c(u, ui), sigma = M)
    pwr4 <- pwr3 / pwr1
    res <- data.frame(delta_a = delta_a, N = N, pwr1 = pwr1, pwr2 = pwr2, pwr3 = pwr3, pwr4 = pwr4)
    pwr_margin <- 1 - pnorm(q = 0, mean = ui, sd = 1)
    pwr_joint <- c()
    for (g in 1:num) {
      pwr_joint_ <- pmvnorm(lower = c(qnorm(1 - alpha), 0), upper = c(Inf, Inf), mean = c(u, ui[g]), sigma = M[c(1, g + 1), c(1, g + 1)])
      pwr_joint <- c(pwr_joint, pwr_joint_)
    }
    pwr_condition <- pwr_joint / pwr1
    L <- list(overall = res, pwr_margin = pwr_margin, pwr_condition = pwr_condition, pwr_joint = pwr_joint)
  }
  if (sim) {
    num <- length(delta_i)
    delta_a <- sum(delta_i * fi)
    lambda0_a <- sum(lambda0_i * fi)
    if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
      N <- getN_Count_Noninf(delta = delta_a, lambda0 = lambda0_a, t = t, k = k, cut = cut, alpha = alpha, beta = beta, N = NA, r = r, direct = direct)$N
    }
    Ni <- N * fi
    da <- data.frame()
    di <- NULL
    for (j in 1:nsim) {
      set.seed(j + seed)
      Xt <- matrix(NA, nrow = N * r / (1 + r), ncol = num)
      Xc <- matrix(NA, nrow = N / (1 + r), ncol = num)
      coef_i <- c()
      for (g in 1:num) {
        nt_g <- round(Ni[g] * r / (1 + r))
        nc_g <- round(Ni[g] / (1 + r))
        if (k == 0) {
          Xt[1:nt_g, g] <- rpois(n = nt_g, lambda = lambda0_i[g] * exp(delta_i[g]) * t)
          Xc[1:nc_g, g] <- rpois(n = nc_g, lambda = lambda0_i[g] * t)
        }
        if (k > 0) {
          Xt[1:nt_g, g] <- rnbinom(n = nt_g, size = k, mu = lambda0_i[g] * exp(delta_i[g]) * t)
          Xc[1:nc_g, g] <- rnbinom(n = nc_g, size = k, mu = lambda0_i[g] * t)
        }
        dat_i <- data.frame(x = c(Xt[1:nt_g, g], Xc[1:nc_g, g]), trt = c(rep(1, nt_g), rep(0, nc_g)))
        if (k == 0) {
          fit_i <- glm(x ~ trt, dat = dat_i, family = poisson(link = "log"))
        }
        if (k > 0) {
          fit_i <- glm.nb(x ~ trt, dat = dat_i)
        }
        coef_i_ <- coef(fit_i)[2]
        coef_i <- c(coef_i, coef_i_)
      }
      Xt <- as.numeric(Xt)
      Xt <- Xt[!is.na(Xt)]
      Xc <- as.numeric(Xc)
      Xc <- Xc[!is.na(Xc)]
      dat_a <- data.frame(x = c(Xt, Xc), trt = c(rep(1, length(Xt)), rep(0, length(Xc))))
      if (k == 0) {
        fit_a <- glm(x ~ trt, dat = dat_a, family = poisson(link = "log"))
      }
      if (k > 0) {
        fit_a <- glm.nb(x ~ trt, dat = dat_a)
      }
      coef_a <- coef(fit_a)[2]
      za <- if_else(direct == 1, summary(fit_a)$coefficients[2, 3] + cut / summary(fit_a)$coefficients[2, 2], summary(fit_a)$coefficients[2, 3] - cut / summary(fit_a)$coefficients[2, 2])
      zi <- if_else(direct == rep(1, num), coef_i + cut, coef_i - cut)
      if (direct == -1) {
        za <- (-1) * za
        zi <- (-1) * zi
      }
      succ_a <- if_else(za > qnorm(1 - alpha), 1, 0)
      succ_i <- if_else(all(zi > 0), 1, 0)
      da <- bind_rows(da, c(succ_a = succ_a, succ_i = succ_i))
      di <- rbind(di, as.numeric(zi > 0))
    }
    res <- suppressMessages({
      da %>%
        summarise(pwr1 = mean(succ_a, na.rm = TRUE), pwr2 = mean(succ_i, na.rm = TRUE), pwr3 = mean(succ_a & succ_i, na.rm = TRUE), pwr4 = mean(succ_i[succ_a == 1], na.rm = TRUE)) %>%
        mutate(delta_a = delta_a, N = N) %>%
        dplyr::select(delta_a, N, pwr1, pwr2, pwr3, pwr4)
    })
    pwr_margin <- colMeans(di)
    pwr_joint <- c()
    for (g in 1:num) {
      pwr_joint_ <- mean(da$succ_a & di[, g])
      pwr_joint <- c(pwr_joint, pwr_joint_)
    }
    pwr_condition <- colMeans(di[da$succ_a == 1, ])
    L <- list(overall = res, pwr_margin = pwr_margin, pwr_condition = pwr_condition, pwr_joint = pwr_joint)
  }
  return(L)
}

#' @rdname getPwr_Count_JM2
#' @export
getPwr_Count_Equi_JM2 <- function(delta_i, lambda0_i, t, k = 0, fi, cut, alpha = 0.025, beta = NA, N = NA, r = 1, sim = FALSE, nsim = 1000, seed = 0) {
  if (!sim) {
    num <- length(delta_i)
    delta_a <- sum(delta_i * fi)
    lambda0_a <- sum(lambda0_i * fi)
    if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
      N <- getN_Count_Equi(delta = delta_a, lambda0 = lambda0_a, t = t, k = k, cut = cut, alpha = alpha, beta = beta, N = NA, r = r)$N
    }
    Ni <- N * fi
    lambda1_i <- exp(delta_i) * lambda0_i
    lambda1_a <- exp(delta_a) * lambda0_a
    sigma1_i <- sqrt(1 / lambda1_i / t + k)
    sigma0_i <- sqrt(1 / lambda0_i / t + k)
    sigma1_a <- sqrt(1 / lambda1_a / t + k)
    sigma0_a <- sqrt(1 / lambda0_a / t + k)
    var_i <- sigma1_i^2 / (r * Ni / (1 + r)) + sigma0_i^2 / (Ni / (1 + r))
    var_a <- sigma1_a^2 / (r * N / (1 + r)) + sigma0_a^2 / (N / (1 + r))
    sei <- sqrt(var_i)
    ui1 <- (delta_i + cut) / sei
    ui2 <- (delta_i - cut) / sei
    se <- sqrt(var_a)
    u1 <- (delta_a + cut) / se
    u2 <- (delta_a - cut) / se
    M <- diag(num + 1)
    M[1, ] <- c(1, sqrt(fi))
    M[, 1] <- c(1, sqrt(fi))
    M1 <- combine(M, M)
    pwr1 <- pmvnorm(lower = c(qnorm(1 - alpha), -Inf), upper = c(Inf, -qnorm(1 - alpha)), mean = c(u1, u2), corr = matrix(1, nrow = 2, ncol = 2))
    pwr_margin <- c()
    for (g in 1:num) {
      pwr_margin_ <- pmvnorm(lower = c(0, -Inf), upper = c(Inf, 0), mean = c(ui1[g], ui2[g]), corr = matrix(1, nrow = 2, ncol = 2))
      pwr_margin <- c(pwr_margin, pwr_margin_)
    }
    pwr2 <- prod(pwr_margin)
    pwr3 <- pmvnorm(lower = c(qnorm(1 - alpha), -Inf, rep(c(0, -Inf), num)), upper = c(Inf, -qnorm(1 - alpha), rep(c(Inf, 0), num)), mean = c(u1, u2, combine(ui1, ui2)), sigma = M1)
    pwr4 <- pwr3 / pwr1
    res <- data.frame(delta_a = delta_a, N = N, pwr1 = pwr1, pwr2 = pwr2, pwr3 = pwr3, pwr4 = pwr4)
    pwr_joint <- c()
    for (g in 1:num) {
      M <- diag(2)
      M[1, ] <- c(1, sqrt(fi[g]))
      M[, 1] <- c(1, sqrt(fi[g]))
      M2 <- combine(M, M)
      pwr_joint_ <- pmvnorm(lower = c(qnorm(1 - alpha), -Inf, 0, -Inf), upper = c(Inf, -qnorm(1 - alpha), Inf, 0), mean = c(u1, u2, ui1[g], ui2[g]), sigma = M2)
      pwr_joint <- c(pwr_joint, pwr_joint_)
    }
    pwr_condition <- pwr_joint / pwr1
    L <- list(overall = res, pwr_margin = pwr_margin, pwr_condition = pwr_condition, pwr_joint = pwr_joint)
  }
  if (sim) {
    num <- length(delta_i)
    delta_a <- sum(delta_i * fi)
    lambda0_a <- sum(lambda0_i * fi)
    if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
      N <- getN_Count_Equi(delta = delta_a, lambda0 = lambda0_a, t = t, k = k, cut = cut, alpha = alpha, beta = beta, N = NA, r = r)$N
    }
    Ni <- N * fi
    da <- data.frame()
    di <- NULL
    for (j in 1:nsim) {
      set.seed(j + seed)
      Xt <- matrix(NA, nrow = N * r / (1 + r), ncol = num)
      Xc <- matrix(NA, nrow = N / (1 + r), ncol = num)
      coef_i <- c()
      for (g in 1:num) {
        nt_g <- round(Ni[g] * r / (1 + r))
        nc_g <- round(Ni[g] / (1 + r))
        if (k == 0) {
          Xt[1:nt_g, g] <- rpois(n = nt_g, lambda = lambda0_i[g] * exp(delta_i[g]) * t)
          Xc[1:nc_g, g] <- rpois(n = nc_g, lambda = lambda0_i[g] * t)
        }
        if (k > 0) {
          Xt[1:nt_g, g] <- rnbinom(n = nt_g, size = k, mu = lambda0_i[g] * exp(delta_i[g]) * t)
          Xc[1:nc_g, g] <- rnbinom(n = nc_g, size = k, mu = lambda0_i[g] * t)
        }
        dat_i <- data.frame(x = c(Xt[1:nt_g, g], Xc[1:nc_g, g]), trt = c(rep(1, nt_g), rep(0, nc_g)))
        if (k == 0) {
          fit_i <- glm(x ~ trt, dat = dat_i, family = poisson(link = "log"))
        }
        if (k > 0) {
          fit_i <- glm.nb(x ~ trt, dat = dat_i)
        }
        coef_i_ <- coef(fit_i)[2]
        coef_i <- c(coef_i, coef_i_)
      }
      Xt <- as.numeric(Xt)
      Xt <- Xt[!is.na(Xt)]
      Xc <- as.numeric(Xc)
      Xc <- Xc[!is.na(Xc)]
      dat_a <- data.frame(x = c(Xt, Xc), trt = c(rep(1, length(Xt)), rep(0, length(Xc))))
      if (k == 0) {
        fit_a <- glm(x ~ trt, dat = dat_a, family = poisson(link = "log"))
      }
      if (k > 0) {
        fit_a <- glm.nb(x ~ trt, dat = dat_a)
      }
      coef_a <- coef(fit_a)[2]
      za1 <- summary(fit_a)$coefficients[2, 3] + cut / summary(fit_a)$coefficients[2, 2]
      za2 <- summary(fit_a)$coefficients[2, 3] - cut / summary(fit_a)$coefficients[2, 2]
      zi1 <- coef_i + cut
      zi2 <- coef_i - cut
      succ_a <- if_else(za1 > qnorm(1 - alpha) & za2 < -qnorm(1 - alpha), 1, 0)
      succ_i <- if_else(all(zi1 > 0) & all(zi2 < 0), 1, 0)
      da <- bind_rows(da, c(succ_a = succ_a, succ_i = succ_i))
      di <- rbind(di, as.numeric(zi1 > 0 & zi2 < 0))
    }
    res <- suppressMessages({
      da %>%
        summarise(pwr1 = mean(succ_a, na.rm = TRUE), pwr2 = mean(succ_i, na.rm = TRUE), pwr3 = mean(succ_a & succ_i, na.rm = TRUE), pwr4 = mean(succ_i[succ_a == 1], na.rm = TRUE)) %>%
        mutate(delta_a = delta_a, N = N) %>%
        dplyr::select(delta_a, N, pwr1, pwr2, pwr3, pwr4)
    })
    pwr_margin <- colMeans(di)
    pwr_joint <- c()
    for (g in 1:num) {
      pwr_joint_ <- mean(da$succ_a & di[, g])
      pwr_joint <- c(pwr_joint, pwr_joint_)
    }
    pwr_condition <- colMeans(di[da$succ_a == 1, ])
    L <- list(overall = res, pwr_margin = pwr_margin, pwr_condition = pwr_condition, pwr_joint = pwr_joint)
  }
  return(L)
}
