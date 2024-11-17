#' Power of mRCT using Japan's Method 2 for continuous endpoints
#'
#' Based on Japan's Method 2, given the global and target region sample sizes, calculate and simulate the marginal probabilities, conditional probabilities, and joint probabilities of global success and efficacy consistency between target region and globally, in clinical trials using superiority, non-inferiority, and equivalence designs with continuous endpoints.
#'
#' @rdname getPwr_Con_JM2
#'
#' @name getPwr_Con_JM2
#'
#' @param delta_i A vector with length equal to number of regions. Mean difference between treatment and control groups in each region.
#' @param sigma Common standard deviation.
#' @param f_i A vector with length equal to number of regions. Proportion of sample size allocated to each region.
#' @param cut A positive value for non-inferiority or equivalence margin in global trial, and this margin is also used in the ith region when the ith element of cut_i is \code{NA}.
#' @param cut_i A vector with length equal to number of regions. Positive value for non-inferiority or equivalence margin in each region. When \code{cut_i = NA} (default), globally margin will be used for each region.
#' @param alpha One-sided type I error rate for global success. Default value is 0.025.
#' @param beta Type II error rate for global success, which is used to calculate global sample size only when \code{N} is \code{NA}.
#' @param N Global sample size. Global sample size. When \code{N} is \code{NA} and \code{beta} is not \code{NA}, \code{N} will be calculated automatically.
#' @param r Ratio of sample sizes of treatment group to control group. Default value is 1.
#' @param direct \code{direct = 1} indicates that a larger mean difference is preferable, while \code{direct = -1} indicates that a smaller mean difference is preferable.
#' @param sim Logical value. When set to \code{FALSE}, theoretical calculation is performed. When set to \code{TRUE}, simulation is used, which is more time-consuming.
#' @param nsim Number of simulations.
#' @param seed Random seed for simulation.
#' @param maxN Maximum possible global sample size (\code{N}) in equivalence design. Default value is 1e+06.
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
#'   \item{\code{cut_i }}{The non-inferiority or equivalence margin in each region (\code{cut_i}).}
#'   \item{\code{pwr_margin }}{The marginal probability that the ith region efficacy is consistent with the global efficacy.}
#'   \item{\code{pwr_condition }}{The conditional probability that the ith region efficacy is consistent with the global efficacy given global success.}
#'   \item{\code{pwr_joint }}{The joint probability of global success and the ith region efficacy being consistent with the global efficacy.}
#' }
#'
#' @details
#' Taking the larger mean difference is preferable as an example. The global success criterion and the efficacy consistency criterion between target region and globally
#'
#' in superiority design:
#' \deqn{Z_a = \frac{\hat \delta_a}{\sqrt{Var(\hat \delta_a)}} > \Phi^{-1}(1 - \alpha)}
#' \deqn{\hat \delta_i > 0 \text{ for i = 1, 2, .., m}}
#'
#' in non-inferiority design:
#' \deqn{Z_a = \frac{\hat \delta_a + \Delta}{\sqrt{Var(\hat \delta_a)}} > \Phi^{-1}(1 - \alpha)}
#' \deqn{\hat \delta_i + \Delta_i > 0 \text{ for i = 1, 2, .., m}}
#'
#' in equivalence design:
#' \deqn{Z_{a_u} = \frac{\hat \delta_a + \Delta}{\sqrt{Var(\hat \delta_a)}} > \Phi^{-1}(1 - \alpha)\text{ and }Z_{a_l} = \frac{\hat \delta_a - \Delta}{\sqrt{Var(\hat \delta_a)}} < \Phi^{-1}(\alpha)}
#' \deqn{\hat \delta_i + \Delta_i > 0 \text{ and } \hat \delta_i - \Delta_i < 0 \text{ for i = 1, 2, .., m}}
#'
#' Where \eqn{\hat \delta} is the mean difference between treatment and control groups, and \eqn{\Delta} and \eqn{\Delta_i} are the non-inferiority or equivalence margins in global trial (\code{cut}) and each region (\code{cut_i}), respectively.
#'
#' @references
#' 1. Quan H, Li M, Chen J, et al. Assessment of Consistency of Treatment Effects in Multiregional Clinical Trials. Drug Information J. 2010;44(5):617-632. doi:10.1177/009286151004400509
#'
#' 2. Liao JJZ, Yu Z, Li Y. Sample size allocation in multiregional equivalence studies. Pharm Stat. 2018;17(5):570-577. doi:10.1002/pst.1871
#'
#' @export
#'
#' @examples
#' getPwr_Con_Super_JM2(
#'   delta_i = c(1, 0.8),
#'   sigma = 4, f_i = c(0.5, 0.5),
#'   alpha = 0.025, beta = NA, N = 200, r = 1, sim = FALSE
#' )
#'
#' # Global mean difference will be calculated based on delta_i and f_i.
#' # Non-inferiority margin in global trial and each region is 2.
#' # Global sample size will be calculated based on alpha and beta.
#' getPwr_Con_Noninf_JM2(
#'   delta_i = c(1, 0),
#'   sigma = 4, f_i = c(0.5, 0.5),
#'   cut = 2, alpha = 0.025, beta = 0.2, N = NA, r = 1, direct = -1, sim = FALSE
#' )
getPwr_Con_Super_JM2 <- function(delta_i, sigma, f_i, alpha = 0.025, beta = NA, N = NA, r = 1, sim = FALSE, nsim = 1000, seed = 0) {
  if (length(delta_i) != length(f_i)) {
    stop("The lengths of delta_i and f_i should be consistent, equal to the number of regions")
  }
  if (length(delta_i) == 1) {
    message("The number of regions in consideration is generally equal to or greater than 2.")
  }
  if (is.na(beta) & is.na(N)) {
    stop("Beta and N cannot both be NA.")
  }
  if (!is.na(beta) & (!is.na(N))) {
    warning("When beta is not NA, N will be automatically calculated.")
  }
  if (!is.logical(sim)) {
    stop("Parameter sim should be one of `TRUE` or `FALSE`.")
  }
  num <- length(delta_i)
  gr <- 2 + r + 1 / r
  delta_a <- sum(delta_i * f_i)
  if (!is.na(beta)) {
    N <- getN_Con_Super(delta = delta_a, sigma = sigma, alpha = alpha, beta = beta, N = NA, r = r)$N
  }
  Ni <- N * f_i
  if (!sim) {
    sei <- sqrt(gr * sigma^2 / Ni)
    ui <- delta_i / sei
    se <- sqrt(gr * sigma^2 / N)
    u <- delta_a / se
    M <- diag(num + 1)
    M[1, ] <- c(1, sqrt(f_i))
    M[, 1] <- c(1, sqrt(f_i))
    if (delta_a < 0) {
      ui <- (-1) * ui
      u <- (-1) * u
    }
    pwr1 <- 1 - pnorm(q = qnorm(1 - alpha), mean = u, sd = 1)
    pwr2 <- prod(1 - pnorm(q = 0, mean = ui, sd = 1))
    pwr3 <- mvtnorm::pmvnorm(lower = c(qnorm(1 - alpha), rep(0, num)), upper = rep(Inf, num + 1), mean = c(u, ui), sigma = M)
    pwr4 <- pwr3 / pwr1
    res <- data.frame(delta_a = delta_a, N = N, pwr1 = pwr1, pwr2 = pwr2, pwr3 = pwr3, pwr4 = pwr4)
    pwr_margin <- 1 - pnorm(q = 0, mean = ui, sd = 1)
    pwr_joint <- c()
    for (k in 1:num) {
      pwr_joint_ <- mvtnorm::pmvnorm(lower = c(qnorm(1 - alpha), 0), upper = c(Inf, Inf), mean = c(u, ui[k]), sigma = M[c(1, k + 1), c(1, k + 1)])
      pwr_joint <- c(pwr_joint, pwr_joint_)
    }
    pwr_condition <- pwr_joint / pwr1
    L <- list(overall = res, pwr_margin = pwr_margin, pwr_condition = pwr_condition, pwr_joint = pwr_joint)
  }
  if (sim) {
    da <- data.frame()
    di <- NULL
    set.seed(seed)
    seed1 <- sample(x = 1:1e8, size = num * nsim, replace = FALSE)
    for (j in 1:nsim) {
      Xt <- matrix(NA, nrow = N * r / (1 + r), ncol = num)
      Xc <- matrix(NA, nrow = N / (1 + r), ncol = num)
      for (k in 1:num) {
        seed2 <- seed1[((j - 1) * num + 1):(j * num)]
        set.seed(seed2[k])
        nt_k <- round(Ni[k] * r / (1 + r))
        nc_k <- round(Ni[k] / (1 + r))
        Xt[1:nt_k, k] <- rnorm(n = nt_k, mean = delta_i[k], sd = sigma)
        Xc[1:nc_k, k] <- rnorm(n = nc_k, mean = 0, sd = sigma)
      }
      ut <- colMeans(Xt, na.rm = TRUE)
      uc <- colMeans(Xc, na.rm = TRUE)
      zi <- ut - uc
      za <- (mean(as.numeric(Xt), na.rm = TRUE) - mean(as.numeric(Xc), na.rm = TRUE)) / sqrt(gr * sigma^2 / N)
      if (delta_a < 0) {
        zi <- (-1) * zi
        za <- (-1) * za
      }
      succ_a <- dplyr::if_else(za > qnorm(1 - alpha), 1, 0)
      succ_i <- dplyr::if_else(all(zi > 0), 1, 0)
      succ_i_ <- dplyr::if_else(zi > 0, 1, 0)
      da <- dplyr::bind_rows(da, c(succ_a = succ_a, succ_i = succ_i))
      di <- rbind(di, succ_i_)
    }
    res <- da %>%
      dplyr::summarise(pwr1 = mean(succ_a, na.rm = TRUE), pwr2 = mean(succ_i, na.rm = TRUE), pwr3 = mean(succ_a & succ_i, na.rm = TRUE), pwr4 = mean(succ_i[succ_a == 1], na.rm = TRUE), .groups = "keep") %>%
      dplyr::mutate(delta_a = delta_a, N = N) %>%
      dplyr::select(delta_a, N, pwr1, pwr2, pwr3, pwr4)
    pwr_margin <- colMeans(di)
    pwr_joint <- c()
    for (k in 1:num) {
      pwr_joint_ <- mean(da$succ_a & di[, k])
      pwr_joint <- c(pwr_joint, pwr_joint_)
    }
    pwr_condition <- colMeans(di[da$succ_a == 1, ])
    L <- list(overall = res, pwr_margin = pwr_margin, pwr_condition = pwr_condition, pwr_joint = pwr_joint)
  }
  return(L)
}


#' @rdname getPwr_Con_JM2
#' @export
getPwr_Con_Noninf_JM2 <- function(delta_i, sigma, f_i, cut, cut_i = NA, alpha = 0.025, beta = NA, N = NA, r = 1, direct = 1, sim = FALSE, nsim = 1000, seed = 0) {
  if (length(delta_i) != length(f_i)) {
    stop("The lengths of delta_i and f_i should be consistent, equal to the number of regions")
  }
  if (length(delta_i) == 1) {
    message("The number of regions in consideration is generally equal to or greater than 2.")
  }
  if (is.na(beta) & cut < 0) {
    warning("Parameter cut should be a positive value.")
  }
  if (!all(is.na(cut_i)) & (length(cut_i) != length(f_i))) {
    stop("The lengths of delta_i, f_i, and cut_i should be consistent, equal to the number of regions")
  }
  if (length(cut_i) == 1 & all(is.na(cut_i))) {
    cut_i <- rep(cut, length(f_i))
  }
  if (length(cut_i) >= 2) {
    cut_i[is.na(cut_i)] <- cut
  }
  if (is.na(beta) & is.na(N)) {
    stop("Beta and N cannot both be NA.")
  }
  if (!is.na(beta) & (!is.na(N))) {
    warning("When beta is not NA, N will be automatically calculated.")
  }
  if (!direct %in% c(-1, 1)) {
    stop("Parameter direct should be one of `1` or `-1`.")
  }
  if (!is.logical(sim)) {
    stop("Parameter sim should be one of `TRUE` or `FALSE`.")
  }
  num <- length(delta_i)
  gr <- 2 + r + 1 / r
  delta_a <- sum(delta_i * f_i)
  if (!is.na(beta)) {
    N <- getN_Con_Noninf(delta = delta_a, sigma = sigma, cut = cut, alpha = alpha, beta = beta, N = NA, r = r)$N
  }
  Ni <- N * f_i
  if (!sim) {
    sei <- sqrt(gr * sigma^2 / Ni)
    ui <- dplyr::if_else(direct == rep(1, num), (delta_i + cut_i) / sei, (delta_i - cut_i) / sei)
    se <- sqrt(gr * sigma^2 / N)
    u <- dplyr::if_else(direct == 1, (delta_a + cut) / se, (delta_a - cut) / se)
    M <- diag(num + 1)
    M[1, ] <- c(1, sqrt(f_i))
    M[, 1] <- c(1, sqrt(f_i))
    if (direct == -1) {
      ui <- (-1) * ui
      u <- (-1) * u
    }
    pwr1 <- 1 - pnorm(q = qnorm(1 - alpha), mean = u, sd = 1)
    pwr2 <- prod(1 - pnorm(q = 0, mean = ui, sd = 1))
    pwr3 <- mvtnorm::pmvnorm(lower = c(qnorm(1 - alpha), rep(0, num)), upper = rep(Inf, num + 1), mean = c(u, ui), sigma = M)
    pwr4 <- pwr3 / pwr1
    res <- data.frame(delta_a = delta_a, cut = cut, N = N, pwr1 = pwr1, pwr2 = pwr2, pwr3 = pwr3, pwr4 = pwr4)
    pwr_margin <- 1 - pnorm(q = 0, mean = ui, sd = 1)
    pwr_joint <- c()
    for (k in 1:num) {
      pwr_joint_ <- mvtnorm::pmvnorm(lower = c(qnorm(1 - alpha), 0), upper = c(Inf, Inf), mean = c(u, ui[k]), sigma = M[c(1, k + 1), c(1, k + 1)])
      pwr_joint <- c(pwr_joint, pwr_joint_)
    }
    pwr_condition <- pwr_joint / pwr1
    L <- list(overall = res, cut_i = cut_i, pwr_margin = pwr_margin, pwr_condition = pwr_condition, pwr_joint = pwr_joint)
  }
  if (sim) {
    da <- data.frame()
    di <- NULL
    set.seed(seed)
    seed1 <- sample(x = 1:1e8, size = num * nsim, replace = FALSE)
    for (j in 1:nsim) {
      Xt <- matrix(NA, nrow = N * r / (1 + r), ncol = num)
      Xc <- matrix(NA, nrow = N / (1 + r), ncol = num)
      for (k in 1:num) {
        seed2 <- seed1[((j - 1) * num + 1):(j * num)]
        set.seed(seed2[k])
        nt_k <- round(Ni[k] * r / (1 + r))
        nc_k <- round(Ni[k] / (1 + r))
        Xt[1:nt_k, k] <- rnorm(n = nt_k, mean = delta_i[k], sd = sigma)
        Xc[1:nc_k, k] <- rnorm(n = nc_k, mean = 0, sd = sigma)
      }
      ut <- colMeans(Xt, na.rm = TRUE)
      uc <- colMeans(Xc, na.rm = TRUE)
      zi <- dplyr::if_else(direct == rep(1, num), ut - uc + cut_i, ut - uc - cut_i)
      za <- dplyr::if_else(direct == 1, (mean(as.numeric(Xt), na.rm = TRUE) - mean(as.numeric(Xc), na.rm = TRUE) + cut) / sqrt(gr * sigma^2 / N), (mean(as.numeric(Xt), na.rm = TRUE) - mean(as.numeric(Xc), na.rm = TRUE) - cut) / sqrt(gr * sigma^2 / N))
      if (direct == -1) {
        zi <- (-1) * zi
        za <- (-1) * za
      }
      succ_a <- dplyr::if_else(za > qnorm(1 - alpha), 1, 0)
      succ_i <- dplyr::if_else(all(zi > 0), 1, 0)
      succ_i_ <- dplyr::if_else(zi > 0, 1, 0)
      da <- dplyr::bind_rows(da, c(succ_a = succ_a, succ_i = succ_i))
      di <- rbind(di, succ_i_)
    }
    res <- da %>%
      dplyr::summarise(pwr1 = mean(succ_a, na.rm = TRUE), pwr2 = mean(succ_i, na.rm = TRUE), pwr3 = mean(succ_a & succ_i, na.rm = TRUE), pwr4 = mean(succ_i[succ_a == 1], na.rm = TRUE), .groups = "keep") %>%
      dplyr::mutate(delta_a = delta_a, cut = cut, N = N) %>%
      dplyr::select(delta_a, cut, N, pwr1, pwr2, pwr3, pwr4)
    pwr_margin <- colMeans(di)
    pwr_joint <- c()
    for (k in 1:num) {
      pwr_joint_ <- mean(da$succ_a & di[, k])
      pwr_joint <- c(pwr_joint, pwr_joint_)
    }
    pwr_condition <- colMeans(di[da$succ_a == 1, ])
    L <- list(overall = res, cut_i = cut_i, pwr_margin = pwr_margin, pwr_condition = pwr_condition, pwr_joint = pwr_joint)
  }
  return(L)
}


#' @rdname getPwr_Con_JM2
#' @export
getPwr_Con_Equi_JM2 <- function(delta_i, sigma, f_i, cut, cut_i = NA, alpha = 0.025, beta = NA, N = NA, r = 1, sim = FALSE, nsim = 1000, seed = 0, maxN = 1e+06) {
  if (length(delta_i) != length(f_i)) {
    stop("The lengths of delta_i and f_i should be consistent, equal to the number of regions")
  }
  if (length(delta_i) == 1) {
    message("The number of regions in consideration is generally equal to or greater than 2.")
  }
  if (is.na(beta) & cut < 0) {
    warning("Parameter cut should be a positive value.")
  }
  if (!all(is.na(cut_i)) & (length(cut_i) != length(f_i))) {
    stop("The lengths of delta_i, f_i, and cut_i should be consistent, equal to the number of regions")
  }
  if (length(cut_i) == 1 & all(is.na(cut_i))) {
    cut_i <- rep(cut, length(f_i))
  }
  if (length(cut_i) >= 2) {
    cut_i[is.na(cut_i)] <- cut
  }
  if (is.na(beta) & is.na(N)) {
    stop("Beta and N cannot both be NA.")
  }
  if (!is.na(beta) & (!is.na(N))) {
    warning("When beta is not NA, N will be automatically calculated.")
  }
  if (!is.logical(sim)) {
    stop("Parameter sim should be one of `TRUE` or `FALSE`.")
  }
  num <- length(delta_i)
  gr <- 2 + r + 1 / r
  delta_a <- sum(delta_i * f_i)
  if (!is.na(beta)) {
    N <- getN_Con_Equi(delta = delta_a, sigma = sigma, cut = cut, alpha = alpha, beta = beta, N = NA, r = r, maxN = maxN)$N
  }
  Ni <- N * f_i
  if (!sim) {
    sei <- sqrt(gr * sigma^2 / Ni)
    ui1 <- (delta_i + cut_i) / sei
    ui2 <- (delta_i - cut_i) / sei
    se <- sqrt(gr * sigma^2 / N)
    u1 <- (delta_a + cut) / se
    u2 <- (delta_a - cut) / se
    M <- diag(num + 1)
    M[1, ] <- c(1, sqrt(f_i))
    M[, 1] <- c(1, sqrt(f_i))
    M1 <- combine(M, M)
    pwr1 <- mvtnorm::pmvnorm(lower = c(qnorm(1 - alpha), -Inf), upper = c(Inf, -qnorm(1 - alpha)), mean = c(u1, u2), corr = matrix(1, nrow = 2, ncol = 2))
    pwr_margin <- c()
    for (k in 1:num) {
      pwr_margin_ <- mvtnorm::pmvnorm(lower = c(0, -Inf), upper = c(Inf, 0), mean = c(ui1[k], ui2[k]), corr = matrix(1, nrow = 2, ncol = 2))
      pwr_margin <- c(pwr_margin, pwr_margin_)
    }
    pwr2 <- prod(pwr_margin)
    pwr3 <- mvtnorm::pmvnorm(lower = c(qnorm(1 - alpha), -Inf, rep(c(0, -Inf), num)), upper = c(Inf, -qnorm(1 - alpha), rep(c(Inf, 0), num)), mean = c(u1, u2, combine(ui1, ui2)), sigma = M1)
    pwr4 <- pwr3 / pwr1
    res <- data.frame(delta_a = delta_a, cut = cut, N = N, pwr1 = pwr1, pwr2 = pwr2, pwr3 = pwr3, pwr4 = pwr4)
    pwr_joint <- c()
    for (k in 1:num) {
      M <- diag(2)
      M[1, ] <- c(1, sqrt(f_i[k]))
      M[, 1] <- c(1, sqrt(f_i[k]))
      M2 <- combine(M, M)
      pwr_joint_ <- mvtnorm::pmvnorm(lower = c(qnorm(1 - alpha), -Inf, 0, -Inf), upper = c(Inf, -qnorm(1 - alpha), Inf, 0), mean = c(u1, u2, ui1[k], ui2[k]), sigma = M2)
      pwr_joint <- c(pwr_joint, pwr_joint_)
    }
    pwr_condition <- pwr_joint / pwr1
    L <- list(overall = res, cut_i = cut_i, pwr_margin = pwr_margin, pwr_condition = pwr_condition, pwr_joint = pwr_joint)
  }
  if (sim) {
    da <- data.frame()
    di <- NULL
    set.seed(seed)
    seed1 <- sample(x = 1:1e8, size = num * nsim, replace = FALSE)
    for (j in 1:nsim) {
      Xt <- matrix(NA, nrow = N * r / (1 + r), ncol = num)
      Xc <- matrix(NA, nrow = N / (1 + r), ncol = num)
      for (k in 1:num) {
        seed2 <- seed1[((j - 1) * num + 1):(j * num)]
        set.seed(seed2[k])
        nt_k <- round(Ni[k] * r / (1 + r))
        nc_k <- round(Ni[k] / (1 + r))
        Xt[1:nt_k, k] <- rnorm(n = nt_k, mean = delta_i[k], sd = sigma)
        Xc[1:nc_k, k] <- rnorm(n = nc_k, mean = 0, sd = sigma)
      }
      ut <- colMeans(Xt, na.rm = TRUE)
      uc <- colMeans(Xc, na.rm = TRUE)
      zi1 <- ut - uc + cut_i
      zi2 <- ut - uc - cut_i
      za1 <- (mean(as.numeric(Xt), na.rm = TRUE) - mean(as.numeric(Xc), na.rm = TRUE) + cut) / sqrt(gr * sigma^2 / N)
      za2 <- (mean(as.numeric(Xt), na.rm = TRUE) - mean(as.numeric(Xc), na.rm = TRUE) - cut) / sqrt(gr * sigma^2 / N)
      succ_a <- dplyr::if_else(za1 > qnorm(1 - alpha) & za2 < -qnorm(1 - alpha), 1, 0)
      succ_i <- dplyr::if_else(all(zi1 > 0) & all(zi2 < 0), 1, 0)
      succ_i_ <- dplyr::if_else(zi1 > 0 & zi2 < 0, 1, 0)
      da <- dplyr::bind_rows(da, c(succ_a = succ_a, succ_i = succ_i))
      di <- rbind(di, succ_i_)
    }
    res <- da %>%
      dplyr::summarise(pwr1 = mean(succ_a, na.rm = TRUE), pwr2 = mean(succ_i, na.rm = TRUE), pwr3 = mean(succ_a & succ_i, na.rm = TRUE), pwr4 = mean(succ_i[succ_a == 1], na.rm = TRUE), .groups = "keep") %>%
      dplyr::mutate(delta_a = delta_a, cut = cut, N = N) %>%
      dplyr::select(delta_a, cut, N, pwr1, pwr2, pwr3, pwr4)
    pwr_margin <- colMeans(di)
    pwr_joint <- c()
    for (k in 1:num) {
      pwr_joint_ <- mean(da$succ_a & di[, k])
      pwr_joint <- c(pwr_joint, pwr_joint_)
    }
    pwr_condition <- colMeans(di[da$succ_a == 1, ])
    L <- list(overall = res, cut_i = cut_i, pwr_margin = pwr_margin, pwr_condition = pwr_condition, pwr_joint = pwr_joint)
  }
  return(L)
}
