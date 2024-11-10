#' Sample size and power for survival endpoints
#'
#' Calculating sample size when given power or power when given sample size for survival endpoints.
#'
#' @rdname getN_Surv
#'
#' @name getN_Surv
#'
#' @param delta log(HR) between treatment and control groups.
#' @param alpha One-sided type I error rate. The default value is 0.025.
#' @param cut A positive value for non-inferiority or equivalence margin. For example, if the non-inferiority margin for HR is 0.6, then \code{cut = -log(0.6)}. If the non-inferiority margin for HR is 1.3, then \code{cut = log(1.3)}.
#' @param beta Type II error rate. When \code{beta} is \code{NA} and \code{N} is not \code{NA}, the power will be returned.
#' @param N Sample size. When \code{N} is \code{NA} and \code{beta} is not \code{NA}, the sample size will be returned.
#' @param r Ratio of the sample sizes of the treatment group to the control group. The default value is 1.
#' @param criterion If \code{criterion = 1}, the success criterion defined on the log(HR) scale will be used. If \code{criterion = 2}, the success criterion defined on the HR scale will be used. See \code{details} for more information.
#' @param direct \code{direct = 1} indicates that a larger HR is preferable, while \code{direct = -1} indicates that a smaller HR is preferable.
#' @param maxN Maximum possible sample size(\code{N}) in equivalence design. Default value is 1e6.

#' @return A data frame containing input parameters and returned sample size or power.
#'
#' @details
#' Taking the larger HR is preferable as an example. Sample size calculation is based on the following Z test, with the success criterion:
#'
#' \code{criterion = 1} in superiority design:
#' \deqn{Z = \frac{\hat \delta}{\sqrt{\frac{1}{N / (r + 1)} + \frac{1}{Nr / (r + 1)}}} > \Phi^{-1}(1 - \alpha)}
#'
#' \code{criterion = 2} in superiority design:
#' \deqn{Z = \frac{1 - e^{\hat \delta}}{e^{\hat \delta}\sqrt{\frac{1}{N / (r + 1)} + \frac{1}{Nr / (r + 1)}}} > \Phi^{-1}(\alpha)}
#'
#' \code{criterion = 1} in non-inferiority design:
#' \deqn{Z = \frac{\hat \delta + \Delta}{\sqrt{\frac{1}{N / (r + 1)} + \frac{1}{Nr / (r + 1)}}} > \Phi^{-1}(1 - \alpha)}
#'
#' \code{criterion = 2} in non-inferiority design:
#' \deqn{Z = \frac{\frac{1}{e^{\Delta}} - e^{\hat \delta}}{e^{\hat \delta}\sqrt{\frac{1}{N / (r + 1)} + \frac{1}{Nr / (r + 1)}}} > \Phi^{-1}(\alpha)}
#'
#' \code{criterion = 1} in equivalence design:
#' \deqn{Z_u = \frac{\hat \delta + \Delta}{\sqrt{\frac{1}{N / (r + 1)} + \frac{1}{Nr / (r + 1)}}} > \Phi^{-1}(1 - \alpha) \text{ and } Z_l = \frac{\hat \delta - \Delta}{\sqrt{\frac{1}{N / (r + 1)} + \frac{1}{Nr / (r + 1)}}} > \Phi^{-1}(\alpha)}
#'
#' \code{criterion = 2} in equivalence design:
#' \deqn{Z_l = \frac{\frac{1}{e^{\Delta}} - e^{\hat \delta}}{e^{\hat \delta}\sqrt{\frac{1}{N / (r + 1)} + \frac{1}{Nr / (r + 1)}}} > \Phi^{-1}(\alpha) \text{ and } Z_u = \frac{e^{\Delta} - e^{\hat \delta}}{e^{\hat \delta}\sqrt{\frac{1}{N / (r + 1)} + \frac{1}{Nr / (r + 1)}}} > \Phi^{-1}(1 - \alpha)}
#'
#' Where \eqn{\hat \delta = log({\hat{HR}})} between treatment and control groups, and \eqn{\Delta} is the non-inferiority or equivalence margin (\code{cut}).
#'
#' For \code{criterion = 2}, by delta method, \eqn{e^{\hat \delta}} approximately follows a distribution of \eqn{N(e^{\delta}, (e^{\delta}\sqrt{\frac{1}{N / (r + 1)} + \frac{1}{Nr / (r + 1)}})^2)}
#'
#' @export
#'
#' @examples
#' (v <- getN_Surv_Super(
#'   delta = log(1.2),
#'   alpha = 0.025, beta = 0.2, N = NA, r = 1, criterion = 1
#' ))
#' getN_Surv_Super(
#'   delta = log(1.2),
#'   alpha = 0.025, beta = NA, N = v$N, r = 1, criterion = 1
#' )
#'
#' (v <- getN_Surv_Noninf(
#'   delta = log(1.1),
#'   cut = log(1.3),
#'   alpha = 0.025, beta = 0.2, N = NA, r = 1, criterion = 2, direct = -1
#' ))
#' getN_Surv_Noninf(
#'   delta = log(1.1),
#'   cut = log(1.3),
#'   alpha = 0.025, beta = NA, N = v$N, r = 1, criterion = 2, direct = -1
#' )
getN_Surv_Super <- function(delta, alpha = 0.025, beta = NA, N = NA, r = 1, criterion = 1) {
  eg <- as.data.frame(expand.grid(delta = delta, alpha = alpha, beta = beta, N = N, r = r, criterion = criterion, stringsAsFactors = FALSE))
  res <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta <- R$delta
    alpha <- R$alpha
    beta <- R$beta
    N <- R$N
    r <- R$r
    criterion <- R$criterion
    if (is.na(beta) & is.na(N)) {
      stop("beta and N cannot be NA simultaneously.")
    }
    if (!is.na(beta) & (!is.na(N))) {
      stop("Set one of beta and N to NA.")
    }
    getPwr <- function(n0) {
      n1 <- r * n0
      if (criterion == 1) {
        z <- delta / sqrt(1 / n1 + 1 / n0)
        pwr <- if_else(delta > 0, 1 - pnorm(q = qnorm(1 - alpha), mean = z, sd = 1), pnorm(q = -qnorm(1 - alpha), mean = z, sd = 1))
      }
      if (criterion == 2) {
        z <- (1 - exp(delta)) / (sqrt(1 / n1 + 1 / n0) * exp(delta))
        pwr <- if_else(delta < 0, 1 - pnorm(q = qnorm(1 - alpha), mean = z, sd = 1), pnorm(q = -qnorm(1 - alpha), mean = z, sd = 1))
      }
      pwr
    }
    if (is.na(N) & (!is.na(beta))) {
      if (criterion == 1) {
        n0 <- (qnorm(1 - alpha) + qnorm(1 - beta))^2 * (1 + 1 / r) / delta^2
      }
      if (criterion == 2) {
        n0 <- (qnorm(1 - alpha) + qnorm(1 - beta))^2 * (1 + 1 / r) * exp(delta)^2 / (1 - exp(delta))^2
      }
      n0 <- ceiling(n0)
      n1 <- r * n0
      N <- n1 + n0
      pwr <- getPwr(n0)
      df <- data.frame(delta, alpha, beta, pwr, r, N, n1, n0, criterion)
    }
    if (is.na(beta) & (!is.na(N))) {
      n0 <- N / (r + 1)
      n1 <- r * n0
      pwr <- getPwr(n0)
      df <- data.frame(delta, alpha, beta, pwr, r, N, n1, n0, criterion)
    }
    df
  }, .options = furrr_options(seed = TRUE))
  return(res)
}

#' @rdname getN_Surv
#' @export
getN_Surv_Noninf <- function(delta, cut, alpha = 0.025, beta = NA, N = NA, r = 1, criterion = 1, direct = 1) {
  eg <- as.data.frame(expand.grid(delta = delta, cut = cut, alpha = alpha, beta = beta, N = N, r = r, criterion = criterion, stringsAsFactors = FALSE))
  res <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta <- R$delta
    cut <- R$cut
    alpha <- R$alpha
    beta <- R$beta
    N <- R$N
    r <- R$r
    criterion <- R$criterion
    if (is.na(beta) & is.na(N)) {
      stop("beta and N cannot be NA simultaneously.")
    }
    if (!is.na(beta) & (!is.na(N))) {
      stop("Set one of beta and N to NA.")
    }
    getPwr <- function(n0) {
      n1 <- r * n0
      if (criterion == 1) {
        z <- if_else(direct == 1, (delta + cut) / sqrt(1 / n1 + 1 / n0), (delta - cut) / sqrt(1 / n1 + 1 / n0))
        pwr <- if_else(direct == 1, 1 - pnorm(q = qnorm(1 - alpha), mean = z, sd = 1), pnorm(q = -qnorm(1 - alpha), mean = z, sd = 1))
      }
      if (criterion == 2) {
        z <- if_else(direct == 1, -(1 / exp(cut) - exp(delta)) / (sqrt(1 / n1 + 1 / n0) * exp(delta)), -(exp(cut) - exp(delta)) / (sqrt(1 / n1 + 1 / n0) * exp(delta)))
        pwr <- if_else(direct == 1, 1 - pnorm(q = qnorm(1 - alpha), mean = z, sd = 1), pnorm(q = -qnorm(1 - alpha), mean = z, sd = 1))
      }
      pwr
    }
    if (is.na(N) & (!is.na(beta))) {
      if (criterion == 1) {
        n0 <- if_else(direct == 1, (qnorm(1 - alpha) + qnorm(1 - beta))^2 * (1 + 1 / r) / (delta + cut)^2, (qnorm(1 - alpha) + qnorm(1 - beta))^2 * (1 + 1 / r) / (delta - cut)^2)
      }
      if (criterion == 2) {
        n0 <- if_else(direct == 1, (qnorm(1 - alpha) + qnorm(1 - beta))^2 * (1 + 1 / r) * exp(delta)^2 / (1 / exp(cut) - exp(delta))^2, (qnorm(1 - alpha) + qnorm(1 - beta))^2 * (1 + 1 / r) * exp(delta)^2 / (exp(cut) - exp(delta))^2)
      }
      n0 <- ceiling(n0)
      n1 <- r * n0
      N <- n1 + n0
      pwr <- getPwr(n0)
      df <- data.frame(delta, cut, alpha, beta, pwr, r, N, n1, n0, criterion, direct)
    }
    if (is.na(beta)) {
      n0 <- N / (r + 1)
      n1 <- r * n0
      pwr <- getPwr(n0)
      df <- data.frame(delta, cut, alpha, beta, pwr, r, N, n1, n0, criterion, direct)
    }
    df
  }, .options = furrr_options(seed = TRUE))
  return(res)
}

#' @rdname getN_Surv
#' @export
getN_Surv_Equi <- function(delta, cut, alpha = 0.025, beta = NA, N = NA, r = 1, criterion = 1, maxN = 1e+06) {
  eg <- as.data.frame(expand.grid(delta = delta, cut = cut, alpha = alpha, beta = beta, N = N, r = r, criterion = criterion, stringsAsFactors = FALSE))
  res <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta <- R$delta
    cut <- R$cut
    alpha <- R$alpha
    beta <- R$beta
    N <- R$N
    r <- R$r
    criterion <- R$criterion
    if (is.na(beta) & is.na(N)) {
      stop("beta and N cannot be NA simultaneously.")
    }
    if (!is.na(beta) & (!is.na(N))) {
      stop("Set one of beta and N to NA.")
    }
    getPwr <- function(n0) {
      n1 <- r * n0
      if (criterion == 1) {
        z1 <- (delta + cut) / sqrt(1 / n1 + 1 / n0)
        z2 <- (delta - cut) / sqrt(1 / n1 + 1 / n0)
        pwr <- pmvnorm(lower = c(qnorm(1 - alpha), -Inf), upper = c(Inf, -qnorm(1 - alpha)), mean = c(z1, z2), sigma = matrix(1, nrow = 2, ncol = 2))
      }
      if (criterion == 2) {
        z1 <- -(1 / exp(cut) - exp(delta)) / (sqrt(1 / n1 + 1 / n0) * exp(delta))
        z2 <- -(exp(cut) - exp(delta)) / (sqrt(1 / n1 + 1 / n0) * exp(delta))
        pwr <- pmvnorm(lower = c(qnorm(1 - alpha), -Inf), upper = c(Inf, -qnorm(1 - alpha)), mean = c(z1, z2), sigma = matrix(1, nrow = 2, ncol = 2))
      }
      pwr
    }
    if (is.na(N) & (!is.na(beta))) {
      n0 <- tryCatch(
        {
          uniroot(f = function(n0) {
            getPwr(n0) - (1 - beta)
          }, interval = c(1e-06, maxN / (r + 1)))$root
        },
        error = function(e) {
          warning(sprintf("The calculated sample size for control group exceeds %s", maxN))
          NA
        }
      )
      n0 <- ceiling(n0)
      n1 <- r * n0
      N <- n1 + n0
      pwr <- getPwr(n0)
      df <- data.frame(delta, cut, alpha, beta, pwr, r, N, n1, n0, criterion)
    }
    if (is.na(beta) & (!is.na(N))) {
      n0 <- N / (r + 1)
      n1 <- r * n0
      pwr <- getPwr(n0)
      df <- data.frame(delta, cut, alpha, beta, pwr, r, N, n1, n0, criterion)
    }
    df
  }, .options = furrr_options(seed = TRUE))
  return(res)
}
