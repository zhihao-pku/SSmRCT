#' Regional sample size allocation using Japan's Method 1 for count endpoints
#'
#' Based on Japan's Method 1, given the global sample size and marginal probability (power) of efficacy consistency between target region and globally, calculate the required sample size allocated to the target region, in clinical trials using superiority, non-inferiority, and equivalence designs with count endpoints.
#'
#' @rdname getN_Count_JM1
#'
#' @name getN_Count_JM1
#'
#' @param delta_a log(RR) between treatment and control groups globally.
#' @param delta_j log(RR) between treatment and control groups for target region.
#' @param lambda0_a Baseline hazard of control group globally.
#' @param lambda0_j Baseline hazard of control group for target region.
#' @param t Average exposure time.
#' @param k The over-dispersion parameter for negative binomial distribution, which is 0 for poisson distribution. The default value is 0.
#' @param pi Proportion of global efficacy to retain. The default value is 0.5, which means retaining half of the efficacy.
#' @param cut A positive value for non-inferiority or equivalence margin. For example, if the non-inferiority margin for RR is 0.6, then \code{cut = -log(0.6)}. If the non-inferiority margin for RR is 1.3, then \code{cut = log(1.3)}.
#' @param alpha One-sided type I error rate for global success, which is used to calculate global sample size only when \code{N} is \code{NA}. The default value is 0.025.
#' @param beta Type II error rate for global success, which is used to calculate global sample size only when \code{N} is \code{NA}.
#' @param beta1 Type II error rate for efficacy consistency between target region and globally. The default value is 0.2.
#' @param N Global sample size. When \code{N} is \code{NA} and \code{alpha} and \code{beta} are not \code{NAs}, \code{N} will be calculated automatically.
#' @param r Ratio of the sample sizes of the treatment group to the control group. The default value is 1.
#' @param direct \code{direct = 1} indicates that a larger RR is preferable, while \code{direct = -1} indicates that a smaller RR is preferable.
#' @param maxN Maximum possible sample size (\code{N}) in equivalence design. Default value is 1e6.
#'
#' @return A data frame where \code{f} is required proportion of sample size allocated to the target region, and \code{Nj} is required sample size for the target region, calculated as \code{Nj = N * f}.
#'
#' @details
#' The global success criterion and the efficacy consistency criterion between target region and globally could be found in \code{\link{getPwr_Count_Super_JM1}}.
#'
#' @references
#' 1. Quan H, Li M, Chen J, et al. Assessment of Consistency of Treatment Effects in Multiregional Clinical Trials. Drug Information J. 2010;44(5):617-632. doi:10.1177/009286151004400509
#'
#' 2. Liao JJZ, Yu Z, Li Y. Sample size allocation in multiregional equivalence studies. Pharm Stat. 2018;17(5):570-577. doi:10.1002/pst.1871
#'
#' @seealso
#' \code{\link{getPwr_Count_Super_JM1}}, \code{\link{getN_Count_Super}}.
#'
#' @export
#'
#' @examples
#' getN_Count_Super_JM1(
#'   delta_a = log(1.4),
#'   delta_j = log(1.3),
#'   lambda0_a = 0.1, lambda0_j = 0.1, t = 5, k = 0, pi = 0.5, beta1 = 0.2,
#'   N = 300, r = 1
#' )
#'
#' getN_Count_Noninf_JM1(
#'   delta_a = log(1.0),
#'   delta_j = log(1.1),
#'   lambda0_a = 0.1, lambda0_j = 0.1, t = 5, k = 0, pi = 0.5, cut = log(1.3),
#'   alpha = 0.025, beta = 0.2, beta1 = 0.2, N = NA, r = 1, direct = -1
#' )
getN_Count_Super_JM1 <- function(delta_a, delta_j, lambda0_a, lambda0_j, t, k = 0, pi = 0.5, alpha = NA, beta = NA, beta1 = 0.2, N = NA, r = 1) {
  eg <- as.data.frame(expand.grid(delta_a = delta_a, delta_j = delta_j, lambda0_a = lambda0_a, lambda0_j = lambda0_j, t = t, k = k, pi = pi, alpha = alpha, beta = beta, beta1 = beta1, N = N, r = r, stringsAsFactors = FALSE))
  res <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta_a <- R$delta_a
    delta_j <- R$delta_j
    lambda0_a <- R$lambda0_a
    lambda0_j <- R$lambda0_j
    t <- R$t
    k <- R$k
    pi <- R$pi
    alpha <- R$alpha
    beta <- R$beta
    beta1 <- R$beta1
    N <- R$N
    r <- R$r
    if ((is.na(alpha) | is.na(beta)) & is.na(N)) {
      stop("alpha and beta, and N cannot be NA simultaneously.")
    }
    if ((!is.na(alpha) | (!is.na(beta))) & (!is.na(N))) {
      stop("Set either alpha and beta, or N to NA.")
    }
    lambda1_j <- exp(delta_j) * lambda0_j
    lambda1_a <- exp(delta_a) * lambda0_a
    sigma1_j <- sqrt(1 / lambda1_j / t + k)
    sigma0_j <- sqrt(1 / lambda0_j / t + k)
    sigma1_a <- sqrt(1 / lambda1_a / t + k)
    sigma0_a <- sqrt(1 / lambda0_a / t + k)
    if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
      N <- getN_Count_Super(delta = delta_a, lambda0 = lambda0_a, t = t, k = k, alpha = alpha, beta = beta, N = NA, r = r)$N
    }
    getPwr <- function(f) {
      Nj <- N * f
      var_j <- sigma1_j^2 / (r * Nj / (1 + r)) + sigma0_j^2 / (Nj / (1 + r))
      var_a <- sigma1_a^2 / (r * N / (1 + r)) + sigma0_a^2 / (N / (1 + r))
      sej <- sqrt(var_j + pi^2 * var_a - 2 * pi * sqrt(f) * sqrt(var_j * var_a))
      uj <- (delta_j - pi * delta_a) / sej
      uj <- if_else(delta_a < 0, (-1) * uj, uj)
      pmvnorm(lower = c(0), upper = c(Inf), mean = c(uj), sigma = 1)
    }
    f <- tryCatch(
      {
        uniroot(f = function(f) getPwr(f) - (1 - beta1), interval = c(1e-06, 1 - 1e-06))$root
      },
      error = function(e) {
        warning("Cannot find an f value that meets the conditions.")
        NA
      }
    )
    delta_nj <- pwr <- NA
    if (!is.na(f)) {
      delta_nj <- (delta_a - f * delta_j) / (1 - f)
      pwr <- getPwr(f)
    }
    data.frame(delta_a, delta_j, delta_nj, lambda0_a, lambda0_j, t, k, pi, alpha, beta, beta1, N, r, pwr, f, Nj = N * f)
  }, .options = furrr_options(seed = TRUE))
  return(res)
}

#' @rdname getN_Count_JM1
#' @export
getN_Count_Noninf_JM1 <- function(delta_a, delta_j, lambda0_a, lambda0_j, t, k = 0, pi = 0.5, cut, alpha = NA, beta = NA, beta1 = 0.2, N = NA, r = 1, direct = 1) {
  eg <- as.data.frame(expand.grid(delta_a = delta_a, delta_j = delta_j, lambda0_a = lambda0_a, lambda0_j = lambda0_j, t = t, k = k, pi = pi, cut = cut, alpha = alpha, beta = beta, beta1 = beta1, N = N, r = r, stringsAsFactors = FALSE))
  res <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta_a <- R$delta_a
    delta_j <- R$delta_j
    lambda0_a <- R$lambda0_a
    lambda0_j <- R$lambda0_j
    t <- R$t
    k <- R$k
    pi <- R$pi
    cut <- R$cut
    alpha <- R$alpha
    beta <- R$beta
    beta1 <- R$beta1
    N <- R$N
    r <- R$r
    if ((is.na(alpha) | is.na(beta)) & is.na(N)) {
      stop("alpha and beta, and N cannot be NA simultaneously.")
    }
    if ((!is.na(alpha) | (!is.na(beta))) & (!is.na(N))) {
      stop("Set either alpha and beta, or N to NA.")
    }
    lambda1_j <- exp(delta_j) * lambda0_j
    lambda1_a <- exp(delta_a) * lambda0_a
    sigma1_j <- sqrt(1 / lambda1_j / t + k)
    sigma0_j <- sqrt(1 / lambda0_j / t + k)
    sigma1_a <- sqrt(1 / lambda1_a / t + k)
    sigma0_a <- sqrt(1 / lambda0_a / t + k)
    if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
      N <- getN_Count_Noninf(delta = delta_a, cut = cut, lambda0 = lambda0_a, t = t, k = k, alpha = alpha, beta = beta, N = NA, r = r)$N
    }
    getPwr <- function(f) {
      Nj <- N * f
      var_j <- sigma1_j^2 / (r * Nj / (1 + r)) + sigma0_j^2 / (Nj / (1 + r))
      var_a <- sigma1_a^2 / (r * N / (1 + r)) + sigma0_a^2 / (N / (1 + r))
      sej <- sqrt(var_j + var_a - 2 * sqrt(f) * sqrt(var_j * var_a))
      uj <- if_else(direct == 1, (delta_j - delta_a + pi * cut) / sej, (delta_j - delta_a - pi * cut) / sej)
      uj <- if_else(direct == -1, (-1) * uj, uj)
      pmvnorm(lower = c(0), upper = c(Inf), mean = c(uj), sigma = 1)
    }
    f <- tryCatch(
      {
        uniroot(f = function(f) getPwr(f) - (1 - beta1), interval = c(1e-06, 1 - 1e-06))$root
      },
      error = function(e) {
        warning("Cannot find an f value that meets the conditions.")
        NA
      }
    )
    delta_nj <- pwr <- NA
    if (!is.na(f)) {
      delta_nj <- (delta_a - f * delta_j) / (1 - f)
      pwr <- getPwr(f)
    }
    data.frame(delta_a, delta_j, delta_nj, lambda0_a, lambda0_j, t, k, pi, cut, alpha, beta, beta1, N, r, direct, pwr, f, Nj = N * f)
  }, .options = furrr_options(seed = TRUE))
  return(res)
}

#' @rdname getN_Count_JM1
#' @export
getN_Count_Equi_JM1 <- function(delta_a, delta_j, lambda0_a, lambda0_j, t, k = 0, pi = 0.5, cut, alpha = NA, beta = NA, beta1 = 0.2, N = NA, r = 1, maxN = 1e6) {
  eg <- as.data.frame(expand.grid(delta_a = delta_a, delta_j = delta_j, lambda0_a = lambda0_a, lambda0_j = lambda0_j, t = t, k = k, pi = pi, cut = cut, alpha = alpha, beta = beta, beta1 = beta1, N = N, r = r, stringsAsFactors = FALSE))
  res <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta_a <- R$delta_a
    delta_j <- R$delta_j
    lambda0_a <- R$lambda0_a
    lambda0_j <- R$lambda0_j
    t <- R$t
    k <- R$k
    pi <- R$pi
    cut <- R$cut
    alpha <- R$alpha
    beta <- R$beta
    beta1 <- R$beta1
    N <- R$N
    r <- R$r
    if ((is.na(alpha) | is.na(beta)) & is.na(N)) {
      stop("alpha and beta, and N cannot be NA simultaneously.")
    }
    if ((!is.na(alpha) | (!is.na(beta))) & (!is.na(N))) {
      stop("Set either alpha and beta, or N to NA.")
    }
    lambda1_j <- exp(delta_j) * lambda0_j
    lambda1_a <- exp(delta_a) * lambda0_a
    sigma1_j <- sqrt(1 / lambda1_j / t + k)
    sigma0_j <- sqrt(1 / lambda0_j / t + k)
    sigma1_a <- sqrt(1 / lambda1_a / t + k)
    sigma0_a <- sqrt(1 / lambda0_a / t + k)
    if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
      N <- getN_Count_Equi(delta = delta_a, cut = cut, lambda0 = lambda0_a, t = t, k = k, alpha = alpha, beta = beta, N = NA, r = r, maxN = maxN)$N
    }
    getPwr <- function(f) {
      Nj <- N * f
      var_j <- sigma1_j^2 / (r * Nj / (1 + r)) + sigma0_j^2 / (Nj / (1 + r))
      var_a <- sigma1_a^2 / (r * N / (1 + r)) + sigma0_a^2 / (N / (1 + r))
      sej <- sqrt(var_j + var_a - 2 * sqrt(f) * sqrt(var_j * var_a))
      uj1 <- (delta_j - delta_a + pi * cut) / sej
      uj2 <- (delta_j - delta_a - pi * cut) / sej
      pmvnorm(lower = c(0, -Inf), upper = c(Inf, 0), mean = c(uj1, uj2), sigma = matrix(1, nrow = 2, ncol = 2))
    }
    f <- tryCatch(
      {
        uniroot(f = function(f) getPwr(f) - (1 - beta1), interval = c(1e-06, 1 - 1e-06))$root
      },
      error = function(e) {
        warning("Cannot find an f value that meets the conditions.")
        NA
      }
    )
    delta_nj <- pwr <- NA
    if (!is.na(f)) {
      delta_nj <- (delta_a - f * delta_j) / (1 - f)
      pwr <- getPwr(f)
    }
    data.frame(delta_a, delta_j, delta_nj, lambda0_a, lambda0_j, t, k, pi, cut, alpha, beta, beta1, N, r, pwr, f, Nj = N * f)
  }, .options = furrr_options(seed = TRUE))
  return(res)
}
