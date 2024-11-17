#' Regional sample size allocation using Japan's Method 1 for count endpoints
#'
#' Based on Japan's Method 1, given the global sample size and marginal probability (power) of efficacy consistency between target region and globally, calculate the required sample size allocated to the target region, in clinical trials using superiority, non-inferiority, and equivalence designs with count endpoints.
#'
#' @rdname getN_Count_JM1
#'
#' @name getN_Count_JM1
#'
#' @param delta_a A vector. log(RR) between treatment and control groups globally.
#' @param delta_j A vector. log(RR) between treatment and control groups for target region.
#' @param lambda0_a A vector. Baseline hazard of control group globally.
#' @param lambda0_j A vector. Baseline hazard of control group for target region.
#' @param t_a A vector. Average exposure time globally.
#' @param t_j A vector. Average exposure time for target region..
#' @param k A vector. The over-dispersion parameter (k > 0) for negative binomial distribution, which is 0 for poisson distribution. Default value is 0.
#' @param pi A vector. Proportion of global efficacy to retain. Default value is 0.5, which means retaining half of the efficacy.
#' @param cut A vector. Positive value for non-inferiority or equivalence margin. For example, if the non-inferiority margin for RR is 0.6, then \code{cut = -log(0.6)}. If the non-inferiority margin for RR is 1.3, then \code{cut = log(1.3)}.
#' @param alpha A vector. One-sided type I error rate for global success, which is used to calculate global sample size only when \code{N} is \code{NA}. Default value is 0.025.
#' @param beta A vector. Type II error rate for global success, which is used to calculate global sample size only when \code{N} is \code{NA}.
#' @param beta1 A vector. Type II error rate for efficacy consistency between target region and globally. Default value is 0.2.
#' @param N A vector. Global sample size. When \code{N} is \code{NA} and \code{alpha} and \code{beta} are not \code{NA}, \code{N} will be calculated automatically.
#' @param r A vector. Ratio of sample sizes of treatment group to control group. Default value is 1.
#' @param direct \code{direct = 1} indicates that a larger RR is preferable, while \code{direct = -1} indicates that a smaller RR is preferable.
#' @param maxN Maximum possible sample size (\code{N}) in equivalence design. Default value is 1e+06.
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
#'   lambda0_a = 0.1, lambda0_j = 0.1, t_a = 5, t_j = 5, k = 0, pi = 0.5, beta1 = 0.2,
#'   N = 300, r = 1
#' )
#'
#' # Global sample size will be calculated based on alpha and beta.
#' getN_Count_Noninf_JM1(
#'   delta_a = log(1.0),
#'   delta_j = log(1.1),
#'   lambda0_a = 0.1, lambda0_j = 0.1, t_a = 5, t_j = 5, k = 0, pi = 0.5, cut = log(1.3),
#'   alpha = 0.025, beta = 0.2, beta1 = 0.2, N = NA, r = 1, direct = -1
#' )
getN_Count_Super_JM1 <- function(delta_a, delta_j, lambda0_a, lambda0_j, t_a, t_j, k = 0, pi = 0.5, alpha = NA, beta = NA, beta1 = 0.2, N = NA, r = 1) {
  eg <- as.data.frame(expand.grid(delta_a = delta_a, delta_j = delta_j, lambda0_a = lambda0_a, lambda0_j = lambda0_j, t_a = t_a, t_j = t_j, k = k, pi = pi, alpha = alpha, beta = beta, beta1 = beta1, N = N, r = r, stringsAsFactors = FALSE))
  res <- furrr::future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta_a <- R$delta_a
    delta_j <- R$delta_j
    lambda0_a <- R$lambda0_a
    lambda0_j <- R$lambda0_j
    t_a <- R$t_a
    t_j <- R$t_j
    k <- R$k
    pi <- R$pi
    alpha <- R$alpha
    beta <- R$beta
    beta1 <- R$beta1
    N <- R$N
    r <- R$r
    if (pi < 0 | pi > 1) {
      warning("Parameter pi generally is between 0 and 1.")
    }
    if ((is.na(alpha) | is.na(beta)) & is.na(N)) {
      stop("The combination of alpha and beta, and N, cannot both be NA.")
    }
    if ((!is.na(alpha) | (!is.na(beta))) & (!is.na(N))) {
      warning("When both alpha and beta are not NA, N will be calculated automatically.")
    }
    lambda1_j <- exp(delta_j) * lambda0_j
    lambda1_a <- exp(delta_a) * lambda0_a
    sigma1_j <- sqrt(1 / lambda1_j / t_j + k)
    sigma0_j <- sqrt(1 / lambda0_j / t_j + k)
    sigma1_a <- sqrt(1 / lambda1_a / t_a + k)
    sigma0_a <- sqrt(1 / lambda0_a / t_a + k)
    if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
      N <- getN_Count_Super(delta = delta_a, lambda0 = lambda0_a, t = t_a, k = k, alpha = alpha, beta = beta, N = NA, r = r)$N
    }
    getPwr <- function(f) {
      Nj <- N * f
      var_j <- sigma1_j^2 / (r * Nj / (1 + r)) + sigma0_j^2 / (Nj / (1 + r))
      var_a <- sigma1_a^2 / (r * N / (1 + r)) + sigma0_a^2 / (N / (1 + r))
      sej <- sqrt(var_j + pi^2 * var_a - 2 * pi * sqrt(f) * sqrt(var_j * var_a))
      uj <- (delta_j - pi * delta_a) / sej
      uj <- dplyr::if_else(delta_a < 0, (-1) * uj, uj)
      mvtnorm::pmvnorm(lower = c(0), upper = c(Inf), mean = c(uj), sigma = 1)
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
    delta_nj <- lambda0_nj <- t_nj <- pwr <- NA
    if (!is.na(f)) {
      delta_nj <- (delta_a - f * delta_j) / (1 - f)
      lambda0_nj <- (lambda0_a - f * lambda0_j) / (1 - f)
      t_nj <- (t_a - f * t_j) / (1 - f)
      pwr <- getPwr(f)
    }
    data.frame(delta_a, delta_j, delta_nj, lambda0_a, lambda0_j, lambda0_nj, t_a, t_j, t_nj, k, pi, alpha, beta, N, r, pwr, beta1, f, Nj = N * f) %>%
      dplyr::do(magrittr::set_rownames(., 1:nrow(.)))
  }, .options = furrr::furrr_options(seed = TRUE))
  return(res)
}

#' @rdname getN_Count_JM1
#' @export
getN_Count_Noninf_JM1 <- function(delta_a, delta_j, lambda0_a, lambda0_j, t_a, t_j, k = 0, pi = 0.5, cut, alpha = NA, beta = NA, beta1 = 0.2, N = NA, r = 1, direct = 1) {
  eg <- as.data.frame(expand.grid(delta_a = delta_a, delta_j = delta_j, lambda0_a = lambda0_a, lambda0_j = lambda0_j, t_a = t_a, t_j = t_j, k = k, pi = pi, cut = cut, alpha = alpha, beta = beta, beta1 = beta1, N = N, r = r, stringsAsFactors = FALSE))
  res <- furrr::future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta_a <- R$delta_a
    delta_j <- R$delta_j
    lambda0_a <- R$lambda0_a
    lambda0_j <- R$lambda0_j
    t_a <- R$t_a
    t_j <- R$t_j
    k <- R$k
    pi <- R$pi
    cut <- R$cut
    alpha <- R$alpha
    beta <- R$beta
    beta1 <- R$beta1
    N <- R$N
    r <- R$r
    if (pi < 0 | pi > 1) {
      warning("Parameter pi generally is between 0 and 1.")
    }
    if (cut < 0) {
      warning("Parameter cut should be a positive value.")
    }
    if ((is.na(alpha) | is.na(beta)) & is.na(N)) {
      stop("The combination of alpha and beta, and N, cannot both be NA.")
    }
    if ((!is.na(alpha) | (!is.na(beta))) & (!is.na(N))) {
      warning("When both alpha and beta are not NA, N will be calculated automatically.")
    }
    if (!direct %in% c(-1, 1)) {
      stop("Parameter direct should be one of `1` or `-1`.")
    }
    lambda1_j <- exp(delta_j) * lambda0_j
    lambda1_a <- exp(delta_a) * lambda0_a
    sigma1_j <- sqrt(1 / lambda1_j / t_j + k)
    sigma0_j <- sqrt(1 / lambda0_j / t_j + k)
    sigma1_a <- sqrt(1 / lambda1_a / t_a + k)
    sigma0_a <- sqrt(1 / lambda0_a / t_a + k)
    if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
      N <- getN_Count_Noninf(delta = delta_a, cut = cut, lambda0 = lambda0_a, t = t_a, k = k, alpha = alpha, beta = beta, N = NA, r = r)$N
    }
    getPwr <- function(f) {
      Nj <- N * f
      var_j <- sigma1_j^2 / (r * Nj / (1 + r)) + sigma0_j^2 / (Nj / (1 + r))
      var_a <- sigma1_a^2 / (r * N / (1 + r)) + sigma0_a^2 / (N / (1 + r))
      sej <- sqrt(var_j + var_a - 2 * sqrt(f) * sqrt(var_j * var_a))
      uj <- dplyr::if_else(direct == 1, (delta_j - delta_a + pi * cut) / sej, (delta_j - delta_a - pi * cut) / sej)
      uj <- dplyr::if_else(direct == -1, (-1) * uj, uj)
      mvtnorm::pmvnorm(lower = c(0), upper = c(Inf), mean = c(uj), sigma = 1)
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
    delta_nj <- lambda0_nj <- t_nj <- pwr <- NA
    if (!is.na(f)) {
      delta_nj <- (delta_a - f * delta_j) / (1 - f)
      lambda0_nj <- (lambda0_a - f * lambda0_j) / (1 - f)
      t_nj <- (t_a - f * t_j) / (1 - f)
      pwr <- getPwr(f)
    }
    data.frame(delta_a, delta_j, delta_nj, lambda0_a, lambda0_j, lambda0_nj, t_a, t_j, t_nj, k, pi, cut, alpha, beta, N, r, direct, pwr, beta1, f, Nj = N * f) %>%
      dplyr::do(magrittr::set_rownames(., 1:nrow(.)))
  }, .options = furrr::furrr_options(seed = TRUE))
  return(res)
}

#' @rdname getN_Count_JM1
#' @export
getN_Count_Equi_JM1 <- function(delta_a, delta_j, lambda0_a, lambda0_j, t_a, t_j, k = 0, pi = 0.5, cut, alpha = NA, beta = NA, beta1 = 0.2, N = NA, r = 1, maxN = 1e+06) {
  eg <- as.data.frame(expand.grid(delta_a = delta_a, delta_j = delta_j, lambda0_a = lambda0_a, lambda0_j = lambda0_j, t_a = t_a, t_j = t_j, k = k, pi = pi, cut = cut, alpha = alpha, beta = beta, beta1 = beta1, N = N, r = r, stringsAsFactors = FALSE))
  res <- furrr::future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta_a <- R$delta_a
    delta_j <- R$delta_j
    lambda0_a <- R$lambda0_a
    lambda0_j <- R$lambda0_j
    t_a <- R$t_a
    t_j <- R$t_j
    k <- R$k
    pi <- R$pi
    cut <- R$cut
    alpha <- R$alpha
    beta <- R$beta
    beta1 <- R$beta1
    N <- R$N
    r <- R$r
    if (pi < 0 | pi > 1) {
      warning("Parameter pi generally is between 0 and 1.")
    }
    if (cut < 0) {
      warning("Parameter cut should be a positive value.")
    }
    if ((is.na(alpha) | is.na(beta)) & is.na(N)) {
      stop("The combination of alpha and beta, and N, cannot both be NA.")
    }
    if ((!is.na(alpha) | (!is.na(beta))) & (!is.na(N))) {
      warning("When both alpha and beta are not NA, N will be calculated automatically.")
    }
    lambda1_j <- exp(delta_j) * lambda0_j
    lambda1_a <- exp(delta_a) * lambda0_a
    sigma1_j <- sqrt(1 / lambda1_j / t_j + k)
    sigma0_j <- sqrt(1 / lambda0_j / t_j + k)
    sigma1_a <- sqrt(1 / lambda1_a / t_a + k)
    sigma0_a <- sqrt(1 / lambda0_a / t_a + k)
    if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
      N <- getN_Count_Equi(delta = delta_a, cut = cut, lambda0 = lambda0_a, t = t_a, k = k, alpha = alpha, beta = beta, N = NA, r = r, maxN = maxN)$N
    }
    getPwr <- function(f) {
      Nj <- N * f
      var_j <- sigma1_j^2 / (r * Nj / (1 + r)) + sigma0_j^2 / (Nj / (1 + r))
      var_a <- sigma1_a^2 / (r * N / (1 + r)) + sigma0_a^2 / (N / (1 + r))
      sej <- sqrt(var_j + var_a - 2 * sqrt(f) * sqrt(var_j * var_a))
      uj1 <- (delta_j - delta_a + pi * cut) / sej
      uj2 <- (delta_j - delta_a - pi * cut) / sej
      mvtnorm::pmvnorm(lower = c(0, -Inf), upper = c(Inf, 0), mean = c(uj1, uj2), sigma = matrix(1, nrow = 2, ncol = 2))
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
    delta_nj <- lambda0_nj <- t_nj <- pwr <- NA
    if (!is.na(f)) {
      delta_nj <- (delta_a - f * delta_j) / (1 - f)
      lambda0_nj <- (lambda0_a - f * lambda0_j) / (1 - f)
      t_nj <- (t_a - f * t_j) / (1 - f)
      pwr <- getPwr(f)
    }
    data.frame(delta_a, delta_j, delta_nj, lambda0_a, lambda0_j, lambda0_nj, t_a, t_j, t_nj, k, pi, cut, alpha, beta, N, r, pwr, beta1, f, Nj = N * f) %>%
      dplyr::do(magrittr::set_rownames(., 1:nrow(.)))
  }, .options = furrr::furrr_options(seed = TRUE))
  return(res)
}
