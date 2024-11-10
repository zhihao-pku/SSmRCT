#' Regional sample size allocation using Japan's Method 1 for continuous endpoints
#'
#' Based on Japan's Method 1, given the global sample size and marginal probability (power) of efficacy consistency between target region and globally, calculate the required sample size allocated to the target region, in clinical trials using superiority, non-inferiority, and  equivalence designs with continuous endpoints.
#'
#' @rdname getN_Con_JM1
#'
#' @name getN_Con_JM1
#'
#' @param delta_a Mean difference between treatment and control groups globally.
#' @param delta_j Mean difference between treatment and control groups in target region.
#' @param sigma Common standard deviation.
#' @param pi Proportion of global efficacy to retain. The default value is 0.5, which means retaining half of the efficacy.
#' @param cut A positive value for non-inferiority or equivalence margin.
#' @param alpha One-sided type I error rate for global success, which is used to calculate global sample size only when \code{N} is \code{NA}. The default value is 0.025.
#' @param beta Type II error rate for global success, which is used to calculate global sample size only when \code{N} is \code{NA}.
#' @param beta1 Type II error rate for efficacy consistency between target region and globally. The default value is 0.2.
#' @param N Global sample size. When \code{N} is \code{NA} and \code{alpha} and \code{beta} are not \code{NAs}, \code{N} will be calculated automatically.
#' @param r Ratio of the sample sizes of the treatment group to the control group. The default value is 1.
#' @param direct \code{direct = 1} indicates that a larger mean difference is preferable, while \code{direct = -1} indicates that a smaller mean difference is preferable.
#' @param maxN Maximum possible sample size (\code{N}) in equivalence design. Default value is 1e6.
#'
#' @return A data frame where \code{f} is required proportion of sample size allocated to the target region, and \code{Nj} is the required sample size for the target region, calculated as \code{Nj = N * f}.
#'
#' @details
#' The global success criterion and the efficacy consistency criterion between target region and globally could be found in \code{\link{getPwr_Con_Super_JM1}}.
#'
#' @references
#' 1. Quan H, Li M, Chen J, et al. Assessment of Consistency of Treatment Effects in Multiregional Clinical Trials. Drug Information J. 2010;44(5):617-632. doi:10.1177/009286151004400509
#'
#' 2. Liao JJZ, Yu Z, Li Y. Sample size allocation in multiregional equivalence studies. Pharm Stat. 2018;17(5):570-577. doi:10.1002/pst.1871
#'
#' @seealso
#' \code{\link{getPwr_Con_Super_JM1}}, \code{\link{getN_Con_Super}}.
#'
#' @export
#'
#' @examples
#' getN_Con_Super_JM1(
#'   delta_a = 0.7, delta_j = 0.5, sigma = 1, pi = 0.5, beta1 = 0.2, N = 100,
#'   r = 1
#' )
#'
#' getN_Con_Noninf_JM1(
#'   delta_a = 0, delta_j = 0.5, sigma = 4, pi = 0.5, cut = 2, alpha = 0.025,
#'   beta = 0.2, beta1 = 0.2, N = NA, r = 1, direct = -1
#' )
getN_Con_Super_JM1 <- function(delta_a, delta_j, sigma, pi = 0.5, alpha = NA, beta = NA, beta1 = 0.2, N = NA, r = 1) {
  eg <- as.data.frame(expand.grid(delta_a = delta_a, delta_j = delta_j, sigma = sigma, pi = pi, alpha = alpha, beta = beta, beta1 = beta1, N = N, r = r, stringsAsFactors = FALSE))
  res <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta_a <- R$delta_a
    delta_j <- R$delta_j
    sigma <- R$sigma
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
    gr <- 2 + r + 1 / r
    if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
      N <- getN_Con_Super(delta = delta_a, sigma = sigma, alpha = alpha, beta = beta, N = NA, r = r)$N
    }
    getPwr <- function(f) {
      Nj <- N * f
      sej <- sqrt(gr * sigma^2 / Nj + pi^2 * gr * sigma^2 / N - 2 * pi * sqrt(f) * sqrt(gr * sigma^2 / Nj * gr * sigma^2 / N))
      uj <- (delta_j - pi * delta_a) / sej
      uj <- if_else(delta_a < 0, (-1) * uj, uj)
      pmvnorm(lower = c(0), upper = c(Inf), mean = c(uj), sigma = 1)
    }
    f <- tryCatch(
      {
        uniroot(f = function(f) {
          getPwr(f) - (1 - beta1)
        }, interval = c(1e-06, 1 - 1e-06))$root
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
    data.frame(delta_a, delta_j, delta_nj, sigma, pi, alpha, beta, beta1, N, r, pwr, f, Nj = N * f)
  }, .options = furrr_options(seed = TRUE))
  return(res)
}

#' @rdname getN_Con_JM1
#' @export
getN_Con_Noninf_JM1 <- function(delta_a, delta_j, sigma, pi = 0.5, cut, alpha = NA, beta = NA, beta1 = 0.2, N = NA, r = 1, direct = 1) {
  eg <- as.data.frame(expand.grid(delta_a = delta_a, delta_j = delta_j, sigma = sigma, pi = pi, cut = cut, alpha = alpha, beta = beta, beta1 = beta1, N = N, r = r, stringsAsFactors = FALSE))
  res <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta_a <- R$delta_a
    delta_j <- R$delta_j
    sigma <- R$sigma
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
    gr <- 2 + r + 1 / r
    if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
      N <- getN_Con_Noninf(delta = delta_a, sigma = sigma, cut = cut, alpha = alpha, beta = beta, N = NA, r = r)$N
    }
    getPwr <- function(f) {
      Nj <- N * f
      sej <- sqrt(gr * sigma^2 / Nj + gr * sigma^2 / N - 2 * sqrt(f) * sqrt(gr * sigma^2 / Nj * gr * sigma^2 / N))
      uj <- if_else(direct == 1, (delta_j - delta_a + pi * cut) / sej, (delta_j - delta_a - pi * cut) / sej)
      uj <- if_else(direct == -1, (-1) * uj, uj)
      pmvnorm(lower = c(0), upper = c(Inf), mean = c(uj), sigma = 1)
    }
    f <- tryCatch(
      {
        uniroot(f = function(f) {
          getPwr(f) - (1 - beta1)
        }, interval = c(1e-06, 1 - 1e-06))$root
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
    data.frame(delta_a, delta_j, delta_nj, sigma, pi, cut, alpha, beta, beta1, N, r, direct, pwr, f, Nj = N * f)
  }, .options = furrr_options(seed = TRUE))
  return(res)
}

#' @rdname getN_Con_JM1
#' @export
getN_Con_Equi_JM1 <- function(delta_a, delta_j, sigma, pi = 0.5, cut, alpha = NA, beta = NA, beta1 = 0.2, N = NA, r = 1, maxN = 1e6) {
  eg <- as.data.frame(expand.grid(delta_a = delta_a, delta_j = delta_j, sigma = sigma, pi = pi, cut = cut, alpha = alpha, beta = beta, beta1 = beta1, N = N, r = r, stringsAsFactors = FALSE))
  res <- future_map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta_a <- R$delta_a
    delta_j <- R$delta_j
    sigma <- R$sigma
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
    gr <- 2 + r + 1 / r
    if (is.na(N) & (!is.na(alpha)) & (!is.na(beta))) {
      N <- getN_Con_Equi(delta = delta_a, sigma = sigma, cut = cut, alpha = alpha, beta = beta, N = NA, r = r, maxN = maxN)$N
    }
    getPwr <- function(f) {
      Nj <- N * f
      sej <- sqrt(gr * sigma^2 / Nj + gr * sigma^2 / N - 2 * sqrt(f) * sqrt(gr * sigma^2 / Nj * gr * sigma^2 / N))
      uj1 <- (delta_j - delta_a + pi * cut) / sej
      uj2 <- (delta_j - delta_a - pi * cut) / sej
      pmvnorm(lower = c(0, -Inf), upper = c(Inf, 0), mean = c(uj1, uj2), sigma = matrix(1, nrow = 2, ncol = 2))
    }
    f <- tryCatch(
      {
        uniroot(f = function(f) {
          getPwr(f) - (1 - beta1)
        }, interval = c(1e-06, 1 - 1e-06))$root
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
    data.frame(delta_a, delta_j, delta_nj, sigma, pi, cut, alpha, beta, beta1, N, r, pwr, f, Nj = N * f)
  }, .options = furrr_options(seed = TRUE))
  return(res)
}
