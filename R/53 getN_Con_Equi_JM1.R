#' Title
#'
#' @param delta_j
#' @param delta_nj
#' @param sigma
#' @param pi
#' @param cut
#' @param beta1
#' @param N
#' @param r
#'
#' @return
#' @export
#'
#' @examples
#' getN_Con_Equi_JM1(
#'   delta_j = -0.1, delta_nj = 0, sigma = 1,
#'   pi = 0.5, cut = 0.3, beta1 = 0.2,
#'   N = seq(100, 400, 100), r = 1
#' )
getN_Con_Equi_JM1 <- function(delta_j, delta_nj, sigma, pi, cut, beta1, N, r) {
  eg <- as.data.frame(expand.grid(
    delta_j = delta_j,
    delta_nj = delta_nj,
    sigma = sigma,
    pi = pi,
    cut = cut,
    beta1 = beta1,
    N = N,
    r = r
  ))
  res <- map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta_j <- R$delta_j
    delta_nj <- R$delta_nj
    sigma <- R$sigma
    pi <- R$pi
    cut <- R$cut
    beta1 <- R$beta1
    N <- R$N
    r <- R$r
    gr <- 2 + r + 1 / r
    getPwr <- function(f) {
      Nj <- N * f
      delta <- delta_j * f + delta_nj * (1 - f)
      sej <- sqrt(gr * sigma^2 / Nj +
        gr * sigma^2 / N -
        2 * sqrt(f) *
          sqrt(gr * sigma^2 / Nj * gr * sigma^2 / N))
      uj1 <- (delta_j - delta + pi * cut) / sej
      uj2 <- (delta_j - delta - pi * cut) / sej
      p2 <- pmvnorm(
        lower = c(0, -Inf),
        upper = c(Inf, 0),
        mean = c(uj1, uj2),
        sigma = matrix(1, nrow = 2, ncol = 2)
      ) - (1 - beta1)
    }
    f <- tryCatch(
      {
        uniroot(f = getPwr, interval = c(1e-6, 1 - 1e-6))$root
      },
      error = function(e) {
        NA
      },
      finally = print("Some N are too small, resulting in NA")
    )
    data.frame(
      delta_j, delta_nj,
      sigma, pi, cut,
      beta1, N, r,
      pwr = 1 - beta1, f, Nj = N * f
    )
  })
  return(res)
}
