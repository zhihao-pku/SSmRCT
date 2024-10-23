#' Title
#'
#' @param delta_j
#' @param delta_nj
#' @param sigma
#' @param pi
#' @param alpha
#' @param beta1
#' @param N
#' @param r
#' @param direct
#'
#' @return
#' @export
#'
#' @examples
#' getN_Con_Super_JM1(
#'   delta_j = 0.5, delta_nj = 0.7, sigma = 1,
#'   pi = 0.5, alpha = 0.025, beta1 = 0.2, N = seq(100, 400, 100), r = 1, direct = 1
#' )
getN_Con_Super_JM1 <- function(delta_j, delta_nj, sigma, pi, alpha, beta1, N, r, direct = 1) {
  eg <- as.data.frame(expand.grid(
    delta_j = delta_j,
    delta_nj = delta_nj,
    sigma = sigma,
    pi = pi,
    alpha = alpha,
    beta1 = beta1,
    N = N,
    r = r,
    direct = direct
  ))
  res <- map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta_j <- R$delta_j
    delta_nj <- R$delta_nj
    sigma <- R$sigma
    pi <- R$pi
    alpha <- R$alpha
    beta1 <- R$beta1
    N <- R$N
    r <- R$r
    direct <- R$direct
    gr <- 2 + r + 1 / r
    getPwr <- function(f) {
      Nj <- N * f
      delta <- delta_j * f + delta_nj * (1 - f)
      sej <- sqrt(gr * sigma^2 / Nj +
        pi^2 * gr * sigma^2 / N -
        2 * pi * sqrt(f) *
          sqrt(gr * sigma^2 / Nj * gr * sigma^2 / N))
      uj <- (delta_j - pi * delta) / sej
      if (direct == -1) {
        uj <- (-1) * uj
        u <- (-1) * u
      }
      pmvnorm(
        lower = c(0),
        upper = c(Inf),
        mean = c(uj),
        sigma = 1
      ) - (1 - beta1)
    }
    f <- uniroot(f = getPwr, interval = c(1e-6, 1 - 1e-6))$root
    data.frame(
      delta_j, delta_nj,
      sigma, pi,
      alpha, beta1, N, r,
      direct,
      pwr = 1 - beta1, f, Nf = N * f
    )
  })
  return(res)
}
