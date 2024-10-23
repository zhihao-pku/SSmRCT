#' Title
#'
#' @param p1_j
#' @param p0_j
#' @param p1_nj
#' @param p0_nj
#' @param pi
#' @param beta1
#' @param N
#' @param r
#' @param direct
#'
#' @return
#' @export
#'
#' @examples
#' getN_Bin_Super_JM1(
#'   p1_j = 0.6, p0_j = 0.4, p1_nj = 0.7, p0_nj = 0.4,
#'   pi = 0.5, beta1 = 0.2, N = seq(100, 400, 100), r = 1, direct = 1
#' )
getN_Bin_Super_JM1 <- function(p1_j, p0_j, p1_nj, p0_nj, pi, beta1, N, r, direct = 1) {
  eg <- as.data.frame(expand.grid(
    p1_j = p1_j,
    p0_j = p0_j,
    p1_nj = p1_nj,
    p0_nj = p0_nj,
    pi = pi,
    beta1 = beta1,
    N = N,
    r = r,
    direct = direct
  ))
  res <- map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    p1_j <- R$p1_j
    p0_j <- R$p0_j
    p1_nj <- R$p1_nj
    p0_nj <- R$p0_nj
    pi <- R$pi
    beta1 <- R$beta1
    N <- R$N
    r <- R$r
    direct <- R$direct
    delta_j <- p1_j - p0_j
    delta_nj <- p1_nj - p0_nj
    p_j <- p0_j / (1 + r) + p1_j * r / (1 + r)
    sigma_j <- sqrt(p_j * (1 - p_j))
    gr <- 2 + r + 1 / r
    getPwr <- function(f) {
      Nj <- N * f
      delta <- delta_j * f + delta_nj * (1 - f)
      p0_a <- p0_j * f + p0_nj * (1 - f)
      p1_a <- p1_j * f + p1_nj * (1 - f)
      p_a <- p0_a / (1 + r) + p1_a * r / (1 + r)
      sigma_a <- sqrt(p_a * (1 - p_a))
      sej <- sqrt(gr * sigma_j^2 / Nj +
        pi^2 * gr * sigma_a^2 / N -
        2 * pi * sqrt(f) *
          sqrt(gr * sigma_j^2 / Nj * gr * sigma_a^2 / N))
      uj <- (delta_j - pi * delta) / sej
      uj <- if_else(direct == -1, (-1) * uj, uj)
      pmvnorm(
        lower = c(0),
        upper = c(Inf),
        mean = c(uj),
        sigma = 1
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
      p0_j, p1_j,
      p0_nj, p1_nj,
      pi,
      beta1, N, r,
      direct,
      pwr = 1 - beta1, f, Nj = N * f
    )
  })
  return(res)
}
