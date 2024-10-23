#' Title
#'
#' @param delta_j
#' @param delta_nj
#' @param pi
#' @param beta1
#' @param N
#' @param r
#' @param criterion
#' @param direct
#'
#' @return
#' @export
#'
#' @examples
#' getN_Surv_Super_JM1(
#'   delta_j = log(0.8), delta_nj = log(0.7),
#'   pi = 0.5, beta1 = 0.2, N = seq(400, 800, 200),
#'   criterion = c(1, 2), r = 1, direct = -1
#' )
getN_Surv_Super_JM1 <- function(delta_j, delta_nj, pi, beta1, N, r, criterion, direct = -1) {
  eg <- as.data.frame(expand.grid(
    delta_j = delta_j,
    delta_nj = delta_nj,
    pi = pi,
    beta1 = beta1,
    N = N,
    r = r,
    criterion = criterion,
    direct = direct
  ))
  res <- map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta_j <- R$delta_j
    delta_nj <- R$delta_nj
    pi <- R$pi
    beta1 <- R$beta1
    N <- R$N
    r <- R$r
    criterion <- R$criterion
    direct <- R$direct
    gr <- 2 + r + 1 / r
    getPwr <- function(f) {
      Nj <- N * f
      delta <- delta_j * f + delta_nj * (1 - f)
      if (criterion == 1) {
        sej <- sqrt(gr / Nj +
          pi^2 * gr / N -
          2 * pi * sqrt(f) *
            sqrt(gr / Nj * gr / N))
        uj <- (delta_j - pi * delta) / sej
      }
      if (criterion == 2) {
        sej <- sqrt(gr / Nj * exp(delta_j)^2 +
          pi^2 * gr / N * exp(delta)^2 -
          2 * pi * gr / N * exp(delta_j) * exp(delta))
        uj <- -(1 - exp(delta_j) - pi * (1 - exp(delta))) / sej
      }
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
      delta_j, delta_nj,
      pi, beta1, N, r, criterion,
      direct,
      pwr = 1 - beta1, f, Nj = N * f
    )
  })
  return(res)
}
