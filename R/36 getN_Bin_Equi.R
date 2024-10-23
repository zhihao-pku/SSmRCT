#' Title
#'
#' @param p1
#' @param p0
#' @param cut
#' @param alpha
#' @param beta
#' @param N
#' @param r
#' @param direct
#'
#' @return
#' @export
#'
#' @examples
#' (v <- getN_Bin_Noninf(
#'   p1 = 0.4, p0 = 0.45, cut = 0.2, alpha = 0.025,
#'   beta = 0.2, N = NA, r = 1
#' ))
#' getN_Bin_Noninf(
#'   p1 = 0.4, p0 = 0.45, cut = 0.2, alpha = 0.025,
#'   beta = NA, N = v$N, r = 1
#' )
getN_Bin_Equi <- function(p1, p0, cut, alpha, beta, N, r, direct = 1) {
  eg <- as.data.frame(expand.grid(
    p1 = p1,
    p0 = p0,
    cut = cut,
    alpha = alpha,
    beta = beta,
    N = N,
    r = r,
    direct = direct
  ))
  res <- map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    p1 <- R$p1
    p0 <- R$p0
    cut <- R$cut
    alpha <- R$alpha
    beta <- R$beta
    N <- R$N
    r <- R$r
    direct <- R$direct
    delta <- p1 - p0
    p_pool <- p1 * r / (1 + r) + p0 / (1 + r)
    sigma <- sqrt(p_pool * (1 - p_pool))
    if (is.na(N)) {
      n2 <- if_else(delta == 0,
        (qnorm(1 - alpha) + qnorm(1 - beta / 2))^2 *
          sigma^2 * (1 + 1 / r) / (abs(delta) - cut)^2,
        (qnorm(1 - alpha) + qnorm(1 - beta))^2 *
          sigma^2 * (1 + 1 / r) / (abs(delta) - cut)^2
      )
      n1 <- r * n2
      n1 <- ceiling(n1)
      n2 <- ceiling(n2)
      N <- n1 + n2
      pwr <- 1 - beta
      df <- data.frame(p1 = p1, p0 = p0, delta, cut, alpha, beta, pwr, r, N, n1, n2)
    }
    if (is.na(beta)) {
      n2 <- N / (r + 1)
      n1 <- r * n2
      z1 <- (delta + cut) / (sigma * sqrt(1 / n1 + 1 / n2))
      z2 <- (delta - cut) / (sigma * sqrt(1 / n1 + 1 / n2))
      pwr <- pmvnorm(
        lower = c(qnorm(1 - alpha), -Inf),
        upper = c(Inf, -qnorm(1 - alpha)),
        mean = c(z1, z2),
        sigma = matrix(1, nrow = 2, ncol = 2)
      )
      df <- data.frame(p1 = p1, p0 = p0, delta, cut, alpha, beta, pwr, r, N, n1, n2)
    }
    df
  })
  return(res)
}
