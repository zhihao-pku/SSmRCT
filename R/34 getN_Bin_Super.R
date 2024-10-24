#' Title
#'
#' @param p1
#' @param p0
#' @param alpha
#' @param beta
#' @param N
#' @param r
#'
#' @return
#' @export
#'
#' @examples
#' (v <- getN_Bin_Super(
#'   p1 = 0.4, p0 = 0.2, alpha = 0.025,
#'   beta = 0.2, N = NA, r = 1
#' ))
#' getN_Bin_Super(
#'   p1 = 0.4, p0 = 0.2, alpha = 0.025,
#'   beta = NA, N = v$N, r = 1
#' )
getN_Bin_Super <- function(p1, p0, alpha, beta, N, r) {
  eg <- as.data.frame(expand.grid(
    p1 = p1,
    p0 = p0,
    alpha = alpha,
    beta = beta,
    N = N,
    r = r
  ))
  res <- map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    p1 <- R$p1
    p0 <- R$p0
    alpha <- R$alpha
    beta <- R$beta
    N <- R$N
    r <- R$r
    delta <- p1 - p0
    p_pool <- p1 * r / (1 + r) + p0 / (1 + r)
    sigma <- sqrt(p_pool * (1 - p_pool))
    if (is.na(N)) {
      n2 <- (qnorm(1 - alpha) + qnorm(1 - beta))^2 *
        sigma^2 * (1 + 1 / r) / delta^2
      n1 <- r * n2
      n1 <- ceiling(n1)
      n2 <- ceiling(n2)
      N <- n1 + n2
      pwr <- 1 - beta
      df <- data.frame(p1, p0, delta, alpha, beta, pwr, r, N, n1, n2)
    }
    if (is.na(beta)) {
      n2 <- N / (r + 1)
      n1 <- r * n2
      z <- delta / (sigma * sqrt(1 / n1 + 1 / n2))
      pwr <- if_else(delta > 0,
        1 - pnorm(q = qnorm(1 - alpha), mean = z, sd = 1),
        pnorm(q = -qnorm(1 - alpha), mean = z, sd = 1)
      )
      df <- data.frame(p1, p0, delta, alpha, beta, pwr, r, N, n1, n2)
    }
    df
  })
  return(res)
}
