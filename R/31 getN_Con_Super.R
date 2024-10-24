#' Title
#'
#' @param delta
#' @param sigma
#' @param alpha
#' @param beta
#' @param N
#' @param r
#'
#' @return
#' @export
#'
#' @examples
#' (v <- getN_Con_Super(
#'   delta = 1, sigma = 4, alpha = 0.025,
#'   beta = 0.2, N = NA, r = 1
#' ))
#' getN_Con_Super(
#'   delta = 1, sigma = 4, alpha = 0.025,
#'   beta = NA, N = v$N, r = 1
#' )
getN_Con_Super <- function(delta, sigma, alpha, beta, N, r) {
  eg <- as.data.frame(expand.grid(
    delta = delta,
    sigma = sigma,
    alpha = alpha,
    beta = beta,
    N = N,
    r = r
  ))
  res <- map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta <- R$delta
    sigma <- R$sigma
    alpha <- R$alpha
    beta <- R$beta
    N <- R$N
    r <- R$r
    if (is.na(N)) {
      n2 <- (qnorm(1 - alpha) + qnorm(1 - beta))^2 *
        sigma^2 * (1 + 1 / r) / delta^2
      n1 <- r * n2
      n1 <- ceiling(n1)
      n2 <- ceiling(n2)
      N <- n1 + n2
      pwr <- 1 - beta
      df <- data.frame(delta, sigma, alpha, beta, pwr, r, N, n1, n2)
    }
    if (is.na(beta)) {
      n2 <- N / (r + 1)
      n1 <- r * n2
      z <- delta / (sigma * sqrt(1 / n1 + 1 / n2))
      pwr <- if_else(delta > 0,
        1 - pnorm(q = qnorm(1 - alpha), mean = z, sd = 1),
        pnorm(q = -qnorm(1 - alpha), mean = z, sd = 1)
      )
      df <- data.frame(delta, sigma, alpha, beta, pwr, r, N, n1, n2)
    }
    df
  })
  return(res)
}
