#' Title
#'
#' @param delta
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
#' (v <- getN_Surv_Super(
#'   delta = log(0.8), alpha = 0.025,
#'   beta = 0.2, N = NA, r = 1
#' ))
#' getN_Surv_Super(
#'   delta = log(0.8), alpha = 0.025,
#'   beta = NA, N = v$N, r = 1
#' )
getN_Surv_Super <- function(delta, alpha, beta, N, r, direct = -1) {
  eg <- as.data.frame(expand.grid(
    delta = delta,
    alpha = alpha,
    beta = beta,
    N = N,
    r = r,
    direct = direct
  ))
  res <- map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta <- R$delta
    alpha <- R$alpha
    beta <- R$beta
    N <- R$N
    r <- R$r
    direct <- R$direct
    if (is.na(N)) {
      n2 <- (qnorm(1 - alpha) + qnorm(1 - beta))^2 *
        (1 + 1 / r) / delta^2
      n1 <- r * n2
      n1 <- ceiling(n1)
      n2 <- ceiling(n2)
      N <- n1 + n2
      pwr <- 1 - beta
      df <- data.frame(delta, alpha, beta, pwr, r, N, n1, n2)
    }
    if (is.na(beta)) {
      n2 <- N / (r + 1)
      n1 <- r * n2
      z <- delta / sqrt(1 / n1 + 1 / n2)
      pwr <- if_else(direct == 1,
        1 - pnorm(q = qnorm(1 - alpha), mean = z, sd = 1),
        pnorm(q = -qnorm(1 - alpha), mean = z, sd = 1)
      )
      df <- data.frame(delta, alpha, beta, pwr, r, N, n1, n2)
    }
    df
  })
  return(res)
}
