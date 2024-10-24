#' Title
#'
#' @param delta
#' @param sigma
#' @param cut
#' @param alpha
#' @param beta
#' @param N
#' @param r
#' @param maxN
#'
#' @return
#' @export
#'
#' @examples
#' (v <- getN_Con_Equi(
#'   delta = 0, sigma = 2, cut = 1, alpha = 0.025,
#'   beta = 0.2, N = NA, r = 1
#' ))
#' getN_Con_Equi(
#'   delta = 0, sigma = 2, cut = 1, alpha = 0.025,
#'   beta = NA, N = v$N, r = 1
#' )
getN_Con_Equi <- function(delta, sigma, cut, alpha, beta, N, r, maxN = 1e6) {
  eg <- as.data.frame(expand.grid(
    delta = delta,
    sigma = sigma,
    cut = cut,
    alpha = alpha,
    beta = beta,
    N = N,
    r = r
  ))
  res <- map_dfr(.x = 1:nrow(eg), .f = function(i) {
    R <- eg[i, ]
    delta <- R$delta
    sigma <- R$sigma
    cut <- R$cut
    alpha <- R$alpha
    beta <- R$beta
    N <- R$N
    r <- R$r
    if (is.na(N)) {
      getN <- function(N) {
        n2 <- N / (r + 1)
        n1 <- r * n2
        z1 <- (delta + cut) / (sigma * sqrt(1 / n1 + 1 / n2))
        z2 <- (delta - cut) / (sigma * sqrt(1 / n1 + 1 / n2))
        pwr <- pmvnorm(
          lower = c(qnorm(1 - alpha), -Inf),
          upper = c(Inf, -qnorm(1 - alpha)),
          mean = c(z1, z2),
          sigma = matrix(1, nrow = 2, ncol = 2)
        ) - (1 - beta)
      }
      N <- uniroot(f = getN, interval = c(1e-6, maxN))$root
      n2 <- ceiling(N / (r + 1))
      n1 <- ceiling(r * n2)
      N <- n1 + n2
      pwr <- 1 - beta
      df <- data.frame(delta, sigma, cut, alpha, beta, pwr, r, N, n1, n2)
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
      df <- data.frame(delta, sigma, cut, alpha, beta, pwr, r, N, n1, n2)
    }
    df
  })
  return(res)
}
