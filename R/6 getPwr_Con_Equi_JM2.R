#' Title
#'
#' @param delta_i
#' @param sigma
#' @param fi
#' @param cut
#' @param alpha
#' @param beta
#' @param N
#' @param r
#' @param sim
#' @param nsim
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
#' f_set <- seq(0.1, 0.9, 0.1)
#' map_dfr(.x = 1:length(f_set), .f = function(i) {
#'   f <- f_set[i]
#'   res <- getPwr_Con_Equi_JM2(
#'     delta_i = c(-0.5, 0), sigma = 4,
#'     fi = c(f, 1 - f), cut = 2,
#'     alpha = 0.025, beta = NA, N = 200, r = 1, sim = FALSE
#'   )$overall
#'   res$M <- "calc"
#'   res$f <- f
#'   res
#' })
#'
#' f_set <- seq(0.1, 0.9, 0.1)
#' map_dfr(.x = 1:length(f_set), .f = function(i) {
#'   f <- f_set[i]
#'   res <- getPwr_Con_Equi_JM2(
#'     delta_i = c(-0.5, 0), sigma = 4,
#'     fi = c(f, 1 - f), cut = 2,
#'     alpha = 0.025, beta = 0.2, N = NA, r = 1, sim = FALSE
#'   )$overall
#'   res$M <- "calc"
#'   res$f <- f
#'   res
#' })
getPwr_Con_Equi_JM2 <- function(delta_i, sigma, fi, cut, alpha = 0.025, beta = NA, N, r = 1, sim = FALSE, nsim = 1000, seed = 0) {
  if (!sim) {
    gr <- 2 + r + 1 / r
    delta <- sum(delta_i * fi)
    if (is.na(N)) {
      N <- getN_Con_Equi(
        delta, sigma, cut, alpha,
        beta, NA, r
      )$N
    }
    Ni <- N * fi
    num <- length(delta_i)
    sei <- sqrt(gr * sigma^2 / Ni)
    ui1 <- (delta_i + cut) / sei
    ui2 <- (delta_i - cut) / sei
    se <- sqrt(gr * sigma^2 / N)
    u1 <- (delta + cut) / se
    u2 <- (delta - cut) / se
    M <- diag(num + 1)
    M[1, ] <- c(1, sqrt(fi))
    M[, 1] <- c(1, sqrt(fi))
    M1 <- combine(M, M)
    p1 <- pmvnorm(
      lower = c(qnorm(1 - alpha), -Inf),
      upper = c(Inf, -qnorm(1 - alpha)),
      mean = c(u1, u2),
      corr = matrix(1, nrow = 2, ncol = 2)
    )
    p_margin <- c()
    for (k in 1:num) {
      p_margin_ <- pmvnorm(
        lower = c(0, -Inf),
        upper = c(Inf, 0),
        mean = c(ui1[k], ui2[k]),
        corr = matrix(1, nrow = 2, ncol = 2)
      )
      p_margin <- c(p_margin, p_margin_)
    }
    p2 <- prod(p_margin)
    p3 <- pmvnorm(
      lower = c(qnorm(1 - alpha), -Inf, rep(c(0, -Inf), num)),
      upper = c(Inf, -qnorm(1 - alpha), rep(c(Inf, 0), num)),
      mean = c(u1, u2, combine(ui1, ui2)),
      sigma = M1
    )
    p4 <- p3 / p1
    res <- data.frame(
      delta = delta, N = N,
      p1 = p1, p2 = p2, p3 = p3, p4 = p4
    )
    p_joint <- c()
    for (k in 1:num) {
      M <- diag(2)
      M[1, ] <- c(1, sqrt(fi[k]))
      M[, 1] <- c(1, sqrt(fi[k]))
      M2 <- combine(M, M)
      p_joint_ <- pmvnorm(
        lower = c(qnorm(1 - alpha), -Inf, 0, -Inf),
        upper = c(Inf, -qnorm(1 - alpha), Inf, 0),
        mean = c(u1, u2, ui1[k], ui2[k]),
        sigma = M2
      )
      p_joint <- c(p_joint, p_joint_)
    }
    p_conditional <- p_joint / p1
    L <- list(
      overall = res,
      p_margin = p_margin,
      p_joint = p_joint,
      p_conditional = p_conditional
    )
  }
  if (sim) {
    gr <- 2 + r + 1 / r
    delta <- sum(delta_i * fi)
    if (is.na(N)) {
      N <- getN_Con_Equi(
        delta, sigma, cut, alpha,
        beta, NA, r
      )$N
    }
    Ni <- N * fi
    num <- length(delta_i)
    da <- data.frame()
    di <- NULL
    for (j in 1:nsim) {
      set.seed(j + seed)
      Xt <- matrix(NA, nrow = N * r / (1 + r), ncol = num)
      Xc <- matrix(NA, nrow = N / (1 + r), ncol = num)
      for (k in 1:num) {
        nt_k <- round(Ni[k] * r / (1 + r))
        nc_k <- round(Ni[k] / (1 + r))
        Xt[1:nt_k, k] <- rnorm(n = nt_k, mean = delta_i[k], sd = sigma)
        Xc[1:nc_k, k] <- rnorm(n = nc_k, mean = 0, sd = sigma)
      }
      ut <- colMeans(Xt, na.rm = T)
      uc <- colMeans(Xc, na.rm = T)
      zi1 <- ut - uc + cut
      zi2 <- ut - uc - cut
      za1 <- (mean(as.numeric(Xt), na.rm = T) -
        mean(as.numeric(Xc), na.rm = T) + cut) / sqrt(gr * sigma^2 / N)
      za2 <- (mean(as.numeric(Xt), na.rm = T) -
        mean(as.numeric(Xc), na.rm = T) - cut) / sqrt(gr * sigma^2 / N)
      succ_a <- if_else(za1 > qnorm(1 - alpha) & za2 < -qnorm(1 - alpha), 1, 0)
      succ_i <- if_else(all(zi1 > 0) & all(zi2 < 0), 1, 0)
      succ_i_ <- if_else(zi1 > 0 & zi2 < 0, 1, 0)
      da <- bind_rows(da, c(succ_a = succ_a, succ_i = succ_i))
      di <- rbind(di, succ_i_)
    }
    res <- da %>%
      summarise(
        p1 = mean(succ_a, na.rm = T),
        p2 = mean(succ_i, na.rm = T),
        p3 = mean(succ_a & succ_i, na.rm = T),
        p4 = mean(succ_i[succ_a == 1], na.rm = T)
      ) %>%
      mutate(
        delta = delta,
        N = N
      ) %>%
      select(delta, N, p1, p2, p3, p4)
    p_margin <- colMeans(di)
    p_joint <- c()
    for (k in 1:num) {
      p_joint_ <- mean(da$succ_a & di[, k])
      p_joint <- c(p_joint, p_joint_)
    }
    p_conditional <- colMeans(di[da$succ_a == 1, ])
    L <- list(
      overall = res,
      p_margin = p_margin,
      p_joint = p_joint,
      p_conditional = p_conditional
    )
  }
  return(L)
}
