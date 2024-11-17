library(SSmRCT)
rm(list = ls())

pplot <- function(a, b) {
  a <- tidyr::pivot_longer(data = a, cols = dplyr::starts_with("pwr"))
  b <- tidyr::pivot_longer(data = b, cols = dplyr::starts_with("pwr"))
  ggplot2::ggplot() +
    ggplot2::geom_point(dat = a, ggplot2::aes(x = f, y = value, color = "calc")) +
    ggplot2::geom_point(dat = b, ggplot2::aes(x = f, y = value, color = "sim")) +
    ggplot2::geom_line(dat = a, ggplot2::aes(x = f, y = value, color = "calc")) +
    ggplot2::geom_line(dat = b, ggplot2::aes(x = f, y = value, color = "sim")) +
    ggplot2::facet_wrap(ggplot2::vars(name)) +
    ggplot2::scale_x_continuous(breaks = seq(0.1, 0.9, 0.2), limits = c(0.1, 0.9)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
    ggplot2::labs(x = "allocation ratio", y = "power", color = "")
}

# -------------------------------------------------------------------------
# getPwr_Con_Super_JM1、getPwr_Con_Noninf_JM1、getPwr_Con_Equi_JM1
# -------------------------------------------------------------------------

a <- getPwr_Con_Super_JM1(delta_j = 0.5, delta_a = 0.7, sigma = 1, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, N = 100, r = 1, sim = FALSE)
b <- getPwr_Con_Super_JM1(delta_j = 0.5, delta_a = 0.7, sigma = 1, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, N = 100, r = 1, sim = TRUE)
(p1_1 <- pplot(a, b))

c <- getPwr_Con_Super_JM1(delta_j = -0.5, delta_nj = -0.7, sigma = 1, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = 0.2, N = NA, r = 1, sim = FALSE)
d <- getPwr_Con_Super_JM1(delta_j = -0.5, delta_nj = -0.7, sigma = 1, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = 0.2, N = NA, r = 1, sim = TRUE)
(p1_2 <- pplot(c, d))

# -------------------------------------------------------------------------

a <- getPwr_Con_Noninf_JM1(delta_j = -0.2, delta_a = -0.1, sigma = 1, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 0.4, alpha = 0.025, beta = NA, N = 400, r = 1, direct = 1, sim = FALSE)
b <- getPwr_Con_Noninf_JM1(delta_j = -0.2, delta_a = -0.1, sigma = 1, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 0.4, alpha = 0.025, beta = NA, N = 400, r = 1, direct = 1, sim = TRUE)
(p1_3 <- pplot(a, b))

c <- getPwr_Con_Noninf_JM1(delta_j = 0.2, delta_nj = 0.1, sigma = 1, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 0.4, alpha = 0.025, beta = 0.2, N = NA, r = 1, direct = -1, sim = FALSE)
d <- getPwr_Con_Noninf_JM1(delta_j = 0.2, delta_nj = 0.1, sigma = 1, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 0.4, alpha = 0.025, beta = 0.2, N = NA, r = 1, direct = -1, sim = TRUE)
(p1_4 <- pplot(c, d))

# -------------------------------------------------------------------------

a <- getPwr_Con_Equi_JM1(delta_j = -0.2, delta_a = -0.1, sigma = 1, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 0.4, alpha = 0.025, beta = NA, N = 400, r = 1, sim = FALSE)
b <- getPwr_Con_Equi_JM1(delta_j = -0.2, delta_a = -0.1, sigma = 1, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 0.4, alpha = 0.025, beta = NA, N = 400, r = 1, sim = TRUE)
(p1_5 <- pplot(a, b))

c <- getPwr_Con_Equi_JM1(delta_j = 0.2, delta_nj = 0.1, sigma = 1, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 0.4, alpha = 0.025, beta = 0.2, N = NA, r = 1, sim = FALSE)
d <- getPwr_Con_Equi_JM1(delta_j = 0.2, delta_nj = 0.1, sigma = 1, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 0.4, alpha = 0.025, beta = 0.2, N = NA, r = 1, sim = TRUE)
(p1_6 <- pplot(c, d))

# -------------------------------------------------------------------------
# getPwr_Con_Super_JM2、getPwr_Con_Noninf_JM2、getPwr_Con_Equi_JM2
# -------------------------------------------------------------------------

f_set <- seq(0.1, 0.9, 0.1)
a <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Con_Super_JM2(delta_i = c(1, 0.8), sigma = 4, f_i = c(f, 1 - f), alpha = 0.025, beta = NA, N = 200, r = 1, sim = FALSE)$overall
  res$f <- f
  res
})
b <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Con_Super_JM2(delta_i = c(1, 0.8), sigma = 4, f_i = c(f, 1 - f), alpha = 0.025, beta = NA, N = 200, r = 1, sim = TRUE)$overall
  res$f <- f
  res
})
(p1_7 <- pplot(a, b))

f_set <- seq(0.1, 0.9, 0.1)
c <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Con_Super_JM2(delta_i = c(-1, -0.8), sigma = 4, f_i = c(f, 1 - f), alpha = 0.025, beta = 0.2, N = NA, r = 1, sim = FALSE)$overall
  res$f <- f
  res
})
d <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Con_Super_JM2(delta_i = c(-1, -0.8), sigma = 4, f_i = c(f, 1 - f), alpha = 0.025, beta = 0.2, N = NA, r = 1, sim = TRUE)$overall
  res$f <- f
  res
})
(p1_8 <- pplot(c, d))

# -------------------------------------------------------------------------

f_set <- seq(0.1, 0.9, 0.1)
a <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Con_Noninf_JM2(delta_i = c(-0.5, 0), sigma = 4, f_i = c(f, 1 - f), cut = 2, cut_i = c(1.8, 2.0), alpha = 0.025, beta = NA, N = 200, r = 1, direct = 1, sim = FALSE)$overall
  res$f <- f
  res
})
b <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Con_Noninf_JM2(delta_i = c(-0.5, 0), sigma = 4, f_i = c(f, 1 - f), cut = 2, cut_i = c(1.8, 2.0), alpha = 0.025, beta = NA, N = 200, r = 1, direct = 1, sim = TRUE)$overall
  res$f <- f
  res
})
(p1_9 <- pplot(a, b))

f_set <- seq(0.1, 0.9, 0.1)
c <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Con_Noninf_JM2(delta_i = c(1, 0), sigma = 4, f_i = c(f, 1 - f), cut = 2, alpha = 0.025, beta = 0.2, N = NA, r = 1, direct = -1, sim = FALSE)$overall
  res$f <- f
  res
})
d <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Con_Noninf_JM2(delta_i = c(1, 0), sigma = 4, f_i = c(f, 1 - f), cut = 2, alpha = 0.025, beta = 0.2, N = NA, r = 1, direct = -1, sim = TRUE)$overall
  res$f <- f
  res
})
(p1_10 <- pplot(c, d))

# -------------------------------------------------------------------------

f_set <- seq(0.1, 0.9, 0.1)
a <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Con_Equi_JM2(delta_i = c(-0.5, 0), sigma = 4, f_i = c(f, 1 - f), cut = 2, cut_i = c(1.8, 2.0), alpha = 0.025, beta = NA, N = 200, r = 1, sim = FALSE)$overall
  res$f <- f
  res
})
b <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Con_Equi_JM2(delta_i = c(-0.5, 0), sigma = 4, f_i = c(f, 1 - f), cut = 2, cut_i = c(1.8, 2.0), alpha = 0.025, beta = NA, N = 200, r = 1, sim = TRUE)$overall
  res$f <- f
  res
})
(p1_11 <- pplot(a, b))

f_set <- seq(0.1, 0.9, 0.1)
c <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Con_Equi_JM2(delta_i = c(1, 0), sigma = 4, f_i = c(f, 1 - f), cut = 2, alpha = 0.025, beta = 0.2, N = NA, r = 1, sim = FALSE)$overall
  res$f <- f
  res
})
d <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Con_Equi_JM2(delta_i = c(1, 0), sigma = 4, f_i = c(f, 1 - f), cut = 2, alpha = 0.025, beta = 0.2, N = NA, r = 1, sim = TRUE)$overall
  res$f <- f
  res
})
(p1_12 <- pplot(c, d))

# -------------------------------------------------------------------------
# getPwr_Bin_Super_JM1、getPwr_Bin_Noninf_JM1、getPwr_Bin_Equi_JM1
# -------------------------------------------------------------------------

a1 <- getPwr_Bin_Super_JM1(p1_j = 0.7, p0_j = 0.5, p1_a = 0.75, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, N = 200, r = 1, scale = "RD", sim = FALSE)
b1 <- getPwr_Bin_Super_JM1(p1_j = 0.7, p0_j = 0.5, p1_a = 0.75, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, N = 200, r = 1, scale = "RD", sim = TRUE)
(p2_1 <- pplot(a1, b1))

a2 <- getPwr_Bin_Super_JM1(p1_j = 0.7, p0_j = 0.5, p1_a = 0.75, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, N = 200, r = 1, scale = "RR", sim = FALSE)
b2 <- getPwr_Bin_Super_JM1(p1_j = 0.7, p0_j = 0.5, p1_a = 0.75, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, N = 200, r = 1, scale = "RR", sim = TRUE)
(p2_2 <- pplot(a2, b2))

a3 <- getPwr_Bin_Super_JM1(p1_j = 0.7, p0_j = 0.5, p1_a = 0.75, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, N = 200, r = 1, scale = "OR", sim = FALSE)
b3 <- getPwr_Bin_Super_JM1(p1_j = 0.7, p0_j = 0.5, p1_a = 0.75, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, N = 200, r = 1, scale = "OR", sim = TRUE)
(p2_3 <- pplot(a3, b3))

c1 <- getPwr_Bin_Super_JM1(p1_j = 0.3, p0_j = 0.5, p1_nj = 0.25, p0_nj = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RD", sim = FALSE)
d1 <- getPwr_Bin_Super_JM1(p1_j = 0.3, p0_j = 0.5, p1_nj = 0.25, p0_nj = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RD", sim = TRUE)
(p2_4 <- pplot(c1, d1))

c2 <- getPwr_Bin_Super_JM1(p1_j = 0.3, p0_j = 0.5, p1_nj = 0.25, p0_nj = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RR", sim = FALSE)
d2 <- getPwr_Bin_Super_JM1(p1_j = 0.3, p0_j = 0.5, p1_nj = 0.25, p0_nj = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RR", sim = TRUE)
(p2_5 <- pplot(c2, d2))

c3 <- getPwr_Bin_Super_JM1(p1_j = 0.3, p0_j = 0.5, p1_nj = 0.25, p0_nj = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "OR", sim = FALSE)
d3 <- getPwr_Bin_Super_JM1(p1_j = 0.3, p0_j = 0.5, p1_nj = 0.25, p0_nj = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "OR", sim = TRUE)
(p2_6 <- pplot(c3, d3))

# -------------------------------------------------------------------------

a1 <- getPwr_Bin_Noninf_JM1(p1_j = 0.4, p0_j = 0.5, p1_a = 0.5, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 0.3, alpha = 0.025, beta = NA, N = 100, r = 1, scale = "RD", direct = 1, sim = FALSE)
b1 <- getPwr_Bin_Noninf_JM1(p1_j = 0.4, p0_j = 0.5, p1_a = 0.5, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 0.3, alpha = 0.025, beta = NA, N = 100, r = 1, scale = "RD", direct = 1, sim = TRUE)
(p2_7 <- pplot(a1, b1))

a2 <- getPwr_Bin_Noninf_JM1(p1_j = 0.4, p0_j = 0.5, p1_a = 0.5, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = -log(0.6), alpha = 0.025, beta = NA, N = 100, r = 1, scale = "RR", direct = 1, sim = FALSE)
b2 <- getPwr_Bin_Noninf_JM1(p1_j = 0.4, p0_j = 0.5, p1_a = 0.5, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = -log(0.6), alpha = 0.025, beta = NA, N = 100, r = 1, scale = "RR", direct = 1, sim = TRUE)
(p2_8 <- pplot(a2, b2))

a3 <- getPwr_Bin_Noninf_JM1(p1_j = 0.4, p0_j = 0.5, p1_a = 0.5, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = -log(0.5), alpha = 0.025, beta = NA, N = 200, r = 1, scale = "OR", direct = 1, sim = FALSE)
b3 <- getPwr_Bin_Noninf_JM1(p1_j = 0.4, p0_j = 0.5, p1_a = 0.5, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = -log(0.5), alpha = 0.025, beta = NA, N = 200, r = 1, scale = "OR", direct = 1, sim = TRUE)
(p2_9 <- pplot(a3, b3))

c1 <- getPwr_Bin_Noninf_JM1(p1_j = 0.6, p0_j = 0.5, p1_nj = 0.5, p0_nj = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 0.3, alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RD", direct = -1, sim = FALSE)
d1 <- getPwr_Bin_Noninf_JM1(p1_j = 0.6, p0_j = 0.5, p1_nj = 0.5, p0_nj = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 0.3, alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RD", direct = -1, sim = TRUE)
(p2_10 <- pplot(c1, d1))

c2 <- getPwr_Bin_Noninf_JM1(p1_j = 0.6, p0_j = 0.5, p1_nj = 0.5, p0_nj = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = log(1.4), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RR", direct = -1, sim = FALSE)
d2 <- getPwr_Bin_Noninf_JM1(p1_j = 0.6, p0_j = 0.5, p1_nj = 0.5, p0_nj = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = log(1.4), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RR", direct = -1, sim = TRUE)
(p2_11 <- pplot(c2, d2))

c3 <- getPwr_Bin_Noninf_JM1(p1_j = 0.6, p0_j = 0.5, p1_nj = 0.5, p0_nj = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = log(1.7), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "OR", direct = -1, sim = FALSE)
d3 <- getPwr_Bin_Noninf_JM1(p1_j = 0.6, p0_j = 0.5, p1_nj = 0.5, p0_nj = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = log(1.7), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "OR", direct = -1, sim = TRUE)
(p2_12 <- pplot(c3, d3))

# -------------------------------------------------------------------------

a1 <- getPwr_Bin_Equi_JM1(p1_j = 0.4, p0_j = 0.5, p1_a = 0.5, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 0.3, alpha = 0.025, beta = NA, N = 100, r = 1, scale = "RD", sim = FALSE)
b1 <- getPwr_Bin_Equi_JM1(p1_j = 0.4, p0_j = 0.5, p1_a = 0.5, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 0.3, alpha = 0.025, beta = NA, N = 100, r = 1, scale = "RD", sim = TRUE)
(p2_13 <- pplot(a1, b1))

a2 <- getPwr_Bin_Equi_JM1(p1_j = 0.4, p0_j = 0.5, p1_a = 0.5, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = -log(0.6), alpha = 0.025, beta = NA, N = 100, r = 1, scale = "RR", sim = FALSE)
b2 <- getPwr_Bin_Equi_JM1(p1_j = 0.4, p0_j = 0.5, p1_a = 0.5, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = -log(0.6), alpha = 0.025, beta = NA, N = 100, r = 1, scale = "RR", sim = TRUE)
(p2_14 <- pplot(a2, b2))

a3 <- getPwr_Bin_Equi_JM1(p1_j = 0.4, p0_j = 0.5, p1_a = 0.5, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = -log(0.5), alpha = 0.025, beta = NA, N = 200, r = 1, scale = "OR", sim = FALSE)
b3 <- getPwr_Bin_Equi_JM1(p1_j = 0.4, p0_j = 0.5, p1_a = 0.5, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = -log(0.5), alpha = 0.025, beta = NA, N = 200, r = 1, scale = "OR", sim = TRUE)
(p2_15 <- pplot(a3, b3))

c1 <- getPwr_Bin_Equi_JM1(p1_j = 0.6, p0_j = 0.5, p1_nj = 0.5, p0_nj = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 0.3, alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RD", sim = FALSE)
d1 <- getPwr_Bin_Equi_JM1(p1_j = 0.6, p0_j = 0.5, p1_nj = 0.5, p0_nj = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 0.3, alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RD", sim = TRUE)
(p2_16 <- pplot(c1, d1))

c2 <- getPwr_Bin_Equi_JM1(p1_j = 0.6, p0_j = 0.5, p1_nj = 0.5, p0_nj = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = log(1.4), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RR", sim = FALSE)
d2 <- getPwr_Bin_Equi_JM1(p1_j = 0.6, p0_j = 0.5, p1_nj = 0.5, p0_nj = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = log(1.4), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RR", sim = TRUE)
(p2_17 <- pplot(c2, d2))

c3 <- getPwr_Bin_Equi_JM1(p1_j = 0.6, p0_j = 0.5, p1_nj = 0.5, p0_nj = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = log(1.7), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "OR", sim = FALSE)
d3 <- getPwr_Bin_Equi_JM1(p1_j = 0.6, p0_j = 0.5, p1_nj = 0.5, p0_nj = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = log(1.7), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "OR", sim = TRUE)
(p2_18 <- pplot(c3, d3))

# -------------------------------------------------------------------------
# getPwr_Bin_Super_JM2、getPwr_Bin_Noninf_JM2、getPwr_Bin_Equi_JM2
# -------------------------------------------------------------------------

f_set <- seq(0.1, 0.9, 0.1)
a1 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Super_JM2(p1_i = c(0.7, 0.75), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), alpha = 0.025, beta = NA, N = 100, r = 1, scale = "RD", sim = FALSE)$overall
  res$f <- f
  res
})
b1 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Super_JM2(p1_i = c(0.7, 0.75), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), alpha = 0.025, beta = NA, N = 100, r = 1, scale = "RD", sim = TRUE)$overall
  res$f <- f
  res
})
(p2_19 <- pplot(a1, b1))

f_set <- seq(0.1, 0.9, 0.1)
a2 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Super_JM2(p1_i = c(0.7, 0.75), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), alpha = 0.025, beta = NA, N = 100, r = 1, scale = "RR", sim = FALSE)$overall
  res$f <- f
  res
})
b2 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Super_JM2(p1_i = c(0.7, 0.75), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), alpha = 0.025, beta = NA, N = 100, r = 1, scale = "RR", sim = TRUE)$overall
  res$f <- f
  res
})
(p2_20 <- pplot(a2, b2))

f_set <- seq(0.1, 0.9, 0.1)
a3 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Super_JM2(p1_i = c(0.7, 0.75), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), alpha = 0.025, beta = NA, N = 100, r = 1, scale = "OR", sim = FALSE)$overall
  res$f <- f
  res
})
b3 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Super_JM2(p1_i = c(0.7, 0.75), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), alpha = 0.025, beta = NA, N = 100, r = 1, scale = "OR", sim = TRUE)$overall
  res$f <- f
  res
})
(p2_21 <- pplot(a3, b3))

f_set <- seq(0.1, 0.9, 0.1)
c1 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Super_JM2(p1_i = c(0.3, 0.25), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RD", sim = FALSE)$overall
  res$f <- f
  res
})
d1 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Super_JM2(p1_i = c(0.3, 0.25), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RD", sim = TRUE)$overall
  res$f <- f
  res
})
(p2_22 <- pplot(c1, d1))

f_set <- seq(0.1, 0.9, 0.1)
c2 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Super_JM2(p1_i = c(0.3, 0.25), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RR", sim = FALSE)$overall
  res$f <- f
  res
})
d2 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Super_JM2(p1_i = c(0.3, 0.25), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RR", sim = TRUE)$overall
  res$f <- f
  res
})
(p2_23 <- pplot(c2, d2))


f_set <- seq(0.1, 0.9, 0.1)
c3 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Super_JM2(p1_i = c(0.3, 0.25), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "OR", sim = FALSE)$overall
  res$f <- f
  res
})
d3 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Super_JM2(p1_i = c(0.3, 0.25), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "OR", sim = TRUE)$overall
  res$f <- f
  res
})
(p2_24 <- pplot(c3, d3))

# -------------------------------------------------------------------------

f_set <- seq(0.1, 0.9, 0.1)
a1 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Noninf_JM2(p1_i = c(0.4, 0.5), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), cut = 0.25, cut_i = c(0.25, 0.27), alpha = 0.025, beta = NA, N = 100, r = 1, scale = "RD", direct = 1, sim = FALSE)$overall
  res$f <- f
  res
})
b1 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Noninf_JM2(p1_i = c(0.4, 0.5), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), cut = 0.25, cut_i = c(0.25, 0.27), alpha = 0.025, beta = NA, N = 100, r = 1, scale = "RD", direct = 1, sim = TRUE)$overall
  res$f <- f
  res
})
(p2_25 <- pplot(a1, b1))

f_set <- seq(0.1, 0.9, 0.1)
a2 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Noninf_JM2(p1_i = c(0.4, 0.5), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), cut = -log(0.6), cut_i = c(-log(0.6), -log(0.58)), alpha = 0.025, beta = NA, N = 100, r = 1, scale = "RR", direct = 1, sim = FALSE)$overall
  res$f <- f
  res
})
b2 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Noninf_JM2(p1_i = c(0.4, 0.5), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), cut = -log(0.6), cut_i = c(-log(0.6), -log(0.58)), alpha = 0.025, beta = NA, N = 100, r = 1, scale = "RR", direct = 1, sim = TRUE)$overall
  res$f <- f
  res
})
(p2_26 <- pplot(a2, b2))

f_set <- seq(0.1, 0.9, 0.1)
a3 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Noninf_JM2(p1_i = c(0.4, 0.5), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), cut = -log(0.5), cut_i = c(-log(0.6), -log(0.48)), alpha = 0.025, beta = NA, N = 100, r = 1, scale = "OR", direct = 1, sim = FALSE)$overall
  res$f <- f
  res
})
b3 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Noninf_JM2(p1_i = c(0.4, 0.5), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), cut = -log(0.5), cut_i = c(-log(0.6), -log(0.48)), alpha = 0.025, beta = NA, N = 100, r = 1, scale = "OR", direct = 1, sim = TRUE)$overall
  res$f <- f
  res
})
(p2_27 <- pplot(a3, b3))

f_set <- seq(0.1, 0.9, 0.1)
c1 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Noninf_JM2(p1_i = c(0.6, 0.5), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), cut = 0.25, alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RD", direct = -1, sim = FALSE)$overall
  res$f <- f
  res
})
d1 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Noninf_JM2(p1_i = c(0.6, 0.5), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), cut = 0.25, alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RD", direct = -1, sim = TRUE)$overall
  res$f <- f
  res
})
(p2_28 <- pplot(c1, d1))

f_set <- seq(0.1, 0.9, 0.1)
c2 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Noninf_JM2(p1_i = c(0.6, 0.5), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), cut = log(1.4), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RR", direct = -1, sim = FALSE)$overall
  res$f <- f
  res
})
d2 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Noninf_JM2(p1_i = c(0.6, 0.5), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), cut = log(1.4), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RR", direct = -1, sim = TRUE)$overall
  res$f <- f
  res
})
(p2_29 <- pplot(c2, d2))

f_set <- seq(0.1, 0.9, 0.1)
c3 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Noninf_JM2(p1_i = c(0.6, 0.5), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), cut = log(1.7), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "OR", direct = -1, sim = FALSE)$overall
  res$f <- f
  res
})
d3 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Noninf_JM2(p1_i = c(0.6, 0.5), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), cut = log(1.7), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "OR", direct = -1, sim = TRUE)$overall
  res$f <- f
  res
})
(p2_30 <- pplot(c3, d3))

# -------------------------------------------------------------------------

f_set <- seq(0.1, 0.9, 0.1)
a1 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Equi_JM2(p1_i = c(0.4, 0.5), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), cut = 0.25, cut_i = c(0.25, 0.27), alpha = 0.025, beta = NA, N = 100, r = 1, scale = "RD", sim = FALSE)$overall
  res$f <- f
  res
})
b1 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Equi_JM2(p1_i = c(0.4, 0.5), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), cut = 0.25, cut_i = c(0.25, 0.27), alpha = 0.025, beta = NA, N = 100, r = 1, scale = "RD", sim = TRUE)$overall
  res$f <- f
  res
})
(p2_31 <- pplot(a1, b1))

f_set <- seq(0.1, 0.9, 0.1)
a2 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Equi_JM2(p1_i = c(0.4, 0.5), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), cut = -log(0.6), cut_i = c(-log(0.6), -log(0.58)), alpha = 0.025, beta = NA, N = 100, r = 1, scale = "RR", sim = FALSE)$overall
  res$f <- f
  res
})
b2 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Equi_JM2(p1_i = c(0.4, 0.5), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), cut = -log(0.6), cut_i = c(-log(0.6), -log(0.58)), alpha = 0.025, beta = NA, N = 100, r = 1, scale = "RR", sim = TRUE)$overall
  res$f <- f
  res
})
(p2_32 <- pplot(a2, b2))

f_set <- seq(0.1, 0.9, 0.1)
a3 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Equi_JM2(p1_i = c(0.4, 0.5), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), cut = -log(0.5), cut_i = c(-log(0.5), -log(0.48)), alpha = 0.025, beta = NA, N = 200, r = 1, scale = "OR", sim = FALSE)$overall
  res$f <- f
  res
})
b3 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Equi_JM2(p1_i = c(0.4, 0.5), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), cut = -log(0.5), cut_i = c(-log(0.5), -log(0.48)), alpha = 0.025, beta = NA, N = 200, r = 1, scale = "OR", sim = TRUE)$overall
  res$f <- f
  res
})
(p2_33 <- pplot(a3, b3))

f_set <- seq(0.1, 0.9, 0.1)
c1 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Equi_JM2(p1_i = c(0.6, 0.5), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), cut = 0.25, alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RD", sim = FALSE)$overall
  res$f <- f
  res
})
d1 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Equi_JM2(p1_i = c(0.6, 0.5), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), cut = 0.25, alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RD", sim = TRUE)$overall
  res$f <- f
  res
})
(p2_34 <- pplot(c1, d1))

f_set <- seq(0.1, 0.9, 0.1)
c2 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Equi_JM2(p1_i = c(0.6, 0.5), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), cut = log(1.4), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RR", sim = FALSE)$overall
  res$f <- f
  res
})
d2 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Equi_JM2(p1_i = c(0.6, 0.5), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), cut = log(1.4), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RR", sim = TRUE)$overall
  res$f <- f
  res
})
(p2_35 <- pplot(c2, d2))

f_set <- seq(0.1, 0.9, 0.1)
c3 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Equi_JM2(p1_i = c(0.6, 0.5), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), cut = log(1.7), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "OR", sim = FALSE)$overall
  res$f <- f
  res
})
d3 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Equi_JM2(p1_i = c(0.6, 0.5), p0_i = c(0.5, 0.5), f_i = c(f, 1 - f), cut = log(1.7), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "OR", sim = TRUE)$overall
  res$f <- f
  res
})
(p2_36 <- pplot(c3, d3))

# -------------------------------------------------------------------------
# getPwr_Surv_Super_JM1、getPwr_Surv_Noninf_JM1、getPwr_Surv_Equi_JM1
# -------------------------------------------------------------------------

a1 <- getPwr_Surv_Super_JM1(delta_j = log(1.3), delta_nj = log(1.4), f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, Ne = 200, r = 1, criterion = 1, sim = FALSE)
b1 <- getPwr_Surv_Super_JM1(delta_j = log(1.3), delta_nj = log(1.4), f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, Ne = 200, r = 1, criterion = 1, sim = TRUE)
(p3_1 <- pplot(a1, b1))

a2 <- getPwr_Surv_Super_JM1(delta_j = log(1.3), delta_a = log(1.4), f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, Ne = 200, r = 1, criterion = 2, sim = FALSE)
b2 <- getPwr_Surv_Super_JM1(delta_j = log(1.3), delta_a = log(1.4), f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, Ne = 200, r = 1, criterion = 2, sim = TRUE)
(p3_2 <- pplot(a2, b2))

c1 <- getPwr_Surv_Super_JM1(delta_j = log(0.8), delta_a = log(0.7), f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = 0.2, Ne = NA, r = 1, criterion = 1, sim = FALSE)
d1 <- getPwr_Surv_Super_JM1(delta_j = log(0.8), delta_a = log(0.7), f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = 0.2, Ne = NA, r = 1, criterion = 1, sim = TRUE)
(p3_3 <- pplot(c1, d1))

c2 <- getPwr_Surv_Super_JM1(delta_j = log(0.8), delta_nj = log(0.7), f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = 0.2, Ne = NA, r = 1, criterion = 2, sim = FALSE)
d2 <- getPwr_Surv_Super_JM1(delta_j = log(0.8), delta_nj = log(0.7), f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = 0.2, Ne = NA, r = 1, criterion = 2, sim = TRUE)
(p3_4 <- pplot(c2, d2))

# -------------------------------------------------------------------------

a1 <- getPwr_Surv_Noninf_JM1(delta_j = log(0.9), delta_nj = log(1.0), f = seq(0.1, 0.9, 0.1), cut = -log(0.7), pi = 0.5, alpha = 0.025, beta = NA, Ne = 400, r = 1, criterion = 1, direct = 1, sim = FALSE)
b1 <- getPwr_Surv_Noninf_JM1(delta_j = log(0.9), delta_nj = log(1.0), f = seq(0.1, 0.9, 0.1), cut = -log(0.7), pi = 0.5, alpha = 0.025, beta = NA, Ne = 400, r = 1, criterion = 1, direct = 1, sim = TRUE)
(p3_5 <- pplot(a1, b1))

a2 <- getPwr_Surv_Noninf_JM1(delta_j = log(0.9), delta_a = log(1.0), f = seq(0.1, 0.9, 0.1), cut = -log(0.7), pi = 0.5, alpha = 0.025, beta = NA, Ne = 400, r = 1, criterion = 2, direct = 1, sim = FALSE)
b2 <- getPwr_Surv_Noninf_JM1(delta_j = log(0.9), delta_a = log(1.0), f = seq(0.1, 0.9, 0.1), cut = -log(0.7), pi = 0.5, alpha = 0.025, beta = NA, Ne = 400, r = 1, criterion = 2, direct = 1, sim = TRUE)
(p3_6 <- pplot(a2, b2))

c1 <- getPwr_Surv_Noninf_JM1(delta_j = log(1.1), delta_a = log(1.0), f = seq(0.1, 0.9, 0.1), cut = log(1.3), pi = 0.5, alpha = 0.025, beta = 0.2, Ne = NA, r = 1, criterion = 1, direct = -1, sim = FALSE)
d1 <- getPwr_Surv_Noninf_JM1(delta_j = log(1.1), delta_a = log(1.0), f = seq(0.1, 0.9, 0.1), cut = log(1.3), pi = 0.5, alpha = 0.025, beta = 0.2, Ne = NA, r = 1, criterion = 1, direct = -1, sim = TRUE)
(p3_7 <- pplot(c1, d1))

c2 <- getPwr_Surv_Noninf_JM1(delta_j = log(1.1), delta_nj = log(1.0), f = seq(0.1, 0.9, 0.1), cut = log(1.3), pi = 0.5, alpha = 0.025, beta = 0.2, Ne = NA, r = 1, criterion = 2, direct = -1, sim = FALSE)
d2 <- getPwr_Surv_Noninf_JM1(delta_j = log(1.1), delta_nj = log(1.0), f = seq(0.1, 0.9, 0.1), cut = log(1.3), pi = 0.5, alpha = 0.025, beta = 0.2, Ne = NA, r = 1, criterion = 2, direct = -1, sim = TRUE)
(p3_8 <- pplot(c2, d2))

# -------------------------------------------------------------------------

a1 <- getPwr_Surv_Equi_JM1(delta_j = log(0.9), delta_nj = log(1.0), f = seq(0.1, 0.9, 0.1), cut = -log(0.7), pi = 0.5, alpha = 0.025, beta = NA, Ne = 400, r = 1, criterion = 1, sim = FALSE)
b1 <- getPwr_Surv_Equi_JM1(delta_j = log(0.9), delta_nj = log(1.0), f = seq(0.1, 0.9, 0.1), cut = -log(0.7), pi = 0.5, alpha = 0.025, beta = NA, Ne = 400, r = 1, criterion = 1, sim = TRUE)
(p3_9 <- pplot(a1, b1))

a2 <- getPwr_Surv_Equi_JM1(delta_j = log(0.9), delta_a = log(1.0), f = seq(0.1, 0.9, 0.1), cut = -log(0.7), pi = 0.5, alpha = 0.025, beta = NA, Ne = 400, r = 1, criterion = 2, sim = FALSE)
b2 <- getPwr_Surv_Equi_JM1(delta_j = log(0.9), delta_a = log(1.0), f = seq(0.1, 0.9, 0.1), cut = -log(0.7), pi = 0.5, alpha = 0.025, beta = NA, Ne = 400, r = 1, criterion = 2, sim = TRUE)
(p3_10 <- pplot(a2, b2))

c1 <- getPwr_Surv_Equi_JM1(delta_j = log(1.1), delta_a = log(1.0), f = seq(0.1, 0.9, 0.1), cut = log(1.3), pi = 0.5, alpha = 0.025, beta = 0.2, Ne = NA, r = 1, criterion = 1, sim = FALSE)
d1 <- getPwr_Surv_Equi_JM1(delta_j = log(1.1), delta_a = log(1.0), f = seq(0.1, 0.9, 0.1), cut = log(1.3), pi = 0.5, alpha = 0.025, beta = 0.2, Ne = NA, r = 1, criterion = 1, sim = TRUE)
(p3_11 <- pplot(c1, d1))

c2 <- getPwr_Surv_Equi_JM1(delta_j = log(1.1), delta_nj = log(1.0), f = seq(0.1, 0.9, 0.1), cut = log(1.3), pi = 0.5, alpha = 0.025, beta = 0.2, Ne = NA, r = 1, criterion = 2, sim = FALSE)
d2 <- getPwr_Surv_Equi_JM1(delta_j = log(1.1), delta_nj = log(1.0), f = seq(0.1, 0.9, 0.1), cut = log(1.3), pi = 0.5, alpha = 0.025, beta = 0.2, Ne = NA, r = 1, criterion = 2, sim = TRUE)
(p3_12 <- pplot(c2, d2))

# -------------------------------------------------------------------------
# getPwr_Surv_Super_JM2、getPwr_Surv_Noninf_JM2、getPwr_Surv_Equi_JM2
# -------------------------------------------------------------------------

f_set <- seq(0.1, 0.9, 0.1)
a <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Surv_Super_JM2(delta_i = c(log(1.2), log(1.4)), f_i = c(f, 1 - f), alpha = 0.025, beta = NA, Ne = 300, r = 1, sim = FALSE)$overall
  res$f <- f
  res
})
b <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Surv_Super_JM2(delta_i = c(log(1.2), log(1.4)), f_i = c(f, 1 - f), alpha = 0.025, beta = NA, Ne = 300, r = 1, sim = TRUE)$overall
  res$f <- f
  res
})
(p3_13 <- pplot(a, b))

f_set <- seq(0.1, 0.9, 0.1)
c <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Surv_Super_JM2(delta_i = c(log(0.7), log(0.8)), f_i = c(f, 1 - f), alpha = 0.025, beta = 0.2, Ne = NA, r = 1, sim = FALSE)$overall
  res$f <- f
  res
})
d <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Surv_Super_JM2(delta_i = c(log(0.7), log(0.8)), f_i = c(f, 1 - f), alpha = 0.025, beta = 0.2, Ne = NA, r = 1, sim = TRUE)$overall
  res$f <- f
  res
})
(p3_14 <- pplot(c, d))

# -------------------------------------------------------------------------

f_set <- seq(0.1, 0.9, 0.1)
a <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Surv_Noninf_JM2(delta_i = c(log(0.9), log(1.0)), f_i = c(f, 1 - f), cut = -log(0.7), cut_i = c(-log(0.7), -log(0.68)), alpha = 0.025, beta = NA, Ne = 300, r = 1, direct = 1, sim = FALSE)$overall
  res$f <- f
  res
})
b <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Surv_Noninf_JM2(delta_i = c(log(0.9), log(1.0)), f_i = c(f, 1 - f), cut = -log(0.7), cut_i = c(-log(0.7), -log(0.68)), alpha = 0.025, beta = NA, Ne = 300, r = 1, direct = 1, sim = TRUE)$overall
  res$f <- f
  res
})
(p3_15 <- pplot(a, b))


f_set <- seq(0.1, 0.9, 0.1)
c <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Surv_Noninf_JM2(delta_i = c(log(1.1), log(1.0)), f_i = c(f, 1 - f), cut = log(1.3), alpha = 0.025, beta = 0.2, Ne = NA, r = 1, direct = -1, sim = FALSE)$overall
  res$f <- f
  res
})
d <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Surv_Noninf_JM2(delta_i = c(log(1.1), log(1.0)), f_i = c(f, 1 - f), cut = log(1.3), alpha = 0.025, beta = 0.2, Ne = NA, r = 1, direct = -1, sim = TRUE)$overall
  res$f <- f
  res
})
(p3_16 <- pplot(c, d))

# -------------------------------------------------------------------------

f_set <- seq(0.1, 0.9, 0.1)
a <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Surv_Equi_JM2(delta_i = c(log(0.9), log(1.0)), f_i = c(f, 1 - f), cut = -log(0.7), cut_i = c(-log(0.7), -log(0.68)), alpha = 0.025, beta = NA, Ne = 300, r = 1, sim = FALSE)$overall
  res$f <- f
  res
})
b <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Surv_Equi_JM2(delta_i = c(log(0.9), log(1.0)), f_i = c(f, 1 - f), cut = -log(0.7), cut_i = c(-log(0.7), -log(0.68)), alpha = 0.025, beta = NA, Ne = 300, r = 1, sim = TRUE)$overall
  res$f <- f
  res
})
(p3_17 <- pplot(a, b))

f_set <- seq(0.1, 0.9, 0.1)
c <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Surv_Equi_JM2(delta_i = c(log(1.1), log(1.0)), f_i = c(f, 1 - f), cut = log(1.3), alpha = 0.025, beta = 0.2, Ne = NA, r = 1, sim = FALSE)$overall
  res$f <- f
  res
})
d <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Surv_Equi_JM2(delta_i = c(log(1.1), log(1.0)), f_i = c(f, 1 - f), cut = log(1.3), alpha = 0.025, beta = 0.2, Ne = NA, r = 1, sim = TRUE)$overall
  res$f <- f
  res
})
(p3_18 <- pplot(c, d))

# -------------------------------------------------------------------------
# getPwr_Count_Super_JM1、getPwr_Count_Noninf_JM1、getPwr_Count_Equi_JM1
# -------------------------------------------------------------------------

a1 <- getPwr_Count_Super_JM1(delta_j = log(1.2), delta_a = log(1.3), lambda0_j = 0.1, lambda0_a = 0.1, t_j = 5, t_a = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, N = 300, r = 1, sim = FALSE)
b1 <- getPwr_Count_Super_JM1(delta_j = log(1.2), delta_a = log(1.3), lambda0_j = 0.1, lambda0_a = 0.1, t_j = 5, t_a = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, N = 300, r = 1, sim = TRUE)
(p4_1 <- pplot(a1, b1))

a2 <- getPwr_Count_Super_JM1(delta_j = log(1.2), delta_nj = log(1.3), lambda0_j = 0.1, lambda0_nj = 0.1, t_j = 5, t_nj = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, N = 300, r = 1, sim = FALSE)
b2 <- getPwr_Count_Super_JM1(delta_j = log(1.2), delta_nj = log(1.3), lambda0_j = 0.1, lambda0_nj = 0.1, t_j = 5, t_nj = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, N = 300, r = 1, sim = TRUE)
(p4_2 <- pplot(a2, b2))

c1 <- getPwr_Count_Super_JM1(delta_j = log(0.8), delta_a = log(0.7), lambda0_j = 0.1, lambda0_a = 0.1, t_j = 5, t_a = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = 0.2, N = NA, r = 1, sim = FALSE)
d1 <- getPwr_Count_Super_JM1(delta_j = log(0.8), delta_a = log(0.7), lambda0_j = 0.1, lambda0_a = 0.1, t_j = 5, t_a = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = 0.2, N = NA, r = 1, sim = TRUE)
(p4_3 <- pplot(c1, d1))

c2 <- getPwr_Count_Super_JM1(delta_j = log(0.8), delta_nj = log(0.7), lambda0_j = 0.1, lambda0_nj = 0.1, t_j = 5, t_nj = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = 0.2, N = NA, r = 1, sim = FALSE)
d2 <- getPwr_Count_Super_JM1(delta_j = log(0.8), delta_nj = log(0.7), lambda0_j = 0.1, lambda0_nj = 0.1, t_j = 5, t_nj = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = 0.2, N = NA, r = 1, sim = TRUE)
(p4_4 <- pplot(c2, d2))

# -------------------------------------------------------------------------

a1 <- getPwr_Count_Noninf_JM1(delta_j = log(0.9), delta_a = log(1.0), lambda0_j = 0.1, lambda0_a = 0.1, t_j = 5, t_a = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = -log(0.7), alpha = 0.025, beta = NA, N = 400, r = 1, direct = 1, sim = FALSE)
b1 <- getPwr_Count_Noninf_JM1(delta_j = log(0.9), delta_a = log(1.0), lambda0_j = 0.1, lambda0_a = 0.1, t_j = 5, t_a = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = -log(0.7), alpha = 0.025, beta = NA, N = 400, r = 1, direct = 1, sim = TRUE)
(p4_5 <- pplot(a1, b1))

a2 <- getPwr_Count_Noninf_JM1(delta_j = log(0.9), delta_nj = log(1.0), lambda0_j = 0.1, lambda0_nj = 0.1, t_j = 5, t_nj = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = -log(0.7), alpha = 0.025, beta = NA, N = 400, r = 1, direct = 1, sim = FALSE)
b2 <- getPwr_Count_Noninf_JM1(delta_j = log(0.9), delta_nj = log(1.0), lambda0_j = 0.1, lambda0_nj = 0.1, t_j = 5, t_nj = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = -log(0.7), alpha = 0.025, beta = NA, N = 400, r = 1, direct = 1, sim = TRUE)
(p4_6 <- pplot(a2, b2))

c1 <- getPwr_Count_Noninf_JM1(delta_j = log(1.1), delta_a = log(1.0), lambda0_j = 0.1, lambda0_a = 0.1, t_j = 5, t_a = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = log(1.3), alpha = 0.025, beta = 0.2, N = NA, r = 1, direct = -1, sim = FALSE)
d1 <- getPwr_Count_Noninf_JM1(delta_j = log(1.1), delta_a = log(1.0), lambda0_j = 0.1, lambda0_a = 0.1, t_j = 5, t_a = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = log(1.3), alpha = 0.025, beta = 0.2, N = NA, r = 1, direct = -1, sim = TRUE)
(p4_7 <- pplot(c1, d1))

c1 <- getPwr_Count_Noninf_JM1(delta_j = log(1.1), delta_nj = log(1.0), lambda0_j = 0.1, lambda0_nj = 0.1, t_j = 5, t_nj = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = log(1.3), alpha = 0.025, beta = 0.2, N = NA, r = 1, direct = -1, sim = FALSE)
d1 <- getPwr_Count_Noninf_JM1(delta_j = log(1.1), delta_nj = log(1.0), lambda0_j = 0.1, lambda0_nj = 0.1, t_j = 5, t_nj = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = log(1.3), alpha = 0.025, beta = 0.2, N = NA, r = 1, direct = -1, sim = TRUE)
(p4_8 <- pplot(c1, d1))

# -------------------------------------------------------------------------

a1 <- getPwr_Count_Equi_JM1(delta_j = log(0.9), delta_a = log(1.0), lambda0_j = 0.1, lambda0_a = 0.1, t_j = 5, t_a = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = -log(0.7), alpha = 0.025, beta = NA, N = 400, r = 1, sim = FALSE)
b1 <- getPwr_Count_Equi_JM1(delta_j = log(0.9), delta_a = log(1.0), lambda0_j = 0.1, lambda0_a = 0.1, t_j = 5, t_a = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = -log(0.7), alpha = 0.025, beta = NA, N = 400, r = 1, sim = TRUE)
(p4_9 <- pplot(a1, b1))

a2 <- getPwr_Count_Equi_JM1(delta_j = log(0.9), delta_nj = log(1.0), lambda0_j = 0.1, lambda0_nj = 0.1, t_j = 5, t_nj = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = -log(0.7), alpha = 0.025, beta = NA, N = 400, r = 1, sim = FALSE)
b2 <- getPwr_Count_Equi_JM1(delta_j = log(0.9), delta_nj = log(1.0), lambda0_j = 0.1, lambda0_nj = 0.1, t_j = 5, t_nj = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = -log(0.7), alpha = 0.025, beta = NA, N = 400, r = 1, sim = TRUE)
(p4_10 <- pplot(a2, b2))

c1 <- getPwr_Count_Equi_JM1(delta_j = log(1.1), delta_a = log(1.0), lambda0_j = 0.1, lambda0_a = 0.1, t_j = 5, t_a = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = log(1.3), alpha = 0.025, beta = 0.2, N = NA, r = 1, sim = FALSE)
d1 <- getPwr_Count_Equi_JM1(delta_j = log(1.1), delta_a = log(1.0), lambda0_j = 0.1, lambda0_a = 0.1, t_j = 5, t_a = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = log(1.3), alpha = 0.025, beta = 0.2, N = NA, r = 1, sim = TRUE)
(p4_11 <- pplot(c1, d1))

c1 <- getPwr_Count_Equi_JM1(delta_j = log(1.1), delta_nj = log(1.0), lambda0_j = 0.1, lambda0_nj = 0.1, t_j = 5, t_nj = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = log(1.3), alpha = 0.025, beta = 0.2, N = NA, r = 1, sim = FALSE)
d1 <- getPwr_Count_Equi_JM1(delta_j = log(1.1), delta_nj = log(1.0), lambda0_j = 0.1, lambda0_nj = 0.1, t_j = 5, t_nj = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = log(1.3), alpha = 0.025, beta = 0.2, N = NA, r = 1, sim = TRUE)
(p4_12 <- pplot(c1, d1))

# -------------------------------------------------------------------------
# getPwr_Count_Super_JM2、getPwr_Count_Noninf_JM2、getPwr_Count_Equi_JM2
# -------------------------------------------------------------------------

f_set <- seq(0.1, 0.9, 0.1)
a <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Count_Super_JM2(delta_i = c(log(1.2), log(1.4)), lambda0_i = c(0.1, 0.1), t_i = c(5, 5), k = 0, f_i = c(f, 1 - f), alpha = 0.025, beta = NA, N = 300, r = 1, sim = FALSE)$overall
  res$f <- f
  res
})
b <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Count_Super_JM2(delta_i = c(log(1.2), log(1.4)), lambda0_i = c(0.1, 0.1), t_i = c(5, 5), k = 0, f_i = c(f, 1 - f), alpha = 0.025, beta = NA, N = 300, r = 1, sim = TRUE)$overall
  res$f <- f
  res
})
(p4_13 <- pplot(a, b))

f_set <- seq(0.1, 0.9, 0.1)
c <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Count_Super_JM2(delta_i = c(log(0.8), log(0.6)), lambda0_i = c(0.1, 0.1), t_i = c(5, 5), k = 0, f_i = c(f, 1 - f), alpha = 0.025, beta = 0.2, N = NA, r = 1, sim = FALSE)$overall
  res$f <- f
  res
})
d <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Count_Super_JM2(delta_i = c(log(0.8), log(0.6)), lambda0_i = c(0.1, 0.1), t_i = c(5, 5), k = 0, f_i = c(f, 1 - f), alpha = 0.025, beta = 0.2, N = NA, r = 1, sim = TRUE)$overall
  res$f <- f
  res
})
(p4_14 <- pplot(c, d))

# -------------------------------------------------------------------------

f_set <- seq(0.1, 0.9, 0.1)
a <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Count_Noninf_JM2(delta_i = c(log(0.9), log(1.0)), lambda0_i = c(0.1, 0.1), t_i = c(5, 5), k = 0, f_i = c(f, 1 - f), cut = -log(0.7), cut_i = c(cut = -log(0.7), cut = -log(0.68)), alpha = 0.025, beta = NA, N = 300, r = 1, direct = 1, sim = FALSE)$overall
  res$f <- f
  res
})
b <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Count_Noninf_JM2(delta_i = c(log(0.9), log(1.0)), lambda0_i = c(0.1, 0.1), t_i = c(5, 5), k = 0, f_i = c(f, 1 - f), cut = -log(0.7), cut_i = c(cut = -log(0.7), cut = -log(0.68)), alpha = 0.025, beta = NA, N = 300, r = 1, direct = 1, sim = TRUE)$overall
  res$f <- f
  res
})
(p4_15 <- pplot(a, b))

f_set <- seq(0.1, 0.9, 0.1)
c <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Count_Noninf_JM2(delta_i = c(log(1.1), log(1.0)), lambda0_i = c(0.1, 0.1), t_i = c(5, 5), k = 0, f_i = c(f, 1 - f), cut = log(1.3), alpha = 0.025, beta = 0.2, N = NA, r = 1, direct = -1, sim = FALSE)$overall
  res$f <- f
  res
})
d <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Count_Noninf_JM2(delta_i = c(log(1.1), log(1.0)), lambda0_i = c(0.1, 0.1), t_i = c(5, 5), k = 0, f_i = c(f, 1 - f), cut = log(1.3), alpha = 0.025, beta = 0.2, N = NA, r = 1, direct = -1, sim = TRUE)$overall
  res$f <- f
  res
})
(p4_16 <- pplot(c, d))

# -------------------------------------------------------------------------

f_set <- seq(0.1, 0.9, 0.1)
a <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Count_Equi_JM2(delta_i = c(log(0.9), log(1.0)), lambda0_i = c(0.1, 0.1), t_i = c(5, 5), k = 0, f_i = c(f, 1 - f), cut = -log(0.7), cut_i = c(cut = -log(0.7), cut = -log(0.68)), alpha = 0.025, beta = NA, N = 300, r = 1, sim = FALSE)$overall
  res$f <- f
  res
})
b <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Count_Equi_JM2(delta_i = c(log(0.9), log(1.0)), lambda0_i = c(0.1, 0.1), t_i = c(5, 5), k = 0, f_i = c(f, 1 - f), cut = -log(0.7), cut_i = c(cut = -log(0.7), cut = -log(0.68)), alpha = 0.025, beta = NA, N = 300, r = 1, sim = TRUE)$overall
  res$f <- f
  res
})
(p4_17 <- pplot(a, b))

f_set <- seq(0.1, 0.9, 0.1)
c <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Count_Equi_JM2(delta_i = c(log(1.1), log(1.0)), lambda0_i = c(0.1, 0.1), t_i = c(5, 5), k = 0, f_i = c(f, 1 - f), cut = log(1.3), alpha = 0.025, beta = 0.2, N = NA, r = 1, sim = FALSE)$overall
  res$f <- f
  res
})
d <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Count_Equi_JM2(delta_i = c(log(1.1), log(1.0)), lambda0_i = c(0.1, 0.1), t_i = c(5, 5), k = 0, f_i = c(f, 1 - f), cut = log(1.3), alpha = 0.025, beta = 0.2, N = NA, r = 1, sim = TRUE)$overall
  res$f <- f
  res
})
(p4_18 <- pplot(c, d))

# -------------------------------------------------------------------------
# getN_Con_Super_JM1、getN_Con_Noninf_JM1、getN_Con_Equi_JM1
# -------------------------------------------------------------------------

v <- getPwr_Con_Super_JM1(delta_j = 0.5, delta_a = 0.7, sigma = 1, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, N = 100, r = 1, sim = FALSE)
getN_Con_Super_JM1(delta_a = 0.7, delta_j = 0.5, sigma = 1, pi = 0.5, beta1 = 1 - v$pwr2, N = 100, r = 1)

v <- getPwr_Con_Super_JM1(delta_j = -0.5, delta_a = -0.7, sigma = 1, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, N = 100, r = 1, sim = FALSE)
getN_Con_Super_JM1(delta_a = -0.7, delta_j = -0.5, sigma = 1, pi = 0.5, beta1 = 1 - v$pwr2, N = 100, r = 1)

# -------------------------------------------------------------------------

v <- getPwr_Con_Noninf_JM1(delta_j = -0.5, delta_a = 0, sigma = 1, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 2, alpha = 0.025, beta = NA, N = 100, r = 1, direct = 1, sim = FALSE)
getN_Con_Noninf_JM1(delta_a = 0, delta_j = -0.5, sigma = 1, pi = 0.5, cut = 2, beta1 = 1 - v$pwr2, N = 100, r = 1, direct = 1)

v <- getPwr_Con_Noninf_JM1(delta_j = 0.5, delta_a = 0, sigma = 1, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 2, alpha = 0.025, beta = NA, N = 100, r = 1, direct = -1, sim = FALSE)
getN_Con_Noninf_JM1(delta_a = 0, delta_j = 0.5, sigma = 1, pi = 0.5, cut = 2, beta1 = 1 - v$pwr2, N = 100, r = 1, direct = -1)

# -------------------------------------------------------------------------

v <- getPwr_Con_Equi_JM1(delta_j = -0.5, delta_a = 0, sigma = 1, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 2, alpha = 0.025, beta = NA, N = 100, r = 1, sim = FALSE)
getN_Con_Equi_JM1(delta_a = 0, delta_j = -0.5, sigma = 1, pi = 0.5, cut = 2, beta1 = 1 - v$pwr2, N = 100, r = 1)

v <- getPwr_Con_Equi_JM1(delta_j = 0.5, delta_a = 0, sigma = 1, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 2, alpha = 0.025, beta = NA, N = 100, r = 1, sim = FALSE)
getN_Con_Equi_JM1(delta_a = 0, delta_j = 0.5, sigma = 1, pi = 0.5, cut = 2, beta1 = 1 - v$pwr2, N = 100, r = 1)

# -------------------------------------------------------------------------
# getN_Bin_Super_JM1、getN_Bin_Noninf_JM1、getN_Bin_Equi_JM1
# -------------------------------------------------------------------------

v <- getPwr_Bin_Super_JM1(p1_j = 0.7, p0_j = 0.5, p1_a = 0.75, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, N = 200, r = 1, scale = "RD", sim = FALSE)
getN_Bin_Super_JM1(p1_a = 0.75, p0_a = 0.5, p1_j = 0.7, p0_j = 0.5, pi = 0.5, alpha = NA, beta = NA, beta1 = 1 - v$pwr2, N = 200, r = 1, scale = "RD")

v <- getPwr_Bin_Super_JM1(p1_j = 0.7, p0_j = 0.5, p1_a = 0.75, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, N = 200, r = 1, scale = "RR", sim = FALSE)
getN_Bin_Super_JM1(p1_a = 0.75, p0_a = 0.5, p1_j = 0.7, p0_j = 0.5, pi = 0.5, alpha = NA, beta = NA, beta1 = 1 - v$pwr2, N = 200, r = 1, scale = "RR")

v <- getPwr_Bin_Super_JM1(p1_j = 0.7, p0_j = 0.5, p1_a = 0.75, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, N = 200, r = 1, scale = "OR", sim = FALSE)
getN_Bin_Super_JM1(p1_a = 0.75, p0_a = 0.5, p1_j = 0.7, p0_j = 0.5, pi = 0.5, alpha = NA, beta = NA, beta1 = 1 - v$pwr2, N = 200, r = 1, scale = "OR")


v <- getPwr_Bin_Super_JM1(p1_j = 0.3, p0_j = 0.5, p1_a = 0.25, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, N = 200, r = 1, scale = "RD", sim = FALSE)
getN_Bin_Super_JM1(p1_a = 0.25, p0_a = 0.5, p1_j = 0.3, p0_j = 0.5, pi = 0.5, alpha = NA, beta = NA, beta1 = 1 - v$pwr2, N = 200, r = 1, scale = "RD")

v <- getPwr_Bin_Super_JM1(p1_j = 0.3, p0_j = 0.5, p1_a = 0.25, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, N = 200, r = 1, scale = "RR", sim = FALSE)
getN_Bin_Super_JM1(p1_a = 0.25, p0_a = 0.5, p1_j = 0.3, p0_j = 0.5, pi = 0.5, alpha = NA, beta = NA, beta1 = 1 - v$pwr2, N = 200, r = 1, scale = "RR")

v <- getPwr_Bin_Super_JM1(p1_j = 0.3, p0_j = 0.5, p1_a = 0.25, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, N = 200, r = 1, scale = "OR", sim = FALSE)
getN_Bin_Super_JM1(p1_a = 0.25, p0_a = 0.5, p1_j = 0.3, p0_j = 0.5, pi = 0.5, alpha = NA, beta = NA, beta1 = 1 - v$pwr2, N = 200, r = 1, scale = "OR")

# -------------------------------------------------------------------------

v <- getPwr_Bin_Noninf_JM1(p1_j = 0.4, p0_j = 0.5, p1_a = 0.5, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 0.3, alpha = 0.025, beta = NA, N = 200, r = 1, scale = "RD", direct = 1, sim = FALSE)
getN_Bin_Noninf_JM1(p1_a = 0.5, p0_a = 0.5, p1_j = 0.4, p0_j = 0.5, pi = 0.5, cut = 0.3, alpha = NA, beta = NA, beta1 = 1 - v$pwr2, N = 200, r = 1, scale = "RD", direct = 1)

v <- getPwr_Bin_Noninf_JM1(p1_j = 0.4, p0_j = 0.5, p1_a = 0.5, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = -log(0.6), alpha = 0.025, beta = NA, N = 200, r = 1, scale = "RR", direct = 1, sim = FALSE)
getN_Bin_Noninf_JM1(p1_a = 0.5, p0_a = 0.5, p1_j = 0.4, p0_j = 0.5, pi = 0.5, cut = -log(0.6), alpha = NA, beta = NA, beta1 = 1 - v$pwr2, N = 200, r = 1, scale = "RR", direct = 1)

v <- getPwr_Bin_Noninf_JM1(p1_j = 0.4, p0_j = 0.5, p1_a = 0.5, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = -log(0.5), alpha = 0.025, beta = NA, N = 200, r = 1, scale = "OR", direct = 1, sim = FALSE)
getN_Bin_Noninf_JM1(p1_a = 0.5, p0_a = 0.5, p1_j = 0.4, p0_j = 0.5, pi = 0.5, cut = -log(0.5), alpha = NA, beta = NA, beta1 = 1 - v$pwr2, N = 200, r = 1, scale = "OR", direct = 1)


v <- getPwr_Bin_Noninf_JM1(p1_j = 0.6, p0_j = 0.5, p1_a = 0.5, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 0.3, alpha = 0.025, beta = NA, N = 200, r = 1, scale = "RD", direct = -1, sim = FALSE)
getN_Bin_Noninf_JM1(p1_a = 0.5, p0_a = 0.5, p1_j = 0.6, p0_j = 0.5, pi = 0.5, cut = 0.3, alpha = NA, beta = NA, beta1 = 1 - v$pwr2, N = 200, r = 1, scale = "RD", direct = -1)

v <- getPwr_Bin_Noninf_JM1(p1_j = 0.6, p0_j = 0.5, p1_a = 0.5, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = log(1.4), alpha = 0.025, beta = NA, N = 200, r = 1, scale = "RR", direct = -1, sim = FALSE)
getN_Bin_Noninf_JM1(p1_a = 0.5, p0_a = 0.5, p1_j = 0.6, p0_j = 0.5, pi = 0.5, cut = log(1.4), alpha = NA, beta = NA, beta1 = 1 - v$pwr2, N = 200, r = 1, scale = "RR", direct = -1)

v <- getPwr_Bin_Noninf_JM1(p1_j = 0.6, p0_j = 0.5, p1_a = 0.5, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = log(1.8), alpha = 0.025, beta = NA, N = 200, r = 1, scale = "OR", direct = -1, sim = FALSE)
getN_Bin_Noninf_JM1(p1_a = 0.5, p0_a = 0.5, p1_j = 0.6, p0_j = 0.5, pi = 0.5, cut = log(1.8), alpha = NA, beta = NA, beta1 = 1 - v$pwr2, N = 200, r = 1, scale = "OR", direct = -1)

# -------------------------------------------------------------------------

v <- getPwr_Bin_Equi_JM1(p1_j = 0.4, p0_j = 0.5, p1_a = 0.5, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 0.3, alpha = 0.025, beta = NA, N = 200, r = 1, scale = "RD", sim = FALSE)
getN_Bin_Equi_JM1(p1_a = 0.5, p0_a = 0.5, p1_j = 0.4, p0_j = 0.5, pi = 0.5, cut = 0.3, alpha = NA, beta = NA, beta1 = 1 - v$pwr2, N = 200, r = 1, scale = "RD")

v <- getPwr_Bin_Equi_JM1(p1_j = 0.4, p0_j = 0.5, p1_a = 0.5, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = -log(0.6), alpha = 0.025, beta = NA, N = 200, r = 1, scale = "RR", sim = FALSE)
getN_Bin_Equi_JM1(p1_a = 0.5, p0_a = 0.5, p1_j = 0.4, p0_j = 0.5, pi = 0.5, cut = -log(0.6), alpha = NA, beta = NA, beta1 = 1 - v$pwr2, N = 200, r = 1, scale = "RR")

v <- getPwr_Bin_Equi_JM1(p1_j = 0.4, p0_j = 0.5, p1_a = 0.5, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = -log(0.4), alpha = 0.025, beta = NA, N = 200, r = 1, scale = "OR", sim = FALSE)
getN_Bin_Equi_JM1(p1_a = 0.5, p0_a = 0.5, p1_j = 0.4, p0_j = 0.5, pi = 0.5, cut = -log(0.4), alpha = NA, beta = NA, beta1 = 1 - v$pwr2, N = 200, r = 1, scale = "OR")


v <- getPwr_Bin_Equi_JM1(p1_j = 0.6, p0_j = 0.5, p1_a = 0.5, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 0.3, alpha = 0.025, beta = NA, N = 200, r = 1, scale = "RD", sim = FALSE)
getN_Bin_Equi_JM1(p1_a = 0.5, p0_a = 0.5, p1_j = 0.6, p0_j = 0.5, pi = 0.5, cut = 0.3, alpha = NA, beta = NA, beta1 = 1 - v$pwr2, N = 200, r = 1, scale = "RD")

v <- getPwr_Bin_Equi_JM1(p1_j = 0.6, p0_j = 0.5, p1_a = 0.5, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = log(1.5), alpha = 0.025, beta = NA, N = 200, r = 1, scale = "RR", sim = FALSE)
getN_Bin_Equi_JM1(p1_a = 0.5, p0_a = 0.5, p1_j = 0.6, p0_j = 0.5, pi = 0.5, cut = log(1.5), alpha = NA, beta = NA, beta1 = 1 - v$pwr2, N = 200, r = 1, scale = "RR")

v <- getPwr_Bin_Equi_JM1(p1_j = 0.6, p0_j = 0.5, p1_a = 0.5, p0_a = 0.5, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = log(2.3), alpha = 0.025, beta = NA, N = 200, r = 1, scale = "OR", sim = FALSE)
getN_Bin_Equi_JM1(p1_a = 0.5, p0_a = 0.5, p1_j = 0.6, p0_j = 0.5, pi = 0.5, cut = log(2.3), alpha = NA, beta = NA, beta1 = 1 - v$pwr2, N = 200, r = 1, scale = "OR")

# -------------------------------------------------------------------------
# getNe_Surv_Super_JM1、getNe_Surv_Noninf_JM1、getN_Surv_Equi_JM1
# -------------------------------------------------------------------------

v <- getPwr_Surv_Super_JM1(delta_j = log(1.3), delta_a = log(1.4), f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, Ne = 200, r = 1, sim = FALSE, criterion = 1)
getNe_Surv_Super_JM1(delta_a = log(1.4), delta_j = log(1.3), pi = 0.5, beta1 = 1 - v$pwr2, Ne = 200, r = 1, criterion = 1)

v <- getPwr_Surv_Super_JM1(delta_j = log(1.3), delta_a = log(1.4), f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, Ne = 200, r = 1, sim = FALSE, criterion = 2)
getNe_Surv_Super_JM1(delta_a = log(1.4), delta_j = log(1.3), pi = 0.5, beta1 = 1 - v$pwr2, Ne = 200, r = 1, criterion = 2)


v <- getPwr_Surv_Super_JM1(delta_j = log(0.8), delta_a = log(0.7), f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, Ne = 200, r = 1, sim = FALSE, criterion = 1)
getNe_Surv_Super_JM1(delta_a = log(0.7), delta_j = log(0.8), pi = 0.5, beta1 = 1 - v$pwr2, Ne = 200, r = 1, criterion = 1)

v <- getPwr_Surv_Super_JM1(delta_j = log(0.8), delta_a = log(0.7), f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, Ne = 200, r = 1, sim = FALSE, criterion = 2)
getNe_Surv_Super_JM1(delta_a = log(0.7), delta_j = log(0.8), pi = 0.5, beta1 = 1 - v$pwr2, Ne = 200, r = 1, criterion = 2)

# -------------------------------------------------------------------------

v <- getPwr_Surv_Noninf_JM1(delta_j = log(0.9), delta_a = log(1.0), f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = -log(0.7), alpha = 0.025, beta = NA, Ne = 200, r = 1, criterion = 1, direct = 1, sim = FALSE)
getNe_Surv_Noninf_JM1(delta_a = log(1.0), delta_j = log(0.9), pi = 0.5, cut = -log(0.7), beta1 = 1 - v$pwr2, Ne = 200, r = 1, criterion = 1, direct = 1)

v <- getPwr_Surv_Noninf_JM1(delta_j = log(0.9), delta_a = log(1.0), f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = -log(0.7), alpha = 0.025, beta = NA, Ne = 200, r = 1, criterion = 2, direct = 1, sim = FALSE)
getNe_Surv_Noninf_JM1(delta_a = log(1.0), delta_j = log(0.9), pi = 0.5, cut = -log(0.7), beta1 = 1 - v$pwr2, Ne = 200, r = 1, criterion = 2, direct = 1)

v <- getPwr_Surv_Noninf_JM1(delta_j = log(1.1), delta_a = log(1.0), f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = log(1.3), alpha = 0.025, beta = NA, Ne = 200, r = 1, criterion = 1, direct = -1, sim = FALSE)
getNe_Surv_Noninf_JM1(delta_a = log(1.0), delta_j = log(1.1), pi = 0.5, cut = log(1.3), beta1 = 1 - v$pwr2, Ne = 200, r = 1, criterion = 1, direct = -1)

v <- getPwr_Surv_Noninf_JM1(delta_j = log(1.1), delta_a = log(1.0), f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = log(1.3), alpha = 0.025, beta = NA, Ne = 200, r = 1, criterion = 2, direct = -1, sim = FALSE)
getNe_Surv_Noninf_JM1(delta_a = log(1.0), delta_j = log(1.1), pi = 0.5, cut = log(1.3), beta1 = 1 - v$pwr2, Ne = 200, r = 1, criterion = 2, direct = -1)

# -------------------------------------------------------------------------

v <- getPwr_Surv_Equi_JM1(delta_j = log(0.9), delta_a = log(1.0), f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = -log(0.7), alpha = 0.025, beta = NA, Ne = 200, r = 1, criterion = 1, sim = FALSE)
getNe_Surv_Equi_JM1(delta_a = log(1.0), delta_j = log(0.9), pi = 0.5, cut = -log(0.7), beta1 = 1 - v$pwr2, Ne = 200, r = 1, criterion = 1)

v <- getPwr_Surv_Equi_JM1(delta_j = log(0.9), delta_a = log(1.0), f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = -log(0.7), alpha = 0.025, beta = NA, Ne = 200, r = 1, criterion = 2, sim = FALSE)
getNe_Surv_Equi_JM1(delta_a = log(1.0), delta_j = log(0.9), pi = 0.5, cut = -log(0.7), beta1 = 1 - v$pwr2, Ne = 200, r = 1, criterion = 2)

v <- getPwr_Surv_Equi_JM1(delta_j = log(1.1), delta_a = log(1.0), f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = log(1.3), alpha = 0.025, beta = NA, Ne = 200, r = 1, criterion = 1, sim = FALSE)
getNe_Surv_Equi_JM1(delta_a = log(1.0), delta_j = log(1.1), pi = 0.5, cut = log(1.3), beta1 = 1 - v$pwr2, Ne = 200, r = 1, criterion = 1)

v <- getPwr_Surv_Equi_JM1(delta_j = log(1.1), delta_a = log(1.0), f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = log(1.3), alpha = 0.025, beta = NA, Ne = 200, r = 1, criterion = 2, sim = FALSE)
getNe_Surv_Equi_JM1(delta_a = log(1.0), delta_j = log(1.1), pi = 0.5, cut = log(1.3), beta1 = 1 - v$pwr2, Ne = 200, r = 1, criterion = 2)

# -------------------------------------------------------------------------
# getPwr_Count_Super_JM1、getPwr_Count_Noninf_JM1、getPwr_Count_Equi_JM1
# -------------------------------------------------------------------------

v <- getPwr_Count_Super_JM1(delta_j = log(1.2), delta_a = log(1.3), lambda0_j = 0.1, lambda0_a = 0.1, t_j = 5, t_a = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, N = 300, r = 1, sim = FALSE)
getN_Count_Super_JM1(delta_a = log(1.3), delta_j = log(1.2), lambda0_a = 0.1, lambda0_j = 0.1, t_a = 5, t_j = 5, k = 0, pi = 0.5, beta1 = 1 - v$pwr2, N = 300, r = 1)

v <- getPwr_Count_Super_JM1(delta_j = log(0.8), delta_a = log(0.7), lambda0_j = 0.1, lambda0_a = 0.1, t_j = 5, t_a = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, alpha = 0.025, beta = NA, N = 300, r = 1, sim = FALSE)
getN_Count_Super_JM1(delta_a = log(0.7), delta_j = log(0.8), lambda0_a = 0.1, lambda0_j = 0.1, t_j = 5, t_a = 5, k = 0, pi = 0.5, beta1 = 1 - v$pwr2, N = 300, r = 1)

# -------------------------------------------------------------------------

v <- getPwr_Count_Noninf_JM1(delta_j = log(0.9), delta_a = log(1.0), lambda0_j = 0.1, lambda0_a = 0.1, t_j = 5, t_a = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = -log(0.7), alpha = 0.025, beta = NA, N = 300, r = 1, direct = 1, sim = FALSE)
getN_Count_Noninf_JM1(delta_a = log(1.0), delta_j = log(0.9), lambda0_a = 0.1, lambda0_j = 0.1, t_j = 5, t_a = 5, k = 0, pi = 0.5, cut = -log(0.7), beta1 = 1 - v$pwr2, N = 300, r = 1, direct = 1)

v <- getPwr_Count_Noninf_JM1(delta_j = log(1.1), delta_a = log(1.0), lambda0_j = 0.1, lambda0_a = 0.1, t_j = 5, t_a = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = log(1.3), alpha = 0.025, beta = NA, N = 300, r = 1, direct = -1, sim = FALSE)
getN_Count_Noninf_JM1(delta_a = log(1.0), delta_j = log(1.1), lambda0_a = 0.1, lambda0_j = 0.1, t_j = 5, t_a = 5, k = 0, pi = 0.5, cut = log(1.3), beta1 = 1 - v$pwr2, N = 300, r = 1, direct = -1)

# -------------------------------------------------------------------------

v <- getPwr_Count_Equi_JM1(delta_j = log(0.9), delta_a = log(1.0), lambda0_j = 0.1, lambda0_a = 0.1, t_j = 5, t_a = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = -log(0.7), alpha = 0.025, beta = NA, N = 300, r = 1, sim = FALSE)
getN_Count_Equi_JM1(delta_a = log(1.0), delta_j = log(0.9), lambda0_a = 0.1, lambda0_j = 0.1, t_j = 5, t_a = 5, k = 0, pi = 0.5, cut = -log(0.7), beta1 = 1 - v$pwr2, N = 300, r = 1)

v <- getPwr_Count_Equi_JM1(delta_j = log(1.1), delta_a = log(1.0), lambda0_j = 0.1, lambda0_a = 0.1, t_j = 5, t_a = 5, k = 0, f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = log(1.3), alpha = 0.025, beta = NA, N = 300, r = 1, sim = FALSE)
getN_Count_Equi_JM1(delta_a = log(1.0), delta_j = log(1.1), lambda0_a = 0.1, lambda0_j = 0.1, t_j = 5, t_a = 5, k = 0, pi = 0.5, cut = log(1.3), beta1 = 1 - v$pwr2, N = 300, r = 1)

# -------------------------------------------------------------------------
# getN_Con_Super、getN_Con_Noninf、getN_Con_Equi
# -------------------------------------------------------------------------

(v <- getN_Con_Super(delta = seq(0.5, 1.5, 0.5), sigma = 4, alpha = 0.025, beta = 0.2, N = NA, r = 1))
getN_Con_Super(delta = seq(0.5, 1.5, 0.5), sigma = 4, alpha = 0.025, beta = NA, N = 300, r = 1)

# -------------------------------------------------------------------------

(v <- getN_Con_Noninf(delta = seq(0, -1.5, -0.5), sigma = 4, cut = 2, alpha = 0.025, beta = 0.2, N = NA, r = 1))
getN_Con_Noninf(delta = seq(0, -1.5, -0.5), sigma = 4, cut = 2, alpha = 0.025, beta = NA, N = 200, r = 1)

# -------------------------------------------------------------------------

(v <- getN_Con_Equi(delta = seq(0, -1.5, -0.5), sigma = 4, cut = 2, alpha = 0.025, beta = 0.2, N = NA, r = 1))
getN_Con_Equi(delta = seq(0, -1.5, -0.5), sigma = 4, cut = 2, alpha = 0.025, beta = NA, N = 200, r = 1)

# -------------------------------------------------------------------------
# getN_Bin_Super、getN_Bin_Noninf、getN_Bin_Equi
# -------------------------------------------------------------------------

(v <- getN_Bin_Super(p0 = 0.4, p1 = seq(0.5, 0.6, 0.05), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RD"))
getN_Bin_Super(p0 = 0.4, p1 = seq(0.5, 0.6, 0.05), alpha = 0.025, beta = NA, N = 300, r = 1, scale = "RD")

(v <- getN_Bin_Super(p0 = 0.4, p1 = seq(0.5, 0.6, 0.05), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RR"))
getN_Bin_Super(p0 = 0.4, p1 = seq(0.5, 0.6, 0.05), alpha = 0.025, beta = NA, N = 300, r = 1, scale = "RR")

(v <- getN_Bin_Super(p0 = 0.4, p1 = seq(0.5, 0.6, 0.05), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "OR"))
getN_Bin_Super(p0 = 0.4, p1 = seq(0.5, 0.6, 0.05), alpha = 0.025, beta = NA, N = 300, r = 1, scale = "OR")

# -------------------------------------------------------------------------

(v <- getN_Bin_Noninf(p0 = 0.5, p1 = seq(0.4, 0.5, 0.05), cut = 0.2, alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RD", direct = 1))
getN_Bin_Noninf(p0 = 0.5, p1 = seq(0.4, 0.5, 0.05), cut = 0.2, alpha = 0.025, beta = NA, N = 200, r = 1, scale = "RD", direct = 1)

(v <- getN_Bin_Noninf(p0 = 0.5, p1 = seq(0.4, 0.5, 0.05), cut = -log(0.6), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RR", direct = 1))
getN_Bin_Noninf(p0 = 0.5, p1 = seq(0.4, 0.5, 0.05), cut = -log(0.6), alpha = 0.025, beta = NA, N = 200, r = 1, scale = "RR", direct = 1)

(v <- getN_Bin_Noninf(p0 = 0.5, p1 = seq(0.4, 0.5, 0.05), cut = -log(0.5), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "OR", direct = 1))
getN_Bin_Noninf(p0 = 0.5, p1 = seq(0.4, 0.5, 0.05), cut = -log(0.5), alpha = 0.025, beta = NA, N = 200, r = 1, scale = "OR", direct = 1)

# -------------------------------------------------------------------------

(v <- getN_Bin_Equi(p0 = 0.5, p1 = seq(0.4, 0.5, 0.05), cut = 0.2, alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RD"))
getN_Bin_Equi(p0 = 0.5, p1 = seq(0.4, 0.5, 0.05), cut = 0.2, alpha = 0.025, beta = NA, N = 200, r = 1, scale = "RD")

(v <- getN_Bin_Equi(p0 = 0.5, p1 = seq(0.4, 0.5, 0.05), cut = -log(0.6), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RR"))
getN_Bin_Equi(p0 = 0.5, p1 = seq(0.4, 0.5, 0.05), cut = -log(0.6), alpha = 0.025, beta = NA, N = 200, r = 1, scale = "RR")

(v <- getN_Bin_Equi(p0 = 0.5, p1 = seq(0.4, 0.5, 0.05), cut = -log(0.5), alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "OR"))
getN_Bin_Equi(p0 = 0.5, p1 = seq(0.4, 0.5, 0.05), cut = -log(0.5), alpha = 0.025, beta = NA, N = 200, r = 1, scale = "OR")

# -------------------------------------------------------------------------
# getNe_Surv_Super、getNe_Surv_Noninf、getNe_Surv_Equi
# -------------------------------------------------------------------------

(v <- getNe_Surv_Super(delta = log(seq(0.7, 0.8, 0.05)), alpha = 0.025, beta = 0.2, Ne = NA, r = 1))
getNe_Surv_Super(delta = log(seq(0.7, 0.8, 0.05)), alpha = 0.025, beta = NA, Ne = 300, r = 1)

(v <- getNe_Surv_Noninf(delta = log(seq(1, 1.1, 0.05)), cut = log(1.3), alpha = 0.025, beta = 0.2, Ne = NA, r = 1, direct = -1))
getNe_Surv_Noninf(delta = log(seq(1, 1.1, 0.05)), cut = log(1.3), alpha = 0.025, beta = NA, Ne = 200, r = 1, direct = -1)

(v <- getNe_Surv_Equi(delta = log(seq(1, 1.1, 0.05)), cut = log(1.5), alpha = 0.025, beta = 0.2, Ne = NA, r = 1))
getNe_Surv_Equi(delta = log(seq(1, 1.1, 0.05)), cut = log(1.5), alpha = 0.025, beta = NA, Ne = 400, r = 1)

# -------------------------------------------------------------------------
# getN_Count_Super、getN_Count_Noninf、getN_Count_Equi
# -------------------------------------------------------------------------

(v <- getN_Count_Super(delta = log(seq(0.7, 0.8, 0.05)), lambda0 = 0.1, t = 5, k = 0, alpha = 0.025, beta = 0.2, N = NA, r = 1))
getN_Count_Super(delta = log(seq(0.7, 0.8, 0.05)), lambda0 = 0.1, t = 5, k = 0, alpha = 0.025, beta = NA, N = 1000, r = 1)

(v <- getN_Count_Noninf(delta = log(seq(0.9, 1, 0.05)), cut = -log(0.7), lambda0 = 0.1, t = 5, k = 1, alpha = 0.025, beta = 0.2, N = NA, r = 1))
getN_Count_Noninf(delta = log(seq(0.9, 1, 0.05)), cut = -log(0.7), lambda0 = 0.1, t = 5, k = 1, alpha = 0.025, beta = NA, N = 1000, r = 1)

(v <- getN_Count_Equi(delta = log(seq(0.9, 1, 0.05)), cut = -log(0.7), lambda0 = 0.1, t = 5, k = 1, alpha = 0.025, beta = 0.2, N = NA, r = 1))
getN_Count_Equi(delta = log(seq(0.9, 1, 0.05)), cut = -log(0.7), lambda0 = 0.1, t = 5, k = 1, alpha = 0.025, beta = NA, N = 1000, r = 1)

# -------------------------------------------------------------------------
# Ne_to_N
# -------------------------------------------------------------------------

Ne_to_N(
  Ne = 100, r = 1, lambda0 = log(2) / 20, lambda1 = log(2) / 20 * 0.75,
  dropoutRate = 0.05, dropoutTime = 12,
  a = 18, f = 18, follow_up = "until_end"
)

Ne_to_N(
  Ne = 100, r = 1, lambda0 = log(2) / 20, lambda1 = log(2) / 20 * 0.75,
  dropoutRate = 0.05, dropoutTime = 12,
  l = 18, follow_up = "fixed_period"
)
