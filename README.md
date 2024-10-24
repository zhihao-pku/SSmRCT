---
output:
  pdf_document: default
  html_document: default
---

# SSmRCT

Sample size allocation for mRCT

```R
library(SSmRCT)
library(ggplot2)
library(tidyr)
rm(list = ls())
```

## 1 Superiority Design for Continuous Endpoint Using Japanese Method 1

### Success Criteria

_Using the High Priority Metric as an Example_

- Global Success Criteria

$$
Z_a = \frac{\hat \delta_a}{\sqrt{Var(\hat \delta_a)}} \gt \Phi^{-1}(1 - \alpha)
$$

- Target Region Consistency Criteria

$$
\hat \delta_j - \pi\hat \delta_a \gt 0
$$

### Case

- Mean difference between experimental and placebo groups: in the target region: `delta_j = 0.5`; in the non-target region: `delta_nj = 0.7`; Common standard deviation: `sigma = 1`.

- Global sample size: `N = 100`; Sample size ratio between experimental and placebo groups: `r = 1`.

- Proportion of sample size allocated to the target region: `seq(0.1, 0.9, 0.1)`.

- Retaining 50% of global efficacy: `pi = 0.5`.

- One sided significance level: `alpha = 0.025.`

- High priority metric (where a larger value is better): `direct = 1`.

- Using theoretical calculation or data simulation method: `sim.`

```R
dat1 <- getPwr_Con_Super_JM1(
  delta_j = 0.5, delta_nj = 0.7, sigma = 1,
  f = seq(0.1, 0.9, 0.1),
  pi = 0.5, alpha = 0.025, N = 100, r = 1, direct = 1, sim = FALSE
)
dat1$M <- "calc"
dat2 <- getPwr_Con_Super_JM1(
  delta_j = 0.5, delta_nj = 0.7, sigma = 1,
  f = seq(0.1, 0.9, 0.1),
  pi = 0.5, alpha = 0.025, N = 100, r = 1, direct = 1, sim = TRUE
)
dat2$M <- "sim"

dat <- bind_rows(dat1, dat2)
dat <- gather(data = dat, key = "pwr class", value = "pwr", p1, p2, p3, p4)

p1 <- ggplot(data = dat, aes(x = f, y = pwr, linetype = M)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(`pwr class`), nrow = 2) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1.0, 0.2)) +
  labs(x = "allocation ratio", linetype = "")
p1
```

![](https://zhihao-pku.oss-cn-beijing.aliyuncs.com/202410212239026.png)

## 2 Non-inferiority Design for Continuous Endpoint Using Japanese Method 1

### Success Criteria

_Using the High Priority Metric as an Example_

Global Success Criteria

$$
Z_a = \frac{\hat \delta_a + \Delta}{\sqrt{Var(\hat \delta_a)}} \gt \Phi^{-1}(1 - \alpha)
$$

Target Region Consistency Criteria

$$
\hat \delta_j - \hat \delta_a + \pi\Delta \gt 0
$$

### Case

- Mean difference between experimental and positive-control groups: in the target region: `delta_j = -0.2`; in the non-target region: `delta_nj = -0.1`; Common standard deviation: `sigma = 1`.
- Global sample size: `N = 400`; Sample size ratio between experimental and positive-control groups: `r = 1`.
- Proportion of sample size allocated to the target region: `seq(0.1, 0.9, 0.1)`.
- Non-inferiority margin: `cut = 0.4`
- Retaining 50% of global efficacy: `pi = 0.5`.
- One sided significance level: `alpha = 0.025`.
- High priority metric (where a larger value is better): `direct = 1`.
- Using theoretical calculation or data simulation method: `sim.`

```R
dat1 <- getPwr_Con_Noninf_JM1(
  delta_j = -0.2, delta_nj = -0.1, sigma = 1,
  f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 0.4, alpha = 0.025,
  N = 400, r = 1, direct = 1, sim = FALSE
)
dat1$M <- "calc"
dat2 <- getPwr_Con_Noninf_JM1(
  delta_j = -0.2, delta_nj = -0.1, sigma = 1,
  f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 0.4, alpha = 0.025,
  N = 400, r = 1, direct = 1, sim = TRUE
)
dat2$M <- "sim"

dat <- bind_rows(dat1, dat2)
dat <- gather(
  data = dat, key = "pwr class", value = "pwr",
  p1, p2, p3, p4
)

p2 <- ggplot(data = dat, aes(x = f, y = pwr, linetype = M)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(`pwr class`), nrow = 2) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1.0, 0.2)) +
  labs(x = "allocation ratio", color = "", linetype = "")
p2
```

![](https://zhihao-pku.oss-cn-beijing.aliyuncs.com/202410212239244.png)

## 3 Equivalence Design for Continuous Endpoint Using Japanese Method 1

### Success Criteria

Global Success Criteria

$$
Z_{a_u} = \frac{\hat \delta_a + \Delta}{\sqrt{Var(\hat \delta_a)}} \gt \Phi^{-1}(1 - \alpha)\text{ and }Z_{a_l} = \frac{\hat \delta_a - \Delta}{\sqrt{Var(\hat \delta_a)}} \lt \Phi^{-1}(\alpha)
$$

Target Region Consistency Criteria

$$
\hat \delta_j - \hat \delta_a + \pi\Delta \gt 0\text{ and }\hat \delta_j - \hat \delta_a - \pi\Delta \lt 0
$$

### Case

- Mean difference between experimental and positive-control groups: in the target region: `delta_j = -0.2`; in the non-target region: `delta_nj = -0.1`; Common standard deviation: `sigma = 1`.
- Global sample size: `N = 400`; Sample size ratio between experimental and positive-control groups: `r = 1`.
- Proportion of sample size allocated to the target region: `seq(0.1, 0.9, 0.1)`.
- Non-inferiority margin: `cut = 0.4`
- Retaining 50% of global efficacy: `pi = 0.5`.
- One sided significance level: `alpha = 0.025`.
- Using theoretical calculation or data simulation method: `sim`.

```R
dat1 <- getPwr_Con_Equi_JM1(
  delta_j = -0.2, delta_nj = -0.1, sigma = 1,
  f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 0.4, alpha = 0.025,
  N = 400, r = 1, sim = FALSE
)
dat1$M <- "calc"
dat2 <- getPwr_Con_Equi_JM1(
  delta_j = -0.2, delta_nj = -0.1, sigma = 1,
  f = seq(0.1, 0.9, 0.1), pi = 0.5, cut = 0.4, alpha = 0.025,
  N = 400, r = 1, sim = TRUE
)
dat2$M <- "sim"

dat <- bind_rows(dat1, dat2)
dat <- gather(
  data = dat, key = "pwr class", value = "pwr",
  p1, p2, p3, p4
)

p3 <- ggplot(data = dat, aes(x = f, y = pwr, linetype = M)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(`pwr class`), nrow = 2) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1.0, 0.2)) +
  labs(x = "allocation ratio", color = "", linetype = "")
p3
```

![](https://zhihao-pku.oss-cn-beijing.aliyuncs.com/202410212240544.png)

## 4 Superiority Design for Continuous Endpoint Using Japanese Method 2

### Success Criteria

_Using the High Priority Metric as an Example_

Global Success Criteria

$$
Z_a = \frac{\hat \delta_a}{\sqrt{Var(\hat \delta_a)}} \gt \Phi^{-1}(1 - \alpha)
$$

Regional Consistency Criteria

$$
\hat \delta_i > 0 \text{ for i = 1, 2, .., m}
$$

### Case

- Mean difference between experimental and positive-control groups for each region: `delta_i = c(1, 0.8)`;

  Common standard deviation: `sigma = 4`.

- Global sample size: `N = 200`; Sample size ratio between experimental and positive-control groups: `r = 1`.

- Proportion of sample size allocated to each region: `c(f, 1 - f)`.

- One sided significance level: `alpha = 0.025`.

- High priority metric (where a larger value is better): `direct = 1`.

- Using theoretical calculation or data simulation method: `sim`.

```R
f_set <- seq(0.1, 0.9, 0.1)
dat1 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Con_Super_JM2(
    delta_i = c(1, 0.8), sigma = 4,
    fi = c(f, 1 - f),
    alpha = 0.025, N = 200, r = 1, direct = 1, sim = FALSE
  )$overall
  res$M <- "calc"
  res$f <- f
  res
})
dat2 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Con_Super_JM2(
    delta_i = c(1, 0.8), sigma = 4,
    fi = c(f, 1 - f),
    alpha = 0.025, N = 200, r = 1, direct = 1, sim = TRUE
  )$overall
  res$M <- "sim"
  res$f <- f
  res
})

dat <- bind_rows(dat1, dat2)
dat <- gather(
  data = dat, key = "pwr class",
  value = "pwr", p1, p2, p3, p4
)

p4 <- ggplot(data = dat, aes(x = f, y = pwr, linetype = M)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(`pwr class`), nrow = 2) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1.0, 0.2)) +
  labs(x = "allocation ratio", linetype = "")
p4
```

![](https://zhihao-pku.oss-cn-beijing.aliyuncs.com/202410212240611.png)

## 5 Non-inferiority Design for Continuous Endpoint Using Japanese Method 2

### Success Criteria

_Using the High Priority Metric as an Example_

Global Success Criteria

$$
Z_a = \frac{\hat \delta_a + \Delta}{\sqrt{Var(\hat \delta_a)}} \gt \Phi^{-1}(1 - \alpha)
$$

Target Region Consistency Criteria

$$
\hat \delta_i + \Delta\gt 0 \text{ for i = 1, 2, .., m}
$$

### Case

- Mean difference between experimental and positive-control groups for each region: `delta_i = c(-0.5, 0)`;

  Common standard deviation: `sigma = 4`.

- Global sample size: `N = 200`; Sample size ratio between experimental and positive-control groups: `r = 1`.

- Proportion of sample size allocated to each region: `c(f, 1 - f)`.

- One sided significance level: `alpha = 0.025`.

- High priority metric (where a larger value is better): `direct = 1`.

- Using theoretical calculation or data simulation method: `sim`.

```R
f_set <- seq(0.1, 0.9, 0.1)
dat1 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Con_Noninf_JM2(
    delta_i = c(-0.5, 0), sigma = 4,
    fi = c(f, 1 - f), cut = 2,
    alpha = 0.025, N = 200, r = 1, direct = 1, sim = FALSE
  )$overall
  res$M <- "calc"
  res$f <- f
  res
})
dat2 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Con_Noninf_JM2(
    delta_i = c(-0.5, 0), sigma = 4,
    fi = c(f, 1 - f), cut = 2,
    alpha = 0.025, N = 200, r = 1, direct = 1, sim = TRUE
  )$overall
  res$M <- "sim"
  res$f <- f
  res
})

dat <- bind_rows(dat1, dat2)
dat <- gather(
  data = dat, key = "pwr class",
  value = "pwr", p1, p2, p3, p4
)

p5 <- ggplot(data = dat, aes(x = f, y = pwr, linetype = M)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(`pwr class`), nrow = 2) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1.0, 0.2)) +
  labs(x = "allocation ratio", linetype = "")
p5
```

![](https://zhihao-pku.oss-cn-beijing.aliyuncs.com/202410212240754.png)

## 6 Equivalence Design for Continuous Endpoint Using Japanese Method 2

### Success Criteria

_Using the High Priority Metric as an Example_

Global Success Criteria

$$
Z_{a_u} = \frac{\hat \delta_a + \Delta}{\sqrt{Var(\hat \delta_a)}} \gt \Phi^{-1}(1 - \alpha)\text{ and }Z_{a_l} = \frac{\hat \delta_a - \Delta}{\sqrt{Var(\hat \delta_a)}} \lt \Phi^{-1}(\alpha)
$$

Regional Consistency Criteria

$$
\hat \delta_i + \Delta\gt 0 \text{ and } \hat \delta_i - \Delta \lt 0 \text{ for i = 1, 2, .., m}
$$

### Case

- Mean difference between experimental and positive-control groups for each region: `delta_i = c(-0.5, 0)`;

  Common standard deviation: `sigma = 4`.

- Global sample size: `N = 200`; Sample size ratio between experimental and positive-control groups: `r = 1`.

- Proportion of sample size allocated to each region: `c(f, 1 - f)`.

- One sided significance level: `alpha = 0.025`.

- Using theoretical calculation or data simulation method: `sim`.

```R
f_set <- seq(0.1, 0.9, 0.1)
dat1 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Con_Equi_JM2(
    delta_i = c(-0.5, 0), sigma = 4,
    fi = c(f, 1 - f), cut = 2,
    alpha = 0.025, N = 200, r = 1, sim = FALSE
  )$overall
  res$M <- "calc"
  res$f <- f
  res
})
dat2 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Con_Equi_JM2(
    delta_i = c(-0.5, 0), sigma = 4,
    fi = c(f, 1 - f), cut = 2,
    alpha = 0.025, N = 200, r = 1, sim = TRUE
  )$overall
  res$M <- "sim"
  res$f <- f
  res
})

dat <- bind_rows(dat1, dat2)
dat <- gather(
  data = dat, key = "pwr class",
  value = "pwr", p1, p2, p3, p4
)

p6 <- ggplot(data = dat, aes(x = f, y = pwr, linetype = M)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(`pwr class`), nrow = 2) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1.0, 0.2)) +
  labs(x = "allocation ratio", linetype = "")
p6
```

![](https://zhihao-pku.oss-cn-beijing.aliyuncs.com/202410212240394.png)

## 7 Superiority Design for Binary Endpoint Using Japanese Method 1

### Success Criteria

_Using the High Priority Metric as an Example_

Global Success Criteria

$$
\begin{gather*}
\hat \delta_a = \hat p_{1a} - \hat p_{0a} \\
Z_a = \frac{\hat \delta_a}{\sqrt{Var(\hat \delta_a)}}> \Phi^{-1}(1 - \alpha)
\end{gather*}
$$

Target Region Consistency Criteria

$$
\begin{gather*}
\hat \delta_j = \hat p_{1j} - \hat p_{0j} \\
\hat \delta_j - \pi\hat \delta_a \gt 0
\end{gather*}
$$

### Case

- Response rates of experimental and placebo groups: in the target region: `p1_j = 0.65, p0_j = 0.5`; in the non-target region: `p1_nj = 0.75, p0_nj = 0.5`.

- Global sample size: `N = 200`; Sample size ratio between experimental and placebo groups: `r = 1`.

- Proportion of sample size allocated to the target region: `seq(0.1, 0.9, 0.1)`.

- Retaining 50% of global efficacy: `pi = 0.5`.

- One sided significance level: `alpha = 0.025.`

- High priority metric (where a larger value is better): `direct = 1`.

- Using theoretical calculation or data simulation method: `sim.`

```R
dat1 <- getPwr_Bin_Super_JM1(
  p1_j = 0.65, p0_j = 0.5, p1_nj = 0.75, p0_nj = 0.5,
  f = seq(0.1, 0.9, 0.1),
  pi = 0.5, alpha = 0.025, N = 200, r = 1, direct = 1, sim = FALSE
)
dat1$M <- "calc"
dat2 <- getPwr_Bin_Super_JM1(
  p1_j = 0.65, p0_j = 0.5, p1_nj = 0.75, p0_nj = 0.5,
  f = seq(0.1, 0.9, 0.1),
  pi = 0.5, alpha = 0.025, N = 200, r = 1, direct = 1, sim = TRUE
)
dat2$M <- "sim"

dat <- bind_rows(dat1, dat2)
dat <- gather(data = dat, key = "pwr class", value = "pwr", p1, p2, p3, p4)

p7 <- ggplot(data = dat, aes(x = f, y = pwr, linetype = M)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(`pwr class`), nrow = 2) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1.0, 0.2)) +
  labs(x = "allocation ratio", linetype = "")
p7
```

![](https://zhihao-pku.oss-cn-beijing.aliyuncs.com/202410212241007.png)

## 8 Non-inferiority Design for Binary Endpoint Using Japanese Method 1

### Success Criteria

_Using the High Priority Metric as an Example_

Global Success Criteria

$$
\begin{gather*}
\hat \delta_a = \hat p_{1a} - \hat p_{0a} \\
Z_a = \frac{\hat \delta_a + \Delta}{\sqrt{Var(\hat \delta_a)}} \gt \Phi^{-1}(1 - \alpha)
\end{gather*}
$$

Target Region Consistency Criteria

$$
\begin{gather*}
\hat \delta_j = \hat p_{1j} - \hat p_{0j} \\
\hat \delta_j - \hat \delta_a + \pi\Delta \gt 0
\end{gather*}
$$

### Case

- Response rates of experimental and positive-control groups: in the target region: `p1_j = 0.55, p0_j = 0.65`; in the non-target region: `p1_nj = 0.65, p0_nj = 0.65`.

- Global sample size: `N = 400`; Sample size ratio between experimental and placebo groups: `r = 1`.

- Proportion of sample size allocated to the target region: `seq(0.1, 0.9, 0.1)`.

- Non-inferiority margin: `cut = 0.2`

- Retaining 50% of global efficacy: `pi = 0.5`.

- One sided significance level: `alpha = 0.025.`

- High priority metric (where a larger value is better): `direct = 1`.

- Using theoretical calculation or data simulation method: `sim.`

```R
dat1 <- getPwr_Bin_Noninf_JM1(
  p1_j = 0.55, p0_j = 0.65, p1_nj = 0.65, p0_nj = 0.65,
  f = seq(0.1, 0.9, 0.1), pi = 0.5,
  cut = 0.2, alpha = 0.025, N = 400, r = 1,
  direct = 1, sim = FALSE
)
dat1$M <- "calc"
dat2 <- getPwr_Bin_Noninf_JM1(
  p1_j = 0.55, p0_j = 0.65, p1_nj = 0.65, p0_nj = 0.65,
  f = seq(0.1, 0.9, 0.1), pi = 0.5,
  cut = 0.2, alpha = 0.025, N = 400, r = 1,
  direct = 1, sim = TRUE
)
dat2$M <- "sim"

dat <- bind_rows(dat1, dat2)
dat <- gather(
  data = dat, key = "pwr class",
  value = "pwr", p1, p2, p3, p4
)

p8 <- ggplot(data = dat, aes(x = f, y = pwr, linetype = M)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(`pwr class`), nrow = 2) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1.0, 0.2)) +
  labs(x = "allocation ratio", color = "", linetype = "")
p8
```

![](https://zhihao-pku.oss-cn-beijing.aliyuncs.com/202410212241758.png)

## 9 Equivalence Design for Binary Endpoint Using Japanese Method 1

### Success Criteria

_Using the High Priority Metric as an Example_

Global Success Criteria

$$
\begin{gather*}
\hat \delta_a = \hat p_{1a} - \hat p_{0a} \\
Z_{a_u} = \frac{\hat \delta_a + \Delta}{\sqrt{Var(\hat \delta_a)}} \gt \Phi^{-1}(1 - \alpha)\text{ and }Z_{a_l} = \frac{\hat \delta_a - \Delta}{\sqrt{Var(\hat \delta_a)}} \lt \Phi^{-1}(\alpha)
\end{gather*}
$$

Target Region Consistency Criteria

$$
\begin{gather*}
\hat \delta_j = \hat p_{1j} - \hat p_{0j} \\
\hat \delta_j - \hat \delta_a + \pi\Delta \gt 0\text{ and }\hat \delta_j - \hat \delta_a - \pi\Delta \lt 0
\end{gather*}
$$

### Case

- Response rates of experimental and positive-control groups: in the target region: `p1_j = 0.55, p0_j = 0.65`; in the non-target region: `p1_nj = 0.65, p0_nj = 0.65`.

- Global sample size: `N = 400`; Sample size ratio between experimental and placebo groups: `r = 1`.

- Proportion of sample size allocated to the target region: `seq(0.1, 0.9, 0.1)`.

- Non-inferiority margin: `cut = 0.2`

- Retaining 50% of global efficacy: `pi = 0.5`.

- One sided significance level: `alpha = 0.025.`

- High priority metric (where a larger value is better): `direct = 1`.

- Using theoretical calculation or data simulation method: `sim.`

```R
dat1 <- getPwr_Bin_Equi_JM1(
  p1_j = 0.55, p0_j = 0.65, p1_nj = 0.65, p0_nj = 0.65,
  f = seq(0.1, 0.9, 0.1), pi = 0.5,
  cut = 0.2, alpha = 0.025, N = 400, r = 1, sim = FALSE
)
dat1$M <- "calc"
dat2 <- getPwr_Bin_Equi_JM1(
  p1_j = 0.55, p0_j = 0.65, p1_nj = 0.65, p0_nj = 0.65,
  f = seq(0.1, 0.9, 0.1), pi = 0.5,
  cut = 0.2, alpha = 0.025, N = 400, r = 1, sim = TRUE
)
dat2$M <- "sim"

dat <- bind_rows(dat1, dat2)
dat <- gather(
  data = dat, key = "pwr class",
  value = "pwr", p1, p2, p3, p4
)

p9 <- ggplot(data = dat, aes(x = f, y = pwr, linetype = M)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(`pwr class`), nrow = 2) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1.0, 0.2)) +
  labs(x = "allocation ratio", color = "", linetype = "")
p9
```

![](https://zhihao-pku.oss-cn-beijing.aliyuncs.com/202410212241052.png)

## 10 Superiority Design for Binary Endpoint Using Japanese Method 2

### Success Criteria

_Using the High Priority Metric as an Example_

Global Success Criteria

$$
\begin{gather*}
\hat \delta_a =  \hat p_{1a} - \hat p_{0a} \\
Z_a = \frac{\hat \delta_a}{\sqrt{Var(\hat \delta_a)}} \gt \Phi^{-1}(1 - \alpha)
\end{gather*}
$$

Regional Consistency Criteria

$$
\hat \delta_i  = \hat p_{1i} - \hat p_{0i} \gt 0 \text{ for i = 1, 2, .., m}
$$

### Case

- Mean difference between experimental and positive-control groups: in the target region: `delta_j = -0.2`; in the non-target region: `delta_nj = -0.1`; Common standard deviation: `sigma = 1`.
- Global sample size: `N = 400`; Sample size ratio between experimental and positive-control groups: `r = 1`.
- Proportion of sample size allocated to the target region: `seq(0.1, 0.9, 0.1)`.
- Non-inferiority margin: `cut = 0.4`
- Retaining 50% of global efficacy: `pi = 0.5`.
- One sided significance level: `alpha = 0.025`.
- Using theoretical calculation or data simulation method: `sim`.

```R
f_set <- seq(0.1, 0.9, 0.1)
dat1 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Super_JM2(
    pt_i = c(0.7, 0.6),
    pc_i = c(0.4, 0.4),
    fi = c(f, 1 - f),
    alpha = 0.025, N = 100, r = 1, direct = 1, sim = FALSE
  )$overall
  res$M <- "calc"
  res$f <- f
  res
})
dat2 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Super_JM2(
    pt_i = c(0.7, 0.6),
    pc_i = c(0.4, 0.4),
    fi = c(f, 1 - f),
    alpha = 0.025, N = 100, r = 1, direct = 1, sim = TRUE
  )$overall
  res$M <- "sim"
  res$f <- f
  res
})

dat <- bind_rows(dat1, dat2)
dat <- gather(
  data = dat, key = "pwr class",
  value = "pwr", p1, p2, p3, p4
)

p10 <- ggplot(data = dat, aes(x = f, y = pwr, linetype = M)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(`pwr class`), nrow = 2) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1.0, 0.2)) +
  labs(x = "allocation ratio", linetype = "")
p10
```

![](https://zhihao-pku.oss-cn-beijing.aliyuncs.com/202410212241602.png)

## 11 Non-inferiority Design for Binary Endpoint Using Japanese Method 2

### Success Criteria

_Using the High Priority Metric as an Example_

Global Success Criteria

$$
\begin{gather*}
\hat \delta_a =  \hat p_{1a} - \hat p_{0a} \\
Z_a = \frac{\hat \delta_a + \Delta}{\sqrt{Var(\hat \delta_a)}} \gt \Phi^{-1}(1 - \alpha)
\end{gather*}
$$

Regional Consistency Criteria

$$
\begin{gather*}
\hat \delta_i  = \hat p_{1i} - \hat p_{0i}\\
\hat \delta_i+ \Delta \gt 0 \text{ for i = 1, 2, .., m}
\end{gather*}
$$

### Case

- Mean difference between experimental and positive-control groups: in the target region: `delta_j = -0.2`; in the non-target region: `delta_nj = -0.1`; Common standard deviation: `sigma = 1`.
- Global sample size: `N = 400`; Sample size ratio between experimental and positive-control groups: `r = 1`.
- Proportion of sample size allocated to the target region: `seq(0.1, 0.9, 0.1)`.
- Non-inferiority margin: `cut = 0.4`
- Retaining 50% of global efficacy: `pi = 0.5`.
- One sided significance level: `alpha = 0.025`.
- Using theoretical calculation or data simulation method: `sim`.

```R
f_set <- seq(0.1, 0.9, 0.1)
dat1 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Noninf_JM2(
    pt_i = c(0.5, 0.6),
    pc_i = c(0.6, 0.6),
    fi = c(f, 1 - f), cut = 0.3,
    alpha = 0.025, N = 100, r = 1, direct = 1, sim = FALSE
  )$overall
  res$M <- "calc"
  res$f <- f
  res
})
dat2 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Noninf_JM2(
    pt_i = c(0.5, 0.6),
    pc_i = c(0.6, 0.6),
    fi = c(f, 1 - f), cut = 0.3,
    alpha = 0.025, N = 100, r = 1, direct = 1, sim = TRUE
  )$overall
  res$M <- "sim"
  res$f <- f
  res
})

dat <- bind_rows(dat1, dat2)
dat <- gather(
  data = dat, key = "pwr class",
  value = "pwr", p1, p2, p3, p4
)

p11 <- ggplot(data = dat, aes(x = f, y = pwr, linetype = M)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(`pwr class`), nrow = 2) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1.0, 0.2)) +
  labs(x = "allocation ratio", linetype = "")
p11
```

![](https://zhihao-pku.oss-cn-beijing.aliyuncs.com/202410212242443.png)

## 12 Equivalence Design for Binary Endpoint Using Japanese Method 2

### Success Criteria

_Using the High Priority Metric as an Example_

Global Success Criteria

$$
\begin{gather*}
\hat \delta_a =  \hat p_{1a} - \hat p_{0a} \\
Z_{a_u} = \frac{\hat \delta_a + \Delta}{\sqrt{Var(\hat \delta_a)}} \gt \Phi^{-1}(1 - \alpha)\text{ and }Z_{a_l} = \frac{\hat \delta_a - \Delta}{\sqrt{Var(\hat \delta_a)}} \lt \Phi^{-1}(\alpha)
\end{gather*}
$$

Regional Consistency Criteria

$$
\begin{gather*}
\hat \delta_i  = \hat p_{1i} - \hat p_{0i}\\
\hat \delta_i + \Delta\gt 0 \text{ and } \hat \delta_i - \Delta \lt 0 \text{ for i = 1, 2, .., m}
\end{gather*}
$$

### Case

- Mean difference between experimental and positive-control groups: in the target region: `delta_j = -0.2`; in the non-target region: `delta_nj = -0.1`; Common standard deviation: `sigma = 1`.
- Global sample size: `N = 400`; Sample size ratio between experimental and positive-control groups: `r = 1`.
- Proportion of sample size allocated to the target region: `seq(0.1, 0.9, 0.1)`.
- Non-inferiority margin: `cut = 0.4`
- Retaining 50% of global efficacy: `pi = 0.5`.
- One sided significance level: `alpha = 0.025`.
- Using theoretical calculation or data simulation method: `sim`.

```R
f_set <- seq(0.1, 0.9, 0.1)
dat1 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Equi_JM2(
    pt_i = c(0.5, 0.6),
    pc_i = c(0.6, 0.6),
    fi = c(f, 1 - f), cut = 0.3,
    alpha = 0.025, N = 100, r = 1, sim = FALSE
  )$overall
  res$M <- "calc"
  res$f <- f
  res
})
dat2 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Bin_Equi_JM2(
    pt_i = c(0.5, 0.6),
    pc_i = c(0.6, 0.6),
    fi = c(f, 1 - f), cut = 0.3,
    alpha = 0.025, N = 100, r = 1, sim = TRUE
  )$overall
  res$M <- "sim"
  res$f <- f
  res
})

dat <- bind_rows(dat1, dat2)
dat <- gather(
  data = dat, key = "pwr class",
  value = "pwr", p1, p2, p3, p4
)

p12 <- ggplot(data = dat, aes(x = f, y = pwr, linetype = M)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(`pwr class`), nrow = 2) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1.0, 0.2)) +
  labs(x = "allocation ratio", linetype = "")
p12
```

![](https://zhihao-pku.oss-cn-beijing.aliyuncs.com/202410212242360.png)

## 13 Superiority Design for Time-to event Endpoint Using Japanese Method 1

### Success Criteria

_Using “the Smaller the HR, the Better” as an Example_

criterion = 1

Global Success Criteria

$$
\begin{gather*}
\hat \delta_a = log(\hat {HR_a}) \\
Z_a = \frac{\hat \delta_a}{\sqrt{Var(\hat \delta_a)}} \lt \Phi^{-1}(\alpha)
\end{gather*}
$$

Target Region Consistency Criteria

$$
\begin{gather*}
\hat \delta_j = log(\hat {HR_j}) \\
\hat \delta_j - \pi\hat \delta_a \lt 0
\end{gather*}
$$

criterion = 2

Global Success Criteria

$$
\begin{gather*}
\hat \delta_a = log(\hat {HR_a}) \\
Z_a = \frac{\hat \delta_a}{\sqrt{Var(\hat \delta_a)}} \gt \Phi^{-1}(\alpha)
\end{gather*}
$$

Target Region Consistency Criteria

$$
1 - e^{\hat \delta_j} - \pi(1 - e^{\hat \delta_a}) \gt 0
$$

### Case

- Mean difference between experimental and positive-control groups: in the target region: `delta_j = -0.2`; in the non-target region: `delta_nj = -0.1`; Common standard deviation: `sigma = 1`.
- Global sample size: `N = 400`; Sample size ratio between experimental and positive-control groups: `r = 1`.
- Proportion of sample size allocated to the target region: `seq(0.1, 0.9, 0.1)`.
- Non-inferiority margin: `cut = 0.4`
- Retaining 50% of global efficacy: `pi = 0.5`.
- One sided significance level: `alpha = 0.025`.
- Using theoretical calculation or data simulation method: `sim`.

```R
dat1 <- getPwr_Surv_Super_JM1(
  delta_j = log(0.8), delta_nj = log(0.6),
  f = seq(0.1, 0.9, 0.1),
  pi = 0.5, alpha = 0.025, N = 200, r = 1,
  criterion = 1, direct = -1, sim = FALSE
)
dat1$M <- "calc"
dat2 <- getPwr_Surv_Super_JM1(
  delta_j = log(0.8), delta_nj = log(0.6),
  f = seq(0.1, 0.9, 0.1),
  pi = 0.5, alpha = 0.025, N = 200, r = 1,
  criterion = 1, direct = -1, sim = TRUE
)
dat2$M <- "sim"
dat3 <- getPwr_Surv_Super_JM1(
  delta_j = log(0.8), delta_nj = log(0.6),
  f = seq(0.1, 0.9, 0.1),
  pi = 0.5, alpha = 0.025, N = 200, r = 1,
  criterion = 2, direct = -1, sim = FALSE
)
dat3$M <- "calc"
dat4 <- getPwr_Surv_Super_JM1(
  delta_j = log(0.8), delta_nj = log(0.6),
  f = seq(0.1, 0.9, 0.1),
  pi = 0.5, alpha = 0.025, N = 200, r = 1,
  criterion = 2, direct = -1, sim = TRUE
)
dat4$M <- "sim"

dat <- bind_rows(dat1, dat2, dat3, dat4)
dat <- gather(data = dat, key = "pwr class", value = "pwr", p1, p2, p3, p4)

p13 <- ggplot(data = dat, aes(
  x = f, y = pwr,
  color = factor(criterion), linetype = M
)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(`pwr class`), nrow = 2) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1.0, 0.2)) +
  labs(x = "allocation ratio", color = "", linetype = "")
p13
```

![](https://zhihao-pku.oss-cn-beijing.aliyuncs.com/202410212242836.png)

## 14 Non-inferiority Design for Time-to event Endpoint Using Japanese Method 1

### Success Criteria

_Using “the Smaller the HR, the Better” as an Example_

Global Success Criteria

$$
\begin{gather*}
\hat \delta_a = log(\hat {HR_a}) \\
Z_a = \frac{\hat \delta_a - \Delta}{\sqrt{Var(\hat \delta_a)}} \lt \Phi^{-1}(\alpha)
\end{gather*}
$$

Target Region Consistency Criteria

$$
\begin{gather*}
\hat \delta_j = log(\hat {HR_j}) \\
\hat \delta_j - \hat \delta_a - \pi\Delta \lt 0
\end{gather*}
$$

### Case

- Mean difference between experimental and positive-control groups: in the target region: `delta_j = -0.2`; in the non-target region: `delta_nj = -0.1`; Common standard deviation: `sigma = 1`.
- Global sample size: `N = 400`; Sample size ratio between experimental and positive-control groups: `r = 1`.
- Proportion of sample size allocated to the target region: `seq(0.1, 0.9, 0.1)`.
- Non-inferiority margin: `cut = 0.4`
- Retaining 50% of global efficacy: `pi = 0.5`.
- One sided significance level: `alpha = 0.025`.
- Using theoretical calculation or data simulation method: `sim`.

```R
dat1 <- getPwr_Surv_Noninf_JM1(
  delta_j = log(1.1), delta_nj = log(1.0),
  f = seq(0.1, 0.9, 0.1), cut = log(1.3),
  pi = 0.5, alpha = 0.025, N = 400, r = 1,
  direct = -1, sim = FALSE
)
dat1$M <- "calc"
dat2 <- getPwr_Surv_Noninf_JM1(
  delta_j = log(1.1), delta_nj = log(1.0),
  f = seq(0.1, 0.9, 0.1), cut = log(1.3),
  pi = 0.5, alpha = 0.025, N = 400, r = 1,
  direct = -1, sim = TRUE
)
dat2$M <- "sim"

dat <- bind_rows(dat1, dat2)
dat <- gather(
  data = dat, key = "pwr class", value = "pwr",
  p1, p2, p3, p4
)

p14 <- ggplot(data = dat, aes(x = f, y = pwr, linetype = M)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(`pwr class`), nrow = 2) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1.0, 0.2)) +
  labs(x = "allocation ratio", color = "", linetype = "")
p14
```

![](https://zhihao-pku.oss-cn-beijing.aliyuncs.com/202410212242224.png)

## 15 Equivalence Design for Time-to event Endpoint Using Japanese Method 1

### Success Criteria

_Using “the Smaller the HR, the Better” as an Example_

Global Success Criteria

$$
\begin{gather*}
\hat \delta_a = log(\hat {HR_a}) \\
Z_{a_u} = \frac{\hat \delta_a + \Delta}{\sqrt{Var(\hat \delta_a)}} \gt \Phi^{-1}(1 - \alpha)\text{ and }Z_{a_l} = \frac{\hat \delta_a - \Delta}{\sqrt{Var(\hat \delta_a)}} \lt \Phi^{-1}(\alpha)
\end{gather*}
$$

Target Region Consistency Criteria

$$
\begin{gather*}
\hat \delta_j = log(\hat {HR_j}) \\
\hat \delta_j - \hat \delta_a + \pi\Delta \gt 0\text{ and }\hat \delta_j - \hat \delta_a - \pi\Delta \lt 0
\end{gather*}
$$

### Case

- Mean difference between experimental and positive-control groups: in the target region: `delta_j = -0.2`; in the non-target region: `delta_nj = -0.1`; Common standard deviation: `sigma = 1`.
- Global sample size: `N = 400`; Sample size ratio between experimental and positive-control groups: `r = 1`.
- Proportion of sample size allocated to the target region: `seq(0.1, 0.9, 0.1)`.
- Non-inferiority margin: `cut = 0.4`
- Retaining 50% of global efficacy: `pi = 0.5`.
- One sided significance level: `alpha = 0.025`.
- Using theoretical calculation or data simulation method: `sim`.

```R
dat1 <- getPwr_Surv_Equi_JM1(
  delta_j = log(1.1), delta_nj = log(1.0),
  f = seq(0.1, 0.9, 0.1), cut = log(1.3),
  pi = 0.5, alpha = 0.025, N = 400, r = 1, sim = FALSE
)
dat1$M <- "calc"
dat2 <- getPwr_Surv_Equi_JM1(
  delta_j = log(1.1), delta_nj = log(1.0),
  f = seq(0.1, 0.9, 0.1), cut = log(1.3),
  pi = 0.5, alpha = 0.025, N = 400, r = 1, sim = TRUE
)
dat2$M <- "sim"

dat <- bind_rows(dat1, dat2)
dat <- gather(
  data = dat, key = "pwr class", value = "pwr",
  p1, p2, p3, p4
)

p15 <- ggplot(data = dat, aes(x = f, y = pwr, linetype = M)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(`pwr class`), nrow = 2) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1.0, 0.2)) +
  labs(x = "allocation ratio", color = "", linetype = "")
p15
```

![](https://zhihao-pku.oss-cn-beijing.aliyuncs.com/202410212243150.png)

## 16 Superiority Design for Time-to event Endpoint Using Japanese Method 2

### Success Criteria

_Using “the Smaller the HR, the Better” as an Example_

Global Success Criteria

$$
\begin{gather*}
\hat \delta_a = log(\hat HR_a) \\
Z_a = \frac{\hat \delta_a}{\sqrt{Var(\hat \delta_a)}} \lt \Phi^{-1}(\alpha)
\end{gather*}
$$

Regional Consistency Criteria

$$
\begin{gather*}
\hat \delta_i = log(\hat {HR_i}) \\
\hat \delta_i \lt 0 \text{ for i = 1, 2, .., m}
\end{gather*}
$$

### Case

- Mean difference between experimental and positive-control groups: in the target region: `delta_j = -0.2`; in the non-target region: `delta_nj = -0.1`; Common standard deviation: `sigma = 1`.
- Global sample size: `N = 400`; Sample size ratio between experimental and positive-control groups: `r = 1`.
- Proportion of sample size allocated to the target region: `seq(0.1, 0.9, 0.1)`.
- Non-inferiority margin: `cut = 0.4`
- Retaining 50% of global efficacy: `pi = 0.5`.
- One sided significance level: `alpha = 0.025`.
- Using theoretical calculation or data simulation method: `sim`.

```R
f_set <- seq(0.1, 0.9, 0.1)
dat1 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Surv_Super_JM2(
    delta_i = c(log(0.8), log(0.6)),
    fi = c(f, 1 - f),
    alpha = 0.025, N = 300, r = 1, direct = -1, sim = FALSE
  )$overall
  res$M <- "calc"
  res$f <- f
  res
})

dat2 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Surv_Super_JM2(
    delta_i = c(log(0.8), log(0.6)),
    fi = c(f, 1 - f),
    alpha = 0.025, N = 300, r = 1, direct = -1, sim = TRUE
  )$overall
  res$M <- "sim"
  res$f <- f
  res
})

dat <- bind_rows(dat1, dat2)
dat <- gather(
  data = dat, key = "pwr class",
  value = "pwr", p1, p2, p3, p4
)

p16 <- ggplot(data = dat, aes(x = f, y = pwr, linetype = M)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(`pwr class`), nrow = 2) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1.0, 0.2)) +
  labs(x = "allocation ratio", linetype = "")
p16
```

![](https://zhihao-pku.oss-cn-beijing.aliyuncs.com/202410212243318.png)

## 17 Non-inferiority Design for Time-to event Endpoint Using Japanese Method 2

### Success Criteria

_Using “the Smaller the HR, the Better” as an Example_

Global Success Criteria

$$
\begin{gather*}
\hat \delta_a = log(\hat {HR_a}) \\
Z_a = \frac{\hat \delta_a - \Delta}{\sqrt{Var(\hat \delta_a)}} \lt \Phi^{-1}(\alpha)
\end{gather*}
$$

Target Region Consistency Criteria

$$
\begin{gather*}
\hat \delta_i = log(\hat {HR_i}) \\
\hat \delta_i - \Delta\lt 0 \text{ for i = 1, 2, .., m}
\end{gather*}
$$

### Case

- Mean difference between experimental and positive-control groups: in the target region: `delta_j = -0.2`; in the non-target region: `delta_nj = -0.1`; Common standard deviation: `sigma = 1`.
- Global sample size: `N = 400`; Sample size ratio between experimental and positive-control groups: `r = 1`.
- Proportion of sample size allocated to the target region: `seq(0.1, 0.9, 0.1)`.
- Non-inferiority margin: `cut = 0.4`
- Retaining 50% of global efficacy: `pi = 0.5`.
- One sided significance level: `alpha = 0.025`.
- Using theoretical calculation or data simulation method: `sim`.

```R
f_set <- seq(0.1, 0.9, 0.1)
dat1 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Surv_Noninf_JM2(
    delta_i = c(log(1.1), log(1.0)),
    fi = c(f, 1 - f), cut = log(1.3),
    alpha = 0.025, N = 300, r = 1, direct = -1, sim = FALSE
  )$overall
  res$M <- "calc"
  res$f <- f
  res
})

dat2 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Surv_Noninf_JM2(
    delta_i = c(log(1.1), log(1.0)),
    fi = c(f, 1 - f), cut = log(1.3),
    alpha = 0.025, N = 300, r = 1, direct = -1, sim = TRUE
  )$overall
  res$M <- "sim"
  res$f <- f
  res
})

dat <- bind_rows(dat1, dat2)
dat <- gather(
  data = dat, key = "pwr class",
  value = "pwr", p1, p2, p3, p4
)

p17 <- ggplot(data = dat, aes(x = f, y = pwr, linetype = M)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(`pwr class`), nrow = 2) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1.0, 0.2)) +
  labs(x = "allocation ratio", linetype = "")
p17
```

![](https://zhihao-pku.oss-cn-beijing.aliyuncs.com/202410212243424.png)

## 18 Equivalence Design for Time-to event Endpoint Using Japanese Method 2

### Success Criteria

_Using “the Smaller the HR, the Better” as an Example_

Global Success Criteria

$$
\begin{gather*}
\hat \delta_a = log(\hat {HR_a}) \\
Z_{a_u} = \frac{\hat \delta_a + \Delta}{\sqrt{Var(\hat \delta_a)}} \gt \Phi^{-1}(1 - \alpha)\text{ and }Z_{a_l} = \frac{\hat \delta_a - \Delta}{\sqrt{Var(\hat \delta_a)}} \lt \Phi^{-1}(\alpha)
\end{gather*}
$$

Regional Consistency Criteria

$$
\begin{gather*}
\hat \delta_i = log(\hat {HR_i})\\
\hat \delta_i + \Delta\gt 0 \text{ and } \hat \delta_i - \Delta \lt 0 \text{ for i = 1, 2, .., m}
\end{gather*}
$$

### Case

- Mean difference between experimental and positive-control groups: in the target region: `delta_j = -0.2`; in the non-target region: `delta_nj = -0.1`; Common standard deviation: `sigma = 1`.
- Global sample size: `N = 400`; Sample size ratio between experimental and positive-control groups: `r = 1`.
- Proportion of sample size allocated to the target region: `seq(0.1, 0.9, 0.1)`.
- Non-inferiority margin: `cut = 0.4`
- Retaining 50% of global efficacy: `pi = 0.5`.
- One sided significance level: `alpha = 0.025`.
- Using theoretical calculation or data simulation method: `sim`.

```R
f_set <- seq(0.1, 0.9, 0.1)
dat1 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Surv_Equi_JM2(
    delta_i = c(log(1.1), log(1.0)),
    fi = c(f, 1 - f), cut = log(1.3),
    alpha = 0.025, N = 300, r = 1, direct = -1, sim = FALSE
  )$overall
  res$M <- "calc"
  res$f <- f
  res
})

dat2 <- map_dfr(.x = 1:length(f_set), .f = function(i) {
  f <- f_set[i]
  res <- getPwr_Surv_Equi_JM2(
    delta_i = c(log(1.1), log(1.0)),
    fi = c(f, 1 - f), cut = log(1.3),
    alpha = 0.025, N = 300, r = 1, direct = -1, sim = TRUE
  )$overall
  res$M <- "sim"
  res$f <- f
  res
})

dat <- bind_rows(dat1, dat2)
dat <- gather(
  data = dat, key = "pwr class",
  value = "pwr", p1, p2, p3, p4
)

p18 <- ggplot(data = dat, aes(x = f, y = pwr, linetype = M)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(`pwr class`), nrow = 2) +
  scale_x_continuous(breaks = seq(0.1, 0.9, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1.0, 0.2)) +
  labs(x = "allocation ratio", linetype = "")
p18
```

![](https://zhihao-pku.oss-cn-beijing.aliyuncs.com/202410212243410.png)
