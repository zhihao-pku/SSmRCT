## Introduction

Package `SSmRCT` is designed to determine the regional sample size allocation strategy in mRCTs (Multi-Regional Clinical Trials). It is applicable to continuous, binary, survival, and count endpoints and is suitable for superiority, non-inferiority, and equivalence designs.

## Installation

Installation requires the `devtools` package to be installed first. Then, run the following command in R to install the `SSmRCT` package:

```R
devtools::install_github("zhihao-pku/SSmRCT")
```

In some countries (e.g., China), installation may fail due to network issues. You can try installing it again after a failure.

## Japan's Methods

### Method 1

Assuming that a larger $\delta$ is more favorable, the criterion of treatment effect consistency between target region and globally is as follows:

Superiority design:
$$
\hat \delta_j - \pi\hat \delta_a > 0
$$
Non-inferiority design:
$$
\hat \delta_j - \hat \delta_a + \pi\Delta > 0
$$
Equivalence design:
$$
\hat \delta_j - \hat \delta_a + \pi\Delta > 0\text{ and }\hat \delta_j - \hat \delta_a - \pi\Delta < 0
$$

- $\delta_j$ and $\delta_a$ are the treatment effects in the target region and globally, respectively.

- $\pi$ is the proportion of the global treatment effect retained in the target region.
- $\Delta$ is the non-inferiority or equivalence margin.

### Method 2

Assuming that a larger $\delta$ is more favorable, the criterion of treatment effect consistency between region $i$ and the globally is as follows:

Superiority design:
$$
\hat \delta_i > 0 \text{ for i = 1, 2, .., m}
$$
Non-inferiority design:
$$
\hat \delta_i + \Delta_i > 0 \text{ for i = 1, 2, .., m}
$$
Equivalence design:
$$
\hat \delta_i + \Delta_i > 0 \text{ and } \hat \delta_i - \Delta_i < 0 \text{ for i = 1, 2, .., m}
$$

- $\Delta_i$ is the non-inferiority or equivalence margin for region $i$.

In above formulas, $\delta$ can represent the mean difference for continuous endpoints, the rate difference, $log(RR)$, or $log(OR)$ for binary endpoints, $log(RR)$ for count endpoints, and $log(HR)$ or $HR - 1$ for survival endpoints.

## SSmRCT Usage Examples

This package mainly includes three series of functions:

Series 1: `getN_xx_xx`

Calculate power given a sample size, or calculate sample size given a power.

Series 2: `getN_xx_xx_JM1`

Calculate the proportion of the target region's sample size to the global sample size given the global sample size and the required consistency probability of the target region's treatment effect with the global treatment effect.

Series 3: `getPwr_xx_xx_JM1` and `getPwr_xx_xx_JM2` 

Calculate the power of the trial given the global sample size and the target region's sample size, including:

- Marginal power: Probability of global trial success
- Marginal power: Probability of consistency of the target region's treatment effect with the global treatment effect
- Joint power: Probability of global trial  success and consistency of the target region's treatment effect with the global treatment effect
- Conditional power: Probability of consistency of the target region's treatment effect with the global treatment effect given global trial success

## Series 1: `getN_xx_xx`

Calculate power given a sample size, or calculate sample size given a power.

### Example 1: Continuous Endpoint, Superiority Design

Mean differences between the treatment group and the control group are 1.2, 1.3, 1.4, 1.5, with a common standard deviation of 4. When `beta` is not `NA` and `N` is `NA`, it will return the sample size.

```R
getN_Con_Super(
  delta = seq(1.2, 1.5, 0.1), sigma = 4, alpha = 0.025, beta = 0.2,
  N = NA, r = 1
)

# delta sigma alpha beta   N  n1  n0 r       pwr
# 1   1.2     4 0.025  0.2 350 175 175 1 0.8013015
# 2   1.3     4 0.025  0.2 298 149 149 1 0.8010063
# 3   1.4     4 0.025  0.2 258 129 129 1 0.8026021
# 4   1.5     4 0.025  0.2 224 112 112 1 0.8013015
```
When `beta` is `NA` and `N` is not `NA`, it will return the power.
```R
getN_Con_Super(
  delta = 1.5, sigma = 4, alpha = 0.025, beta = NA,
  N = seq(140, 200, 20), r = 1
)

# delta sigma alpha beta   N  n1  n0 r       pwr
# 1   1.5     4 0.025   NA 140  70  70 1 0.6020149
# 2   1.5     4 0.025   NA 160  80  80 1 0.6597366
# 3   1.5     4 0.025   NA 180  90  90 1 0.7107621
# 4   1.5     4 0.025   NA 200 100 100 1 0.7554329
```

### Example 2: Continuous Endpoint, Non-Inferiority Design

Mean difference between the treatment group and the control group is 1, with a common standard deviation of 4. Non-inferiority margins are 2, 2.5, 3, and `direct = -1` indicates that a smaller mean difference is more favorable.

```R
getN_Con_Noninf(
  delta = 1, sigma = 4, cut = seq(2, 3, 0.5), alpha = 0.025, beta = 0.2,
  N = NA, r = 1, direct = -1
)

# delta sigma cut alpha beta direct   N  n1  n0 r       pwr
# 1     1     4 2.0 0.025  0.2     -1 504 252 252 1 0.8013015
# 2     1     4 2.5 0.025  0.2     -1 224 112 112 1 0.8013015
# 3     1     4 3.0 0.025  0.2     -1 126  63  63 1 0.8013015
```

### Example 3: Survival Endpoint, Non-Inferiority Design

The hazard ratio (HR) of the treatment group compared to the control group is 1.1, and the non-inferiority margin for HR is 1.3. `direct = -1` indicates that a smaller HR is more favorable. This function returns the number of events rather than the sample size.

```R
getNe_Surv_Noninf(
  delta = log(1.1),
  cut = log(1.3),
  alpha = 0.025, beta = 0.2, Ne = NA, r = 1, direct = -1
)

# delta       cut alpha beta direct   Ne ne1 ne0 r       pwr
# 1 0.09531018 0.2623643 0.025  0.2     -1 1126 563 563 1 0.8003475
```

## Number of Event to Sample Size: `Ne_to_N` 

Given the number of events, calculate the required sample size using the `Ne_to_N` function.

### Example 1:  Followed Until the End of Study

Required event number is 100, sample size ratio of treatment group to control group is 1:1 (`r = 1`), median survival time of control group is 20 months, hazard ratio (HR) of treatment group compared to control group is 0.75, annual dropout rate is 5%, enrollment duration is 18 months, follow-up durations are 16, 18, and 24 months, and participants are followed until the end of the study.

For example, when the follow-up duration is 18 months, the event rates for the treatment group, control group, and total are 0.47, 0.57, and 0.52, respectively. The required sample sizes for the treatment group, control group, and total are 95, 95, and 190, respectively.

```R
Ne_to_N(
  Ne = 100, r = 1, lambda0 = log(2) / 20, lambda1 = log(2) / 20 * 0.75,
  dropoutRate = 0.05, dropoutTime = 12,
  a = 18, f = c(16, 18, 24), follow_up = "until_end"
)

# Ne r    lambda0    lambda1 dropoutRate dropoutTime  a  f  l follow_up eventRate0
# 1 100 1 0.03465736 0.02599302        0.05          12 18 16 NA until_end  0.5469831
# 2 100 1 0.03465736 0.02599302        0.05          12 18 18 NA until_end  0.5727036
# 3 100 1 0.03465736 0.02599302        0.05          12 18 24 NA until_end  0.6388738
# eventRate1 eventRate        N0        N1        N
# 1  0.4508670 0.4989250 100.21545 100.21545 200.4309
# 2  0.4748355 0.5237696  95.46183  95.46183 190.9237
# 3  0.5386219 0.5887478  84.92600  84.92600 169.8520
```

### Example 2: Followed a fixed period

Each participant is followed for fixed 16, 18, and 20 months after enrollment.

```R
Ne_to_N(
  Ne = 100, r = 1, lambda0 = log(2) / 20, lambda1 = log(2) / 20 * 0.75,
  dropoutRate = 0.05, dropoutTime = 12,
  l = c(16, 18, 20), follow_up = "fixed_period"
)

# Ne r    lambda0    lambda1 dropoutRate dropoutTime  a  f  l    follow_up eventRate0
# 1 100 1 0.03465736 0.02599302        0.05          12 NA NA 16 fixed_period  0.4127430
# 2 100 1 0.03465736 0.02599302        0.05          12 NA NA 18 fixed_period  0.4485171
# 3 100 1 0.03465736 0.02599302        0.05          12 NA NA 20 fixed_period  0.4816120
# eventRate1 eventRate       N0       N1        N
# 1  0.3296716 0.3712073 134.6956 134.6956 269.3912
# 2  0.3607570 0.4046370 123.5675 123.5675 247.1351
# 3  0.3900169 0.4358145 114.7277 114.7277 229.4554
```

## Series 2: `getN_xx_xx_JM1`

Given the global sample size, the required consistency probability of the target region's treatment effect with the global treatment effect, and based on Japan's Method 1, calculate the proportion of the target region's sample size to the global sample size.

### Example 1: Binary Endpoint, Superiority Design

The global sample size is 200, with a 1:1 ratio between the treatment group and the control group. The effective rates for the global treatment group and control group are 0.75 and 0.5, respectively. The effective rates for the target region's treatment group and control group are 0.7 and 0.5, respectively. The endpoint is the rate difference (`scale = "RD"`). The target region requires retaining 50% of the global treatment effect (`pi = 0.5`), and the consistency probability of the target region's treatment effect with the global treatment effect is 80% (`beta1 = 0.2`). Under these conditions, the calculated proportion of the target region's sample size to the global sample size is 0.41, with the number of participants being 81.

```R
getN_Bin_Super_JM1(
  p1_a = 0.75, p0_a = 0.5, p1_j = 0.7, p0_j = 0.5, pi = 0.5, beta1 = 0.2,
  N = 200, r = 1, scale = "RD"
)

# p1_a p0_a p1_j p0_j     p1_nj p0_nj  pi alpha beta   N r scale pwr beta1        f       Nj
# 1 0.75  0.5  0.7  0.5 0.7841571   0.5 0.5    NA   NA 200 1    RD 0.8   0.2 0.405873 81.17461
```

### Example 2: Binary Endpoint, Non-Inferiority Design

The effective rates for the global treatment group and control group are 0.5 and 0.5, respectively. The endpoint is the risk ratio (RR), and a smaller RR is better (`direct = -1`). Given the global trial's one-sided alpha = 0.025 and beta = 0.1, the global trial's sample size is calculated to be 192. 

The effective rates for the target region's treatment group and control group are 0.6 and 0.5, respectively. The target region requires retaining 50% of the global treatment effect, and the consistency probability of the target region's treatment effect with the global treatment effect is 80% (`beta1 = 0.2`). Under these conditions, the calculated proportion of the target region's sample size to the global sample size is 0.82, with the number of participants being 158.

In this example, the global trial's sample size is not provided and will be automatically calculated based on the global trial's efficacy parameters, alpha, and beta.

```R
getN_Bin_Noninf_JM1(
  p1_a = 0.5, p0_a = 0.5, p1_j = 0.6, p0_j = 0.5, pi = 0.5, cut = log(1.6),
  alpha = 0.025, beta = 0.1, beta1 = 0.2, N = NA, r = 1, scale = "RR",
  direct = -1
)

# p1_a p0_a p1_j p0_j      p1_nj p0_nj  pi       cut alpha beta   N r scale direct
# 1  0.5  0.5  0.6  0.5 0.03823555   0.5 0.5 0.4700036 0.025  0.1 192 1    RR     -1
# pwr beta1         f      Nj
# 1 0.8000028   0.2 0.8219894 157.822
```

The `getN_xx_xx_JM1` series of functions can calculate multiple parameter values simultaneously, such as when considering non-inferiority margins for RR of 1.5, 1.6, and 1.7.

It is important to note that when the non-inferiority margins are 1.5 and 1.6, although the corresponding `f` can be calculated, the `p1_nj` calculated based on `p1_a` and `p1_j` is negative, so the `f` at this time is also not met the conditions.

```R
getN_Bin_Noninf_JM1(
  p1_a = 0.5, p0_a = 0.5, p1_j = 0.6, p0_j = 0.5, pi = 0.5,
  cut = log(seq(1.5, 1.7, 0.1)), alpha = 0.025, beta = 0.2,
  beta1 = 0.2, N = NA, r = 1, scale = "RR",
  direct = -1
)

# p1_a p0_a p1_j p0_j      p1_nj p0_nj  pi       cut alpha beta   N r scale direct       pwr
# 1  0.5  0.5  0.6  0.5 -3.5374434   0.5 0.5 0.4054651 0.025  0.2 192 1    RR     -1 0.8000002
# 2  0.5  0.5  0.6  0.5 -0.1244289   0.5 0.5 0.4700036 0.025  0.2 144 1    RR     -1 0.8000005
# 3  0.5  0.5  0.6  0.5  0.1851602   0.5 0.5 0.5306283 0.025  0.2 112 1    RR     -1 0.7999981
# beta1         f        Nj
# 1   0.2 0.9758305 187.35945
# 2   0.2 0.8619602 124.12227
# 3   0.2 0.7589431  85.00162
```

## Regional Sample Size under Japan's Method 2

Given the global sample size, the consistency probability of the treatment effect of region $i$ with the global treatment effect, and based on Japan's Method 2, calculate the proportion of region $i$ sample size to the global sample size.

With Japan's Method 2, to prove that the efficacy in region $i$ is consistent with the global efficacy, it is only necessary to ensure that the direction of the efficacy estimate in region $i$ is consistent with the global efficacy direction. The sample size for region $i$ can be calculated using a general sample size formula, setting one-sided alpha = 0.5.

### Example 1: Count Endpoint, Superiority Design

The risk ratio (RR) for the treatment group compared to the control group in region $i$ is 1.2, with a baseline hazard rate of 0.1 per month for the control group. The average exposure duration is 5 months, and the number of events per unit time follows a Poisson distribution (`k = 0` for Poisson distribution and `k > 0` for Negative distribution). The consistency probability of the region $i$ treatment effect with the global treatment effect is 80% (`beta = 0.2`). Under these conditions, the sample size for this region should be at least 158.

```R
getN_Count_Super(
  delta = log(1.2),
  lambda0 = 0.1, t = 5, k = 0, alpha = 0.5, beta = 0.2, N = NA, r = 1
)

# delta lambda0 t k alpha beta   N n1 n0 r       pwr
# 1 0.1823216     0.1 5 0   0.5  0.2 158 79 79 1 0.8013027
```

### Example 2: Count Endpoint, Equivalence Design

The risk ratio (RR) for the treatment group compared to the control group in region $i$ is 1.1, with a equivalence margin for RR being 1.4. The consistency probability of the region $i$ treatment effect with the global treatment effect is 80%. Under these conditions, the sample size for this region should be at least 196.

```R
getN_Count_Equi(
  delta = log(1.1),
  lambda0 = 0.1, t = 5, k = 1, cut = log(1.4),
  alpha = 0.5, beta = 0.2, N = NA, r = 1
)

# delta lambda0 t k       cut alpha beta   N n1 n0 r       pwr
# 1 0.09531018     0.1 5 1 0.3364722   0.5  0.2 196 98 98 1 0.8006631
```

## Series 3: `getPwr_xx_xx_JM1`, `getPwr_xx_xx_JM2`

Given the global and target region sample sizes, calculate and simulate the global success probability and the consistency probability of the target region's treatment effect with the global treatment effect based on Japan's Method 1 and Method 2.

### Example 1: Continuous Endpoint, Superiority Design, Japan's Method 1

The global trial sample size is 100, with a 1:1 ratio between the treatment group and the control group (`r = 1`). The mean differences for the target region and globally are 0.5 and 0.7, respectively, with a common standard deviation of 1. The global trial's one-sided alpha is 0.025, and the target region requires retaining 50% of the global treatment effect (`pi = 0.5`). When the target region's sample size is 50% of the global sample size (`f = 0.5`), based on theoretical calculation (`sim = FALSE`), the results are as follows:

- `pwr1`: The marginal probability of global trial success is 0.94.
- `pwr2`: The marginal probability that the target region's treatment effect is consistent with the global treatment effect is 0.75.
- `pwr3`: The joint probability that the global trial is successful and the target region's treatment effect is consistent with the global treatment effect is 0.72.
- `pwr4`: The conditional probability that the target region's treatment effect is consistent with the global treatment effect, given global trial success, is 0.77.

```R
getPwr_Con_Super_JM1(
  delta_j = 0.5, delta_a = 0.7, sigma = 1, f = 0.5,
  pi = 0.5, alpha = 0.025, beta = NA, N = 100, r = 1, sim = FALSE
)

# delta_a delta_j delta_nj sigma   f  pi alpha beta   N r      pwr1      pwr2      pwr3
# 1     0.7     0.5      0.9     1 0.5 0.5 0.025   NA 100 1 0.9382242 0.7488325 0.7235823
# pwr4
# 1 0.7712253
```

The `getPwr_xx_xx_JM1` series of functions can calculate multiple parameter values simultaneously, such as `f = seq(0.1, 0.9, 0.1)`.

```R
getPwr_Con_Super_JM1(
  delta_j = 0.5, delta_a = 0.7, sigma = 1, f = seq(0.1, 0.9, 0.1),
  pi = 0.5, alpha = 0.025, beta = NA, N = 100, r = 1, sim = FALSE
)

# delta_a delta_j  delta_nj sigma   f  pi alpha beta   N r      pwr1      pwr2
# 1     0.7     0.5 0.7222222     1 0.1 0.5 0.025   NA 100 1 0.9382242 0.5973905
# 2     0.7     0.5 0.7500000     1 0.2 0.5 0.025   NA 100 1 0.9382242 0.6419976
# 3     0.7     0.5 0.7857143     1 0.3 0.5 0.025   NA 100 1 0.9382242 0.6796171
# 4     0.7     0.5 0.8333333     1 0.4 0.5 0.025   NA 100 1 0.9382242 0.7146248
# 5     0.7     0.5 0.9000000     1 0.5 0.5 0.025   NA 100 1 0.9382242 0.7488325
# 6     0.7     0.5 1.0000000     1 0.6 0.5 0.025   NA 100 1 0.9382242 0.7832890
# 7     0.7     0.5 1.1666667     1 0.7 0.5 0.025   NA 100 1 0.9382242 0.8187115
# 8     0.7     0.5 1.5000000     1 0.8 0.5 0.025   NA 100 1 0.9382242 0.8555778
# 9     0.7     0.5 2.5000000     1 0.9 0.5 0.025   NA 100 1 0.9382242 0.8939983
# pwr3      pwr4
# 1 0.5684372 0.6058650
# 2 0.6139979 0.6544255
# 3 0.6524882 0.6954502
# 4 0.6883889 0.7337147
# 5 0.7235823 0.7712253
# 6 0.7592017 0.8091901
# 7 0.7961059 0.8485242
# 8 0.8350790 0.8900634
# 9 0.8771662 0.9349217
```

When the above results are calculated using a simulation method, `sim = FALSE` is changed to `sim = TRUE`, defaulting to 1000 simulations per parameter setting and using 2 CPU cores for computation.

```R
getPwr_Con_Super_JM1(
  delta_j = 0.5, delta_a = 0.7, sigma = 1, f = seq(0.1, 0.9, 0.1),
  pi = 0.5, alpha = 0.025, beta = NA, N = 100, r = 1, sim = TRUE
)

# delta_a delta_j  delta_nj sigma   f  pi alpha beta   N r  pwr1  pwr2  pwr3
# 1     0.7     0.5 0.7222222     1 0.1 0.5 0.025   NA 100 1 0.947 0.589 0.564
# 2     0.7     0.5 0.7500000     1 0.2 0.5 0.025   NA 100 1 0.930 0.637 0.610
# 3     0.7     0.5 0.7857143     1 0.3 0.5 0.025   NA 100 1 0.936 0.664 0.640
# 4     0.7     0.5 0.8333333     1 0.4 0.5 0.025   NA 100 1 0.937 0.714 0.683
# 5     0.7     0.5 0.9000000     1 0.5 0.5 0.025   NA 100 1 0.946 0.747 0.728
# 6     0.7     0.5 1.0000000     1 0.6 0.5 0.025   NA 100 1 0.938 0.788 0.765
# 7     0.7     0.5 1.1666667     1 0.7 0.5 0.025   NA 100 1 0.938 0.810 0.787
# 8     0.7     0.5 1.5000000     1 0.8 0.5 0.025   NA 100 1 0.940 0.874 0.852
# 9     0.7     0.5 2.5000000     1 0.9 0.5 0.025   NA 100 1 0.941 0.890 0.875
# pwr4
# 1 0.5955649
# 2 0.6559140
# 3 0.6837607
# 4 0.7289221
# 5 0.7695560
# 6 0.8155650
# 7 0.8390192
# 8 0.9063830
# 9 0.9298618
```

### Example 2: Binary Endpoint, Non-Inferiority Design, Japan's Method 1

The effective rates for the treatment group and control group in the target region are 0.6 and 0.5, respectively. The effective rates for the treatment group and control group in the non-target region are 0.5 and 0.5, respectively. Given the proportion of the target region's sample size (`f = 0.5`), the effective rates for the global treatment group and control group will be automatically calculated.

The endpoint is the risk ratio (RR), and the non-inferiority margin for RR is 1.4, with `direct = -1` indicating that a smaller RR is better. Given the global trial's one-sided alpha = 0.025 and beta = 0.2, to achieve 80% power, the global trial's sample size is 492. The target region requires retaining 50% of the global treatment effect (`pi = 0.5`). Under these conditions, the power when the target region's sample size is 50% of the global sample size (`f = 0.5`) is calculated.

```R
getPwr_Bin_Noninf_JM1(
  p1_j = 0.6, p0_j = 0.5, p1_nj = 0.5, p0_nj = 0.5, f = 0.5,
  pi = 0.5, cut = log(1.4),
  alpha = 0.025, beta = 0.2, N = NA, r = 1, scale = "RR", direct = -1,
  sim = FALSE
)

# p1_a p0_a p1_j p0_j p1_nj p0_nj    delta_a   delta_j scale   f  pi       cut alpha
# 1 0.55  0.5  0.6  0.5   0.5   0.5 0.09531018 0.1823216    RR 0.5 0.5 0.3364722 0.025
# beta   N r direct      pwr1     pwr2      pwr3      pwr4
# 1  0.2 492 1     -1 0.8009997 0.837892 0.6681649 0.8341637
```

### Example 3: Survival Endpoint, Equivalence Design, Japan's Method 2

A multi-regional clinical trial (mRCT) is conducted in two regions with a total of 400 events, and the  number of events between the treatment group and the control group is 1:1 (a default setting in some sample size calculation software like EAST, i.e., estimating the variance of $log(HR)$ under H0 assumption). Each region is allocated to half of the events (`f_i = c(0.5, 0.5)`). The hazard ratios (HRs) for the treatment group compared to the control group in region 1 and region 2 are 1.1 and 1.0, respectively. The equivalence margin for HR in the global trial is 1.3, and the equivalence margins for HR in region 1 and region 2 are 1.3 and 1.35, respectively. The calculated results are as follows:

- `pwr1`: The probability of global trial success is 0.45.
- `pwr2`: The probability that the treatment effects in region 1 and region 2 are consistent with the global treatment effect is 0.85, with 0.88 for region 1 and 0.97 for region 2.
- pwr3: The probability that the global trial is successful and the treatment effects in region 1 and region 2 are consistent with the global treatment effect is 0.43, with 0.44 for region 1 and 0.44 for region 2.
- pwr4: The probability that the treatment effects in region 1 and region 2 are consistent with the global treatment effect, given global trial success, is 0.97, with 0.97 for region 1 and 0.99 for region 2.

```R
getPwr_Surv_Equi_JM2(
  delta_i = c(log(1.1), log(1.0)),
  f_i = c(0.5, 0.5),
  cut = log(1.3), cut_i = c(log(1.3), log(1.35)),
  alpha = 0.025, beta = NA, Ne = 400, r = 1, sim = FALSE
)

# $overall
# delta_a       cut  Ne      pwr1      pwr2      pwr3      pwr4
# 1 0.04765509 0.2623643 400 0.4471244 0.8459097 0.4334924 0.9695117
# 
# $cut_i
# [1] 0.2623643 0.3001046
# 
# $pwr_margin
# [1] 0.8755313 0.9661673
# 
# $pwr_joint
# [1] 0.4352447 0.4434537
# 
# $pwr_condition
# [1] 0.9734309 0.9917905
```

