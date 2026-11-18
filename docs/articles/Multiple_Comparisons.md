# Multiple Comparisons for MixMod Objects

## Multiple Comparisons with MixMod Objects

In this vignette we illustrate how to correct p-values for multiple
comparisons using the
[**multcomp**](https://cran.r-project.org/package=multcomp) and
[**emmeans**](https://cran.r-project.org/package=emmeans) packages.

### Additive Model

We start by simulating some data for a binary longitudinal outcome:

``` r

set.seed(1234)
n <- 300 # number of subjects
K <- 4 # number of measurements per subject
t_max <- 15 # maximum follow-up time

# we construct a data frame with the design: 
# everyone has a baseline measurement, and then measurements at K time points
DF <- data.frame(id = rep(seq_len(n), each = K),
                 time = gl(K, 1, n*K, labels = paste0("Time", 1:K)),
                 sex = rep(gl(2, n/2, labels = c("male", "female")), each = K))

# design matrices for the fixed and random effects
X <- model.matrix(~ sex * time, data = DF)
Z <- model.matrix(~ 1, data = DF)

betas <- c(-2.13, 1, rep(c(1.2, -1.2), K-1)) # fixed effects coefficients
D11 <- 1 # variance of random intercepts

# we simulate random effects
b <- cbind(rnorm(n, sd = sqrt(D11)))
# linear predictor
eta_y <- as.vector(X %*% betas + rowSums(Z * b[DF$id, ]))
# we simulate binary longitudinal data
DF$y <- rbinom(n * K, 1, plogis(eta_y))
```

We fit a mixed effects logistic regression for `y` assuming random
intercepts for the random-effects part, and the main effects of `sex`
and `time` for the fixed-effects part.

``` r

fm <- mixed_model(fixed = y ~ sex + time, random = ~ 1 | id, data = DF,
                  family = binomial())
```

The uncorrected p-values for the 4 time points are give by the
[`summary()`](https://rdrr.io/r/base/summary.html) method:

``` r

summary(fm)
#> 
#> Call:
#> mixed_model(fixed = y ~ sex + time, random = ~1 | id, data = DF, 
#>     family = binomial())
#> 
#> Data Descriptives:
#> Number of Observations: 1200
#> Number of Groups: 300 
#> 
#> Model:
#>  family: binomial
#>  link: logit 
#> 
#> Fit statistics:
#>    log.Lik      AIC      BIC
#>  -583.4797 1178.959 1201.182
#> 
#> Random effects covariance matrix:
#>               StdDev
#> (Intercept) 1.447826
#> 
#> Fixed effects:
#>             Estimate Std.Err z-value    p-value
#> (Intercept)  -2.2865  0.2496 -9.1601    < 1e-04
#> sexfemale     0.8860  0.2447  3.6200 0.00029455
#> timeTime2     0.3066  0.2266  1.3534 0.17592327
#> timeTime3    -0.5053  0.2463 -2.0516 0.04020481
#> timeTime4     0.5632  0.2237  2.5181 0.01179888
#> 
#> Integration:
#> method: adaptive Gauss-Hermite quadrature rule
#> quadrature points: 11
#> 
#> Optimization:
#> method: hybrid EM and quasi-Newton
#> converged: TRUE
```

To perform the pairwise comparisons and obtain corrected p-values, we
load the **multcomp** package and use the
[`glht()`](https://rdrr.io/pkg/multcomp/man/glht.html) function. Because
no specific methods exist for `MixMod` object returned by
[`mixed_model()`](https://drizopoulos.github.io/GLMMadaptive/reference/mixed_model.md),
we need to specify the `vcov.` and `coef.` arguments of
[`glht()`](https://rdrr.io/pkg/multcomp/man/glht.html), i.e.,

``` r

library("multcomp")
#> Loading required package: mvtnorm
#> Loading required package: survival
#> Loading required package: TH.data
#> Loading required package: MASS
#> 
#> Attaching package: 'MASS'
#> The following object is masked from 'package:GLMMadaptive':
#> 
#>     negative.binomial
#> 
#> Attaching package: 'TH.data'
#> The following object is masked from 'package:MASS':
#> 
#>     geyser
fm_mc <- glht(fm, linfct = mcp(time = "Tukey"),
           vcov. = vcov(fm, "fixed"), coef. = fixef)

summary(fm_mc)
#> 
#>   Simultaneous Tests for General Linear Hypotheses
#> 
#> Multiple Comparisons of Means: Tukey Contrasts
#> 
#> 
#> Fit: mixed_model(fixed = y ~ sex + time, random = ~1 | id, data = DF, 
#>     family = binomial())
#> 
#> Linear Hypotheses:
#>                    Estimate Std. Error z value Pr(>|z|)    
#> Time2 - Time1 == 0   0.3066     0.2266   1.353  0.52807    
#> Time3 - Time1 == 0  -0.5053     0.2463  -2.052  0.16892    
#> Time4 - Time1 == 0   0.5632     0.2237   2.518  0.05695 .  
#> Time3 - Time2 == 0  -0.8119     0.2423  -3.351  0.00463 ** 
#> Time4 - Time2 == 0   0.2566     0.2165   1.185  0.63542    
#> Time4 - Time3 == 0   1.0685     0.2406   4.441  < 0.001 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> (Adjusted p values reported -- single-step method)
```

### Interaction Model

We continue our illustration by including the interaction term between
`sex` and `time`, and we focus on the difference between males and
females for the various time points. We start by fitting the model:

``` r

gm <- mixed_model(fixed = y ~ sex * time, random = ~ 1 | id, data = DF,
                  family = binomial())
```

To compute the estimated log odds for males and females at the different
time points we use the
[`emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.html)
functions from the **emmeans** package:

``` r

library("emmeans")
gm_mc <- emmeans(gm, ~ sex | time)

gm_mc
#> time = Time1:
#>  sex    emmean    SE  df asymp.LCL asymp.UCL
#>  male    -2.74 0.342 Inf     -3.41    -2.067
#>  female  -1.21 0.255 Inf     -1.71    -0.707
#> 
#> time = Time2:
#>  sex    emmean    SE  df asymp.LCL asymp.UCL
#>  male    -1.68 0.275 Inf     -2.21    -1.135
#>  female  -1.44 0.263 Inf     -1.95    -0.920
#> 
#> time = Time3:
#>  sex    emmean    SE  df asymp.LCL asymp.UCL
#>  male    -4.40 0.546 Inf     -5.47    -3.331
#>  female  -1.39 0.261 Inf     -1.90    -0.876
#> 
#> time = Time4:
#>  sex    emmean    SE  df asymp.LCL asymp.UCL
#>  male    -1.34 0.261 Inf     -1.85    -0.825
#>  female  -1.25 0.256 Inf     -1.75    -0.749
#> 
#> Results are given on the logit (not the response) scale. 
#> Confidence level used: 0.95
```

The corresponding pairwise comparisons are performed by the
[`pairs()`](https://rdrr.io/r/graphics/pairs.html) function:

``` r

pairs(gm_mc)
#> time = Time1:
#>  contrast      estimate    SE  df z.ratio p.value
#>  male - female  -1.5305 0.410 Inf  -3.735  0.0002
#> 
#> time = Time2:
#>  contrast      estimate    SE  df z.ratio p.value
#>  male - female  -0.2400 0.366 Inf  -0.655  0.5125
#> 
#> time = Time3:
#>  contrast      estimate    SE  df z.ratio p.value
#>  male - female  -3.0126 0.586 Inf  -5.145 <0.0001
#> 
#> time = Time4:
#>  contrast      estimate    SE  df z.ratio p.value
#>  male - female  -0.0851 0.355 Inf  -0.240  0.8106
#> 
#> Results are given on the log odds ratio (not the response) scale.
```

For additional examples in testing interactions with the **emmeans**
package check the vignette:
[`vignette("interactions", package = "emmeans")`](https://rvlenth.github.io/emmeans/articles/interactions.html).
