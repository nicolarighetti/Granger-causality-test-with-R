Four ways to perform Granger Causality test with R
================

The Granger causal test is widely used in time series analysis, but the
various ways in which it may be performed in R are poorly documented. I
published before a function to perform the Toda-Yamamoto bivariate and
multivariate versions of the Granger causality test (see at
[Toda-Yamamoto-Causality-Test](https://github.com/nicolarighetti/Toda-Yamamoto-Causality-Test/blob/main/README.md)).
In the few lines that follow, I am going to focus on its standard
version, suggesting four different approaches.

For starters, a definition. A time series X is said to Granger cause
another time series Y if past values of X and Y predict Y significantly
better than past values of Y only. (Granger, 1969, click
[here](https://github.com/nicolarighetti/Toda-Yamamoto-Causality-Test/blob/main/README.md)
for other details). We can test that hypothesis using linear regression
modeling and a Wald test for linear restrictions. The Wald test compares
the performance of a restricted model for Y, which excludes X, to a
non-restricted model for Y, which includes X.

Firstly, let’s create two simulated time series. It is important to note
that the data are represented as time series using the base R *ts*
function. Some functions described here do not work well if the data is
not in this format. Also note that the Granger causality test has
several assumptions. I’m not going to go into that. I just add that it’s
based on linear modeling, so the linear regression assumptions apply to
it.

## First approach: lmtest::grangertest

The first way to perform a Granger causality test in R is to use the
*grangertest* function provided by the
[lmtest](https://cran.r-project.org/web/packages/lmtest/index.html)
package. We can test if the series *ts1* Granger-causes the series *t2*
with a simple line of code. It turns out the *ts1* Granger-causes *ts2*.

``` r
lmtest::grangertest(ts1, ts2, order=3, test="F")
```

    ## Granger causality test
    ## 
    ## Model 1: ts2 ~ Lags(ts2, 1:3) + Lags(ts1, 1:3)
    ## Model 2: ts2 ~ Lags(ts2, 1:3)
    ##   Res.Df Df      F  Pr(>F)  
    ## 1     90                    
    ## 2     93 -3 2.7758 0.04585 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

While *grangertest* is a quick way to perform a Granger causal test, it
is quite rough, as it does not rely on any modelling of the series. It
is important to fit an appropriate linear model to the series, also just
to check if the assumptions hold.

## Second approach: dynlm::dynlm + lmtest::waldtest

A second, more flexible approach is based on performing a Wald test on
two fitted linear regression models. As defined in the Granger causality
test, linear regression models include past lags of the variables. This
type of model can be easily fitted with the function
[dynlm](https://cran.r-project.org/web/packages/dynlm/index.html) of the
homonymous R package. Next, if the assumptions of the linear models and
Granger test hold, a Wald test can be performed by using the function
*waldtest* of the already mentioned library *lmtest*.

``` r
library(dynlm)
```

``` r
unrestricted_ts2 <- dynlm(ts2 ~ L(ts2, 1:3) + L(ts1, 1:3))
restricted_ts2 <- dynlm(ts2 ~ L(ts2, 1:3))

lmtest::waldtest(unrestricted_ts2, restricted_ts2, test="F")                  
```

    ## Wald test
    ## 
    ## Model 1: ts2 ~ L(ts2, 1:3) + L(ts1, 1:3)
    ## Model 2: ts2 ~ L(ts2, 1:3)
    ##   Res.Df Df      F  Pr(>F)  
    ## 1     90                    
    ## 2     93 -3 2.7758 0.04585 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

This approach has the advantage that the series can be modeled in a
flexible manner before testing Granger’s causality. For instance, by
adding a seasonal term or other variables when necessary.

Moreover, this approach permits the specification of
heteroscedasticity-consistent (HC) or heteroscedasticity and
autocorrelation consistent standard errors (HAC), also in the Newey West
form, for instance by using the functions provided by the
[sandwich](https://cran.r-project.org/web/packages/sandwich/index.html)
package.

``` r
lmtest::waldtest(unrestricted_ts2, restricted_ts2, test="F", vcov = sandwich::vcovHC) 
```

    ## Wald test
    ## 
    ## Model 1: ts2 ~ L(ts2, 1:3) + L(ts1, 1:3)
    ## Model 2: ts2 ~ L(ts2, 1:3)
    ##   Res.Df Df      F  Pr(>F)  
    ## 1     90                    
    ## 2     93 -3 2.7386 0.04801 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
lmtest::waldtest(unrestricted_ts2, restricted_ts2, test="F", vcov = sandwich::vcovHAC) 
```

    ## Wald test
    ## 
    ## Model 1: ts2 ~ L(ts2, 1:3) + L(ts1, 1:3)
    ## Model 2: ts2 ~ L(ts2, 1:3)
    ##   Res.Df Df      F  Pr(>F)  
    ## 1     90                    
    ## 2     93 -3 3.3839 0.02157 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
lmtest::waldtest(unrestricted_ts2, restricted_ts2, test="F", vcov = sandwich::NeweyWest) 
```

    ## Wald test
    ## 
    ## Model 1: ts2 ~ L(ts2, 1:3) + L(ts1, 1:3)
    ## Model 2: ts2 ~ L(ts2, 1:3)
    ##   Res.Df Df      F  Pr(>F)  
    ## 1     90                    
    ## 2     93 -3 3.9935 0.01017 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Third approach: dynlm::dynlm + aod::wald.test

The third approach is a version of the latter. A Wald test can also be
performed with the wald.test function in the
[aod](https://cran.r-project.org/web/packages/aod/index.html) package.
This is trickier because you have to manually specify the terms to be
tested. However, this extra flexibility may be useful in certain cases.
By default, *wald.test* will calculate a Wald Chi-square test. The
Chi-squared version of the test can be also computed with the *waldtest*
function of *lmtest* by specifying the option test=“Chisq”. Similarly,
the F or Chi-squared test can be performed by using the the function
*linearHypothesis* in the
[car](https://cran.r-project.org/web/packages/car/index.html) package.

``` r
aod::wald.test(b=coef(unrestricted_ts2), 
               Sigma=vcov(unrestricted_ts2),
               Terms = 5:7,
               verbose = T)
```

    ## Wald test:
    ## ----------
    ## 
    ## Coefficients:
    ##  (Intercept) L(ts2, 1:3)1 L(ts2, 1:3)2 L(ts2, 1:3)3 L(ts1, 1:3)1 L(ts1, 1:3)2 
    ##        0.094        0.114       -0.038       -0.132        0.091       -0.139 
    ## L(ts1, 1:3)3 
    ##        0.279 
    ## 
    ## Var-cov matrix of the coefficients:
    ##              (Intercept) L(ts2, 1:3)1 L(ts2, 1:3)2 L(ts2, 1:3)3 L(ts1, 1:3)1
    ## (Intercept)   0.01196    -0.00076     -0.00061     -0.00073      0.00169    
    ## L(ts2, 1:3)1 -0.00076     0.01070     -0.00111      0.00045     -0.00031    
    ## L(ts2, 1:3)2 -0.00061    -0.00111      0.01079     -0.00103      0.00097    
    ## L(ts2, 1:3)3 -0.00073     0.00045     -0.00103      0.01065     -0.00036    
    ## L(ts1, 1:3)1  0.00169    -0.00031      0.00097     -0.00036      0.01146    
    ## L(ts1, 1:3)2  0.00170    -0.00119     -0.00022      0.00076     -0.00133    
    ## L(ts1, 1:3)3  0.00215     0.00084     -0.00129     -0.00030     -0.00023    
    ##              L(ts1, 1:3)2 L(ts1, 1:3)3
    ## (Intercept)   0.00170      0.00215    
    ## L(ts2, 1:3)1 -0.00119      0.00084    
    ## L(ts2, 1:3)2 -0.00022     -0.00129    
    ## L(ts2, 1:3)3  0.00076     -0.00030    
    ## L(ts1, 1:3)1 -0.00133     -0.00023    
    ## L(ts1, 1:3)2  0.01180     -0.00120    
    ## L(ts1, 1:3)3 -0.00120      0.01169    
    ## 
    ## Test-design matrix:
    ##    (Intercept) L(ts2, 1:3)1 L(ts2, 1:3)2 L(ts2, 1:3)3 L(ts1, 1:3)1 L(ts1, 1:3)2
    ## L1           0            0            0            0            1            0
    ## L2           0            0            0            0            0            1
    ## L3           0            0            0            0            0            0
    ##    L(ts1, 1:3)3
    ## L1            0
    ## L2            0
    ## L3            1
    ## 
    ## Positions of tested coefficients in the vector of coefficients: 5, 6, 7 
    ## 
    ## H0:  L(ts1, 1:3)1 = 0; L(ts1, 1:3)2 = 0; L(ts1, 1:3)3 = 0 
    ## 
    ## Chi-squared test:
    ## X2 = 8.3, df = 3, P(> X2) = 0.04

## Fourth approach: vars::VAR + vars::causality

The fourth way to perform a Granger causality test in R is probably the
most commonly used, and consists in using the package
[vars](https://cran.r-project.org/web/packages/vars/index.html) to fit a
VAR model (which is basically a system of linear regression models) and
subsequently, run the Granger test using the function *causality*
provided by the same package.

``` r
tsDat <- ts.union(ts1, ts2) 
tsVAR <- vars::VAR(tsDat, p = 3)
vars::causality(tsVAR, cause = "ts1")$Granger
```

    ## 
    ##  Granger causality H0: ts1 do not Granger-cause ts2
    ## 
    ## data:  VAR object tsVAR
    ## F-Test = 2.7758, df1 = 3, df2 = 180, p-value = 0.04276

The *causation* function offers several interesting options. For
example, it implements a bootstrapping procedure. It is also possible to
specify the heteroscedasticity-consistent estimation of the covariance
matrix of the coefficient estimates in the regression models (HC), but
other types of estimation of the covariance matrix (e.g., HAC) are not
supported.
