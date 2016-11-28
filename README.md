
<!-- README.md is generated from README.Rmd. Please edit that file -->
r2glmm
======

This package computes model and semi partial R squared with confidence limits for the linear and generalized linear mixed model (LMM and GLMM). The R squared measure from Edwards et.al (2008) is extended to the GLMM using penalized quasi-likelihood (PQL) estimation (see Jaeger et al. 2016).

-   Changes: Version 0.1.1

1.  Updated the r2beta function with an optional data argument. Users who wish to use a function (i.e. log(x)) in their model formula should specify the original data frame used by the model when using the r2beta function.
2.  Included support for GLMMs fitted using the glmer function from the lme4 package.
3.  Incorporated a weighted adjustment for binary mixed models for improved model selection.
4.  Added generic plot and print functions for R2 objects.

-   Why use this package?

The R Squared statistic is a well known tool that describes goodness-of-fit for a statistical model. In the linear model, R Squared may be interpreted as the proportion of variance in the data explained by the fixed predictors, which is often of interest to investigators in social and biological sciences. Other popular uses of R squared include meta-analyses and standardized measurement of effect sizes using semi-partial R squared. The semi-partial R squared statistic corresponds to a select subset of fixed predictors in a fitted model and measures the relative increase in association between the dependent and independent variables resulting from the inclusion of the specified subset. Semi-partial R squared allow investigators to select a set of predictors based on both statistical significance and relative importance.

Currently information criteria dominate the applied practice of selecting the most parsimonious mixed model. These criteria provide guidance, but cannot be used to measure goodness-of-fit. Further, they can only be used to compare models fitted to the same data. Lastly, the information criteria do not allow investigators to assess individual fixed predictors. Thus, it is beneficial to apply the information criteria in conjunction with R Squared statistics when conducting statistical inference on mixed models.

-   Instructions for installation:

Currently, the r2glmm package is available at my github site. After installing and loading the devtools package, run this code from the R console:

devtools::install\_github('bcjaeger/r2glmm')

The files should then be downloaded and installed.

-   How to use this package

The main function in this package is called r2beta. A user may fit a mixed model for one of the supported model types and then apply the r2beta function using the specified model as input. Additionally, the investigator may specify whether semi-partial R Squared are computed (they are by default) and what type of method to employ for computation. Three methods of computation are currently provided:

1.  An approach using standardized generalized variance (SGV) that can be used for covariance model selection. The SGV approach chooses the model which explains the greatest amount of generalized variance. Below is an example using the dental data from the nlme package.

``` r

library(lme4)
#> Loading required package: Matrix
library(nlme)
#> 
#> Attaching package: 'nlme'
#> The following object is masked from 'package:lme4':
#> 
#>     lmList
library(r2glmm)

data(Orthodont)

# Linear mixed models

# m1 is a compound symmetric (CS) or random intercept model.
m1 = lme(distance ~ age*Sex,  ~1|Subject, data = Orthodont)

# m2 is an order 1 autoregressive (AR1) model with
# heterogeneity of variance between genders
m2 = lme(distance ~ age*Sex, data=Orthodont, 
         correlation = corAR1(form=~1|Subject),
         weights = varIdent(form=~1|Sex))

# Compare the models
(r2m1 = r2beta(m1,method='sgv'))
#>          Effect   Rsq upper.CL lower.CL
#> 1         Model 0.559    0.669    0.447
#> 2           age 0.392    0.527    0.256
#> 4 age:SexFemale 0.038    0.144    0.000
#> 3     SexFemale 0.004    0.067    0.000
(r2m2 = r2beta(m2,method='sgv'))
#>          Effect   Rsq upper.CL lower.CL
#> 1         Model 0.616    0.713    0.514
#> 2           age 0.454    0.580    0.323
#> 4 age:SexFemale 0.051    0.165    0.001
#> 3     SexFemale 0.006    0.075    0.000
```

The SGV \(R^2\) statistic follows an asymptotic non-central beta distribution of type I and the difference between two SGV \(R^2\) statistics can be tested with the r2dt function.

``` r

# Use cor = TRUE since the models were fit using the same data
(test = r2dt(r2m1, r2m2, cor = TRUE))
#> Warning in stats::qbeta(samps + stats::rnorm(nsims, sd = 0.1), shape1 = yy
#> $v1/2, : NaNs produced
#> $d
#> [1] -0.05663145
#> 
#> $ci
#> [1] -0.09578086 -0.01631509
#> 
#> $p
#> [1] 0.002002002
# d is the difference between model R^2
# ci is a 95% confidence interval for d
# p is the p-value corresponding to H0: d = 0
```

Finally, the plot function can summarize our results. On the top left and right there are panelled plots of model and semi-partial \(R^2\). Below is a plot with the estimated difference between model \(R^2\) with 95% confidence interval. A vertical line at 0 helps illustrate statistical significance.

``` r

plot(x=r2m1, y=r2m2, r2labs = c('CS', 'AR1'),
     txtsize = 8, r2mthd='sgv', cor = TRUE)
#> Warning in stats::qbeta(samps + stats::rnorm(nsims, sd = 0.1), shape1 = yy
#> $v1/2, : NaNs produced
```

![](README-unnamed-chunk-4-1.png)

1.  The Kenward-Roger approach applies the small sample approximation to the F statistic using the pbkrtest package, and is recommended for selecting fixed effects. Rhe Kenward-Roger approach is limited to mermod objects.

``` r

# Compute the R2 statistic using the Kenward-Roger approach.

mermod = lmer(distance ~ age*Sex + (1|Subject), data = Orthodont)

r2beta(mermod, method = 'kr')
#>          Effect   Rsq upper.CL lower.CL
#> 1         Model 0.671    0.771    0.563
#> 2           age 0.565    0.680    0.437
#> 4 age:SexFemale 0.074    0.212    0.004
#> 3     SexFemale 0.004    0.065    0.000
```

1.  The method introduced by Nakagawa and Schielzeth (2013) computes the marginal R squared described by Nakagawa and Schielzeth, which was later modified by Johnson (2014). Additionally, this package computes confidence limits and semi-partial R Squared for the marginal coefficient of determination. The r2glmm package only computes marginal R squared for the LMM and does not generalize the statistic to the GLMM.

``` r

# Compute the R2 statistic using Nakagawa and Schielzeth's approach.
(r2nsj = r2beta(mermod, method = 'nsj', partial = TRUE))
#>          Effect   Rsq upper.CL lower.CL
#> 1         Model 0.410    0.540    0.290
#> 2           age 0.261    0.398    0.137
#> 4 age:SexFemale 0.021    0.105    0.000
#> 3     SexFemale 0.002    0.055    0.000

# Check the result with MuMIn's r.squaredGLMM
r2nsj_mum = MuMIn::r.squaredGLMM(mermod)

all.equal(r2nsj[1,'Rsq'],as.numeric(r2nsj_mum[1]), tolerance = 1e-5)
#> [1] TRUE
```

The semi-partial R squared and confidence limits are two useful and exclusive features of the r2beta function which are used for covariate selection and interpretation of effect size (or multivariate association with the outcome). The plot function may be used to visualize these statistics.

``` r

plot(r2nsj,r2mthd='nsj')
```

![](README-unnamed-chunk-7-1.png)

Update 0.1.1 allows the r2beta function to handle models fitted using the glmer function from the lme4 package. Additionally, using the data argument in the r2beta function allows users to assess goodness-of-fit for models that have a function (i.e. log(x) or cbind(events, trials)) in their model's formula.

``` r

library(lattice)

cbpp$period = as.numeric(cbpp$period)
gm1 <- glmer(
  formula=cbind(incidence, size-incidence) ~ poly(period,2) + (1|herd),
  data = cbpp, family = binomial)

r2beta(model = gm1, method = 'sgv', data = cbpp)
#>             Effect   Rsq upper.CL lower.CL
#> 1            Model 0.387    0.589    0.202
#> 2 poly(period, 2)1 0.383    0.579    0.187
#> 3 poly(period, 2)2 0.011    0.147    0.000
```

Longitudinal analyses sometimes have multiple valid approaches to model the progression of time. For example, the the toenail oncholysis data from Backer et al 1998 can be modelled with a continuous time covariate or a visit number. AIC and BIC favor the latter, but do not provide a direct measure of how well the favored model explains the data.

``` r
library(HSAUR2)
#> Loading required package: tools

gm2_time <- glmer(outcome~treatment*time+(1|patientID),
                 data=toenail,
                 family=binomial,nAGQ=20)

gm2_visit <- glmer(outcome~treatment*visit+(1|patientID),
                 data=toenail,
                 family=binomial,nAGQ=20)

AIC(gm2_time)
#> [1] 1260.75
AIC(gm2_visit)
#> [1] 1252.388

BIC(gm2_time)
#> [1] 1288.519
BIC(gm2_visit)
#> [1] 1280.157
```

The R Squared is in agreement with information criteria, but demonstrates that the favored model is substantially more well suited to the observed toenail data.

``` r

r2beta(gm2_time, method = 'sgv')
#>                      Effect   Rsq upper.CL lower.CL
#> 1                     Model 0.337    0.369    0.306
#> 3                      time 0.119    0.147    0.093
#> 4 treatmentterbinafine:time 0.000    0.004    0.000
#> 2      treatmentterbinafine 0.000    0.003    0.000

r2beta(gm2_visit, method = 'sgv')
#>                       Effect   Rsq upper.CL lower.CL
#> 1                      Model 0.704    0.722    0.686
#> 3                      visit 0.312    0.344    0.280
#> 4 treatmentterbinafine:visit 0.012    0.024    0.004
#> 2       treatmentterbinafine 0.000    0.003    0.000
```

The R squared statistic and semi-partial R squared may be used to measure effect size. Sometimes large samples can have significant estimates of effect size (\(p < .05\)), but the estimated values are practically trivial.

``` r

library(geepack)
#> Warning: package 'geepack' was built under R version 3.3.1

# The ohio dataset is a subset of a six-city study, 
# longitudinal study of the health effects of air pollution.

data(ohio)

fit <- glmer(resp ~ age + smoke + age:smoke + (1|id),
             data=ohio,family=binomial)

summary(fit)
#> Generalized linear mixed model fit by maximum likelihood (Laplace
#>   Approximation) [glmerMod]
#>  Family: binomial  ( logit )
#> Formula: resp ~ age + smoke + age:smoke + (1 | id)
#>    Data: ohio
#> 
#>      AIC      BIC   logLik deviance df.resid 
#>   1599.3   1627.7   -794.7   1589.3     2143 
#> 
#> Scaled residuals: 
#>     Min      1Q  Median      3Q     Max 
#> -1.3995 -0.1778 -0.1589 -0.1276  2.6024 
#> 
#> Random effects:
#>  Groups Name        Variance Std.Dev.
#>  id     (Intercept) 5.502    2.346   
#> Number of obs: 2148, groups:  id, 537
#> 
#> Fixed effects:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -3.40171    0.27882 -12.201   <2e-16 ***
#> age         -0.21704    0.08678  -2.501   0.0124 *  
#> smoke        0.47824    0.29923   1.598   0.1100    
#> age:smoke    0.10465    0.13911   0.752   0.4519    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Correlation of Fixed Effects:
#>           (Intr) age    smoke 
#> age        0.272              
#> smoke     -0.442 -0.193       
#> age:smoke -0.146 -0.621  0.280
```

Here we see that the effect of age has a statistically significant p-value, but this may be due to the large sample size. R squared can help assess whether the effect is substantial or not (it isn't).

``` r

r2beta(fit,method='sgv')
#>      Effect   Rsq upper.CL lower.CL
#> 1     Model 0.042    0.061    0.028
#> 2       age 0.028    0.043    0.016
#> 3     smoke 0.009    0.019    0.003
#> 4 age:smoke 0.003    0.009    0.000
```
