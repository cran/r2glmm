% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r2beta_Function.R
\name{r2beta}
\alias{r2beta}
\title{r2beta Compute R Squared for Mixed Models}
\usage{
r2beta(model, partial = TRUE, method = "sgv", data = NULL)
}
\arguments{
\item{model}{a fitted mermod, lme, or glmmPQL model.}

\item{partial}{if TRUE, semi-partial R squared are calculated for each
fixed effect in the mixed model.}

\item{method}{Specifies the method of computation for R squared beta:
if \code{method} = 'sgv' then the standardized generalized variance
approach is applied. This method is recommended for covariance model
selection.
if \code{method} = 'kr', then the Kenward Roger approach is applied.
This option is only available for \code{\link[nlme]{lme}} models.
if \code{method} = 'nsj',then the Nakagawa and Schielzeth approach
is applied. This option is available for
\code{\link[lme4]{lmer}} and \code{\link[nlme]{lme}} objects.
if \code{method} = 'lm', the classical R squared from the
linear model is computed. This method should only be used
on glm and lm object.}

\item{data}{The data used by the fitted model. This argument is required
for models with special expressions in their formula, such as
offset, log, cbind(sucesses, trials), etc.}
}
\value{
A dataframe containing the model F statistic, numerator
  and denominator degrees of freedom, non-centrality parameter,
  and R squared statistic with 95 percent confidence limits.
  If partial = TRUE, then the dataframe also contains partial R
  squared statistics for all fixed effects in the model.
}
\description{
Computes coefficient of determination (R squared) from
edwards et al., 2008 and the generalized R squared from Jaeger et al., 2016.
Currently implemented for linear mixed models with
\code{\link[lme4]{lmer}} and \code{\link[nlme]{lme}} objects. For
generalized linear mixed models, only \code{\link[MASS]{glmmPQL}} are supported.
}
\examples{
library(nlme)
library(lme4)
data(Orthodont)

# Linear mixed models
mermod = lmer(distance ~ age*Sex + (1|Subject), data = Orthodont)
lmemod = lme(distance ~ age*Sex, random = ~1|Subject, data = Orthodont)

# The Kenward-Roger approach
r2beta(mermod, method = 'kr')

# Standardized Generalized Variance
r2beta(mermod, method = 'sgv')
r2beta(lmemod, method = 'sgv')

# The marginal R squared by Nakagawa and Schielzeth (extended by Johnson)
r2beta(mermod, method = 'nsj')

# linear and generalized linear models

library(datasets)
dis = data.frame(discoveries)
dis$year = 1:nrow(dis)

lmod = lm(discoveries ~ year + I(year^2), data = dis)
glmod = glm(discoveries ~ year + I(year^2), family = 'poisson', data = dis)

# Using an inappropriate link function (normal) leads to
# a poor fit relative to the poisson link function.

r2beta(lmod)
r2beta(glmod)

# PQL models
# Currently only SGV method is supported
library(MASS)
PQL_bac = glmmPQL(y ~ trt + I(week > 2), random = ~ 1 | ID,
                  family = binomial, data = bacteria,
                  verbose = FALSE)

r2beta(PQL_bac, method='sgv')

}
\references{
Edwards, Lloyd J., et al. "An R2 statistic for fixed effects in
   the linear mixed model." Statistics in medicine 27.29 (2008): 6137-6157.

Nakagawa, Shinichi, and Holger Schielzeth. "A general and simple method for
   obtaining R2 from generalized linear mixed effects models." Methods in
   Ecology and Evolution 4.2 (2013): 133-142.

Jaeger, Byron C., et al., "An R Squared Statistic for Fixed Effects in the
   Generalized Linear Mixed Model." Journal of Applied Statistics (2016).
}
