

#'  R Squared Difference Test (R2DT). Test for a statistically significant difference in genealized explained variance between two candidate models.
#'
#' @param x An R2 object from the r2beta function.
#' @param y An R2 object from the r2beta function.
#' @param fancy if TRUE, the output values are rounded and changed to characters.
#' @param cor if TRUE, the R squared statistics are assumed to be positively correlated and a simulation based approach is used. If FALSE, the R squared are assumed independent and the difference of independent beta distributions is used.
#' @param onesided if TRUE, the alternative hypothesis is that one model explains a larger proportion of generalized variance. If false, the alternative is that the amount of generalized variance explained by the two candidate models is not equal.
#' @param nsims number of samples to draw when simulating correlated non-central beta random variables. This parameter is only relevant if cor=TRUE.
#' @param clim Desired confidence level for interval estimates regarding the difference in generalized explained variance.
#' @return A confidence interval for the difference in R Squared statistics and a p-value corresponding to the null hypothesis of no difference.
#' @examples
#' library(nlme)
#' library(lme4)
#' library(r2glmm)
#'
#' data(Orthodont)
#'
#' # Comparing two linear mixed models
#' m1 = lmer(distance ~ age*Sex+(1|Subject), Orthodont)
#' m2 = lmer(distance ~ age*Sex+(1+age|Subject), Orthodont)
#'
#' m1r2 = r2beta(model=m1,partial=FALSE)
#' m2r2 = r2beta(model=m2,partial=FALSE)
#'
#' # Accounting for correlation can make a substantial difference.
#'
#' r2dt(x=m1r2, y = m2r2, cor = TRUE)
#' r2dt(x=m1r2, y = m2r2, cor = FALSE)
#'
#'
#' @export

r2dt = function(x, y, cor = FALSE, fancy=FALSE,
                onesided=TRUE, clim=95,
                nsims=1000){

    alpha = ifelse(clim>1, 1-clim/100, 1-clim)
    limits = c(0+alpha/2, 1-alpha/2)

  pbetas <- function (z,a=1,b=1,c=1,d=1,ncp1=0,ncp2=0){
    n=length(z)
    tt=z
    for (i in 1:length(z)){
      f=function (xx) {
        dbetas(xx,a=a,b=b,c=c,d=d,ncp1=ncp1,ncp2=ncp2)
      }
      tt[i]=stats::integrate(f,lower=-1,upper=z[i])$value
    }
    return(tt)
  }

  dbetas <- function (z,a=1,b=1,c=1,d=1,ncp1=0,ncp2=0){
    n=length(z)
    tt=rep(0,n)
    for (i in 1:length(z))
    {f=function (xx) {
      stats::dbeta(xx,shape1=a,shape2=b,ncp=ncp1)*
        stats::dbeta(xx-z[i],shape1=c,shape2=d,ncp=ncp2)
    }
    if (z[i]>0&z[i]<=1){
      tt[i]=stats::integrate(f,lower=z[i],upper=1)$value
    }
    if (z[i]>=-1&z[i]<=0){
      tt[i]=stats::integrate(f,lower=0,upper=1+z[i])$value}
    }
    return(tt)
  }

  xx = x[1,]; yy = y[1,]; diff = xx$Rsq - yy$Rsq

  if(cor){

    # It is assumed that there is a positive
    # correlation between R squared statistics
    # when they have identical mean models.

    samps = stats::rnorm(nsims, mean = 0.50, sd = 0.08)
    samps[samps > 0.75] = 0.75
    samps[samps < 0.25] = 0.25

    df = data.frame(
      r2x = stats::qbeta(samps+stats::rnorm(nsims, sd=0.1),
                  shape1 = xx$v1/2,
                  shape2 = xx$v2/2,
                  ncp    = xx$ncp),
      r2y = stats::qbeta(samps+stats::rnorm(nsims, sd=0.1),
                  shape1 = yy$v1/2,
                  shape2 = yy$v2/2,
                  ncp    = yy$ncp))

    diffs = df$r2x-df$r2y
    diffs = diffs[!is.na(diffs)]

    pval = ifelse(onesided, 1, 2) *
      mean(sign(stats::median(diffs)) != sign(diffs))

    ci = as.numeric(stats::quantile(diffs,probs=limits))

  } else {

    # Null Distribution
    # The difference distribution centered on 0
    # Using the CDF, calculate P(r2diff > observed)
    # multiply by 2 if test is two sided

    pval <- ifelse(onesided, 1, 2) *
      (1 - pbetas(z = abs(diff),
                  a = xx$v1/2, b = xx$v2/2,
                  c = xx$v1/2, d = xx$v2/2,
                  ncp1 = xx$ncp, ncp2 = xx$ncp))

    r2_diff_ci = function(start){

      lower = upper = start

      # Upper confidence limit

      repeat{

        qq = pbetas(z = upper,
                    a = xx$v1/2, b = xx$v2/2,
                    c = yy$v1/2, d = yy$v2/2,
                    ncp1 = xx$ncp, ncp2 = yy$ncp)

        if (qq >= limits[2] | abs(qq-limits[2])<0.01) break

        upper = upper+0.005

      }

      # Lower confidence limit

      repeat{

        qq = pbetas(z = lower,
                    a = xx$v1/2, b = xx$v2/2,
                    c = yy$v1/2, d = yy$v2/2,
                    ncp1 = xx$ncp, ncp2 = yy$ncp)

        if (qq <= limits[1] | abs(qq-limits[1])<0.01) break

        lower = lower-0.005

      }

      return(c(lower, upper))

    }

    ci = r2_diff_ci(start = diff)

  }

  res = list(d=diff, ci=ci,p=pval)

  if (fancy == T){

    make.ci = function(x, upper=NULL, lower=NULL, dig = 2){

      # make.ci gives nice confidence intervals

      fr = function(x, dig=2){

        # fr formats and rounds numbers.

        res=trimws(paste(format(round(x, dig), nsmall=dig)))
        return(res)
      }

      c1 = is.null(upper)
      c2 = is.null(lower)
      c3 = length(x)==3

      parens <- function(left, right = NULL){

        if(is.null(right)){

          paste0('(',fr(left),")")

        } else {

          paste0('(',fr(left),' ',fr(right),')')

        }

      }

      if(c1&c2&c3){

        ci = paste(fr(x[1]), parens(left=x[2], right = x[3]))

      } else {

        ci = paste(fr(x), parens(left=lower, right=upper))

      }

      return(ci)

    }

    res$p = ifelse(res$p < 0.001,
                   '<0.001',
                   sprintf("%.3f",round(pval,3)))
    res$ci = make.ci(diff, upper = res$ci[2], lower=res$ci[1])

  }

  return(res)

}




