
#' @export

r2beta.lme <- function(model, partial=TRUE, method='sgv',
                       data = NULL){

  if(is.null(data)) data = model$data

  # Get model matrices
  X=stats::model.matrix(stats::formula(model), data = data)

  # Get number of observations
  n <- nrow(X)

  # Get grouping information from the model
  clust.id = names(summary(model)$groups)[1]
  obsperclust = as.numeric(table(data[,clust.id]))
  mobs = mean(obsperclust)
  nclusts = length(obsperclust)

  # The Kenward Roger Approach can't be applied to these models
  if (toupper(method) == 'KR'){
    stop('The Kenward Roger approach is only compatible with lmerMod objects.')
  }

  if(toupper(method) == 'SGV' | toupper(method) == 'NSJ'){

    # Get fixed effects
    beta = nlme::fixef(model)
    p <- length(beta)

    if(p==1) stop('Model must have at least one fixed effect')

    # Get covariance matrix from the model

    mlist = mgcv::extract.lme.cov2(model, data, start.level=1)[['V']]

    if(toupper(method)=='NSJ'){

      # NS approach takes the mean of the trace of SigHat
      SigHat = mean(diag(as.matrix(Matrix::bdiag(mlist))))

    }

    # SGV approach takes the standardized generalized variance of SigHat
    if(toupper(method)=='SGV'){

      SigHat = calc_sgv(nblocks = nclusts,
                        blksizes = obsperclust,
                        vmat = mlist)

    }

    # C matrix defines the Wald Test for Fixed Effects
    C = list()
    assign <- attr(model$fixDF, "assign")
    nTerms <- length(assign)
    nms = c('Model', names(assign)[-1])

    # Define the model Wald statistic for all fixed effects
    C[['Model']] = cbind(rep(0, p-1), diag(p-1))

    # For partial R2 statistics:
    if (partial == T){

      # add the partial contrast matrices to C
      for(i in 2:(nTerms)) {
        C[[nms[i]]] = make.partial.C(rows=p-1, cols = p, index = assign[[i]])
      }

    }

    # Compute the specified R2
    r2=lapply(C, FUN=cmp_R2, x=X, SigHat=SigHat, beta=beta,
              method=method, obsperclust=obsperclust, nclusts=nclusts)

    # initialize a dataframe to hold results
    R2 = data.frame(Effect = names(r2))

    # place results in the dataframe
    for(i in names(r2[[1]])){
      R2[,i] = as.vector(unlist(lapply(r2, function(x) x[i])))
    }
  }

  R2 = within(R2, {
    lower.CL = stats::qbeta(0.025, R2$v1/2, R2$v2/2, R2$ncp)
    upper.CL = stats::qbeta(0.975, R2$v1/2, R2$v2/2, R2$ncp)
  } )

  R2 = R2[order(-R2$Rsq),]

  class(R2) <- c('R2', 'data.frame')

  return(R2)

}
