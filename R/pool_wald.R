#' Pool multivariate Wald test statistic
#'
#' Calculate a pooled multivariate Wald statistic for a set of models fit on multiply imputed data sets.
#'
#'@param model1 The model which has a larger number of variables.
#'@param model0 The model which has a smaller number of variables.
#'
#'@return Multivariate Wald statistic.
#'
#'@references
#' Li, K. H., T. E. Raghunathan, and D. B. Rubin. 1991.
#' Large-Sample Significance Levels from Multiply Imputed Data Using
#' Moment-Based Statistics and an F Reference Distribution.
#' \emph{Journal of the American Statistical Association}, 86(416): 1065–73.
#'
#'
#'
#'@export
pool_wald = function(model1, model0){


  # define function that will pool invariant effect estimates
  pool_invariate = function(model.fit.list, effect.path = "invar_effects_estimates"){
    # extract the time-invariant beta coefficients from the models fit on multiply imputed data sets
    mi.estimates = t(sapply(model.fit.list, function(x)
      x[[effect.path]][,"estimate"]))
    # extract the time-invariant standard errors of the beta coefficients from the models fit on multiply imputed data sets
    mi.se = t(sapply(model.fit.list, function(x)
      x[[effect.path]][,"standard_error"]))
    # create a pooled estimate of the beta coefficients fit on the multiply imputed data sets
    mi.pool = lapply(1:ncol(mi.estimates), function(x)
      mice::pool.scalar(mi.estimates[,x], (mi.se[,x]^2), n = Inf))

    return(mi.pool)
  }

  # define function that will pool variate effect estimates
  pool_variate = function(model.fit.list, effect.path = "grid_fitted_coefficients"){
    effects = names(model.fit.list[[1]][[effect.path]])
    mi.estimates = sapply(model.fit.list, function(x)
      sapply(effects, function(y)
        x[[effect.path]][[y]][["estimate"]]))
    mi.se = sapply(model.fit.list, function(x)
      sapply(effects, function(y)
        x[[effect.path]][[y]][["standard_error"]]))
    mi.pool = lapply(1:nrow(mi.estimates), function(x)
      mice::pool.scalar(mi.estimates[x,], (mi.se[x,]^2), n = Inf))

    return(mi.pool)
  }


  # determine number of imputations
  m = length(model1)


  # extract information about invariant effects (if any were specified)
  ## pools invariant parameter estimates across imputations
  ## extracts the names of the invariant parameters

  # full model
  if(is.null(model1[[1]]$invar_effects_estimates)){
    model1_pooled_invariate = NULL
    model1_invariate = NULL
  }  else {
    model1_pooled_invariate = pool_invariate(model1)
    model1_invariate = 1:length(model1_pooled_invariate)
  }
  # reduced model
  if(is.null(model0[[1]]$invar_effects_estimates)){
    model0_pooled_invariate = NULL
    model0_invariate = NULL
  } else {
    model0_pooled_invariate = pool_invariate(model0)
    model0_invariate = 1:length(model0_pooled_invariate)
  }

  # extract information about the time-varying coefficients (code assumes that at least one was specified)
  ## pool time-varying parameter estimates
  model1_pooled_variate = pool_variate(model1)
  model0_pooled_variate = pool_variate(model0)
  ## extract names time-varying parameters
  model1_variate = names(model1[[1]]$grid_fitted_coefficients)
  model0_variate = names(model0[[1]]$grid_fitted_coefficients)

  # determine which variables distinguish the full from null model
  new_names_variate = which(!c(model1_variate %in% model0_variate))
  new_names_invariate = which(!c(model1_invariate %in% model0_invariate))
  target_variate = c()
  for(i in 1:length(new_names_variate)){
    temp = (((new_names_variate[[i]]*100)-99)+(length(model1_invariate) - length(new_names_invariate))) : ((new_names_variate[[i]]*100) + (length(model1_invariate)) - length(new_names_invariate))
    target_variate = c(target_variate, temp)
  }

  # define Q, which will be used in matrix math to zero out parameter estimates (qbar and ubar) that are common to both the full and reduced model
  dimQ1 = length(model1_pooled_invariate) + length(model1_pooled_variate)
  dimQ2 = length(model0_pooled_invariate) + length(model0_pooled_variate)
  Q = diag(dimQ1)
  names(Q) = c(model1_invariate, rep(model1_variate, each = 100))
  Q = Q[c(new_names_invariate, target_variate), , drop = FALSE]

  # pooled parameter estimates (for coefficients that are in full model but not the reduced model)
  qbar = Q %*% c(
    dplyr::bind_rows(lapply(model1_pooled_invariate, function(x) x[c("qbar")]))$qbar,
    dplyr::bind_rows(lapply(model1_pooled_variate, function(x) x[c("qbar")]))$qbar
    )

  # within imputation variance (for coefficients that are in full model but not the reduced model)
  ubar = Q %*% diag(
    c(
      dplyr::bind_rows(lapply(model1_pooled_invariate, function(x) x[c("ubar")]))$ubar,
      dplyr::bind_rows(lapply(model1_pooled_variate, function(x) x[c("ubar")]))$ubar
      )
    ) %*% (t(Q))

  # between imputation variance
  Bm = Q %*% diag(
    c(
      dplyr::bind_rows(lapply(model1_pooled_invariate, function(x) x[c("b")]))$b,
      dplyr::bind_rows(lapply(model1_pooled_variate, function(x) x[c("b")]))$b
      )
    ) %*% (t(Q))

  # calculate the relative increase in variance due to nonresponse
  rm = (1 + 1 / m) * sum(diag(Bm %*% (solve(ubar)))) / nrow(qbar)

  # calculate multivariate wald statistic
  D1 = (t(qbar) %*% solve((1 + rm)*ubar) %*% qbar) / nrow(qbar)

  # return test statistic as the function output
  ## note that a degrees of freedom and p-value are not returned because this function is set up for permutation testing
  ### empirical degrees of freedom = the number of unique variables in the full model versus the reduced model
  ### empirical p value = (1 - (number of times that the permuted data have larger test statistics than the observed data / # permutations))
  output = as.vector(D1)

  return(output)

}
