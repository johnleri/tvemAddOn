#' Pool invariant coefficient estimates
#'
#'@description
#' Accepts a list of arbitrary length, for which each list element contains a tvem model fit on a multiply imputed dataset. Returns pooled statistics.
#'
#'
#'@param model.fit.list List composed of elements which are 'tvem' objects. Each 'tvem' object should be a model fit using the tvem package on one iteration of a multiply imputed dataset.
#'@param effect.path A string that denotes the subelement of each 'tvem' object where the invariant effects are located. Defaults to 'invar_effects_estimates', which is the default string used by 'tvem' objects.
#'@param return.term A string that determines if t-statistics or beta coefficients and standard errors are returned. Accepts 'tstat' or 'beta' as inputs.
#'
#'@returns Returns a dataframe composed of either t-statistics or beta coefficients and standard errors.
#'
#'
#'@export
pool.invariate.estimates = function(model.fit.list, effect.path = "invar_effects_estimates", return.term = "tstat"){

  # extract the time-invariant beta coefficients from the models fit on multiply imputed data sets
  mi.estimates = t(sapply(model.fit.list, function(x)
    x[[effect.path]][,"estimate"]))

  # extract the time-invariant standard errors of the beta coefficients from the models fit on multiply imputed data sets
  mi.se = t(sapply(model.fit.list, function(x)
    x[[effect.path]][,"standard_error"]))

  # create a pooled estimate of the beta coefficients fit on the multiply imputed data sets
  mi.pool = lapply(1:ncol(mi.estimates), function(x)
    mice::pool.scalar(mi.estimates[,x], (mi.se[,x]^2), n = Inf))

  pooled.estimate = dplyr::bind_rows(lapply(mi.pool, function(x) x[c("qbar", "t")]))

  if(return.term == "tstat"){
    pooled.estimate = pooled.estimate$qbar / sqrt(pooled.estimate$t)
  }
  if(return.term == "beta"){
    pooled.estimate$t = sqrt(pooled.estimate$t)
  }

  return(pooled.estimate)
}
