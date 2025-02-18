#' Pool variant coefficient estimates
#'
#'@description
#' Accepts a list of arbitrary length, for which each list element contains a tvem model fit on a multiply imputed dataset. Returns pooled statistics.
#'
#'
#'@param model.fit.list List composed of elements which are 'tvem' objects. Each 'tvem' object should be a model fit using the tvem package on one iteration of a multiply imputed dataset.
#'@param effect.path A string that denotes the subelement of each 'tvem' object where the variant effects are located. Defaults to 'grid_fitted_coefficients', which is the default string used by 'tvem' objects.
#'@param return.term A string that determines if t-statistics or beta coefficients and standard errors are returned. Accepts 'tstat' or 'beta' as inputs.
#'
#'@returns Returns a list of dataframes (one list element per variant effect) composed of either t-statistics or beta coefficients and standard errors (one value for each point on the time grid).
#'
#'@export
pool.variate.estimates = function(model.fit.list, effect.path = "grid_fitted_coefficients", return.term = "tstat"){

  effects = names(model.fit.list[[1]][[effect.path]])


  mi.estimates = sapply(model.fit.list, function(x)
    sapply(effects, function(y)
      x[[effect.path]][[y]][["estimate"]]
    ))


  mi.se = sapply(model.fit.list, function(x)
    sapply(effects, function(y)
      x[[effect.path]][[y]][["standard_error"]]
    ))


  pooled.estimate = lapply(1:nrow(mi.estimates), function(x)
    mice::pool.scalar(mi.estimates[x,], (mi.se[x,]^2), n = Inf )
  )

  pooled.estimate = dplyr::bind_rows(lapply(pooled.estimate, function(x) x[c("qbar", "t")]))

  # returns the t-statistic (beta / se)
  if(return.term == "tstat"){
    pooled.estimate = pooled.estimate$qbar / sqrt(pooled.estimate$t)

    output = list()
    for(i in 1:length(effects)){
      output[[i]] = pooled.estimate[(i*100-99) : (i*100)]
    }

    names(output) = effects
  }

  # returns the pooled beta value and standard error
  if(return.term == "beta"){
    pooled.estimate$t = sqrt(pooled.estimate$t)

    output = list()
    for(i in 1:length(effects)){
      output[[i]] = pooled.estimate[(i*100-99) : (i*100),]
    }

    names(output) = effects
  }

  # returns a list with a length equal to the number of variate effects
  # each list element will have 100 values; one for each point on the 'time' curve of the tvem
  return(output)



}
