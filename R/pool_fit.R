#' Extract median fit statistics from across multiple imputations

median_fit = function(model.fit.list){

  aic = median(sapply(model.fit.list, function(x)
    x[["back_end_model"]][["aic"]]
    )
  )

  bic = median(sapply(model.fit.list, function(x)
    x[["model_information"]][["pseudo_bic"]]
  )
  )

  deviance = median(sapply(model.fit.list, function(x)
    x[["back_end_model"]][["deviance"]]
  )
  )

  output = list(deviance, aic, bic); names(output) = c("deviance", "aic", "bic")

  return(output)
}
