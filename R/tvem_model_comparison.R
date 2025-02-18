#' Compare nested tvem models
#'
#' tvem_model_comparison pools coefficients and test statistics.
#'
#' @details tvem_model_comparison pools coefficients and test statistics across imputations to allow for nested model comparisons and the derviation of permutation p-values.
#'
#' @param full_model A list of length two. The first list element is assumed to be the output from the observed model. The second list element is assumed to be the output from the permuted models. These models should include more parameters than the reduced (AKA - the null) model.
#' @param reduced_model A list of length two. The first list element is assumed to be the output from the observed model. The second list element is assumed to be the output from the permuted models. These models should include fewer parameters than the full model.
#'
#' @returns Returns a list of length two. The first list element will be the pooled invariate and variate statistics with their associated permuted p-values. The second list element will be the multivariate wald statistic for the nested models and its associated permuted p-value.
#'
#' @export
tvem_model_comparison = function(full_model, reduced_model){

  # extract the pooled coefficient values from the observed model
  pooled.coefficients = pool_tvem(full_model[[1]], full_model[[2]])

  # extract the pooled test statistic (wald) from observed data
  observed_teststatistic = pool_wald(full_model[[1]], reduced_model[[1]])
  # extract the pooled test statistic (wald) from the permuted data (i.e., one test statistic for each permuted dataset)
  permuted_teststatistics = as.vector(do.call(rbind, lapply(1:100, function(x) pool_wald(full_model[[2]][[x]], reduced_model[[2]][[x]]))))
  # calculate p-value for observed test-statistic based on permuted distribution
  teststatistic_p = 1 - (sum(observed_teststatistic > permuted_teststatistics) / (length(permuted_teststatistics) + 1))


  # extract the median fit statistics (deviance, aic, bic)
  ## extracting median because there isn't an operable way to pool these statistics
  observed_fit = median_fit(full_model[[1]])
  permuted_fit = lapply(full_model[[2]], median_fit)

  # compare deviance
  #observed_deviance = observed_fit[[1]]
  #deviance_p = abs(sum(observed_fit[[1]] < as.vector(do.call(rbind, lapply(permuted_fit, function(x) x[[1]])))) - length(permuted_fit)) / (length(permuted_fit) + 1)

  # compare AIC
  #observed_aic = observed_fit[[2]]
  #aic_p = abs(sum(observed_fit[[2]] < as.vector(do.call(rbind, lapply(permuted_fit, function(x) x[[2]])))) - length(permuted_fit)) / (length(permuted_fit) + 1)

  # compare BIC
  #observed_bic = observed_fit[[3]]
  #bic_p = abs(sum(observed_fit[[3]] < as.vector(do.call(rbind, lapply(permuted_fit, function(x) x[[3]])))) - length(permuted_fit)) / (length(permuted_fit) + 1)

  output = list(
    "pooled.coefficients" = pooled.coefficients,
    "wald" = list("wald_statistic" = observed_teststatistic, "wald_p" = teststatistic_p)
    #"deviance" = list("deviance_statistic" = observed_deviance, "deviance_p" = deviance_p),
    #"aic" = list("aic_statistic" = observed_aic, "aic_p" = aic_p),
    #"bic" = list("bic_statistic" = observed_bic, "bic_p" = bic_p)
  )

  # note that only the pooled coefficients and multivariate wald tests are returned
  ## the deviance, aic, and bic appear to be not sensitive
  output = output[1:2]

  return(output)

}
