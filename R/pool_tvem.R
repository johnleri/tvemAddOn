#' Pool coefficients from `tvem` models
#'
#' Utility function used to pool coefficients from `tvem` models fit on multiply imputed data.
#'
#' @details Uses the `pool.variate.estimates` and `pool.invariate.estimates` functions to pool beta values and generate permuted p-values across multiply imputed datasets.
#'
#' @param obs.model List of `tvem` models fit on multiply imputed data.
#' @param perm.model List of `tvem` models fit on permuted data.
#'
#' @returns A list with a length of two. First list element contains the beta and p-values for invariant effects. Second list element contains the beta and p-values for the variate effects.
#'
#' @export
pool_tvem = function(obs.model, perm.model){

  #############################################################################
  # pool and format variate effects

  if(is.null(obs.model[[1]]$grid_fitted_coefficients)){
    obs.variate.output = NULL
  } else {

    pooled.variate = lapply(perm.model, pool.variate.estimates)

    temp.list = list()
    list.names = names(pooled.variate[[1]])
    for(effects in 1:length(pooled.variate[[1]])){
      # formats data so that rows are timepoints and columns are permutations
      loop.data = t(do.call(rbind, lapply(pooled.variate, function(x) x[[effects]])))

      temp.list[[effects]] = loop.data

    }

    pooled.variate = temp.list; names(pooled.variate) = list.names


    obs.variate = pool.variate.estimates(obs.model)
    n.perm = dim(pooled.variate[[1]])[[2]]
    n.time = dim(pooled.variate[[1]])[1]

    obs.variate.p = lapply(1:length(obs.variate), function(x)
      sapply(1:n.time, function(y)
        1 - (sum(abs(obs.variate[[x]][y]) > abs(pooled.variate[[x]][y, ])) / (n.perm + 1))
      )
    )

    obs.variate.output = lapply(1:length(obs.variate), function(x)
      data.frame(obs.variate[[x]], obs.variate.p[[x]])
    )
    obs.variate.output = lapply(obs.variate.output, setNames, c("beta", "p"))
    names(obs.variate.output) = names(obs.model[[1]][["grid_fitted_coefficients"]])

  }

  #############################################################################
  # pool and format invariant effects

  if(is.null(obs.model[[1]]$invar_effects_estimates)){
    obs.invariate_output = NULL
  } else {

    pooled.invariate = t(do.call(rbind, lapply(perm.model, pool.invariate.estimates)))

    # create summary for invariant effects
    obs.invariate = pool.invariate.estimates(obs.model)

    obs.invariate.p = sapply(1:length(obs.invariate), function(x)
      1 - ((sum(abs(obs.invariate[x]) > abs(pooled.invariate[x,]))) / (ncol(pooled.invariate) + 1))
      )

    obs.invariate_output = data.frame(obs.invariate, obs.invariate.p)
    colnames(obs.invariate_output) = c("beta", "p")

  }
  #############################################################################
  # produce function output
  output = list(obs.invariate_output, obs.variate.output); names(output) = c("invariate", "variate")

  return(output)
}
