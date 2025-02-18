#' Fit tvem model(s)
#'
#'  fit_tvem is a wrapper for the 'tvem' function in the 'tvem' package.
#'
#'  @details The primary utility of fit_tvem is to accept lists of data objects and return a list of fitted tvem models.
#'
#'  @param tvem_outcome_variable A string of length 1 which denotes the name of the dependent variable.
#'  @param tvem_time_variable A sting of length 1 which denotes the name of the variable which represents time (e.g., age or minutes).
#'  @param tvem_id_variable A string of length 1 which denotes the name of the variable which represents each particpant.
#'  @param tvem_variate_effects A string of any length (which is reasonable to model) which denotes the names of the variables which should be modeled as time-varying effects.
#'  @param tvem_invariate_effects A string of any length (which is reasonable to model) which denotes theh names of the variables which should be modeled as time-invariant effects.
#'  @param tvem_data Either a dataframe or a list of dataframes
#'  @param tvem_knots The number of knots to use in the tvem model. For an explination of knots, see the documentation for the 'tvem' package.
#'
#'  @returns Returns the output from the 'tvem' function (which is a list containing information about model fit). This function returns the tvem output itself (if the inputed dataset is a dataframe) or a list of tvem output (if the inputed data was a list of dataframes) fit for each inputed dataframe.
#'
#' @export
fit_tvem = function(
    tvem_outcome_variable,
    tvem_time_variable,
    tvem_id_variable,
    tvem_variate_effects = NULL,
    tvem_invariate_effects,
    tvem_data,
    tvem_knots,
    fit_permuted = FALSE,
    permutation_sets = NULL
){

  #####################################################################################################
  # Ensure that the variable names are compatible with TVEM - return an error if names are not compatible
  variable_names_wperiod = grep("\\.", c(tvem_outcome_variable, tvem_variate_effects, tvem_invariate_effects, tvem_time_variable, tvem_id_variable))

  if(length(variable_names_wperiod) > 0){
    stop("Variable names cannot include periods, switch them to underscores _ ")
  }

  #####################################################################################################
  # Specify formula for invariant effects
  if(is.null(tvem_invariate_effects)){
    tvem_invariate_formula = reformulate(termlabels = "1", response = "SUD")
  }else{
    tvem_invariate_formula = reformulate(termlabels = tvem_invariate_effects, response = NULL)
  }

  #####################################################################################################
  # Specify formula for variate effects
  if(is.null(tvem_variate_effects)){
    tvem_variate_formula = reformulate(termlabels = "1", response = tvem_outcome_variable, intercept = T)
  } else {
    tvem_variate_formula = reformulate(termlabels = tvem_variate_effects, response = tvem_outcome_variable, intercept = T)
  }

  #####################################################################################################
  #If a list of dataframes was entered as the data input - then fit tvem for each dataframe in the list.
  #####################################################################################################

  if(class(tvem_data) == "list"){

    colnames(tvem_data[[1]])[colnames(tvem_data[[1]]) == tvem_time_variable] = "time"
    colnames(tvem_data[[1]])[colnames(tvem_data[[1]]) == tvem_id_variable] = "record_id"
    these.colnames = colnames(tvem_data[[1]])

    tvem_data = lapply(tvem_data, setNames, these.colnames)

    # trim the datasets to just the required columns
    target.colnames = which(these.colnames %in% c("record_id", "time", tvem_invariate_effects, tvem_variate_effects, tvem_outcome_variable))
    tvem_data = lapply(tvem_data, function(x)
      x[,target.colnames])

    # fit observed model
    model = lapply(tvem_data, function(x)
      tvem::tvem(
        data = x,
        formula = tvem_variate_formula,
        invar_effects = tvem_invariate_formula,
        family=binomial(),
        id = record_id,
        time = time,
        num_knots = tvem_knots
      )
    )

    if(fit_permuted == TRUE){

      # extract outcome variable as a vector
      outcome = tvem_data[[1]][,tvem_outcome_variable]
      # permute the outcome variable
      outcome = apply(permutation_sets, 1, function(x) outcome[x])

      ####
      # trim the datasets to just the required columns
      target.colnames = which(colnames(tvem_data[[1]]) %in% c("record_id", "time", tvem_invariate_effects, tvem_variate_effects))
      tvem_data = lapply(tvem_data, function(x)
        x[,target.colnames])


      # define helper function to reduce size of the fitted tvem models
      tvem_trim = function(model.fit){

        model.fit$model_information[-which(names(model.fit$model_information) %in% c("pseudo_aic", "pseudo_bic"))] = NULL
        model.fit$back_end_model[-which(names(model.fit$back_end_model) %in% c("deviance", "aic", "null.deviance", "df.null", "df.residual"))] = NULL

        return(model.fit)
      }


      parallel.tvem = function(outcome_data, mi.datasets = tvem_data, outcome.variable = tvem_outcome_variable){

        mi.length = length(mi.datasets)

        outcome.data = data.frame(outcome_data); colnames(outcome.data) = outcome.variable

        mi.data = lapply(mi.datasets, cbind, outcome.data)


        fit = lapply(mi.data, function(index)
          tvem::tvem(
            data = index,
            formula = tvem_variate_formula,
            invar_effects = tvem_invariate_formula,
            family=binomial(),
            id = record_id,
            time = time,
            num_knots = tvem_knots
          )
        )

        fit = lapply(fit, tvem_trim)

        return(fit)
      }


      # set up the parallel processing
      cl <- parallel::makeCluster(8)
      parallel::clusterExport(cl, c("parallel.tvem", "tvem_trim", "outcome", "tvem_data", "tvem_outcome_variable", "tvem_variate_formula", "tvem_invariate_formula", "tvem_knots"), envir = environment())
      parallel::clusterEvalQ(cl, {
        c(library(tvem))
      })

      # run tvem models in parallel
      models_permuted <- parallel::parApply(
        cl = cl,
        outcome,
        2,
        parallel.tvem
      )

      stopCluster(cl)

      output = list(model, models_permuted); names(output) = c("observed_model", "permuted_models")
      return(output)

    }else{
      return(model)
    }


  }

  #####################################################################################################
  # If a single dataframe was entered as the data input - then fit a single tvem.
  #####################################################################################################

  if(class(tvem_data) == "data.frame"){

    # make sure that the time and ID variables can be recognized by the TVEM function
    colnames(tvem_data)[colnames(tvem_data) == tvem_time_variable] = "time"
    colnames(tvem_data)[colnames(tvem_data) == tvem_id_variable] = "record_id"

    # Fit the time-varying effect model (uses the TVEM package)
    model = tvem::tvem(
      data = tvem_data,
      formula = tvem_variate_formula,
      invar_effects = tvem_invariate_formula,
      family=binomial(),
      id = record_id,
      time = time,
      num_knots = tvem_knots
    )
  }




  return(model)

}
