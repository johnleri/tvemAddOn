#' Permute outcome variable
#'
#' Generates permutations.
#'
#' Accepts a tvem model or a list of tvem models, extracts the outcome variable from the first model, and returns the desired number of permuted sets of the outcome variable.
#'
#' @param observed.model.fit A vector of values to be imputed.
#' @param n.perm An integer that designates the number of permuted datasets to create.
#'
#' @returns A matrix of y*x, where y is the length of the outcome variable vector and x is equal to the number of permutations. Each permuted vector of the output variable is stored as a column in the matrix.
#'
#' @export
permute_outcome = function(variable.to.permute, n.perm){

    sample.n = length(variable.to.permute)
    n.perm = n.perm

    sets = permute::shuffleSet(n = sample.n, nset = n.perm)
    permuted.outcome = apply(sets, 1, function(x) variable.to.permute[x])

  return(permuted.outcome)
}

