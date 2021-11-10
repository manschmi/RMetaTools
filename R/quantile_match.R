#' Quantile Matching
#'
#' Selects from list of vectors, items from each vector with same count and quantile distribution
#'
#' @param list_of_values list of sets of values.
#' @param n_quantiles number of quantiles used for matching (default=20 but check details)
#'
#' @details Splits the smallest set values into n_quantiles as defined. For each quantile gets a similar sized sample of values within this quantile from each element in list_of_values.
#'
#' @return list with 2 elements: indices for each element in list_of_values
#'
#' @examples
#'
#' values1 <- rnorm(10000)
#' values2 <- rnorm(100,2,1)
#' values3 <- rnorm(1000,3,2)
#' summary(values1)
#' summary(values2)
#' summary(values3)
#'
#' matched_idx <- match_quantiles(list(values1, values2, values3))
#' summary(values1[matched_idx[[1]]])
#' summary(values2[matched_idx[[2]]])
#' summary(values3[matched_idx[[3]]])
#'
#' @export
match_quantiles <- function(list_of_values, n_quantiles=20) {

  set_lengths <- sapply(list_of_values, length)
  min_set <- which(set_lengths == min(set_lengths))[1]
  qs <- quantile(list_of_values[[min_set]], probs = seq(0,1,1/n_quantiles))

  matched <- lapply(1:(length(qs)-1), function(x) {
    q_sets <- lapply(list_of_values, function(values) which(values > qs[x] & values <= qs[x+1]))
    min_qsize <- min(sapply(q_sets, length))
    q_sets <- lapply(q_sets, function(q_set) sample(q_set, min_qsize))
    q_sets
  })

  res <- lapply(seq_along(list_of_values), function(i) unlist(sapply(matched, function(x) x[[i]])))

  return( res )
}
