# ---- Main consumption funciton ----

#' Create diet samples, based on consumption probabilities
#'
#' @param available vector of available values for each food item, from 
#' 0 to 1. This can be interpreted as the probability of encountering the item.
#' @param preference vector of preferences for each food item, from 0 to 1. 
#' This can be interpreted as the probability of trying to consume the item
#' given an encounter.
#' @param success vector of success rates for consuming each food item, 
#' from 0 to 1. This can be interpreted as the probability of consuming the
#' item given an encounter and an attempt to consume.
#' @param samples number of samples to take, passed as the \code{n} argument
#' of \code{rmultinom()} 
#' @param target vector of total possible encounters for each sample, passed
#' as the \code{size} argument of \code{rmultinom()}. 
#' @param random logical, should the outcome be stochastic, through use of 
#' \code{rmultinom()}.
#' @param prop logical, should the outcome be converted be converted to a
#' proportion of \code{target}?
#' @param items optional vector of names of each item. The default is
#' \code{seq_len(k)} where \code{k} is the number of items. Alternatively,
#' if \code{available} is a named vector (and \code{items} is not provided), 
#' those names will be used.
#' @param fill logical, should information for a compensatory resource be 
#' included? If \code{FALSE}, then the given items will be forced to comprise 
#' the entirety of the diet (the value \code{p = available*prefernce*success} 
#' is scaled to sum to 1). If \code{TRUE}, then an additional group
#' (named by \code{.other}) is created with a probability of (\code{1 - p}).
#' As such, this group is filled by failures to consume any of the specified
#' items.
#' @param .other a string indicating the name to give the additional group. 
#' Defaults to ".other".
#'
#' @seealso [stats::rmultinom()]
#'
#' @returns a consumption matrix, with rows corresponding to food items and
#' columns corresponding to samples. 
#' 
#' @details
#' The expected value for the consumption rate of a resource is
#' \code{p} = \code{available} \times \code{preference} \times 
#' \code{success}. When \code{random = TRUE}, \code{p} is the probability 
#' of a an item being consumed during a single consumption event. Otherwise it 
#' is a fixed proportion of the diet.
#' 
#' When \code{preference = 1} and \code{success = 1} 
#' (the default), then \code{available} is the expected item-specific 
#' consumption proportion.
#' 
#' @export
#'
#' @examples
#' # binomial analog: was our item of interest consumed?
#' consume_diet(.8, samples = 10) # .other=1 means item of interest not consumed
#' 
#' # binary choice between two items
#' consume_diet(c("veggies" = 1, "fruit" = 1), samples = 10)
#' 
#' #' # binary choice between two items, with preference for one
#' consume_diet(
#'   available = c("veggies" = 1, "fruit" = 1), preference = c(0.1, 0.9),
#'   samples = 10
#' )
#' 
#' # allow for "nothing" to be consumed, named items
#' consume_diet(
#'   c(0.5, 0.2), samples = 10, items = c("A", "B"), .other = "nothing"
#' )
#' 
#' # consumed amounts of 3 items, from a total of 100
#' consume_diet(c(0.5, 0.4, 0.1), target = 100)
#' 
#' # consumed prop of 3 items
#' consume_diet(c(0.5, 0.4, 0.1), target = 100, prop = TRUE)
#' 
#' # average diet composition from 1000 samples
#' consume_diet(
#'   c(0.5, 0.4, 0.1), target = 50, samples = 1000, prop = TRUE
#' ) |> rowMeans()
#' 
#' # fixed amounts
#' consume_diet(c(0.5, 0.2), target = 100, random = FALSE)
#' 
#' # fixed props
#' consume_diet(
#'   c(0.4, 0.1, 0.1), target = 100, samples = 10,
#'   random = FALSE, prop = TRUE
#' )
#' 
#' # fixed props *among only the items of interest* (excluding ".other")
#' consume_diet(
#'   c(b = 0.2, a = 0.1), target = 100, fill = FALSE, 
#'   random = FALSE, prop = TRUE
#' ) # sums to 1
#' 
#' # imperfect consumption success, with rare resource
#' consume_diet(
#'   c(0.5, 0.5, 0.1, 0.1), success = c(1, 0.3, 1, 0.3), target = 1000,
#'   items = c("common_easy", "common_hard", "rare_easy", "rare_hard"),
#'   prop = TRUE, .other = ".nothing"
#' )
#' 
#' # complete abundance (totally random)
#' consume_diet(
#'   c(1, 1, 1, 1), preference = 1, success = 1, samples = 3,
#'   target = 100, prop = TRUE
#' )
consume_diet <- function(
    available, preference = 1, success = 1, samples = 1, target = 1, 
    random = TRUE, prop = FALSE, items = names(available),
    fill = TRUE, .other = ".other"
    ){
  
  # TODO: check inputs
  
  p = available*preference*success
  psum = sum(p)
  
  
  if (missing(items)) {
    if (is.null(names(available))) {
      items <- seq_len(length(p))
    } else {
      items <- names(available)
    }
  }
  
  
  if (fill & psum < 1) {
    p = c(p, 1 - psum)
    items = c(items, .other)
  } else {
    p = p / sum(p) # normalize to sum to 1.
  }
  
  k = length(target)
  if (k == 1) {
    if (random) {
      consumed = rmultinom(n = samples, size = target, prob = p)
    } else {
      consumed = replicate(samples, p * target, simplify = TRUE)
    }
    
    rownames(consumed) <- items
    
    if (prop) {
      return(consumed / target)
    } 
    return(consumed)
  
    
  } else if (k == samples){
    # create empty matrix to fill
    consumed <- matrix(
      NA, nrow = length(p), ncol = k, dimnames = list(items, NULL)
    )
    for (i in seq_len(k)) {
      if (random) {
        consumed[, i] <- rmultinom(n = 1, size = target[i], p)
      } else {
        consumed[, i] <- p*target[i]
      }
    }
    
    if (prop) {
      consumed = t(apply(consumed, 1, function(r){r/target}))
    } 
    return(consumed)
  
    
  } else {
    stop("length of target must equal 1 or the value of samples")
  }
}
