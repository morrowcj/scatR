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

# ---- Function to creat an example scat sample table ----

#' Title
#'
#' @param strats a table of probabilities for each item (cols) by group (rows)
#' @param pop_size the size of the total population
#' @param samples the size of the samples drawn from the population
#' @param recap_rate a length-2 vector whose elements are 1) the proportion of 
#' the sample that will be captured a second time and 2) the proportion of
#' recaptures that will be captured a third time.
#' @param target passed as the \code{target} argument in \code{consume_diet()}.
#'
#' @returns a data frame with individual ID, group the individual belongs to,
#' the capture number for the individual, and proportions of each food item
#' consumed.
#' 
#' @examples
example_scat_samples <- function(
    strats, 
    pop_size = 1000, samples = 100, recap_rate = c(0.2, 0.5), target = 100
  ){
  # calculate the number of individuals, recaptures, and re-recaptures
  n_recaps <- round(samples*recap_rate[1])
  n_rerecaps <- round(n_recaps * recap_rate[2])
  n_individuals <- samples - (n_recaps + n_rerecaps)
  
  
  # sample the full population and assign individual IDs
  iid <- sample(samples, size = n_individuals)
  # subsample individuals that will be recaptured
  recaps <- sample(iid, size = n_recaps)
  rerecaps <- sample(recaps, size = n_rerecaps)

  # calculate the number of different groups
  n_groups = nrow(strats)
  
  # assign individuals to the different groups
  groups <- sample(n_groups, n_individuals, replace = TRUE)
  # match groups to recaps
  recap_groups <- groups[match(recaps, iid)]
  rerecap_groups <- groups[match(rerecaps, iid)]
  
  
  df <- tibble(
    iid = c(iid, recaps, rerecaps), 
    group = c(groups, recap_groups, rerecap_groups),
    capture = c(rep(1, n_individuals), rep(2, n_recaps), rep(3, n_rerecaps))
  ) |> arrange(group, iid, capture)
  
  results <- tibble()
  for (i in seq_len(nrow(df))) {
    df_row = df[i, ]
    strat = unlist(strats[df_row$group, ])
    cons <- consume_diet(
      available = strat, target = target, prop = TRUE, .other = "other"
    ) |> t() |> data.frame() |> tibble() 
    
    # stack results
    results <- results |> bind_rows(cons)
  }
    
  observed_average <- colMeans(results, na.rm = TRUE)
  
  if ("other" %in% names(results)) {
    results <- results |> mutate(other = replace_na(other, 0))
  }
  
  # join the results with the ID cols.
  df <- bind_cols(df, results)
  
  # return the results (shuffled rows)
  return(
    df |> arrange(sample(samples))
  )
}

# ---- Example 1 ---- 

set.seed(97)

# create a table of 4 different "strategies" that the population will have
strategies <- tribble(
  ~A, ~B, ~C,
  1/3, 1/3, 1/3,  # equal probability of consuming each resource
  0.5, 0.3, 0.2, # decreasing consumption probability
  0.495, 0.495, 0.01, # one rarely consumed resource
  0.01, 0.05, 0.05 # only rarely consumed items (compensate with an extra)
) |> tibble()

# generate the data
scats_df <- example_scat_samples(strats = strategies, target = 100)

# proportion of scat comprised of a food item.
scats <- scats_df |> select(-c(iid:capture))

# observed average proportion of each food item in scat
colMeans(scats)

# observed averages for each group
scats_df |> group_by(group) |> summarize(across(A:other, ~mean(.x)))
                                         

# indicator matrix for whether the food item was observed in scat
scats_bin <- scats |> mutate(across(everything(), ~as.numeric(.x > 0)))

# observed proportion of scat samples in which the items were detected
colMeans(scats_bin)

# observed proportion of samples containing each item within groups
scats_df |> group_by(group) |> summarize(across(A:other, ~mean(.x > 0)))


## TODO: create a function that applies the scat analysis method on 
## `scats` and/or `scats_bin`. Note that the function should be able to handle
## any number of rows and columns. The next example creates results for 10 food
## items.

# ---- Example 2 ----
# generate consumption probabilities for 10 items
weights = c(2, 1, 1, 1, 0.5, 0.5, 0.2, 0.2, 0.1, 0.01)
probs = weights / sum(weights)
# put them in a table (by row)
prob_tab <- t(tibble(probs))
# name the food items
colnames(prob_tab) <- paste("item", 1:10, sep = ".")

# generate the sample
scat_samples <- example_scat_samples(prob_tab, target = 10)

# remove ID columns
scats2 <- scat_samples |> select(-c(iid:capture))

# average observed proportions (vs expected)
colMeans(scats2) |> 
  unlist() %>% cbind(expected = probs, observed = .) |> 
  round(digits = 3)

# proportions of samples in which eat item was found
scats2 |> summarize(across(everything(), ~mean(.x > 0)))
