# Test parameters
n = 10
k = 1
Sigma <- matrix(data = c(1.0, 0.5, 0.2,
                         0.5, 1.0, -0.7,
                         0.2, -0.7, 1.0),
                ncol = 3, dimnames = list(c("A", "B", "C"), c("A", "B", "C")))

#' Repeat a vector as rows of a matrix
#'
#' @param n the number of rows in the matrix
#' @param v the vector to be replicated for each row
#'
#' @return a matrix of the replcated vector
#'
#' @examples
#' replicate_to_matrix(5, c(1, 2, 3))
replicate_to_matrix <- function(n, v){
  matrix(replicate(n, v), nrow = n, byrow = TRUE)
}

#' Simulate a baseline population diet
#'
#' @param size integer, the size of the population
#' @param means vector, the mean consumption value for each food item
#' @param varcov matrix, the variance-covariance matrix defining consumptive relationships among food items
#'
#' @description
#' Simulates a multivariate normal distribution representing consumption of food items within a population
#'
#' @return a baseline population consumption matrix
sim_baseline_pop_diet <- function(size, means, varcov){
  stopifnot(all.equal(size, length(means), ncol(varcov), nrow(varcov)))
  diet <- mvtnorm::rmvnorm(n = size, mean = means, sigma = varcov)
  colnames(diet) <- names(means)
}

#' Create a null preference factor matrix
#'
#' @param n integer, population size
#' @param items integer, number of diet items
#'
#' @return a matrix of ones
get_null_pop_strat <- function(n = 100, items = 3){
  matrix(1, ncol = items, nrow = n)
}

#' Create a preference factor matrix
#'
#' @param n integer, population size
#' @param strat_table a table of diet strategies. Each column represents a strategy and each row represents a food item.
#' Values represent multipliers to the baseline consumption values.
#' @param strat_props proportion of the population that belongs to each strategy. Any remainders will use the population
#' default (i.e., factor of 1 for all items)
#'
#' @return a preference factor matrix
#'
#' @examples
#' # Example for 3 diet items
#' strat_tab = data.frame(avoider = c(0, 1, 1), # avoids first item
#'                        specialist = c(.2, .2, 1.2) # specializes in third item, reduced consumption of others.
#'                        )
#'
#' define_ind_preferences(10, strat_tab, strat_props = c(.2, .1))
#'
define_ind_preferences <- function(n, strat_table, strat_props){
  stopifnot(sum(strat_props) <= 1)
  stopifnot(length(strat_props) == ncol(strat_table))

  n_strat = length(strat_props)
  n_food = nrow(strat_table)

  # create the null matrix
  pref_factors <- get_null_pop_strat(n, n_food)

  ids <- seq(1, n)
  for (i in seq_len(n_strat)){
    count = round(n * strat_props[i])
    samp <- sample(ids, count)
    pref_factors[samp, ] <- replicate_to_matrix(count, strat_table[,i])
    # update available IDs
    ids <- ids[!ids %in% samp]
  }
  pref_factors
}











sim_population_diet <- function(pop_size = 1000, time_length = 1, diet_means = c(0.6, 0.4, 0.2),
                                diet_probs = NA, diet_varcov = NA, diet_trends = NA, ind_eff_sd = NA,
                                sample_n = 20, AR_burnin = 2000, detect_threshold = NA
                                ){
  # count the food items
  p_food = length(diet_means)

  # get the names of food items
  food_names = names(diet_means)

  # make diet_probs if not given
  if (is.null(diet_probs) | all(is.na(diet_probs))){
    diet_probs = rep(1, p_food)
  } else if (length(diet_probs == 1)){
    diet_probs = rep(diet_probs, p_food)
  }

  # make diet_trends if not given
  if (is.null(diet_trends) | all(is.na(diet_trends))){
    diet_trends = rep(0, p_food)
  } else if (length(diet_trends == 1)){
    diet_trends = rep(diet_trends, p_food)
  }

  # make varcov if not provided
  if (is.null(diet_varcov)){
    diet_varcov = matrix(0, nrow = length(diet_means), ncol = length(diet_means))
    diag(diet_varcov) <- 1
  }

  # check that dimensions all match
  stopifnot(length(diet_probs) == p_food)
  stopifnot(length(diet_trends) == p_food)
  stopifnot(all.equal(ncol(diet_varcov), nrow(diet_varcov), p_food))

  # add names to dimensions
  names(diet_probs) <- names(diet_trends) <- colnames(diet_varcov) <- rownames(diet_varcov) <- food_names


  # simulate the diet starting conditions
  diet_start <- mvtnorm::rmvnorm(n = pop_size, mean = diet_means, sigma = diet_varcov)

  # simulate diet through time
  if (time_length > 1) {

  } else {
    # TODO
  }

  # Calculate the individual random effects, which remain constant through time.
  if (is.na(ind_eff_sd)){
    individual_randefs <- 0
  } else {
    individual_randefs <- matrix(rnorm(n = pop_size*p_food, sd = ind_eff_sd), nrow = pop_size,
                                 dimnames = list(rownames(diet_start), food_names))
  }



}


#' Title
#'
#' @param n integer, number of individuals in the population
#' @param k integer, number of time points to simulate the population for
#' @param items numeric vector
#'
#' @return
#' @export
#'
#' @examples
sim_binary_consumption <- function(n = 10, k = 1, probs = c("berries", "meat", "fish")){

}

probmat = t(replicate(n, c("berries" = .6, "meat" = .4, "fish" = .2)))

#' Title
#'
#' @param probs
#' @param k
#'
#' @return
#' @export
#'
#' @examples
make_binary_cons <- function(probs = probmat, k = 1){
  apply(probs, 2, function(x)rmultinom(n = k, size = 1, prob = x))
}

