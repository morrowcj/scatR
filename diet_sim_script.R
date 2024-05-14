library(ggplot2)
library(dplyr)

set.seed(888)

# logit functions
logit <- function(p) {log(p / (1 - p))}
inverse_logit <- function(l) {1 / (1 + exp(-l))}

# ---- Initial Parameters ----
# simulation dimensions
pop_size = 10000 # population size # n
time_points = 5 # time points # k

# population consumption parameters
food_items = c("berries", "preyA", "preyB") # p = 3
n_food = length(food_items)
consume_probs = c(0.75, 0.25, 0.10) # probability that each food is consumed by a random individual
consume_means = c(0.5, 0.0, -1.0) # given that an item is consumed, how much was consumed?
names(consume_means) = names(consume_probs) = food_items # add names to the probs and means

# Variance-covariance of food consumption
# Note, if all food items have variance=1, then the varriance-covariance matrix is equivalent to the correlation matrix.
# Otherwise, additional steps are needed to calculate covariance from correlations.
food_varcov = matrix(data = c(1.0,  0.5,  0.2,
                              0.5,  1.0, -0.7,
                              0.2, -0.7,  1.0), # variance-covariance of food consumption
                     ncol = 3, dimnames = list(food_items, food_items))



# ---- Baseline Population Consumption ----
# Generate random population baseline consumption (Multivariate-normal) from the covariance matrix
pop_consumption = mvtnorm::rmvnorm(n=pop_size, mean=consume_means, sigma=food_varcov)

# Estimates of population statistics (these are highly accurate with large pop size, e.g., 2000)
(diet_mean = apply(X = pop_consumption, MARGIN = 2, FUN = function(x){round(mean(x), 1)})) # mean item consumption
(diet_cor = round(cov(pop_consumption), 1)) # covariance of item consumption
(diet_sd = round(sqrt(diag(diet_cor)), 1)) # sd of item consumption

# Visualize the diet distribution
pop_consumption %>%
  reshape2::melt() %>%
  transmute(ID = Var1, food = Var2, value = value) %>%
  ggplot() +
  geom_density(aes(x = value, fill = food, col = food), alpha = .2)

# ---- Baseline Individual Consumption ----
individual_id = seq_len(pop_size) # numeric ID for each individual

# Individual preferences
# For simplicity, 3 preference groups
strategy_table = data.frame(opportunist = c(1, 1, 1), # equally likely to eat all items
                            avoider = c(0, 1, 1), # complete avoidance of first item
                            specialist = c(.2, .2, 1.2) # reduced consumption of items 1 and 2, increased consumption of item 3
)
strat_n = nrow(strategy_table) # number of distinct strategies
strategy_probs = c(.7, .05, .25) # proportion of the population that belongs to each preference group
stopifnot(sum(strategy_probs) == 1) # ensure rates add up to 1

# randomly sample from the strategy index, based on the specified rates
strat_indx = sample(x = seq_len(strat_n), size = pop_size, replace = TRUE, prob = strategy_probs)
# Create preference table, using corresponding strategies
individual_strategies = as.matrix(strategy_table[strat_indx, ])
dimnames(individual_strategies) <- list(individual_id, food_items) # set names

# create baseline preference matrix for each individual, from consumption probabilities
consumption_base = matrix(replicate(n = pop_size, consume_probs), nrow = pop_size, byrow = TRUE)

# full individual preferences
individual_consumption_probs = individual_strategies * consumption_base


# ---- Calculate binary consumption over time ----
binary_array = array(NA, dim = c(pop_size, time_points, n_food), dimnames = list(NULL, NULL, food_items))
for(j in seq_len(pop_size)) for (i in seq_len(n_food)) {
  multinom_timeseries = rmultinom(n = time_points, size = 1, prob = individual_consumption_probs[j, ])
  binary_array[j, ,i] = multinom_timeseries[i, ]
}

# ---- Autoregression ----
AR_array = array(NA, dim = c(pop_size, time_points, n_food), dimnames = list(NULL, NULL, food_items))

# burn in for AR generation
burn_in = 2000

# AR coefficients for each food item, in the population.
food_ARcoefs = list("berries" = c(0.1, 0.05), "preyA" = 0.1, "preyB" = 0.1)
stopifnot(length(food_ARcoefs) == n_food)

# Generate an AR time series for each food type within each individual
for (i in seq_len(n_food)) for (j in seq_len(pop_size)) {
  AR_timeseries = arima.sim(n = time_points, model = list(ar = food_ARcoefs[[i]]), n.start = burn_in)
  AR_array[j, ,i] = AR_timeseries
}

# ---- Individual Quantitative Effects ----
randef_sd = 1
# generate random effect of each individual for each food that will remain consistent over time.
individual_randef = matrix(rnorm(n = pop_size*n_food, sd = randef_sd), nrow = pop_size,
                           dimnames = list(individual_id, food_items))


# ---- Calculate Quantitative consumption over time ----
consumption_array = array(NA, dim = c(pop_size, time_points, n_food), dimnames = list(NULL, NULL, food_items))
for (i in seq_len(n_food)) {
  base_mat = replicate(time_points, pop_consumption[, i]) # population consumption, repeated for each time point
  randef_mat = replicate(time_points, individual_randef[, i]) # individual effect, repeated for each time point
  food_consumed = base_mat + AR_array[, , i] + randef_mat # regression
  consumption_array[, , i] = food_consumed # add value into array
}
# Remove value, if item not consumed by individual
consumption_array[binary_array == 0] = NA

print(round(consumption_array, 2))

## --- Calculate raw statistics ----
# Average population binary consumption across time
cons_prop_avg = apply(consumption_array, 3, function(x)1 - mean(is.na(x)))

(prop_avg_table = data.frame(true = consume_probs, expected = round(apply(individual_consumption_probs, 2, mean), 2),
                             observed = round(cons_prop_avg, 2)))

# Average population quantitative consumption across time
item_pop_avg = apply(consumption_array, 3, function(x)mean(x, na.rm=TRUE))
(pop_avg_table = data.frame(true = food_means, observed = round(item_pop_avg, 2)) %>% mutate(diff = observed - true))

# Population varcovariance across time (if possible)
consumption_permute = aperm(consumption_array, c(1, 3, 2)) # rearrange the array
cov_list = list(true = food_varcov) # create a list of covariances
for (t in seq_len(time_points)){
  name = paste0("time", t)
  cov_list[[name]] = round(cor(consumption_permute[, , t], use = "na.or.complete"), 2) # add to covariance list
}
cov_list
