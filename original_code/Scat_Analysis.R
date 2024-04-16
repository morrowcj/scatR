
library(reshape2)
library(dplyr)
library(tidyr)
library(purrr)
library(utils)
library(ggplot2)
library(jagsUI)

# Individuals (ind)
# Sampling occasions (surv)
# Proportion of Item in Diet (prop)
# Detection probability of item (det)
# Still need to include missing data somehow?

simhist <- function(ind, surv, prop, det)
{
  hist <- matrix(NA, ind, surv)
  z <- rbinom(ind, 1, prop)
  for (i in 1:ind){
    hist[i,]<-rbinom(surv, 1, det*z[i])
  }
  return(hist)
}

sim.1 <- simhist(ind = 10000, surv = 5, prop = .5, det = .5)

colMeans(sim.1)
rows <- rowSums(sim.1)
sum(replace(rows, rows >0, 1))/nrow(sim.1) # individuals that consumed item

# Traditional Frequency of Occurrence

sum(sim.1)/(nrow(sim.1)*ncol(sim.1)) # frequency of item occurance

# Bootstrap estimate of all scats
sim.1.df <- data.frame(sim.1)
sim.1.melt <- melt(sim.1.df)
sim.boot <- data.frame(matrix(sample(sim.1.melt$value,
                                       size = nrow(sim.1.melt)*1000,
                                       replace = T),
                                nrow = nrow(sim.1.melt)))
mean(colMeans(sim.boot))
sd(colMeans(sim.boot))

# Bootstrap estimate using separate surveys and a binomial distribution

sim.1.df$ind <- c(1:nrow(sim.1.df))

sim.data.melt <- melt(sim.1.df, id.var = "ind")

cluster_boot_function <- function(x){

clusted_boot <- sim.data.melt %>%
    group_by(ind) %>%
    nest() %>%
    mutate(samps = map(data, ~ sample(.$value, size = 1, replace = T))) %>%
    dplyr::select(ind, samps) %>%
    unnest(cols = samps)

results <- clusted_boot %>%
    group_by(ind)

results
}

out.binom <- function(x){

  out <- purrr::map_df(.x = 1, cluster_boot_function, .id = "iteration")
  y = mean(rbinom(n = 100, size = 1, prob = mean(out$samps)))
  return(y)

}

binom.data <- data.frame(sapply(rep(1, 100), out.binom))
names(binom.data) <- "binom.estimate"

### Plot it
hist(binom.data$binom.estimate, xlab = "Proportion in Diet")

### Mean and SD
mean(binom.data$binom.estimate)
sd(binom.data$binom.estimate)




########################
## Scat model in JAGS ##
########################

### (1) bundle and summarize data set ###
jags.data <- list( y = sim.1.df[,1:5],                              # detection/non-detection data
                   M = nrow(sim.1.df), J = ncol(sim.1.df[, 1:5])   # create objects describing data dimensions
)

str(jags.data) # look at list summary

### Step 2: set initial values for MCMC procedure ###
# Initial values: must give for same quantities as priors given!
zst <- apply(sim.1.df[, 1:5], 1, max, na.rm=TRUE)  # Avoid data/model/inits conflict, also need na.rm=TRUE if data contains NAs

inits <- function(){
  list(
    z = zst
  )
}

### Step 3: tell JAGS what parameters you want it to "monitor" at each MCMC iteration ###
# Parameters monitored
params <- c("detection", "proportion",
             "N.Ind", "Prop.Ind") # note these are specified in the JAGS model file

### Step 4: tell JAGS specifications for MCMC analysis ###
# MCMC settings
n.interations <- 10000 # number of iterations
thin.rate <- 10    # thin rate
burn.in <- 2000  # burn-in length
chains <- 3     # number of chains

# Call JAGS from R
jags.output <- jags(data = jags.data,              # Step 1
                inits = inits,                     # Step 2
                parameters.to.save = params,       # Step 3
                n.thin = thin.rate,                # Step 4
                n.chains= chains,                  # ''
                n.burnin = burn.in,                # ''
                n.iter = n.interations,            # ''
                model.file = "scat_JAGS.R") # tell JAGS the model file

print(jags.output, dig = 3)


# model diagnostics
summary(jags.output)
traceplot(jags.output)

# Estimates
hist(plogis(jags.output$sims.list$detection), xlab = "Detection", main = NULL)
hist(plogis(jags.output$sims.list$proportion), xlab = "Proportion", main = NULL)

logit_detect = plogis(jags.output$sims.list$detection)
logit_prop = plogis(jags.output$sims.list$proportion)

hist(logit_detect * logit_prop)


