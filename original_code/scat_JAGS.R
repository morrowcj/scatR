
# JAGS model for scat

model {

  ### Priors ###
  detection ~ dnorm(0, 1)  # detection intercept
  proportion ~ dnorm(0, 1)   # psi intercept

  ### Likelihood ###
  ## Ecological model for the partially observed true state, psi
  for (i in 1:M) {

    #psi at site level
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- proportion

    ## Observation model for p
    for (j in 1:J) {
      y[i,j] ~ dbern(z[i] * p[i,j])  # detection-nondetection at i and j
      logit(p[i,j]) <- detection

    } #j
  } #i

  ### Derived quantities ###
  N.Ind <- sum(z[])   # number of M sites occupied
  Prop.Ind <- N.Ind/M   # proportion of M sites occupied

}
