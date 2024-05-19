
## Estimate the SD with N=6
# Allen Downey used grid approximation: https://allendowney.github.io/DataQnA/gauss_bayes.html
# We will try and different approach

# We have 6 samples of blood potassium levels. Is there a way
# to estimate from these 6 observations what a range of probable
# values would be for blood potassium if we had 100 or 1000 samples?

## Six observations
d <- c(4.0, 3.9, 4.0, 4.0, 4.7, 4.2)

# place observation info into their own objects
obs_n <- length(d)
obs_df <- obs_n - 1
obs_mu <- mean(d)
obs_sd <- sd(d)

obs_mu
obs_sd

# plot
hist(d, main = "blood samples")

#### Prior Mu
## Need a prior for mu and it's standard error
# The normal rage for Potassium (usually reported as the
# 5th and 95th percentiles) in the population is 3.5 - 5.4.
# We can use this to determine a prior mu

low95 <- 3.5
high95 <- 5.4

prior_mu <- (low95 + high95) / 2
prior_se <- 0.6

# plot prior mu
plot(x = seq(from = 2, to = 7, length.out = 50),
     y = dnorm(x = seq(from = 2, to = 7, length.out = 50), mean = prior_mu, sd = prior_se),
     type = "l",
     main = "Prior Mu Density",
     ylab = "PDF",
     xlab = "Blood Potassium")


## Prior Sigma
# Downey's example uses a gamma distribution so let's plot that
alpha <- 2
beta <- 0.5

plot(x = seq(from = 0.01, to = 5, length.out = 100),
     y = dgamma(x = seq(from = 0.01, to = 5, length.out = 100), shape = alpha, scale = beta),
     type = "l",
     main = "Prior Sigma Density (Gamma Prior)",
     ylab = "PDF",
     xlab = "Sigma for Blood Potassium")


## Calculate a posterior SD using Jeffrey's Prior
# calculate the sum of squared error for each observation relative to the prior_mu
sse <- sum((d - prior_mu)^2)
sse

set.seed(333)
posterior_sigma <- mean(sqrt(sse *  1/rchisq(n = 1000, df = obs_n)))
posterior_sigma

# This value is slightly larger than what we see with the SD for the observations
obs_sd

### Calculate a posterior mean
posterior_mu <- ((1 / prior_se^2 * prior_mu) + (1 / obs_sd^2 * obs_mu)) / ((1/obs_sd^2) + (1/prior_se^2))
posterior_se <- sqrt(1 / (1/obs_sd^2 + 1/prior_se^2))

posterior_mu
posterior_se

## What would we expect the standard deviation to be if we had 100 samples?
set.seed(456)
posterior_mu100 <- rnorm(n = 100, mean = posterior_mu, sd = posterior_se)
posterior_sd100 <- sqrt(sse *  1/rchisq(n = 100, df = obs_n))

df100 <- data.frame(mu = posterior_mu100, sigma = posterior_sd100)

head(df100)

# Create 100 simulations (with replacement) from the mu and sigma values by randomly sampling 
# from df100. Repeat the process 10 times
sim_storage100 <- matrix(data = NA, nrow = 100, ncol = 10)

for(i in 1:10){
  
  row_id <- sample(1:nrow(df100), size = 1, replace = TRUE)
  sim_storage100[, i] <- rnorm(n = 100, mean = df100[row_id, "mu"], sd = df100[row_id, "sigma"])
}

head(sim_storage100)

# plot the ten 1000 simulations
plot(density(sim_storage100[, 1]), 
     col = 'red', 
     xlim = c(1, 7), 
     ylim = c(0, 1.6), 
     lwd = 2,
     main = "100 Simulations repeated 10 times")
lines(density(sim_storage100[, 2]), col = 'green', lwd = 2)
lines(density(sim_storage100[, 3]), col = 'blue', lwd = 2)
lines(density(sim_storage100[, 4]), col = 'brown', lwd = 2)
lines(density(sim_storage100[, 5]), col = 'black', lwd = 2)
lines(density(sim_storage100[, 6]), col = 'orange', lwd = 2)
lines(density(sim_storage100[, 7]), col = 'yellow', lwd = 2)
lines(density(sim_storage100[, 8]), col = 'pink', lwd = 2)
lines(density(sim_storage100[, 9]), col = 'grey', lwd = 2)
lines(density(sim_storage100[, 10]), col = 'palegreen', lwd = 2)

## get mean and SD of each of the columns
apply(X = sim_storage100, MARGIN = 2, FUN = mean)
apply(X = sim_storage100, MARGIN = 2, FUN = sd)

# summarize each simulation with the median ± 90% credible intervals
apply(X = sim_storage100, MARGIN = 2, FUN = quantile, probs = c(0.5, 0.05, 0.95))

# Calculate the median ± 90% CI for the mean of each column
quantile(apply(X = sim_storage100, MARGIN = 2, FUN = mean), probs = c(0.5, 0.05, 0.95))

