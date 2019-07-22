#Example
load("exampleData.RData")

params = rep(NA, 2)
names(params) = c("K", "sigma")
params[1] = 1000
params[2] = 100

# perform maximum likelihood estimation using the example data provided in exampleData.RData
mle <- optim(params, trapLogLikelihood, method = "L-BFGS-B", lower = c(1, 10), upper = c(1E20, 1000), control = list(fnscale = -1),  trapLocations = list(trapLocations), releaseLocations = list(releaseLocations), trapData = list(trapData), distributionType = "poisson", includeDrift = FALSE, includeDeath = FALSE)

