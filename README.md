## MRRk


## Estimation of Insect Dispersal Parameters via Diffusion and Survival Models
**Authors**: Dan Pagendam





### Package installation

To install the package from GitHub, you will first need to install the devtools package in R using the command:

```install.packages("devtools")```

Once installed, you will need to load the devtools R package and install the GEIC R package using:

```
library(devtools)
install_github("dpagendam/GEIC")
```

### Using this package

To use MRRk with the packaged example data, try:

```
library(MRRk)
data(exampleData)
params = rep(NA, 2)
names(params) = c("K", "sigma")
params[1] = 1000
params[2] = 100
mle <- optim(params, trapLogLikelihood, method = "L-BFGS-B", lower = c(1, 10), upper = c(1E20, 1000), control = list(fnscale = -1),  trapLocations = list(trapLocations), releaseLocations = list(releaseLocations), trapData = list(trapData), distributionType = "poisson", includeDrift = FALSE, includeDeath = FALSE)
```
