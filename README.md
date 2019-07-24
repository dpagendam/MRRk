## MRRk


## Estimation of Insect Dispersal Parameters via Diffusion and Survival Models
**Authors**: Dan Pagendam, CSIRO Data61 (dan.pagendam@csiro.au)
**Contributors**: Brendan Trewin, CSIRO Health & Biosecurity (Brendan.Trewin@csiro.au)

### About
This package was developed as a simple alternative to using temporally-static dispersion kernels for the analysis of Mark-Release-Recapture experiments.  The approach uses a Gaussian dispersion kernel that evolves through time in a manner analogous to isotropic diffusion in 2D space.  The approach also allows for the inclusion of reduced trap catches over time as a result of death and for the specification of different distributional models for the trap counts.  The mean of the trap count distribution is modelled as being proportional to the time-integral of the diffusion kernel at the trap location over the time interval that the trap was catching for.



### Package installation

To install the package from GitHub, you will first need to install the devtools package in R using the command:

```install.packages("devtools")```

Once installed, you will need to load the devtools R package and install the MRRk R package using:

```
library(devtools)
install_github("dpagendam/MRRk")
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

### Need help?
The package currently only contains two functions and help can be accessed for these within R by typing:
```
?trapLogLikelihood
?densityIntegral
```

For further support, please email :-)