# This function computes the log-likelihood of the observed trappind data given a set of parameters, release locations releasenumbers and trap locations
# Trap locations, release locations and trapData are all lists that contain an element for each 
# mark-release-recapture experiment that was performed (i.e. independent experiments to use for parameter estimation).
# Each of these lists contains the following object types:

# trapLocations is a list containing matrices or dataframes whose rows contain coordinates of individual traps.  The matrix must have three columns: 
# the first column should be a unique trapID (either numeric or character) for each trap
# the second column should contain the Easting of the trap
# the third column should contain the Northing of the trap.

# releaseLocations is a list containing matrices or dataframes whose rows contain coordinates and times for the releases used.
# Because individual releases can sometimes use different types of marking, the releaseLocations matrix should contain
# five columns.  
# The first column should be a markerID which can either be numeric IDs like 1,2 and 3 or a character such as "red", "blue", "green"
# The second column should contain the Easting of the release locations.
# The third column should contain the Northing of the release locations.
# The fourth column should contain the number of days since the beginning of the MRR study at which the releases occured.
# Typically the release event(s) at the beginning of the experiment will be have a zero for in the fourth column
# (i.e. releases occured at the start of the study).
# The fifth column should contain the numbers of individuals released at each release event.

# trapData is a list containing matrices or dataframes that pertain to the numbers of individuals caught in traps, the markings
# on these individuals and the time intervals in which these individuals were caught.
# the matrices/dataframes in trapData should be 5 columns.
# The first column should contain the trapID at which the trap counts were collected
# The second column should contain the number of insects caught in the trap for a particular marking type
# The third column should contain the markerID that identifies whether there are different markers (e.g. colours) used
# The fourth column should contain the time (days since the start of the study) at which this trap started collecting the individuals in this trap count record
# The fifth column should contain the time (days since the start of the study) at which this trap stopped collecting the individuals in this trap count record

# The parameter distributionType can be one of: "poisson", "zero_inflated_poisson", "negative_binomial" or "zero_inflated_negative_binomial"
# These are the models that are used to model the actual count data, where the means of these distributions are derived from the integral of the
# time-varying Gaussian kernel.

# includeDrift is a boolean that allows you to include a constant drift applied to the dispersal kernel over time.  If TRUE, then
# the drift parameters (mu_Easting and muNorthing) will be estimated as part of the parameter vector.

# includeDeath is a boolean that allows you to apply a constant death rate to the individuals in your study.  This creates the 
# realistic situation where the numbers of individuals caught in traps decreases over time.  The use of a constant death rate corresponds to 
# the assumption of an exponentially distributed lifetime distribution on the individuals in the study.

# params must be a named numeric vector containing the parameter values at which to compute the
# log-likelihood.  At a minimum, the vector should contain parameters named "sigma" and "K".  If includeDeath
# is set to TRUE, then it should also contain a parameter named "deathRate".  If includeDrift is TRUE
# then params should also contain names parameters "muX" and "muY".
# If distributionType is equal to "zero_inflated_poisson" or "zero_inflated_negative_binomial" then params should
# also include a value for a parameter named pZero that is greater that is in the interval [0, 1].
# If distributionType is equal to "zero_inflated_negative_binomial" of "negative_binomial" then params
# should also include a parameter named "phi".

trapLogLikelihood = function(params, trapLocations, releaseLocations, trapData, distributionType = "poisson", includeDrift = FALSE, includeDeath = TRUE)
{
  #extract the model parameters from the params vector
  muX = params["muX"]
  muY = params["muY"]
  sigma = params["sigma"]
  deathRate = params["deathRate"]
  K = params["K"]
  pZero = params["pZero"]
  phi = params["phi"]
  
  if(!includeDrift)
  {
    muX = 0
    muY = 0
  }
  if(!includeDeath)
  {
    deathRate = 0
  }
  
  if(is.na(muX) | is.na(muY) | is.na(sigma) | is.na(deathRate) | is.na(K))
  {
    stop("One or more of your parameters was NA.  This can occur if your parameter vector is not correctly named.")
  }
  
  if(is.infinite(muX) | is.infinite(muY) | is.infinite(sigma) | is.infinite(deathRate) | is.infinite(K))
  {
    stop("One or more of your parameters was infinite. Please ensure parameter values are finite.")
  }
  
  if(is.nan(muX) | is.nan(muY) | is.nan(sigma) | is.nan(deathRate) | is.nan(K))
  {
    stop("One or more of your parameters was NaN. Please ensure parameter values are properly defined.")
  }
  
  if(sigma <= 0)
  {
    stop("sigma must be positive.")
  }
  
  if(deathRate < 0)
  {
    stop("deathRate must be non-negative")
  }
  
  if(K <= 0)
  {
    stop("K must be positive.")
  }
  
  if((distributionType %in% c("zero_inflated_negative_binomial", "negative_binomial")) & (phi <= 0))
  {
    stop("Phi must be positive")
  }
  
  if((distributionType %in% c("zero_inflated_negative_binomial", "zero_inflated_poisson")) & (pZero < 0 | pZero > 1))
  {
    stop("Phi must be positive")
  }
  
  test = c(class(trapLocations), class(releaseLocations), class(trapData))
  if(!all(test == "list"))
  {
    stop("One or more of the inputs named trapLocations, releaseLocations or trapData is not a list.")
  }
  
  if((length(trapLocations) != length(releaseLocations)) | (length(trapLocations) != length(trapData)) | (length(releaseLocations) != length(trapData)))
  {
    stop("One or more of the lists named trapLocations, releaseLocations or trapData is not the same length as the others.")
  }
  
    
  for(l in 1:length(trapData))
  {
    ll = 0
    thisTrapData = trapData[[l]]
    thisTrapLocations = trapLocations[[l]]
    thisReleaseLocations = releaseLocations[[l]]
    
    test = all(thisTrapData[, 1] %in% thisTrapLocations[, 1])
    if(!test)
    {
      stop("One or more of the trapIDs specified in trapData is not present in trapLocations.")
    }
    
    test = all(thisTrapData[, 3] %in% thisReleaseLocations[, 1])
    if(!test)
    {
      stop("One or more of the markerIDs specified in the trapData is not present in trapLocations.")
    }
    
    
    # set the start of the study to time zero and set all times to be relative to this
    minReleaseTime <- min(thisReleaseLocations[, 4])
    thisReleaseLocations[, 4] = thisReleaseLocations[, 4] - minReleaseTime
    thisTrapData[, 4] <- thisTrapData[, 4] - minReleaseTime
    thisTrapData[, 5] <- thisTrapData[, 5] - minReleaseTime
    
    for(i in 1:nrow(thisTrapData))
    {
      startTime <- thisTrapData[i, 4]
      endTime <- thisTrapData[i, 5]
      ind <- which(thisTrapLocations[, 1] == thisTrapData[i, 1])
      trapLocation <- matrix(c(thisTrapLocations[ind, 2], thisTrapLocations[ind, 3]), 1, 2)
      marker <- thisTrapData[i, 3]
      thisCount <- thisTrapData[i, 2]
      
      #which releases are relevant
      ind <- which(thisReleaseLocations[, 4] < endTime & thisReleaseLocations[, 1] == marker)
      if(length(ind) > 0)
      {
        releaseXY <- matrix(thisReleaseLocations[ind, 2:3], 1, 2)
        releaseTimes <- thisReleaseLocations[ind, 4]
        releaseNumbers <- thisReleaseLocations[ind, 5]
        Lambda <- densityIntegral(startTime, endTime, trapLocation, releaseXY, releaseTimes, releaseNumbers, muX, muY, sigma, deathRate)
        if(distributionType == "poisson")
        {
          ll = ll + dpois(thisCount, Lambda*K, log = TRUE)
        }
        else if(distributionType == "zero_inflated_poisson")
        {
          if(trapCounts[i] == 0)
          {
            ll = ll + log(pZero + (1 - pZero)*dpois(0, Lambda*K))
          }
          else
          {
            ll = ll + log(1 - pZero) + dpois(thisCount, Lambda*K, log = TRUE)
          }
        }
        else if(distributionType == "negative_binomial")
        {
          ll = ll + dnbinom(trapCounts[i], mu = Lambda*K, size = phi, log = TRUE)
        }
        else if(distributionType == "zero_inflated_negative_binomial")
        {
          if(trapCounts[i] == 0)
          {
            ll = ll + log(pZero + (1 - pZero)*dnbinom(thisCount, mu = Lambda*K, size = phi))
          }
          else
          {
            ll = ll + log(1 - pZero) + dnbinom(thisCount, mu = Lambda*K, size = phi, log = TRUE)
          }
        }
      }
    }
  }
  
  return(as.numeric(ll))
}


