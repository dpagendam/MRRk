#' @title Compute the time integral of the diffusion kernel at a single trap location. 
#' @description \code{densityIntegral} approximates the time integral (via quadrature) of a gaussian pdf that is evolving via a process of diffusion. We model the time-varying standard deviation of the gaussian pdf as sqrt(t)*sigma in both the x and y directions.  We also allow for a time-homogeneous drift in the mean which is governed by parameters muX and muY.  We also allow for the dispersing particles (e.g. insects) to die over time according to a time-homogeneous death rate (deathRate).  The integral is computed over the interval [startTime, endTime] at locations in the rows of the trapLocations matrix.  We assume that particles were released from the locations specified in the rows of the matrix releaseLocations with the numbers released and the times of release specified in the vectors releaseNumbers and releaseTimes respectively (one element for each row of releaseLocations).  
#'
#' @param startTime is a numeric value and is the time at which to start integrating the diffusion kernel.
#' @param endTime is a numeric value and is the time at which to stop integrating the diffusion kernel.
#' @param trapLocation is a matrix or vector containing two elements, the first containing the trap's Easting and the second containing the trap's Northing.
#' @param releaseXY is a matrix with two columns containing the locations of relevant releases.  The first column should contain the Eastings of the release locations and the second column should contain Northings.
#' @param releaseTimes is a vector containing the time (in Days) into the experiment that the releases occurred.  There should be one element of releaseTimes for each row of releaseXY.
#' @param releaseNumbers is a vector containing the number of individuals released at each release location in releaseXY.  There should be one element of releaseNumbers for each row of releaseXY.
#' @param muX is a numeric value that corresponds to the drift in the orientation of the east-west axis that should be applied to the diffusion kernel over time.  muX should be in units of distance per day, where distance is in units of Easting and Northings used to specify trap locations and released locations.
#' @param muY is a numeric value that corresponds to the drift in the orientation of the north-south axis that should be applied to the diffusion kernel over time.  muY should be in units of distance per day, where distance is in units of Easting and Northings used to specify trap locations and released locations.
#' @param sigma is a numeric value that corresponds to the standard deviation of the diffusion kernel in the east-west and north-south directions at one day post release.
#' @param deathRate is a numeric value that corresponds to the decay parameter of an exponential decay model.  This model is used to model the death of individuals over time.  The percent daily mortality rate can be calculated as 100*(1 - exp(-lambda))
#' 
#'
#' @return the function returns a single numeric value equal to the integrals of the diffusion kernels (one from each release location) at the trap location over the time interval [startTime, endTime].
#' @export
#'


# This function approximates the time integral of a gaussian pdf that is evolving via a process of diffusion
# We model the time-varying standard deviation of the gaussian pdf as sqrt(t)*sigma in both the x and y directions
# We also allow for a time-homogeneous drift in the mean which is governed by parameters muX and muY.
# We also allow for the dispersing particles (e.g. insects) to die over time according to a time-homogeneous death rate (deathRate)
# The integral is computed over the interval [startTime, endTime] at locations in the rows of the trapLocations matrix.
# We assume that particles were released from the locations specified in the rows of the matrix releaseLocations
# with the numbers released and the times of release specified in the vectors releaseNumbers and releaseTimes respectively (one element for each row of releaseLocations).

densityIntegral = function(startTime, endTime, trapLocation, releaseXY, releaseTimes,  releaseNumbers, muX, muY, sigma, deathRate)
{
  f <- function(timeT)
  {
    n = length(releaseTimes)
    elapsedTimes <- rep(timeT, n) - releaseTimes
    s <- sigma*sqrt(elapsedTimes)
    xDiffs <- rep(trapLocation[1], n) - releaseXY[, 1] - muX*elapsedTimes
    yDiffs <- rep(trapLocation[2], n) - releaseXY[, 2] - muY*elapsedTimes
    numberAlive <- releaseNumbers*exp(-deathRate*elapsedTimes)
    return(sum(numberAlive*dnorm(xDiffs, 0, s)*dnorm(yDiffs, 0, s)))
  }

  integratedConcentration <- integrate(Vectorize(f), startTime, endTime, subdivisions = 1000)
  integratedConcentration <- integratedConcentration$value
  return(integratedConcentration)
}
