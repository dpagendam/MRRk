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
    return(sum(releaseNumbers*dnorm(xDiffs, 0, s)*dnorm(yDiffs, 0, s)))
  }

  integratedConcentration <- integrate(Vectorize(f), startTime, endTime, subdivisions = 1000)
  integratedConcentration <- integratedConcentration$value
  return(integratedConcentration)
}
