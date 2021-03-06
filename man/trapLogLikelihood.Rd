% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trapLogLikelihood.R
\name{trapLogLikelihood}
\alias{trapLogLikelihood}
\title{Compute the log-likelihood of MRR datasets for given set of parameters.}
\usage{
trapLogLikelihood(params, trapLocations, releaseLocations, trapData,
  distributionType = "poisson", includeDrift = FALSE,
  includeDeath = TRUE)
}
\arguments{
\item{params}{is a named numeric vector containing the parameter values at which to compute the
log-likelihood.  At a minimum, the vector should contain parameters named "sigma" and "K".  If includeDeath
is set to TRUE, then it should also contain a parameter named "deathRate".  If includeDrift is TRUE
then params should also contain names parameters "muX" and "muY".
If distributionType is equal to "zero_inflated_poisson" or "zero_inflated_negative_binomial" then params should
also include a value for a parameter named pZero that is greater that is in the interval [0, 1].
If distributionType is equal to "zero_inflated_negative_binomial" of "negative_binomial" then params
should also include a parameter named "phi".}

\item{trapLocations}{is a list containing matrices or dataframes whose rows contain coordinates of individual traps.  The matrix must have three columns:
(i) the first column should be a unique trapID (either numeric or character) for each trap
(ii) the second column should contain the Easting of the trap.
(iii) the third column should contain the Northing of the trap.}

\item{releaseLocations}{is a list containing matrices or dataframes whose rows contain coordinates and times for the releases used.
Because individual releases can sometimes use different types of marking, the releaseLocations matrix should contain five columns:  
(i) the first column should be a markerID which can either be numeric IDs like 1,2 and 3 or a character such as "red", "blue", "green".
(ii) the second column should contain the Easting of the release locations.
(iii) the third column should contain the Northing of the release locations.
(iv) the fourth column should contain the number of days since the beginning of the MRR study at which the releases occured.  Typically the release event(s) at the beginning of the experiment will have a zero for in the fourth column (i.e. releases occured at the start of the study).
(v) the fifth column should contain the numbers of individuals released at each release event.}

\item{trapData}{is a list containing matrices or dataframes that pertain to the numbers of individuals caught in traps, the markings
on these individuals and the time intervals in which these individuals were caught.
the matrices/dataframes in trapData should be 5 columns:
(i) the first column should contain the trapID at which the trap counts were collected.
(ii) the second column should contain the number of insects caught in the trap for a particular marking type.
(iii) the third column should contain the markerID that identifies whether there are different markers (e.g. colours) used.
(iv) the fourth column should contain the time (days since the start of the study) at which this trap started collecting the individuals in this trap count record.
(v) the fifth column should contain the time (days since the start of the study) at which this trap stopped collecting the individuals in this trap count record.}

\item{distributionType}{can be one of: "poisson", "zero_inflated_poisson", "negative_binomial" or "zero_inflated_negative_binomial"}

\item{includeDrift}{is a boolean that allows you to include a constant drift applied to the dispersal kernel over time.  If TRUE, then
the drift parameters (mu_Easting and muNorthing) will be estimated as part of the parameter vector.}

\item{includeDeath}{is a boolean that allows you to apply a constant death rate to the individuals in your study.  This creates the 
realistic situation where the numbers of individuals caught in traps decreases over time.  The use of a constant death rate corresponds to 
the assumption of an exponentially distributed lifetime distribution on the individuals in the study.}
}
\value{
the function returns the log-likelihood of the mark-release-recapture data given the parameters under the chosen probability distribution for the counts.
}
\description{
\code{trapLogLikelihood} computes the log-likelihood of the observed trappind data given a set of parameters, release locations releasenumbers and trap locations.  Trap locations, release locations and trapData are all lists that contain an element for each mark-release-recapture experiment that was performed (i.e. independent experiments to use for parameter estimation).  Each of these lists contains a dataframe or matrix and are outlines as follows.
}
