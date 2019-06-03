\encoding{utf8}
\name{plotCIN1Incidence}
\alias{plotCIN1Incidence}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculates and plots the CIN1 incidence.
}
\description{
  Calculates and plots the CIN1 incidence for one or several prevention strategies.
}
\usage{
  plotCIN1Incidence(..., current=NULL, labels=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{
  one or several microsimulated cohort corresponding to one or several microsimulated cohorts.
}
  \item{current}{
  real CIN 1 incidence in the population of interest.
  }
  \item{labels}{
  labels to be used in the plot.
  }
}

\value{
Returns a list with CIN 1 incidence for each age group.
}

\references{
Georgalis L, de Sanjosé S, Esnaola M, Bosch F X, Diaz M. Present and future of cervical cancer prevention in Spain: a cost-effectiveness analysis. European Journal of Cancer Prevention 2016;\bold{25}(5):430-439.

Moriña D, de Sanjosé S, Diaz M. Impact of model calibration on cost-effectiveness analysis of cervical cancer prevention 2017;\bold{7}.
}

\author{
David Moriña, Pedro Puig and Mireia Diaz
}

\seealso{
\code{\link{mSimCC-package}}, \code{\link{microsim}}, \code{\link{costs}}, \code{\link{le}},
\code{\link{bCohort}}, \code{\link{plotCIN2Incidence}}, \code{\link{plotCIN3Incidence}},
\code{\link{plotIncidence}}, \code{\link{plotMortality}}, \code{\link{plotPrevalence}},
\code{\link{qalys}}, \code{\link{yls}}
}

\examples{
data(probs)
nsim       <- 3
p.men      <- 0
size       <- 20
min.age    <- 10
max.age    <- 84

#### Natural history
hn <- microsim(seed=1234, nsim, probs, abs_states=c(10, 11), sympt_states=c(5, 6, 7, 8), 
               prob_sympt=c(0.11, 0.23, 0.66, 0.9), 
                size, p.men, min.age, max.age, 
                utilityCoefs = c(1, 1, 0.987, 0.87, 0.87, 0.76, 0.67, 0.67, 0.67, 0.938, 0, 0),
                costCoefs.md = c(0, 0, 254.1, 1495.9, 1495.9, 5546.8, 12426.4, 23123.4, 
                                 34016.6, 0, 0, 0),
                costCoefs.nmd = c(0, 0, 81.4, 194.1, 194.1, 219.1, 219.1, 219.1, 219.1, 0, 0, 0),
                costCoefs.i = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), disc=3, 
                treatProbs=c(0,0,1,1,1,0.9894,0.9422,0.8262,0.5507,0,0,0),
                nCores=1) ### individual level
hn_c <- bCohort(hn)
plotCIN1Incidence(hn_c) ### Aggregated level
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ mSimCC }
