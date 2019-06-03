\encoding{utf8}
\name{microsim}
\alias{microsim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate microsimulated cohorts
}
\description{
  Generates several microsimulated cohorts with desired specifications.
}
\usage{
  microsim(seed=1234, nsim, transition, abs_states, sympt_states, prob_sympt, size, 
           p_men, min_age, max_age, utilityCoefs, costCoefs.md, costCoefs.nmd, 
           costCoefs.i, disc=3, vacc=FALSE, vacc.age=NULL, ndoses=NULL, vacc.cov=NULL, 
           vacc.eff=NULL, vacc.type=NULL, vacc.prop=NULL, vaccprice.md=NULL, 
           vaccprice.nmd=NULL, vaccprice.i=NULL, screening=FALSE, screenType=0, 
           scrSchema=0, screenPeriod=NULL, cytoType=NULL, screenPrice.md=NULL, 
           screenPrice.nmd=NULL, screenPrice.i=NULL, colpoPrice.md=NULL, 
           colpoPrice.nmd=NULL, colpoPrice.i=NULL, hpvTestPrice.md=NULL, 
           hpvTestPrice.nmd=NULL, hpvTestPrice.i=NULL, cytoHpvPrice.md=NULL, 
           cytoHpvPrice.nmd=NULL, cytoHpvPrice.i=NULL, biopsPrice.md=NULL, 
           biopsPrice.nmd=NULL, biopsPrice.i=NULL, screenCoverage=NULL, screenSensi=NULL,
           screenSensi2=NULL, screenSensi3=NULL, colpoSensi=NULL, biopSensi=NULL, 
           hpvTestSensi=NULL, treatProbs, nAnnualVisits=0, nAnnualVisitsLSIL=0, 
           nAnnualVisitsHSIL=0, cytoHPVPeriod=0, cytoHPVPostColpo=0, 
           cytoHPVPostBiop=NULL, cytoLSILperiod=0, cytoHSILperiod=0, switchAge=0, 
           C_period=NULL, hpvPeriod=0, nCores=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seed}{
  seed to be used in the simulation. Default value is \code{1234}.
}
  \item{nsim}{
  number of cohorts to be simulated.
}
  \item{transition}{
  transition probabilities matrix.
}
  \item{abs_states}{
  vector with the absorbing states.
}
  \item{sympt_states}{
  vector with the health states that might present symptoms.
}
  \item{prob_sympt}{
  vector with the probability of presenting symptoms for each health state that might present symptoms. Should have the same length of \code{sympt_states}.
}
  \item{size}{
  number of individuals on each simulated cohort.
}
  \item{p_men}{
  proportion of men in the simulated cohorts.
}
  \item{min_age}{
  lowest age in the cohort.
}
  \item{max_age}{
  largest age in the cohort.
}
  \item{utilityCoefs}{
  vector with the utilities for each health state.
}
  \item{costCoefs.md}{
  vector with the direct medical costs for each health state.
}
  \item{costCoefs.nmd}{
  vector with the direct non medical costs for each health state.
}
  \item{costCoefs.i}{
  vector with the indirect costs for each health state.
}
  \item{disc}{
  discount rate in percentage. Default value is \code{3}.
}
  \item{vacc}{
  boolean value specifying if the considered scenario includes vaccination. Default value is \code{FALSE}.
}
  \item{vacc.age}{
  vector with ages at vaccination if the considered scenario includes vaccination. Default value is \code{NULL}.
}
  \item{ndoses}{
  number of doses of vaccine if the considered scenario includes vaccination. Default value is \code{NULL}.
}
  \item{vacc.cov}{
  vaccine coverage if the considered scenario includes vaccination. Default value is \code{NULL}.
}
  \item{vacc.eff}{
  vaccine effectivity if the considered scenario includes vaccination. Default value is \code{NULL}.
}
  \item{vacc.type}{
  type of vaccine if the considered scenario includes vaccination, character with values \code{biv} for bivalent, \code{quad} for quadrivalent and \code{nona} for nonavalent vaccines. Default value is \code{NULL}.
}
  \item{vacc.prop}{
  proportion of vaccinated women on each age group if the considered scenario includes vaccination. Default value is \code{NULL}.
}
  \item{vaccprice.md}{
  vaccine direct medical costs if the considered scenario includes vaccination. Default value is \code{NULL}.
}
  \item{vaccprice.nmd}{
  vaccine direct non medical costs if the considered scenario includes vaccination. Default value is \code{NULL}.
}
  \item{vaccprice.i}{
  vaccine indirect if the considered scenario includes vaccination. Default value is \code{NULL}.
}
  \item{screening}{
  boolean specifying if the considered scenario includes screening of any type. Default value is \code{FALSE}.
}
  \item{screenType}{
  type of screening. \code{1} stands for organized screening, \code{2} stands for opportunistic screening. Default value is \code{0} (no screening).
}
  \item{scrSchema}{
  screening schema. \code{1} stands for cytology alone with repeat cytology for triage, \code{2} stands for cytology with HPV triage, \code{3} stands for HPV with cytology triage and \code{4} stands for HPV genotyping with cytology triage. Default value is \code{0} (no screening).
}
  \item{screenPeriod}{
  screening period (in years). Default value is \code{NULL} (no screening).
}
  \item{cytoType}{
  type of cytology. \code{0} stands for conventional cytology, \code{1} stands for Liquid Based Cytology (LBC). Default value is \code{NULL} (no cytology).
}
  \item{screenPrice.md}{
  medical direct cost of cytology. Default value is \code{NULL}.
}
  \item{screenPrice.nmd}{
  non-medical direct cost of cytology. Default value is \code{NULL}.
}
  \item{screenPrice.i}{
  indirect cost of cytology. Default value is \code{NULL}.
}
  \item{colpoPrice.md}{
  medical direct cost of colposcopy. Default value is \code{NULL}.
}
  \item{colpoPrice.nmd}{
  non-medical direct cost of colposcopy. Default value is \code{NULL}.
}
  \item{colpoPrice.i}{
  indirect cost of colposcopy. Default value is \code{NULL}.
}
  \item{hpvTestPrice.md}{
  medical direct cost of HPV test. Default value is \code{NULL}.
}
  \item{hpvTestPrice.nmd}{
  non-medical direct cost of HPV test. Default value is \code{NULL}.
}
  \item{hpvTestPrice.i}{
  indirect cost of HPV test. Default value is \code{NULL}.
}
  \item{cytoHpvPrice.md}{
  medical direct cost of HPV reflex test, in case \code{cytoType=1}. Default value is \code{NULL}.
}
  \item{cytoHpvPrice.nmd}{
  non-medical direct cost of HPV reflex test, in case \code{cytoType=1}. Default value is \code{NULL}.
}
  \item{cytoHpvPrice.i}{
  indirect cost of HPV reflex test, in case \code{cytoType=1}. Default value is \code{NULL}.
}
  \item{biopsPrice.md}{
  medical direct cost of biopsy. Default value is \code{NULL}.
}
  \item{biopsPrice.nmd}{
  non-medical direct cost of biopsy. Default value is \code{NULL}.
}
  \item{biopsPrice.i}{
  indirect cost of biopsy. Default value is \code{NULL}.
}
  \item{screenCoverage}{
  cytology coverage for each age group. Default value is \code{NULL}.
}
  \item{screenSensi}{
  cytology sensitivity for each age group. Default value is \code{NULL}.
}
  \item{screenSensi2}{
  cytology sensitivity after cytology for each age group. Default value is \code{NULL}.
}
  \item{screenSensi3}{
  cytology sensitivity after HPV test for each age group. Default value is \code{NULL}.
}
  \item{colpoSensi}{
  colposcopy sensitivity for each age group. Default value is \code{NULL}.
}
  \item{biopSensi}{
  biopsy sensitivity for each age group. Default value is \code{NULL}.
}
  \item{hpvTestSensi}{
  HPV test sensitivity for each age group. Default value is \code{NULL}.
}
  \item{treatProbs}{
  probability of recuperation after treatment for each FIGO I - FIGO IV states.
}
  \item{nAnnualVisits}{
  number of annual visits after colposcopy for screening schema \code{1}. Default value is \code{0}.
}
  \item{nAnnualVisitsLSIL}{
  number of annual visits after LSIL for screening schema \code{2}. Default value is \code{0}.
}
  \item{nAnnualVisitsHSIL}{
  number of annual visits after HSIL for screening schema \code{2}. Default value is \code{0}.
}
  \item{cytoHPVPeriod}{
  cytology and HPV test protocol period for screening schemas \code{3} and \code{4}. Default value is \code{0}.
}
  \item{cytoHPVPostColpo}{
  cytology and HPV test protocol period after colposcopy protocol for screening schemas \code{3} and \code{4}. Default value is \code{0}.
}
  \item{cytoHPVPostBiop}{
  cytology and HPV test protocol period after biopsy protocol for screening schemas \code{2}. Default value is \code{NULL}.
}
  \item{cytoLSILperiod}{
  period for cytology after LSIL detection for screening schame \code{2}. Default value is \code{0}.
}
  \item{cytoHSILperiod}{
  period for cytology after HSIL detection for screening schame \code{2}. Default value is \code{0}.
}
  \item{switchAge}{
  age at which screening protocol changes for screening schemas \code{3} and \code{4}. Default value is \code{0}.
}
  \item{C_period}{
  vector with screening periods (in years) before and after switch age for screening schemas \code{3} and \code{4}. Default value is \code{NULL}.
}
  \item{hpvPeriod}{
  period for HPV test in screening schema \code{2}. Default value is \code{0}.
}
  \item{nCores}{
  number of cores of the computer. Default value is \code{1}.
}
}

\value{
  Data frame containing the simulated cohorts and the individual history for each person in each simulated cohort.
}

\references{
Georgalis L, de Sanjosé S, Esnaola M, Bosch F X, Diaz M. Present and future of cervical cancer prevention in Spain: a cost-effectiveness analysis. European Journal of Cancer Prevention 2016;\bold{25}(5):430-439.

Moriña D, de Sanjosé S, Diaz M. Impact of model calibration on cost-effectiveness analysis of cervical cancer prevention 2017;\bold{7}.
}

\author{
David Moriña, Pedro Puig and Mireia Diaz
}

\seealso{
\code{\link{mSimCC-package}}, \code{\link{bCohort}}, \code{\link{costs}}, \code{\link{le}},
\code{\link{plotCIN1Incidence}}, \code{\link{plotCIN2Incidence}}, \code{\link{plotCIN3Incidence}},
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ mSimCC }
