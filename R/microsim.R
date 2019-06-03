microsim <- function(seed=1234, nsim, transition, abs_states, sympt_states, prob_sympt, size, p_men, min_age, max_age,
                      utilityCoefs, costCoefs.md, costCoefs.nmd, costCoefs.i, disc=3, vacc=FALSE, vacc.age=NULL, ndoses=NULL,
                      vacc.cov=NULL, vacc.eff=NULL, vacc.type=NULL, vacc.prop=NULL, vaccprice.md=NULL, vaccprice.nmd=NULL, vaccprice.i=NULL,
                      screening=FALSE, screenType=0, scrSchema=0, screenPeriod=NULL, cytoType=NULL, screenPrice.md=NULL, screenPrice.nmd=NULL, screenPrice.i=NULL,
                      colpoPrice.md=NULL, colpoPrice.nmd=NULL, colpoPrice.i=NULL, hpvTestPrice.md=NULL, hpvTestPrice.nmd=NULL, hpvTestPrice.i=NULL,
                      cytoHpvPrice.md=NULL, cytoHpvPrice.nmd=NULL, cytoHpvPrice.i=NULL,
                      biopsPrice.md=NULL, biopsPrice.nmd=NULL, biopsPrice.i=NULL,
                      screenCoverage=NULL, screenSensi=NULL, screenSensi2=NULL, screenSensi3=NULL, colpoSensi=NULL, biopSensi=NULL, hpvTestSensi=NULL,
                      treatProbs, nAnnualVisits=0, nAnnualVisitsLSIL=0, nAnnualVisitsHSIL=0, cytoHPVPeriod=0, cytoHPVPostColpo=0, cytoHPVPostBiop=NULL,
                      cytoLSILperiod=0, cytoHSILperiod=0, switchAge=0, C_period=NULL, hpvPeriod=0, nCores=1)
{
  ### vaccination inputs checking
  if (length(vacc.age) != length(vacc.cov) | length(vacc.age) != length(vacc.eff) | length(vacc.cov) != length(vacc.eff) |
      length(vacc.age) != length(ndoses)) stop ("Wrong vaccination inputs")
  if (vacc==FALSE & (!is.null(vacc.age) | !is.null(ndoses) | !is.null(vacc.cov) | !is.null(vacc.eff) | !is.null(vacc.prop) |
                     !is.null(vaccprice.md) | !is.null(vaccprice.nmd) | !is.null(vaccprice.i))) stop ("Wrong vaccination inputs")
  if (vacc==TRUE & (is.null(vacc.age) | is.null(ndoses) | is.null(vacc.cov) | is.null(vacc.eff) | is.null(vacc.prop) |
                    is.null(vaccprice.md) | is.null(vaccprice.nmd) | is.null(vaccprice.i))) stop ("Wrong vaccination inputs")
  if (vacc==TRUE & (is.null(vacc.type))) stop ("Wrong vaccination inputs")
  Call               <- match.call()
  transition$min.age <- as.numeric(substr(transition[, 1], 1, 2))
  transition$max.age <- as.numeric(substr(transition[, 1], 4, 5))
  dat.fin            <- data.frame(id=seq(1, size))
  steps              <- (max_age-min_age+1)*2 ### Step unit is semester and not year
  newHPV             <- matrix(nrow=nsim, ncol=steps)
  newCC              <- matrix(nrow=nsim, ncol=steps)
  newCIN1            <- matrix(nrow=nsim, ncol=steps)
  newCIN2            <- matrix(nrow=nsim, ncol=steps)
  newCIN3            <- matrix(nrow=nsim, ncol=steps)
  newCCDeath         <- matrix(nrow=nsim, ncol=steps)
  newODeath          <- matrix(nrow=nsim, ncol=steps)
  newCytos           <- matrix(nrow=nsim, ncol=steps)
  newColpos          <- matrix(nrow=nsim, ncol=steps)
  newHPVTests        <- matrix(nrow=nsim, ncol=steps)
  newBiops           <- matrix(nrow=nsim, ncol=steps)
  abs_states         <- as.list(abs_states)
  abs_states         <- lapply(abs_states, as.integer)
  sympt_states       <- as.list(sympt_states)
  sympt_states       <- lapply(sympt_states, as.numeric)
  prob_sympt         <- as.list(prob_sympt)
  prob_sympt         <- lapply(prob_sympt, as.numeric)
  dage               <- transition$max.age-transition$min.age+1
  dage               <- as.list(dage)
  dage               <- lapply(dage, as.integer)
  transition2        <- transition[, 2:(dim(transition)[2])]
  utilityCoefs       <- as.list(utilityCoefs)
  utilityCoefs       <- lapply(utilityCoefs, as.numeric)
  costCoefs.md       <- as.list(costCoefs.md)
  costCoefs.md       <- lapply(costCoefs.md, as.numeric)
  costCoefs.nmd      <- as.list(costCoefs.nmd)
  costCoefs.nmd      <- lapply(costCoefs.nmd, as.numeric)
  costCoefs.i        <- as.list(costCoefs.i)
  costCoefs.i        <- lapply(costCoefs.i, as.numeric)
  treatProbs         <- as.list(treatProbs)
  treatProbs         <- lapply(treatProbs, as.numeric)

  if (screening==FALSE)
  {
    screenCoverage  <- rep(0, dim(transition)[2]-3)
    screenSensi     <- rep(0, dim(transition)[2]-3)
    screenSensi2    <- rep(0, dim(transition)[2]-3)
    screenSensi3    <- rep(0, dim(transition)[2]-3)
    screenPrice.md  <- 0
    screenPrice.nmd <- 0
    screenPrice.i   <- 0
    colpoPrice.md   <- 0
    colpoPrice.nmd  <- 0
    colpoPrice.i    <- 0
    hpvTestPrice.md <- 0
    hpvTestPrice.nmd<- 0
    hpvTestPrice.i  <- 0
    biopsPrice.md   <- 0
    biopsPrice.nmd  <- 0
    biopsPrice.i    <- 0
    cytoHpvPrice.md <- 0
    cytoHpvPrice.nmd <- 0
    cytoHpvPrice.i  <- 0
    colpoSensi      <- rep(0, dim(transition)[2]-3)
    hpvTestSensi    <- rep(0, dim(transition)[2]-3)
    biopSensi       <- rep(0, dim(transition)[2]-3)
    screenPeriod    <- 0
    C_period        <- rep(0, 2)
    cytoHPVPostBiop <- rep(0, 2)
  }
  if ((screening==TRUE & scrSchema==1))
  {
    screenCoverage     <- as.list(screenCoverage)
    screenCoverage     <- lapply(screenCoverage, as.numeric)
    screenSensi        <- as.list(screenSensi)
    screenSensi        <- lapply(screenSensi, as.numeric)
    screenSensi2       <- as.list(screenSensi)
    screenSensi2       <- lapply(screenSensi, as.numeric)
    screenSensi3       <- rep(0, dim(transition)[2]-3)
    colpoSensi         <- as.list(colpoSensi)
    colpoSensi         <- lapply(colpoSensi, as.numeric)
    hpvTestSensi       <- rep(0, dim(transition)[2]-3)
    biopSensi          <- rep(0, dim(transition)[2]-3)
    hpvTestPrice.md    <- 0
    hpvTestPrice.nmd   <- 0
    hpvTestPrice.i     <- 0
    biopsPrice.md      <- 0
    biopsPrice.nmd     <- 0
    biopsPrice.i       <- 0
    cytoHpvPrice.md    <- 0
    cytoHpvPrice.nmd   <- 0
    cytoHpvPrice.i     <- 0
    C_period           <- rep(0, 2)
    cytoHPVPostBiop    <- rep(0, 2)
  }
  if (screening==TRUE & scrSchema==2)
  {
     screenCoverage     <- as.list(screenCoverage)
     screenCoverage     <- lapply(screenCoverage, as.numeric)
     screenSensi        <- as.list(screenSensi)
     screenSensi        <- lapply(screenSensi, as.numeric)
     screenSensi2       <- as.list(screenSensi2)
     screenSensi2       <- lapply(screenSensi2, as.numeric)
     screenSensi3       <- as.list(screenSensi3)
     screenSensi3       <- lapply(screenSensi3, as.numeric)
     biopSensi          <- as.list(biopSensi)
     biopSensi          <- lapply(biopSensi, as.numeric)
     colpoSensi         <- as.list(colpoSensi)
     colpoSensi         <- lapply(colpoSensi, as.numeric)
     hpvTestSensi       <- as.list(hpvTestSensi)
     hpvTestSensi       <- lapply(hpvTestSensi, as.numeric)
     C_period           <- rep(0, 2)
     cytoHPVPostBiop    <- as.list(cytoHPVPostBiop)
     cytoHPVPostBiop    <- lapply(cytoHPVPostBiop, as.integer)
     if (cytoType==0)
     {
       cytoHpvPrice.md    <- 0
       cytoHpvPrice.nmd   <- 0
       cytoHpvPrice.i     <- 0
     }
  }
  if (screening==TRUE & ((scrSchema==3) | (scrSchema==4)))
  {
    screenCoverage     <- as.list(screenCoverage)
    screenCoverage     <- lapply(screenCoverage, as.numeric)
    screenSensi        <- as.list(screenSensi)
    screenSensi        <- lapply(screenSensi, as.numeric)
    screenSensi2       <- as.list(screenSensi2)
    screenSensi2       <- lapply(screenSensi2, as.numeric)
    screenSensi3       <- as.list(screenSensi3)
    screenSensi3       <- lapply(screenSensi3, as.numeric)
    biopSensi          <- as.list(biopSensi)
    biopSensi          <- lapply(biopSensi, as.numeric)
    colpoSensi         <- as.list(colpoSensi)
    colpoSensi         <- lapply(colpoSensi, as.numeric)
    hpvTestSensi       <- as.list(hpvTestSensi)
    hpvTestSensi       <- lapply(hpvTestSensi, as.numeric)
    C_period           <- as.list(C_period)
    C_period           <- lapply(C_period, as.integer)
    screenPeriod       <- 0
    cytoHPVPostBiop    <- as.list(cytoHPVPostBiop)
    cytoHPVPostBiop    <- lapply(cytoHPVPostBiop, as.integer)
  }
  if (vacc==FALSE)
  {
    ndoses        <- 0
    vacc.cov      <- 0
    vacc.eff      <- 0
    vaccprice.md  <- 0
    vaccprice.nmd <- 0
    vaccprice.i   <- 0
    vacc.age      <- 0
    vacc.prop     <- 0
  }
  pCC <- 0
  if (vacc==TRUE)
  {
    if (vacc.type=="biv" | vacc.type=="quad") pCC <- 0.7
    if (vacc.type=="nona") pCC <- 0.94
}
  vacc.cov     <- vacc.cov*vacc.prop
  vacc.age     <- as.list(vacc.age)
  vacc.age     <- lapply(vacc.age, as.numeric)
  vacc.cov     <- as.list(vacc.cov)
  vacc.cov     <- lapply(vacc.cov, as.numeric)
  vacc.eff     <- vacc.eff*pCC
  vacc.eff     <- as.list(vacc.eff)
  vacc.eff     <- lapply(vacc.eff, as.numeric)
  nVaccPeriods <- ifelse(vacc==TRUE, length(vacc.age), 0)
  nabs_states  <- length(abs_states)
  ndoses       <- as.list(ndoses)
  ndoses       <- lapply(ndoses, as.numeric)

  genCohort <- function(k)
  {
    set.seed(seed+k)
    res2  <- .Call("Cmicrosim", transition2, as.integer(dim(transition)[2]-3), abs_states, as.integer(nabs_states), sympt_states,
                   as.integer(length(sympt_states)), prob_sympt, as.integer(size), p_men, as.integer(steps), as.integer(min_age),
                   as.integer(dage)[1], utilityCoefs, costCoefs.md, costCoefs.nmd, costCoefs.i, as.numeric(disc), as.integer(vacc),
                   vacc.age, ndoses, vacc.cov, vacc.eff, vaccprice.md, vaccprice.nmd, vaccprice.i, as.integer(nVaccPeriods),
                   as.integer(screening), as.integer(screenType), as.integer(scrSchema), as.integer(screenPeriod), as.integer(cytoType),
                   screenPrice.md, screenPrice.nmd, screenPrice.i, colpoPrice.md, colpoPrice.nmd, colpoPrice.i, biopsPrice.md, biopsPrice.nmd,
                   biopsPrice.i, screenCoverage, screenSensi, screenSensi2, screenSensi3, colpoSensi, biopSensi,
                   cytoHpvPrice.md, cytoHpvPrice.nmd, cytoHpvPrice.i,
                   hpvTestSensi, hpvTestPrice.md, hpvTestPrice.nmd, hpvTestPrice.i, treatProbs, as.integer(nAnnualVisits),
                   as.integer(nAnnualVisitsLSIL), as.integer(nAnnualVisitsHSIL), as.integer(cytoHPVPeriod), as.integer(cytoHPVPostColpo),
                   cytoHPVPostBiop, as.integer(cytoLSILperiod), as.integer(cytoHSILperiod), as.integer(switchAge), C_period, as.integer(hpvPeriod))
    dat.fin  <- data.frame(cbind(res2, sim=rep(k, size+11)))
    return(dat.fin)
  }

  # Initiate cluster
  k  <- NULL
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl))
  res1 <- foreach(k=1:nsim, .combine=rbind) %dopar% genCohort(k)
  dat.fin <- as.data.frame(res1)
  idx <- vector()
  for (i in 1:nsim)
  {
    dat.fin[i*(size+10)+i-10, (steps+1):(ncol(dat.fin)-1)] <- NA
    dat.fin[i*(size+10)+i-9, (steps+1):(ncol(dat.fin)-1)]  <- NA
    dat.fin[i*(size+10)+i-8, (steps+1):(ncol(dat.fin)-1)]  <- NA
    dat.fin[i*(size+10)+i-7, (steps+1):(ncol(dat.fin)-1)]  <- NA
    dat.fin[i*(size+10)+i-6, (steps+1):(ncol(dat.fin)-1)]  <- NA
    dat.fin[i*(size+10)+i-5, (steps+1):(ncol(dat.fin)-1)]  <- NA
    dat.fin[i*(size+10)+i-4, (steps+1):(ncol(dat.fin)-1)]  <- NA
    dat.fin[i*(size+10)+i-3, (steps+1):(ncol(dat.fin)-1)]  <- NA
    dat.fin[i*(size+10)+i-2, (steps+1):(ncol(dat.fin)-1)]  <- NA
    dat.fin[i*(size+10)+i-1, (steps+1):(ncol(dat.fin)-1)]  <- NA
    dat.fin[i*(size+10)+i, (steps+1):(ncol(dat.fin)-1)]    <- NA
    newHPV[i, ]     <- as.matrix(dat.fin[i*(size+10)+i-10, 1:(steps)])
    newCC[i, ]      <- as.matrix(dat.fin[i*(size+10)+i-9, 1:(steps)])
    newCIN1[i, ]    <- as.matrix(dat.fin[i*(size+10)+i-8, 1:(steps)])
    newCIN2[i, ]    <- as.matrix(dat.fin[i*(size+10)+i-7, 1:(steps)])
    newCIN3[i, ]    <- as.matrix(dat.fin[i*(size+10)+i-6, 1:(steps)])
    newCCDeath[i, ] <- as.matrix(dat.fin[i*(size+10)+i-5, 1:(steps)])
    newODeath[i, ]  <- as.matrix(dat.fin[i*(size+10)+i-4, 1:(steps)])
    newCytos[i, ]   <- as.matrix(dat.fin[i*(size+10)+i-3, 1:(steps)])
    newColpos[i, ]  <- as.matrix(dat.fin[i*(size+10)+i-2, 1:(steps)])
    newHPVTests[i, ]<- as.matrix(dat.fin[i*(size+10)+i-1, 1:(steps)])
    newBiops[i, ]   <- as.matrix(dat.fin[i*(size+10)+i, 1:(steps)])
  }
  for (i in 1:nsim)
  {
    idx <- c(idx,i*(size+10)+i-10, i*(size+10)+i-9, i*(size+10)+i-8, i*(size+10)+i-7, i*(size+10)+i-6,
             i*(size+10)+i-5, i*(size+10)+i-4, i*(size+10)+i-3, i*(size+10)+i-2, i*(size+10)+i-1,
             i*(size+10)+i)
  }
  dat.fin <- dat.fin[-idx, ]

  for (i in (steps):2)
  {
    newHPV[, i]     <- newHPV[, i-1]
    newCC[, i]      <- newCC[, i-1]
    newCIN1[, i]    <- newCIN1[, i-1]
    newCIN2[, i]    <- newCIN2[, i-1]
    newCIN3[, i]    <- newCIN3[, i-1]
    newCCDeath[, i] <- newCCDeath[, i-1]
    newODeath[, i]  <- newODeath[, i-1]
    newCytos[, i]   <- newCytos[, i-1]
    newColpos[, i]  <- newColpos[, i-1]
    newHPVTests[, i]<- newHPVTests[, i-1]
    newBiops[, i]   <- newBiops[, i-1]
  }
  newHPV[, 1]     <- 0
  newCC[, 1]      <- 0
  newCIN1[, 1]    <- 0
  newCIN2[, 1]    <- 0
  newCIN3[, 1]    <- 0
  newCCDeath[, 1] <- 0
  newODeath[, 1]  <- 0
  newCytos[, 1]   <- 0
  newColpos[, 1]  <- 0
  newHPVTests[, 1]<- 0
  newBiops[, 1]   <- 0
  attr(dat.fin, "Call")       <- as.list(Call)[-1]
  attr(dat.fin, "size")       <- size
  attr(dat.fin, "min.age")    <- min_age
  attr(dat.fin, "max.age")    <- max_age
  attr(dat.fin, "hstates")    <- colnames(transition)[2:(ncol(transition)-2)]
  attr(dat.fin, "nsim")       <- nsim
  attr(dat.fin, "sex")        <- dat.fin[, steps+1]
  attr(dat.fin, "utility")    <- sum(dat.fin[, steps+2])/nsim
  attr(dat.fin, "md_cost")    <- sum(dat.fin[, steps+3])/nsim
  attr(dat.fin, "nmd_cost")   <- sum(dat.fin[, steps+4])/nsim
  attr(dat.fin, "i_cost")     <- sum(dat.fin[, steps+5])/nsim
  attr(dat.fin, "utilityD")   <- sum(dat.fin[, steps+6])/nsim
  attr(dat.fin, "md_costD")   <- sum(dat.fin[, steps+7])/nsim
  attr(dat.fin, "nmd_costD")  <- sum(dat.fin[, steps+8])/nsim
  attr(dat.fin, "i_costD")    <- sum(dat.fin[, steps+9])/nsim
  attr(dat.fin, "le")         <- sum(dat.fin[, steps+10])/nsim
  attr(dat.fin, "leD")        <- sum(dat.fin[, steps+11])/nsim
  attr(dat.fin, "nVaccs")     <- sum(dat.fin[, steps+12])/nsim
  attr(dat.fin, "nCyto")      <- newCytos
  attr(dat.fin, "nColpo")     <- newColpos
  attr(dat.fin, "nHPVTests")  <- newHPVTests
  attr(dat.fin, "nBiops")     <- newBiops
  attr(dat.fin, "newHPV")     <- newHPV
  attr(dat.fin, "newCases1")  <- newCC
  attr(dat.fin, "newCIN1")    <- newCIN1
  attr(dat.fin, "newCIN2")    <- newCIN2
  attr(dat.fin, "newCIN3")    <- newCIN3
  attr(dat.fin, "newCCDeath") <- newCCDeath
  attr(dat.fin, "newODeath")  <- newODeath
  class(dat.fin)              <- c("data.frame", "micro.sim.individual")
  colnames(dat.fin) <- c(paste0("Sem_", seq(1:steps)), "sex", "utilities", "md_cost", "nmd_cost", "i_cost",
                         "utilitiesD", "md_costD", "nmd_costD", "i_costD", "le", "leD", "Vaccinated", "Treated", "Symptoms",
                         "Cytology", "HPVTest", "Colposcopy", "Biopsy", "sim")
  return(dat.fin)
}
