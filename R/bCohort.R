bCohort <- function(ind)
{
  if (class(ind)[2]!="micro.sim.individual") stop("Object is not a micro simulated cohort")
  min.age  <- attr(ind, "min.age")
  max.age  <- attr(ind, "max.age")
  steps    <- 2*(attr(ind, "max.age") - attr(ind, "min.age") + 1)
  age.vars <- colnames(ind)[1:steps]
  ages     <- seq(min.age, max.age+0.5, 0.5)

  dat.cohort <- data.frame(ages)
  for (j in 1:length(attr(ind, "hstates")))
  {
    var.name <- attr(ind, "hstates")[j]
    var.val  <- j-1
    if (!exists(eval(parse(text="var.name")))) assign(var.name, vector())
    for (i in 1:length(age.vars))
    {
      assign(var.name, c(eval(parse(text=var.name)),
                         eval(parse(text=paste0("length(ind$", age.vars[i], "[which(ind$", age.vars[i], "==j-1)])/attr(ind, 'nsim')")))))
    }
    dat.cohort <- cbind(dat.cohort, eval(parse(text=var.name)))
  }
  dat.cohort <- as.data.frame(dat.cohort)
  colnames(dat.cohort) <- c("Ages", attr(ind, "hstates"))
  attr(dat.cohort, "size")       <- attr(ind, "size")
  attr(dat.cohort, "newHPV")     <- attr(ind, "newHPV")
  attr(dat.cohort, "newCases1")  <- attr(ind, "newCases1")
  attr(dat.cohort, "newCIN1")    <- attr(ind, "newCIN1")
  attr(dat.cohort, "newCIN2")    <- attr(ind, "newCIN2")
  attr(dat.cohort, "newCIN3")    <- attr(ind, "newCIN3")
  attr(dat.cohort, "newCCDeath") <- attr(ind, "newCCDeath")
  attr(dat.cohort, "newODeath")  <- attr(ind, "newODeath")
  attr(dat.cohort, "nCyto")      <- attr(ind, "nCyto")
  class(dat.cohort)              <- c("data.frame", "micro.sim.cohort")
  return(dat.cohort)
}
