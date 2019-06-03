plotPrevalence <- function(..., current=NULL, labels=NULL)
{
  if (dev.cur() != 1) dev.off()
  dots <- list(...)
  if (length(dots) < 1) stop("At least one scenario should be defined")
  
  if (class(dots[[1]])[2]!="micro.sim.cohort") stop("Not a proper cohort provided")
  prevalence1 <- function(x){
    num <- x[3] 
    den <- sum(x[2:11])
    return(num/den)
  }
  newCases1 <- attr(dots[[1]], "newHPV")
  newCasesMean1 <- apply(newCases1, 2, mean)
  newCasesMedian1 <- apply(newCases1, 2, median)
  newCasesCentralMean1 <- apply(newCases1, 2, function(x){
    xS <- sort(x)
    mean(xS[(round(length(xS)/2)-2):(round(length(xS)/2)+2)])
  })
  tmp <- lapply(as.list(unique(dots[[1]]$Ages)), function(i, dat){
    aux <- dots[[1]][dots[[1]]$Ages == i, ]
    apply(aux[,1:(ncol(aux)-2)], 2, mean)
  }, dat = dots[[1]])
  meanDat1 <- as.data.frame(t(as.matrix(as.data.frame(tmp))))
  names(newCasesMean1) <- 1:length(newCasesMean1)
  names(newCasesCentralMean1) <- 1:length(newCasesCentralMean1)
  totPrev1 <- apply(meanDat1, 1, prevalence1)
  
  age <- seq(min(dots[[1]]$Ages), max(dots[[1]]$Ages), 0.5)
  par(xpd=TRUE)
  colorets <- "#000000"
  plot(age, totPrev1, pch = 20, xlab = "Age", ylab = "Prevalence",
       main = "HPV Prevalence by age", axes = FALSE, type = "l", col=colorets)
  axis(2)
  axis(1)
  if (length(dots)>1)
  {
    colorets <- c(colorets, rainbow(length(dots)-1))
    for (j in 2:length(dots))
    {
      newCases1 <- attr(dots[[j]], "newHPV")
      newCasesMean1 <- apply(newCases1, 2, mean)
      newCasesMedian1 <- apply(newCases1, 2, median)
      newCasesCentralMean1 <- apply(newCases1, 2, function(x){
        xS <- sort(x)
        mean(xS[(round(length(xS)/2)-2):(round(length(xS)/2)+2)])
      })
      tmp <- lapply(as.list(unique(dots[[j]]$Ages)), function(i, dat){
        aux <- dots[[j]][dots[[j]]$Ages == i, ]
        apply(aux[,1:(ncol(aux)-2)], 2, mean)
      }, dat = dots[[j]])
      meanDat1 <- as.data.frame(t(as.matrix(as.data.frame(tmp))))
      names(newCasesMean1) <- 1:length(newCasesMean1)
      names(newCasesCentralMean1) <- 1:length(newCasesCentralMean1)
      totPrev1 <- apply(meanDat1, 1, prevalence1)
      
      age <- min(dots[[j]]$Ages):max(dots[[j]]$Ages)
      lines(age, totPrev1, col=colorets[j])
    }
  }
  if (!is.null(current)) 
  {
    colorets <- c(colorets, "#00ff00")
    h        <- as.numeric(substr(current[1, 1], 4, 5))[1]-as.numeric(substr(current[1, 1], 1, 2))[1]
    age      <- as.numeric(substr(current[1, 1], 1, 2)):as.numeric(substr(current[dim(current)[1], 1], 4, 5))
    incid    <- rep(current[, 2], each=h+1)
    lines(age, incid, col=colorets[length(colorets)])
  }
  if (is.null(labels)) labels <- ""
  if (!is.null(current)) labels <- c(labels, "Current prev.")
  legend("topright", labels, lty=1, col=colorets)
  return(list(data.frame(age=age, val=cbind(as.numeric(apply(newCases1, 2, mean)))), 
              data.frame(age=age, val=cbind(as.numeric(totPrev1)))))
}