plotMortality <- function(..., current=NULL, labels=NULL)
{
  if (dev.cur() != 1) dev.off()
  dots <- list(...)
  if (length(dots) < 1) stop("At least one scenario should be defined")
  
  mort1 <- function(i, x, CCdeathsMean1)
  {
    num <- mean(CCdeathsMean1[i])
    surv <- meanDat1[i, 1:9]
    den <- mean(apply(surv, 1, sum))
    num/den
  }
  
  CCdeaths1 <- attr(dots[[1]], "newCCDeath")
  CCdeathsMean1 <- apply(CCdeaths1, 2, mean)
  
  tmp <- lapply(as.list(unique(dots[[1]]$Ages)), function(i, dat){
    aux <- dots[[1]][dots[[1]]$Ages == i, ]
    apply(aux[,1:(ncol(aux)-2)], 2, mean)
  }, dat = dots[[1]])
  meanDat1 <- as.data.frame(t(as.matrix(as.data.frame(tmp))))
  names(CCdeathsMean1) <- 1:length(CCdeathsMean1)
  totMort1 <- unlist(lapply(as.list(1:nrow(meanDat1)), mort1,
                           x = meanDat1, CCdeathsMean1 = CCdeathsMean1))
  
  age <- seq(min(dots[[1]]$Ages), max(dots[[1]]$Ages), 0.5)
  par(xpd=TRUE)
  colorets <- "#000000"
  plot(age, totMort1*attr(dots[[1]], "size"), pch = 20, xlab = "Age", ylab = paste0("Mortality x ", attr(dots[[1]], "size")),
       main = "CC mortality by age", axes = FALSE, type = "l", col=colorets)
  axis(2)
  axis(1)
  if (length(dots)>1)
  {
    colorets <- c(colorets, rainbow(length(dots)-1))
    for (j in 2:length(dots))
    {
      CCdeaths1 <- attr(dots[[j]], "CCdeaths")
      CCdeathsMean1 <- apply(CCdeaths1, 2, mean)
      
      tmp <- lapply(as.list(unique(dots[[j]]$Ages)), function(i, dat){
        aux <- dots[[j]][dots[[j]]$Ages == i, ]
        apply(aux[,1:(ncol(aux)-2)], 2, mean)
      }, dat = dots[[j]])
      meanDat1 <- as.data.frame(t(as.matrix(as.data.frame(tmp))))
      names(CCdeathsMean1) <- 1:length(CCdeathsMean1)
      totMort1 <- unlist(lapply(as.list(1:nrow(meanDat1)), mort1,
                                x = meanDat1, CCdeathsMean1 = CCdeathsMean1))
      
      age <- min(dots[[j]]$Ages):max(dots[[j]]$Ages)
      lines(age, totMort1*attr(dots[[1]], "size"), col=colorets[j])
    }
  }
  if (!is.null(current)) 
  {
    colorets <- c(colorets, "#00ff00")
    h        <- as.numeric(substr(current[1, 1], 4, 5))[1]-as.numeric(substr(current[1, 1], 1, 2))[1]
    age      <- as.numeric(substr(current[1, 1], 1, 2)):as.numeric(substr(current[dim(current)[1], 1], 4, 5))
    mort     <- rep(current[, 2], each=h+1)
    lines(age, mort, col=colorets[length(colorets)])
  }
  if (is.null(labels)) labels <- ""
  if (!is.null(current)) labels <- c(labels, "Current mort.")
  legend("topright", labels, lty=1, col=colorets)
  return(list(data.frame(age=age, val=cbind(as.numeric(apply(CCdeaths1, 2, mean)))), 
              data.frame(age=age, val=cbind(as.numeric(totMort1*attr(dots[[1]], "size"))))))
}