plotCIN1Incidence <- function(..., current=NULL, labels=NULL)
{
  if (dev.cur() != 1) dev.off()
  dots <- list(...)
  if (length(dots) < 1) stop("At least one scenario should be defined")
  
  CIN1Incidence1 <- function(i, x, aux){
    num <- aux[i]
    den <- sum(x[i, 1:2])
    num/den
  }
  newCIN1Cases1     <- attr(dots[[1]], "newCIN1")
  newCIN1CasesMean1 <- apply(newCIN1Cases1, 2, mean)
  
  tmp <- lapply(as.list(unique(dots[[1]]$Ages)), function(i, dat){
    aux <- dots[[1]][dots[[1]]$Ages == i, ]
    apply(aux[,1:(ncol(aux)-2)], 2, mean)
  }, dat = dots[[1]])
  meanDat1 <- as.data.frame(t(as.matrix(as.data.frame(tmp))))
  names(newCIN1CasesMean1) <- 1:length(newCIN1CasesMean1)
  totCIN1Inc1 <- unlist(lapply(as.list(1:nrow(meanDat1)), CIN1Incidence1,
                           x = meanDat1, aux = newCIN1CasesMean1))
  age <- seq(min(dots[[1]]$Ages), max(dots[[1]]$Ages), 0.5)
  par(xpd=TRUE)
  colorets <- "#000000"
  plot(age, totCIN1Inc1*attr(dots[[1]], "size"), pch = 20, xlab = "Age", ylab = paste0("CIN 1 incidence x ", attr(dots[[1]], "size")),
       main = "CIN 1 Incidence by age", axes = FALSE, type = "l", col=colorets)
  axis(2)
  axis(1)
  if (length(dots)>1)
  {
    colorets <- c(colorets, rainbow(length(dots)-1))
    for (j in 2:length(dots))
    {
      newCIN1Cases1     <- attr(dots[[j]], "newCin1Cases")
      newCIN1CasesMean1 <- apply(newCIN1Cases1, 2, mean)
      
      tmp <- lapply(as.list(unique(dots[[j]]$Ages)), function(i, dat){
        aux <- dots[[j]][dots[[j]]$Ages == i, ]
        apply(aux[,1:(ncol(aux)-2)], 2, mean)
      }, dat = dots[[j]])
      meanDat1 <- as.data.frame(t(as.matrix(as.data.frame(tmp))))
      names(newCIN1CasesMean1) <- 1:length(newCIN1CasesMean1)
      totCIN1Inc1 <- unlist(lapply(as.list(1:nrow(meanDat1)), CIN1Incidence1,
                                   x = meanDat1, aux = newCIN1CasesMean1))
      age <- min(dots[[j]]$Ages):max(dots[[j]]$Ages)
      lines(age, totCIN1Inc1*attr(dots[[1]], "size"), col=colorets[j])
    }
  }
  if (!is.null(current)) 
  {
    colorets  <- c(colorets, "#00ff00")
    h         <- as.numeric(substr(current[1, 1], 4, 5))[1]-as.numeric(substr(current[1, 1], 1, 2))[1]
    age       <- as.numeric(substr(current[1, 1], 1, 2)):as.numeric(substr(current[dim(current)[1], 1], 4, 5))
    CIN1incid <- rep(current[, 2], each=h+1)
    lines(age, CIN1incid, col=colorets[length(colorets)])
  }
  if (is.null(labels)) labels <- ""
  if (!is.null(current)) labels <- c(labels, "Current CIN1 inc.")
  legend("topright", labels, lty=1, col=colorets)
  return(list(data.frame(age=age, val=cbind(as.numeric(apply(newCIN1Cases1, 2, mean)))), 
              data.frame(age=age, val=cbind(as.numeric(totCIN1Inc1*attr(dots[[1]], "size"))))))
}