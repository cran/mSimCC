qalys <- function(scenario, disc=FALSE)
{
  if (disc==FALSE)
  {
    res <- attr(scenario, "utility")
  }else{
    res <- attr(scenario, "utilityD")
  }
  return(c(res, res/attr(scenario, "size")))
}