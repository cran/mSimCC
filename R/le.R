le <- function(scenario, disc=FALSE)
{
  if (disc==FALSE)
  {
    res <- attr(scenario, "le")
  }else{
    res <- attr(scenario, "leD")
  }
  return(res/attr(scenario, "size"))
}