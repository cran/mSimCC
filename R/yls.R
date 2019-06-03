yls <- function(scenario1, scenario2, disc = FALSE)
{
  if (scenario1[1, 1] != scenario2[1, 1]) stop("Different cohort sizes")
  res1 <- le(scenario1, disc)
  res2 <- le(scenario2, disc)
  
  res <- res2 - res1
  return(res)
}