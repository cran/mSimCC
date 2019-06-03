costs <- function(scenario, disc=FALSE)
{
  if (disc==FALSE)
  {
    mdres  <- attr(scenario, "md_cost")
    nmdres <- attr(scenario, "nmd_cost")
    ires   <- attr(scenario, "i_cost")
  }else{
    mdres  <- attr(scenario, "md_costD")
    nmdres <- attr(scenario, "nmd_costD")
    ires   <- attr(scenario, "i_costD")
  }
  res    <- mdres + nmdres + ires

  return(c(mdres, nmdres, ires, res, mdres/attr(scenario, "size"), nmdres/attr(scenario, "size"), ires/attr(scenario, "size"),
         res/attr(scenario, "size")))
}