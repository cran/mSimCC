byAgeGroup2 <- function(data)
{
  data <- as.data.frame(data)
  if (dim(data)[2]!=2) stop("Wrong length")
  data$group = rep(1:(nrow(data)/10), each=10)
  res <- data.frame(Age.groups = c("10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44",
                                   "45-49", "50-54", "55-59", "60-64", "65-69", "70-74",
                                   "75-79", "80-84"))
  res$value <- tapply(data[, 2], data$group, FUN = sum)
  return(res)
}
