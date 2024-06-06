RMSE <- function(model = NULL, obs = NULL, pred = NULL, na.rm = TRUE, rm.dup = FALSE) {
  # version 1.0 (24 Dec 2022)

  obspred <- inputMunch(model, obs, pred, na.rm = na.rm, rm.dup = rm.dup)
  obs <- obspred[ , "obs"]
  pred <- obspred[ , "pred"]

  squared_error <- (obs - pred) ^ 2
  # return(sd(obs - pred))
  return(sqrt(mean(squared_error)))
}
