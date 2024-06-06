mod2obspred <- function(model, obs.only = FALSE) {
  # version 1.1 (22 Dez 2022)

  if (is(model, "glm") || is(model, "gam")) {
    obs <- model$y
    if (!obs.only) pred <- model$fitted.values
  } else if (is(model, "gbm")) {
    obs <- model$data$y
    # logit <- function(x) exp(x) / (1 + exp(x))
    if (!obs.only) pred <- suppressMessages(predict(model, type = "response"))  # logit(model$fit)
  } else if (is(model, "randomForest")) {
    obs <- as.integer(as.character(model$y))
    if (!obs.only) pred <- predict(model, type = "prob")[ , "1"]
  } else if (is(model, "bart")) {
    if (is.null(model$fit$data)) stop("'$fit$data' section not available in 'model'. Try running 'bart' with 'keeptrees=TRUE'.")
    obs <- model$fit$data@y  # requires model ran with keeptrees=TRUE
    if (!obs.only) pred <- fitted(model, type = "response")
  } else stop("'model' is of a non-implemented class.")

  if (obs.only) return(data.frame(obs = obs))
  data.frame(obs = obs, pred = pred)
}
