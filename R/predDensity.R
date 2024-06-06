predDensity <- function(model = NULL, obs = NULL, pred = NULL, separate = TRUE, type = "both", ci = NA, legend.pos = "topright", main = "Density of predicted values", na.rm = TRUE, rm.dup = FALSE, xlim = NULL, ...) {
  # version 1.7 (2 Fev 2024)

  stopifnot(is.null(xlim) || (is.numeric(xlim) && is.finite(xlim) && length(xlim) == 2))

  if (!is.na(ci)) {
    # if (separate) {
      # message("Because 'ci' is not NULL, 'separate' was changed to FALSE.")
      separate <- FALSE
    # }
    if (!is.numeric(ci) || length(ci) != 1 || ci < 0 || ci > 1) stop("'ci' must be either NULL or a numeric value between 0 and 1.")  # new
  }

  if (is.null(obs) && is.null(model)) {
    if (is.null(pred)) stop ("You must provide either 'model' or 'pred'.")
    separate <- FALSE
    obs <- sample(c(0, 1), length(pred), replace = TRUE)  # dummy variable
  }

  obspred <- inputMunch(model, obs, pred, na.rm = na.rm, rm.dup = rm.dup)
  obs <- obspred[ , "obs"]
  pred <- obspred[ , "pred"]

  # if (is.null(obs)) {
  #   if (is.null(pred)) stop ("You must provide either 'model' or 'pred'.")
  #   # if (separate) message("'obs' not provided, so 'separate' automatically set to FALSE.")
  #   separate <- FALSE
  #   obs <- sample(c(0, 1), length(pred), replace = TRUE)  # dummy variable
  # } else {
  #   if (length(obs) != length(pred)) stop("'obs' and 'pred' must have the same length.")
  # }

  pred0 <- pred[obs == 0]
  pred1 <- pred[obs == 1]

  type <- match.arg(type, c("histogram", "density", "both"))

  #yy <- NULL
  #if ("histogram" %in% type) yy <- c(yy, pmax(hist0$density, hist1$density))
  #if ("density" %in% type) yy <- c(yy, pmax(hist0$density, hist1$density))
  #plot(x = c(0, 1), y = pmax(hist0$density, hist1$density), type = "n", xlab = "Predicted value", ylab = "Density")

  rslt <- vector("list")

  if (type %in% c("density", "both")) {  # "density" %in% type
    if (!separate) {
      dens <- density(pred)
      xrange <- xlim
      if (is.null(xlim)) xrange <- range(dens$x, finite = TRUE)
      yrange <- range(dens$y, finite = TRUE)
      rslt[["density"]] <- dens
    } else {
      dens0 <- density(pred0)
      dens1 <- density(pred1)
      xrange <- xlim
      if (is.null(xlim)) xrange <- range(c(dens0$x, dens1$x), finite = TRUE)
      yrange <- range(dens0$y, dens1$y, finite = TRUE)
      rslt[["density_obs1"]] <- dens1
      rslt[["density_obs0"]] <- dens0
    }
    plot(x = xrange, y = yrange, xlab = "Predicted value", ylab = "Density", type = "n", main = main)
  }

  if (type %in% c("histogram", "both")) {  # "histogram" %in% type
    #brks <- ifelse(is.null(breaks), "Sturges", breaks)
    hist0 <- hist(pred0, plot = FALSE, ...)
    hist1 <- hist(pred1, plot = FALSE, ...)
    if (type == "histogram") {  # !("density" %in% type)
      yrange <- range(hist0$density, hist1$density, finite = TRUE)
      plot(x = c(0, 1), y = yrange, type = "n", xlab = "Predicted value", ylab = "Density", main = main)
    }
    if (!separate) {
      histogram <- hist(c(pred0, pred1), freq = FALSE, col = "grey40", add = TRUE, ...)
      rslt[["histogram"]] <- histogram
      } else {
      hist(pred1, freq = FALSE, col = "grey40", add = TRUE, ...)
      hist(pred0, freq = FALSE, col = "darkgrey", density = 40, angle = 45, add = TRUE, ...)
      rslt[["histogram_obs1"]] <- hist1
      rslt[["histogram_obs0"]] <- hist0
      if (legend.pos != "n" && type == "histogram") legend(legend.pos, legend = c("absences", "presences"), fill = c("darkgrey", "grey40"), border = NA, density = c(40, NA), bty = "n")
    }
  }

  if (type %in% c("density", "both")) {  # "density" %in% type
    if (!separate) {
      lines(dens, col = "grey10", lwd = 2)
    } else {
      lines(dens1, col = "grey10", lwd = 2)
      lines(dens0, col = "darkgrey", lty = 5, lwd = 2)
      if (!is.na(legend.pos) && legend.pos != "n" && type == "density") legend(legend.pos, legend = c("absences", "presences"), col = c("darkgrey", "grey10"), lty = c(5, 1), bty = "n")
      if (!is.na(legend.pos) && legend.pos != "n" && type == "both") legend(legend.pos, legend = c("absences", "presences"), fill = c("darkgrey", "grey40"), border = NA, lty = c(5, 1), col = c("darkgrey", "grey20"), density = c(40, NA), bty = "n")
    }
  }

  if (!is.na(ci)) {
    quants <- quantile(pred, probs = c(1 - ci, ci))
    maxs <- vector("numeric")
    if ("density" %in% names(rslt)) maxs <- c(maxs, max(rslt$density$y))
    if ("histogram" %in% names(rslt)) maxs <- c(maxs, max(rslt$histogram$density))
    # abline(v = quants, col = "darkgreen")
    rect(quants[1], 0, quants[2], max(maxs), col = adjustcolor("darkgreen", alpha.f = 0.3), border = NA)
    abline(v = mean(pred), lwd = 2, col = "darkblue")
  }

  return(rslt)
}
