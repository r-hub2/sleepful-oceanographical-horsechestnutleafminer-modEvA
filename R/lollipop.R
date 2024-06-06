#' Title Lollipop chart
#'
#' @param x a numeric vector.
#' @param names a vector of the same length as 'x' with the names to be plotted below the lollipops. If this argument is left NULL and 'x' has names, then these will be used.
#' @param ymin numeric value for the lower limit of the y axis. The default is zero. If set to NA, the minimum of 'x' will be used.
#' @param sticks logical value indicating whether the sticks of the lollipops should be drawn. The default is TRUE.
#' @param col colour for the lollipops.
#' @param grid logical, whether or not to add a grid to the plot. The default is TRUE.
#' @param cex numeric value indicating the size of the lollipops. Will be passed as 'cex' to 'points' and as 'lwd' to 'arrows' (the lines or lollipop sticks).
#' @param cex.axis numeric value indicating the size of the x and y axis labels.
#' @param las argument to pass to 'plot' indicating the orientation of the axis labels.
#' @param ... additional arguments that can be used for the plot, e.g. 'main'.
#'
#' @return This function produces a lollipop chart of the values in 'x'.
#' @export
#'
#' @examples lollipop(mtcars[,1], names = rownames(mtcars), las = 2, ylab = names(mtcars)[1], cex.axis = 0.6, main = "Lollipop chart")

lollipop <- function(x, names = NULL, ymin = 0, sticks = TRUE, col = "royalblue", grid = TRUE, cex = 1, cex.axis = 1, las = 2, ...) {
  # version 1.1 (13 Jan 2023)

  if (is.na(ymin))  ymin <- min(x, na.rm = TRUE)
  plot(c(ymin, x), axes = FALSE, type = "n", xlab = "", ...)
  if (grid) grid()
  if (is.null(names)) names <- names(x)
  axis(1, at = 1:length(x), labels = names, las = las, cex.axis = cex.axis)
  axis(2, ylim = c(ymin, max(x, na.rm = TRUE)), cex.axis = cex.axis)
  points(x, pch = 20, col = col, cex = cex)
  if (sticks)  arrows(x0 = 1:length(x), x1 = 1:length(x), y0 = 0, y1 = x, length = 0, col = col, lwd = cex)
}
