\name{lollipop}
\alias{lollipop}
\title{Lollipop chart}
\description{
This function creates a lollipop chart from a (optionally named) numeric vector.
}
\usage{
  lollipop(x, names = NULL, ymin = 0, sticks = TRUE, col = "royalblue",
  grid = TRUE, cex = 1, cex.axis = 1, las = 2, ...)
}
\arguments{
  \item{x}{a numeric vector.}
  \item{names}{a vector of the same length as 'x' with the names to be plotted below the lollipops. If this argument is left NULL and 'x' has names, then these will be used.}
  \item{ymin}{numeric value for the lower limit of the y axis. The default is zero. If set to NA, the minimum of 'x' will be used.}
  \item{sticks}{logical value indicating whether the sticks of the lollipops should be drawn. The default is TRUE.}
  \item{col}{colour for the lollipops.}
  \item{grid}{logical, whether or not to add a grid to the plot. The default is TRUE.}
  \item{cex}{numeric value indicating the size of the lollipops. Will be passed as 'cex' to 'points' and as 'lwd' to 'arrows' (the lines or lollipop sticks).}
  \item{cex.axis}{numeric value indicating the size of the x and y axis labels.}
  \item{las}{argument to pass to \code{\link{par}} indicating the orientation of the axis labels.}
  \item{\dots}{additional arguments that can be used for the plot, e.g. 'main'.}
}
\details{
  According to modern data viz recommendations, lollipop charts are generally a better alternative to bar charts, as they reduce the visual distortion caused by the length of the bars, making it easier to compare the values.
}
\value{
  This function produces a lollipop chart of the values in 'x'.
}
\author{
  A. Marcia Barbosa
}
\seealso{
  \code{\link{barplot}}
}
\examples{
lollipop(mtcars[,1], names = rownames(mtcars), las = 2, ylab = names(mtcars)[1],
cex.axis = 0.6, main = "Lollipop chart")

lollipop(mtcars[,1], names = rownames(mtcars), las = 2, ylab = names(mtcars)[1],
cex.axis = 0.6, main = "Lollipop chart", ymin = NA)
}
