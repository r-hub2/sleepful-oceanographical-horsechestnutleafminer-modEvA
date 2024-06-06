\name{mod2obspred}
\alias{mod2obspred}
\title{
Extract observed and predicted values from a model object.
}
\description{
This function takes a model object and returns the observed and (optionally) the fitted values in that model.
}
\usage{
mod2obspred(model, obs.only = FALSE)
}
\arguments{
  \item{model}{a model object of class "glm", "gam", "gbm", "randomForest" or "bart" from which the response variable and fitted (predicted) values can be extracted. Note that, for "randomForest" models, only the out-of-bag prediction is available from the model object (see \code{?predict.randomForest} if you have that package installed), so here you'll get different results if you provide 'model' or the modelled 'obs' (and corresponding 'pred') values.}
  \item{obs.only}{logical value indicating whether only 'obs' should be obtained (saves computing time when 'pred' not needed -- used e.g. by \code{\link{prevalence}}). Defaults to FALSE.}
}
\value{
A data frame with one column containing the observed and (if obs.only=FALSE, the default) another column containing the predicted values from 'model'.
}
\author{
A. Marcia Barbosa
}
\seealso{
\code{\link{prevalence}}
}
\examples{
data(rotif.mods)
mod <- rotif.mods$models[[1]]
obspred <- mod2obspred(mod)
head(obspred)
}
