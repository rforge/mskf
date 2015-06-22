\name{mskf-package}
\alias{mskf-package}
\docType{package}
\title{
\packageTitle{mskf}
}
\description{
\packageDescription{mskf}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{mskf}
\packageIndices{mskf}
The mskf package is used for fitting state space models to multivariate time series in which the dynamics changes according to switches between unobserved regimes. The package supports multiple state indicator (endogenous) timeseries, exogenous timeseries, path structure specification, equality and box-constraints, and handles missing data.

The main functions for fitting a model to a multivariate time series (\code{mts}) object are:

\code{mskf.model.skeleton} which for an \code{mts} object returns a model skeleton that contains the necessary model parameter matrices for the indicated dimension of the state space process, the number of regimes, and the (optional) exogenous time series

\code{mskf.model} which turns the model skeleton into a fully specified model by providing influence path structures, fixed values, starting values and box constraints (upper and lower bounds on estimated parameters).

\code{mskf} does the actual fitting of the model.

\code{plot} draws a graphical representation of the state space model under the different regimes

\code{summary} gives model fit summary in terms of likelihood, goodness of fit (AIC), parameter estimates and confidence intervals

\code{simulate} simulates data from a model
}
\author{
\packageAuthor{mskf}

Maintainer: \packageMaintainer{mskf}
}
\references{
Kim, C.J. (1994). Dynamic linear models with Markov-switching. \emph{Journal of Econometrics, 60}, 1–22.

Kim, C.-J., and Nelson, C. R. (1999). \emph{State-Space Models with Regime Switching:
 Classical and Gibbs-Sampling Approaches with Applications.} Cambridge, MIT press.

Hamaker, E. L. & Grasman, R. P. P. P. (2012) Regime Switching State-Space Models Applied To Psychological Processes: Handling Missing Data And Making Inferences. \emph{Psychometrika, 77}:2, 400–422. URL \href{http://doi.org/10.1007/S11336-012-9254-8}{http://doi.org/10.1007/S11336-012-9254-8}
}
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
~~ simple examples of the most important functions ~~
}
