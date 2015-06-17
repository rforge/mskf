\name{mskf-package}
\alias{mskf-package}
\docType{package}
\title{
Markov Switching Kalman filter
}
\description{
This package contains functions to fit Markov governed regime switching time series models.
}
\details{
\tabular{ll}{
Package: \tab mskf\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2008-04-17\cr
License: \tab GPL v2 or higher\cr
}

\section{
The mskf package is used for fitting state space models to multivariate time series in which the dynamics changes according to switches between unobserved regimes. The package supports multiple state indicator (endogenous) timeseries, exogenous timeseries, path structure specification, equality and box-constraints, and handles missing data.

The main functions for fitting a model to a multivariate time series (\code{mts}) object are:

\item{\code{mskf.model.skeleton}} which for an \code{mts} object returns a model skeleton that contains the necessary model parameter matrices for the indicated dimension of the state space process, the number of regimes, and the (optional) exogenous time series

\item{\code{mskf.model}} which turns the model skeleton into a fully specified model by providing influence path structures, fixed values, starting values and box constraints (upper and lower bounds on estimated parameters).

\item \code{mskf} does the actual fitting of the model.

\item \code{plot} draws a graphical representation of the state space model under the different regimes

\item \code{summary} gives model fit summary in terms of likelihood, goodness of fit (AIC), parameter estimates and confidence intervals

\item \code{simulate} simulates data from a model
}
\author{
Ellen L. Hamaker, Raoul P. P. P. Grasman

Maintainer: Ellen L. Hamaker <E.L.Hamaker@uu.nl>
}
\references{
Hamaker, E. L., & Grasman, R. P. P. P. (2012). Regime Switching State-Space Models Applied To Psychological Processes: Handling Missing Data And Making Inferences. Psychometrika, 77:2, 400--422.
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
~~ the R documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
~~ simple examples of the most important functions ~~
}
