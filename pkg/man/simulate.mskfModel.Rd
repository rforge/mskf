\name{simulate.mskfModel}
\alias{simulate.mskfModel}
\alias{sim.msar}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate Responses from \code{mskf} Models
}
\description{
Simulate one or more responses from the distribution corresponding to a fitted model object.
Using the fitted values of the parameters from the \code{mskfModel} object, nsim new data from this model are simulated.
}
\usage{
\method{simulate}{mskfModel}(object, nsim = 1, seed = NA, nt = x$nt, use.start = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object from class "\code{\link{mskfModel}}" representing a Markov switching state space model.
}
  \item{nsim}{
an optional positive integer specifying the number of simulations to perform. Deafults to 1.
}
  \item{seed}{
an optional positive integer that is passed to \code{set.seed}. If set, the returned object will have a "\code{seed}" attribute with this value.
}
  \item{nt}{
an optional positive integer that specifies the length (number of samples) of the simulated time series.
}
  \item{use.start}{
Boolean. Should starting values be used in stead of estimated parameter values?
}
  \item{\dots}{
optional additional arguments. None are used.
}
}
\details{
Currently only (univariate) AR(p) processes are implemented.
}
\value{
List
\item{y }{simulated time series}
\item{a}{object returned by \code{sim.msar} from which the time series is constructed}
\item{regime}{vector of length \code{length(y)} with the states/regimes corresponding to the elements in \code{y}}

The object has a \code{seed} attribute which is either \code{NA} or a positive integer (see \code{seed} parameter).
}
\references{
Hamaker, E. L. & Grasman, R. P. P. P. (2012) Regime Switching State-Space Models Applied To Psychological Processes: Handling Missing Data And Making Inferences. \emph{Psychometrika, 77}:2, 400???422. URL \href{http://doi.org/10.1007/S11336-012-9254-8}{http://doi.org/10.1007/S11336-012-9254-8}
}
\author{
Raoul P. P. P. Grasman
}


\seealso{
	\code{\link{mskf}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
