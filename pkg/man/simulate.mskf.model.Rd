\name{simulate.mskf.model}
\alias{simulate.mskf.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate Responses from \code{mskf} Models
}
\description{
Simulate one or more responses from the distribution corresponding to a fitted model object.
Using the fitted values of the parameters from the \code{mskf.model} object, nsim new data from this model are simulated.
}
\usage{
simulate.mskf.model(object, nsim = 1, seed = NA, nt = x$nt, use.start = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object from class "\code{\link{mskf.mode}}" representing a Markov switching state space model.
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
%% ~put references to the literature/web site here ~
}
\author{
Raoul P. P. P. Grasman
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, nsim = 1, seed = NA, nt = x$nt, use.start = FALSE,
    ...)
{
    getStart = function(x, pa) na.omit(start(x)[as.character(pa)])
    x = object
    c = mac(x)
    H = maH(x)
    p = map(x)
    G = maG(x)
    W = maW(x)
    B = maB(x)
    R = maR(x)
    if (use.start) {
        c[pac(x) != 0] = getStart(x, pac(x))
        H[paH(x) != 0] = getStart(x, paH(x))
        p[pap(x) != 0] = getStart(x, pap(x))
        G[paG(x) != 0] = getStart(x, paG(x))
        W[paW(x) != 0] = getStart(x, paW(x))
        B[paB(x) != 0] = getStart(x, paB(x))
        R[paR(x) != 0] = getStart(x, paR(x))
    }
    X = x$x
    sqrtR = R = if (is.null(R))
        array(diag(x$ny), c(x$ny, x$ny, x$nm))
    else R
    for (i in 1:x$ne) {
        sqrtR[, , i] = chol(R[, , i])
    }
    if (!is.na(seed)) {
        set.seed(seed)
    }
    r = replicate(nsim, {
        a = sim.msar(nt, p = p, H = H, m = c, G = G)
        y = array(NA, dim = c(nt, x$ny))
        e = matrix(rnorm(nt * x$ny), nt)
        for (i in 1:x$nm) {
            y[a$regime == i, ] = a$y[a$regime == i, , drop = FALSE] \%*\%
                t(array(W[, , i], dim(W)[-3])) + e[a$regime ==
                i, , drop = FALSE] \%*\% sqrtR[, , i]
            if (!is.null(X) && !is.null(B)) {
                ind = a$regime == i
                y[ind, ] = y[ind, ] + X[ind, ] \%*\% t(array(B[,
                  , i], dim(B)[-3]))
            }
        }
        list(y = as.ts(y), a = a$y, regime = a$regime)
    })
    attr(r, "seed") = seed
    drop(r)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
