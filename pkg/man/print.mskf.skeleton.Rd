\name{print.mskf.skeleton}
\alias{print.mskf.skeleton}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{print}{mskf.skeleton}(obj, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ ~~Describe \code{obj} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~

 %~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(obj, ...)
{
    with(obj, {
    cat("\nMarkov switching state space model skeleton\n\n", ...);
    cat("\nNumber of endogenous variables : ", ny, ...);
    cat("\nNumber of exogenous variables  : ", nx, ...);
    cat("\nNumber of latent variables     : ", ne, ...);
    cat("\nNumber of states               : ", nm, ...);
    cat("\nNumber of time points          : ", nt, ...);
    tmp = lapply(const,dim)
    tmp = sapply(tmp, paste, collapse=" x ")
    tmp = paste(c("\n","",""), names(tmp), " (", tmp, ")", sep="", collapse=", ")
    cat("\n\nMatrices used                : ", tmp, ...);
    tmp = lapply(pattern,dim)
    tmp = sapply(tmp, paste, collapse=" x ")
    tmp = paste(c("\n","",""), names(tmp), " (", tmp, ")", sep="", collapse=", ")
    cat("\n\nMatrices with free parameters used : ", tmp, ...);
    })
    invisible(obj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
