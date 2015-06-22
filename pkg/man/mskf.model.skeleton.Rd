\name{mskf.model.skeleton}
\alias{mskf.model.skeleton}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create a model skeleton for the specification of a state space model with
regime switches for time series}
\description{
  This function generates a \code{mskf.skeleton} object to be used with \code{mskf.model}
  for the specification of a state space model with markovian switching between regimes.
}
\usage{ mskf.model.skeleton(y, ne = 1, nm = 2, x = NA, ipat = c(W = 1,\\
    B = nx, R = 0, c = 2, H = 2, G = 1, K = 2, p = 2))
}
\arguments{
  \item{y}{a univariate or multivariate time series (\code{\link{ts}})
  object or vector or matrix object (with rows as time points), possibly with missing values.}
  \item{ne}{number of (underlying/latent) state variables.}
  \item{nm}{number of regimes. Should be 2 or larger.}
  \item{x}{an (optional) uni- or multivariate time series object of predictors for the measurement equation.}
  \item{ipat}{a named numeric vector with entry values 0, 1 or 2. Specifies the model matrices to be used (see Details).}
}
\details{
See \code{\link{mskf}} for complete definition of the regime switching state-space model.

The values in \code{ipat} indicate which model matrices (\code{W},
\code{B}, \code{R}, \code{c}, \code{H}, \code{G}, \code{K}, and
\code{p}; see `Value' below for a description) are used, and whether
or not they contain free parameters: For each named element \code{0}
implies the model matrix is not use; \code{1} implies the
corresponding model matrix contains fixed values only; \code{2}
implies the model matrix contains free parameters.

}
\value{
An object of class \code{mskf.skeleton} with components
 \item{const}{A list wich contains all of the components below for which \code{ipat} is not equal to zero.
                \itemize{
                \item \code{maW} {numeric \code{ny} by \code{ne} by \code{nm} \code{\link{array}}
 for fixed values in the \bold{W} matrices which relate
 the observations to the states.}
                 \item \code{maB} {numeric \code{ny} by \code{nx} by \code{nm} array
 for fixed values in the \bold{B} matrices which relate
 the observations to the predictors.}
                 \item \code{maR} {numeric  \code{ny} by \code{ny} by \code{nm} array
 for fixed values in the \bold{R} matrices that specify the
 (co)variances of residuals in the measurement equation.}
                 \item \code{mac} {numeric   \code{ne} by 1 by \code{nm} array
 for fixed values in the state intercepts \bold{c}.}
                \item \code{maH} {numeric  of \code{ne} by \code{ne} by \code{nm} array
 for fixed values in the \bold{H} matrices which relate the current
 state to the previous state.}
                \item \code{maG} {numeric \code{ne} by \code{ne} by \code{nm} array
 for fixed values in the \bold{G} matrices.}
                \item \code{maK} {numeric   \code{ne} by \code{ne} by \code{nm} array
 for fixed values in the \bold{K} matrices that specify the (co)variances
 of the state innovations.}
                \item \code{map} {numeric \code{nm} by \code{nm} array
 for fixed values in the regime transition matrix \bold{p}.}
                }
 }
 \item{pattern}{A list wich contains all of the components below for which \code{ipat} is equal to two.
                Numeric arrays \code{paW} (\code{ny} by \code{ne} by \code{nm}), \code{paB}
 (\code{ny} by \code{nx} by \code{nm}), \code{paR} (\code{ny} by \code{ny} by \code{nm}),
 \code{pac} (\code{ne} by \code{1} by \code{nm}), \code{paH}
 (\code{ne} by \code{ne} by \code{nm}), \code{paG} (\code{ne}
 by \code{ne} by \code{nm}), \code{paK} (\code{ne} by \code{ne}
 by \code{nm}), and \code{pap} (\code{nm} by \code{nm})
 for the free parameters in the model
 matrices \bold{W}, \bold{B}, \bold{R}, \bold{pac}, \bold{H}, \bold{G}, \bold{K}, and \bold{p} respectively.
                }
  \item{y}{time series \code{y} to be modeled as entered as the first argument.}
  \item{x}{(optional) predictor time series \code{x}.
  If none was present, a \code{1} by \code{1} array with value \code{NA}.}
  \item{ny}{number of time series in \code{y}.}
  \item{nx}{number of time series in \code{x}.}
  \item{nt}{length of the time series in \code{y}.}
  \item{ne}{number of (underlying/latent) state variables.}
  \item{nm}{number of regimes.}
  \item{ipat}{\code{ipat} as entered or its default value.}
}
\references{
Kim, C.-J. (1994). Dynamic linear models with Markov-switching. \emph{Journal of Econometrics}, 60, 1-22.

Kim, C.-J., and Nelson, C. R. (1999). \emph{State-Space Models with Regime Switching:
 Classical and Gibbs-Sampling Approaches with Applications.} Cambridge, MIT press.
}

\author{Ellen L. Hamaker and Raoul P. P. P. Grasman}

\note{ The
object returned by this function is to be used with \code{\link{mskf.model}} for full model specification. }

\seealso{See also \code{\link{mskf}}, \code{\link{mskf.model}} }

\examples{ ##---- Should be DIRECTLY executable !! ---- ##-- ==>
Define data, use random, ##--	or do  help(data=index)  for the
standard data sets.

## The function is currently defined as
function(y, ne=1, nm=2, x = NA, ipat = c(W=1,B=nx,R=0,c=2,H=2,G=1,K=2,p=2))
{
    y = as.matrix(y)
    ny = ncol(y)
    nt = nrow(y)
    nx = if(is.na(x)) 0 else NCOL(x)
    x  = as.matrix(x)

    val.ipat.nms = c("W", "B", "R", "c", "H", "G", "K", "p")
    if(!all(names(ipat) \%in\% val.ipat.nms)||!all(val.ipat.nms \%in\% names(ipat)))
        stop('ipat should have elements named "W", "B", "R", "c", "H", "G", "K", "p"')

	dims = list(
		W = c(ny, ne, nm),
		B = c(ny, nx, nm),
		R = c(ny, ny, nm),
		c = c(ne,  1, nm),
		H = c(ne, ne, nm),
		G = c(ne, ne, nm),
		K = c(ne, ne, nm),
		p = c(nm, nm)
	)
	const = list();
	pattr = list();
	for(nam in names(ipat)){
		pa = paste("pa",nam,sep="")
		ma = paste("ma",nam,sep="")
		if(ipat[nam]>=1)
			const[[ma]] = array(NA, dims[[nam]]);
		if(ipat[nam]==2)
			pattr[[pa]] = array(NA, dims[[nam]]);
	}
	mdl <- list(const = const, pattern = pattr, y=y, x=x, ny=ny, nx=nx,
				nt=nt, ne=ne, nm=nm, ipat=ipat)
	class(mdl) <- "mskf.skeleton"
	mdl
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
