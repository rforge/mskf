\name{mskfModel.mskfSkeleton}
\alias{mskfModel.mskfSkeleton}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generate an \code{mskfModel} from an \code{mskfSkeleton} }
\description{
Generate an \code{mskfModel} from an \code{mskfSkeleton} specifying values for the different model parameters, starting values, and upper and or lower bounds for the parameter values.
}
\usage{
\method{mskfModel}{mskfSkeleton}(object, ..., start = NA, lower = NA, upper = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object returned by \code{mskfSkeleton}}
  \item{\dots}{Extra arguments which can be zero or more of
  the arrays \code{maW}, \code{maB}, \code{maR}, \code{mac},
   \code{maH}, \code{maG}, and \code{maK}, \code{map}, \code{paW}, \code{paB}, \code{paR}, \code{pac},
   \code{paH}, \code{paG}, and \code{paK}, and/or the matrix \code{pap} (see Details).
    }
  \item{start}{a named numeric vector with starting values for the freely estimated parameters; see details}
  \item{lower}{a named vector, compatible with \code{start}, of lower bounds for the freely estimated parameters}
  \item{upper}{a named vector, compatible with \code{start}, of upper bounds for the freely estimated parameters}
}
\details{
 See \code{\link{mskf}} for precise definition of the regime switching state-space model.
Arrays and matrix may be specified as strings.
Can be used to name the free parameters using numbers, and to
constrain free parameters to be equal (by using the same
number). Parameter estimates will be presented in the same
order as they were numbered. Numbering does not need to be
consecutively.

 The vector \code{start} can be used to initialize the parameters at specific values. The names in this
 vector are the parameter names as given for the model matrices. As an example, if \code{paW} is specified as
 \code{paW=diag(1:2)}, and \code{paH} is specified as \code{paH=2}, then if \code{start=c('1'=0.25, '2'=0.8)},
 then \code{paW} will be initialized to \code{diag(c(0.25, 0.8))} and \code{paH} will be initialized to
 \code{matrix(0.8,1,1)}.

}
\value{
An object of class \code{mskfSkeleton} with components
 \item{const}{A list wich contains all of the components below for which \code{ipat} is not equal to zero.
                \itemize{
                \item \code{maW} {numeric \code{ny} by \code{ne} by \code{nm} \code{\link{array}}
 with fixed values in the \bold{W} matrices which relate
 the observations to the states.}
                 \item \code{maB} {numeric \code{ny} by \code{nx} by \code{nm} array
 with fixed values in the \bold{B} matrices which relate
 the observations to the predictors.}
                 \item \code{maR} {numeric  \code{ny} by \code{ny} by \code{nm} array
 with fixed values in the \bold{R} matrices that specify the
 (co)variances of residuals in the measurement equation.}
                 \item \code{mac} {numeric   \code{ne} by 1 by \code{nm} array
 with fixed values in the state intercepts \bold{c}.}
                \item \code{maH} {numeric  of \code{ne} by \code{ne} by \code{nm} array
 with fixed values in the \bold{H} matrices which relate the current
 state to the previous state.}
                \item \code{maG} {numeric \code{ne} by \code{ne} by \code{nm} array
 with fixed values in the \bold{G} matrices.}
                \item \code{maK} {numeric   \code{ne} by \code{ne} by \code{nm} array
 with fixed values in the \bold{K} matrices that specify the (co)variances
 of the state innovations.}
                \item \code{map} {numeric \code{nm} by \code{nm} array
 with fixed values in the regime transition matrix \bold{p}.}
                }
 }
 \item{pattern}{A list wich contains all of the components below for which \code{ipat} is equal to two.
                Numeric arrays \code{paW} (\code{ny} by \code{ne} by \code{nm}), \code{paB}
                 (\code{ny} by \code{nx} by \code{nm}), \code{paR} (\code{ny} by \code{ny} by \code{nm}),
                  \code{pac} (\code{ne} by \code{1} by \code{nm}), \code{paH}
                    (\code{ne} by \code{ne} by \code{nm}), \code{paG}
                     (\code{ne} by \code{ne} by \code{nm}), \code{paK}
                      (\code{ne} by \code{ne} by \code{nm}), \code{pap}
                       (\code{nm} by \code{nm}) that indicate which parameters
                        are to be estimated in the model matrices
                        \bold{W}, \bold{B}, \bold{R}, \bold{pac}, \bold{H}, \bold{G}, \bold{K},
                        and \bold{p} respectively.
                }
  \item{y}{time series \code{y} to be modeled as entered as the first argument.}
  \item{x}{(optional) predictor time series \code{x}.
  If none was present, a \code{1} by \code{1} array with value \code{NA}.} # MOET JE DIT NOG STEEDS ZO ZEGGEN? HET IS IN DE VORIGE FUNCTIE NAMELIJK AL GEMAAKT
  \item{ny}{number of time series in \code{y}.}
  \item{nx}{number of time series in \code{x}.}
  \item{nt}{length of the time series in \code{y}.}
  \item{ne}{number of state variables.}
  \item{nm}{number of regimes.}
  \item{ipat}{\code{ipat} as entered or its default value.} # MOET DAT HIER WEER ZO WORDEN GEZEGD? HIJ MAAKT HIER GEEN DEFAULTS MEER NAMELIJK
  \item{theta}{a named vector with starting values with the same names as \code{start}}
  \item{lobo}{a vector with lower bound values}
  \item{upbo}{a vector with upper bound values}
  \item{a0}{a vector with ne zeros as the initial state vector}
  \item{P0}{a matrix of ne by ne with 100 on the diagonal, which serves as the initial mean squared error matrix.}
}

\references{
Kim, C.-J. (1994). Dynamic linear models with Markov-switching. \emph{Journal of Econometrics}, 60, 1-22.

Kim, C.-J., and Nelson, C. R. (1999). \emph{State-Space Models with Regime Switching:
 Classical and Gibbs-Sampling Approaches with Applications.} Cambridge, MIT press.
}

\author{Ellen L. Hamaker and Raoul P. P. P. Grasman} \note{ The
object returned by this function is to be used with
\code{\link{mskfModel}} for full model specification. }

\seealso{See also \code{\link{mskf}}, \code{\link{mskfModel}} }

\examples{
1:10
}
