\name{mskf.model}
\alias{mskf.model}
\alias{print.mskf.model}

\title{ Specify a Markov regime switching state-space model }
\description{
  Specify a Markov regime switching state-space model for uni- or multivariate time series.
}

\usage{
    \method{mskf.model}{mskf.skeleton}(x, ...)
}

\arguments{
  \item{x}{ A \code{\link{mskf.skeleton}} object. }
  \item{\dots}{ Specifications of the model vectors and matrices. See Details below. }
}
\details{
    The \code{x} object is a basic skeleton for a Markov regime switching state-space model for uni- or multivariate time series, the details of which are specified through this function. The time serie(s) is (are) assumed to be (an) indicator(s) for an underlying, (latent) stochastic state process that is described by a state space (vector autoregressive; VAR) model. Both uni- and multi-dimensional latent state spaces are supported. The dynamics of the latent series is assumed to switch between a discrete set of unobserved regimes, each associated with their own set of parameters. Switching between regimes is modeled with a Hidden Markov Model.

    Let \code{ny} be the number of (endogenous) observed time series, \code{ne} the dimension of the (latent) state process, \code{nm} the number of regimes, and \code{nx} be the number of exogonous time series. The regime switching state-space model consists of the measurement equation and the state transition equation, and a discrete Markov Chain.

    The measurement equation is given by
    \deqn{
        \mathbf{y}_t = \mathbf{W}_{S_t}\mathbf{a}_t + \mathbf{B}_{S_t}\mathbf{x}_t+\mathbf{e}_t
    }{%
        y[t] = W(S[t]) a[t] + B(S[t]) x[t] + e[t]
    }
    where \eqn{\mathbf{y}_t}{y[t]} is an \code{ny}-vector measurement made at occasion \eqn{t}; \eqn{\mathbf{a}_t}{a[t]} is an \code{ne}-vector (latent) state vector for occasion \eqn{t}; \eqn{\mathbf{x}_t}{x[t]} is an \code{nx}-vector with exogenous measurements; \eqn{W_{St}}{W(S[t])} is the \code{ny} by \code{ne} matrix with regression coefficients for predicting the observations from the states; \eqn{\mathbf{B}_{S_t}}{B(S[t])} is an \code{ny} by \code{nx} matrix with regression coefficients for \eqn{\amthbf{y}_t}{y[t]} from the exogenous measurements; and \eqn{\mathbf{e}_t \sim N(\mathbf{0}, \mathbf{R}_{S_t})}{e[t] ~ N(0, R(S[t]))}, \eqn{t = \ldots, -1, 0, 1, \ldots}, are \code{ny}-dimensional i.i.d. vectors.

    The transition equation for the states is given by
    \deqn{
      \mathbf{a}_t = \mathbf{c}_{S_t}  +  \mathbf{H}_{S_t} \mathbf{a}_{t-1}  +  \mathbf{G}_{S_t} \mathbf{u}_t
    }{%
            a[t] = c(S[t]) + H(S[t]) a[t-1] + G(S[t]) u[t]
    }
    where \eqn{\mathbf{c}_{S_t}}{c(S[t])} is a vector with \code{ne} intercepts; \eqn{\mathbf{H}_{S_t}}{H(S[t])} is a \code{ne} by \code{ne} matrix of vector-autoregression coefficients; \eqn{\mathbf{G}_{S_t}}{G(S[t])} is an \code{ne} by \code{ne} design matrix; and \eqn{\mathbf{u}_t \sim N(\mathbf{0},\mathbf{K}_{S_t})}{u[t] ~ N(0, K(S_t))}, \eqn{t = \ldots, -1, 0, 1, \ldots}, are i.i.d. innovations.

    The model parameters \eqn{\mathbf{W}_{S_t}}{W(S[t])}, \eqn{\mathbf{B}_{S_t}}{B(S[t])}, \eqn{\mathbf{H}_{S_t}}{H(S[t])}, \eqn{\mathbf{R}_{S_t}}{R(S[t])}, \eqn{\mathbf{G}_{S_t}}{G(S[t])}, \eqn{\mathbf{K}_{S_t}}{K(S[t])}, and \eqn{\mathbf{c}_{S_t}}{c(S[t])}, all depend on the current regime \eqn{S_t}{S[t]}. The regimes swith from one to the other according to
    \deqn{
        P(S[t] = s) = \sum_{s'} \pi_{s,s'} P(S[t-1] = s').
    }{%
        P(S[t] = s) = \sum \pi(s,s') P(S[t-1] = s').
    }
    The regimes are labeled with integers from 1 to \code{nm}.

    The parameters \eqn{\mathbf{B}}{B}, \eqn{\mathbf{c}}{c}, and \eqn{\mathbf{G}}{G} are optional. For the model to be identifiable, the parameter matrices/vectors need to be structured/constraint. All parameters can be constraint to a fixed value, or can be set to be estimated freely. Individual parameters can be fixed or freed through the \dQuote{patterns} and \dQuote{matrix} specifications of the model matrices/vectors. To fix model matrix elements to a specific value, the \code{maW}, \code{mac}, \code{maH}, \code{maG}, \code{maK}, \code{map}, arguments are passed with that value. To free model matrix elements to be freely estimated, the \code{paW}, \code{pac}, \code{paH}, \code{paG}, \code{paK}, \code{pap}, arguments are passed specifying (unique) names for the elements. For example, in a two regime model, if \eqn{W} is to be the identity matrix (e.g., to equate the states with the measurements; see the example below), this is specified by
    \verb{
        mdl = mskf.model(x, maW = diag(ny)) # set H to fixed value diag(ny) for each regime
    }
    If in addition, \eqn{c} and \eqn{H} are to be estimated, this is specified by
    \verb{
        mdl = mskf.model(x, maW = diag(ny),
            pac = 1:4, paH = array(5:12,c(2,2,2))) # elements of c and H freely estimated
    }
    Names in \code{paW}, \code{pac}, \code{paH}, \code{paG}, \code{paK}, \code{pap} argument values must be integers. Entries in these arguments that are equal to 0 are fixed and not estimated (their values are obtained from \code{maW}, \code{mac}, \code{maH}, \code{maG}, \code{maK}, \code{map}, respectively). Non-unique names in these arguments are treated as equality constraints. For instance, if in the above example the \eqn{c}'s should be equal for the two regimes, this is specified by
    \verb{
        mdl = mskf.model(x, maW = diag(ny),
            pac = 1:2, paH = array(5:12,c(2,2,2))) # intercept c same across regimes
    }
    where the names provided through \code{pac} (1 and 2) are automatically recycled untill al the \eqn{c} elements have a name. If there shouldn't be any cross-laged influences between the state variables in the example above, this can be specified with
    \verb{
        mdl = mskf.model(x, maW = diag(ny),
            pac = 1:2, paH = c(diag(5:6),diag(7:8))) # no cross-lagged influences
    }

    Non-zero entries in pattern specifications always supersede matrix value specifications of the corresponding entries.

    Parameters can be further constraint with box constraints by passing the \code{lower} and/or \code{upper} arguments a named vector of lower and upperbounds, where the names are the names is in the pattern specifications.

    Patterns left unspecified are automatically given all unique names.

    Starting values can be set by passing the \code{start} argument a vector in the same way as \code{lower} and \code{upper}.
}
\value{
  Object of class \code{mskf.model}.
}
\references{ \cite{Hamaker, E.L., & Grasman, R.P.P.P. (2011) ..., Psychometrika, ?, ??--??.} }
\author{ Ellen Hamaker & Raoul Grasman }
%\note{ ~~further notes~~
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{ \code{mskf}, \code{\link{paW}}, \code{\link{pac}}, \code{\link{paH}}, \code{\link{paG}}, \code{\link{paK}}, \code{\link{pap}}, \code{\link{maW}}, \code{\link{mac}}, \code{\link{maH}}, \code{\link{maG}}, \code{\link{maK}}, \code{\link{map}}, \code{\link[=start.mskf.model]{start}}, \code{\link[=plot.mskf.model]{plot}} }
\examples{
    # simulate a time series
    y = ts(matrix(rnorm(200), ncol = 2));

    ##
    ## build the model
    ##

    ## first a skeleton
    skel = mskf.model.skeleton(y, ne = 2, nm = 2) # 2 latent state variables, 3 regimes

    ## specify the details
    mdl = mskf.model(skel)  # default: y[t] = a[t], same intercept across regimes,
    print(mdl)              #          no cross-lagged influences, unit variance innovations,
                            #          equal switching probabilities

    # freely estimate c across regimes
    mdl = mskf.model(skel, pac = c(1, 2, 3, 4))
    print(mdl)

    # c free across regimes, a cross-lagged influence in one direction
    mdl = mskf.model(skel, pac = c(1, 2, 3, 4), paH = c(5,6,0,7,8,9,0,10))
    print(mdl)
    plot(mdl)

    # set (some) starting values
    mdl = mskf.model(skel, pac = c(1, 2, 3, 4), paH = c(5,6,0,7,8,9,0,10),
            start = c('1'=1, '2'=1, '3'=2, '4'=2, '5'=0.7, '8'=0.7))

    ##
    ## changing the model
    ##

    # restrict auto-regressions to be same across regimes
    paH(mdl) = c(5,6,0,7)
    print(mdl)

    # no cross-lagged influence in the 2nd regime
    paH(mdl)[1,2, "regime 2"] = 0
    plot(mdl)

    # change (some) upper and lower bounds
    upper(mdl) = c('1' = 1, '3'= 1)
    upper(mdl, "paH") = array(1, dim(paH(mdl)))
    lower(mdl, "paH") = array(-1, dim(paH(mdl)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
