\name{mskf}
\alias{mskf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit A Markov Switching State-Space Model} % DIT LIJKT ME GEEN GOEDE TITEL, MAAR WEET NIET WAT ANDERS: ZONDER DAT ESTIMATING A?
\description{
  Fit a regime switching state-space model to a (multivariate) time series.
}
\usage{
mskf(model, method = "L-BFGS-B", hessian = TRUE, control = list(trace = 3, maxit = 500),
    debug = FALSE)
}
\arguments{
  \item{model}{Object returned by \code{\link{mskfModel}}}
  \item{method}{Optimization method to be used by \code{\link{optim}}}
  \item{hessian}{Logical. Should a numerically differentiated Hessian matrix be returned?}
  \item{control}{A list of control parameters. See Details in \code{\link{optim}}}
  \item{debug}{ Logical. For debugging purposes. }
}
\details{
The \code{model} object describes a Markov regime switching state-space model for uni- or multivariate time series. The multivariate time series are assumed to be indicators for an underlying stochastic time series that is described by a state space model. The dynamics of the underlying series is assumed to switch between a discrete set of states, or regimes, each associated with their own set of parameters. Switching between regimes is modeled as a discrete time Markov chain. See \code{\link{mskfModel}} for a detail model specification. The (unconstraint) parameters and switching probabilities are estimated by means of maximum likelihood.

Missing values are treated as discussed in the last reference.
}

\value{
An object of class \code{mskf} with components returned by \code{\link{optim}}.
}

\references{
\cite{Kim, C.-J. (1994). Dynamic linear models with Markov-switching. \emph{Journal of Econometrics}, 60, 1-22.}

\cite{Kim, C.-J., and Nelson, C. R. (1999). \emph{State-Space Models with Regime Switching:
 Classical and Gibbs-Sampling Approaches with Applications.} Cambridge, MIT press.}

\cite{Hamaker, E. L., & Grasman, R. P. P. P. (2012). Regime Switching State-Space Models Applied To Psychological Processes: Handling Missing Data And Making Inferences. Psychometrika, 77:2, 400--422.}
}

\author{Ellen L. Hamaker and Raoul P. P. P. Grasman}

%\note{ ~~further notes~~
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}

\seealso{ See Also as \code{\link{mskfModel}}, \code{\link{mskfModelSkeleton}}, \code{\link{optim}}  }% HIER VERWIJZEN NAAR OPTIM?

%    # Univariate AR(1) with two regimes
%
%    # parameters
%    phi01 = 0       # constant in regime 1
%    phi02 = 3       # constant in regime 2
%    phi11 = 0.5     # AR(1) parameter in regime 1
%    phi12 = 0.3     # AR(1) parameter in regime 2
%    rsd1  = 1       # residual standard deviation in regime 1
%    rsd2  = 1       # residual standard deviation in regime 2
%
%    # number of time point and burn in
%    nt = 150
%    burn.in = 100
%
%    # Markov transition matrix
%    p = matrix(c(.7,.3,.2,.8), 2, 2)
%    lp = eigen(p)$vectors[,1]
%    lp = lp/sum(lp) # equilibrium regime distribution
%
%    true.par = c(phi01, phi02, phi11, phi12, rsd1, rsd2, log(p[1,1]/p[2,1]),
%    	log(p[2,2]/p[1,2]))
%
%    # Create a discrete Markov series
%    r = sample(1:2, 1, p = lp);
%    for (t in 2:(nt+burn.in))
%    	r[t] = sample(1:2, 1, prob = p[,r[t-1]]);
%    y = rnorm(1)
%    for (t in 2:(nt+burn.in))
%    	y[t]  = if (r[t]==1)  phi01 + phi11 * y[t-1] + rnorm(1,0,rsd1) else
%    		phi02 + phi12 * y[t-1] + rnorm(1,0,rsd2)
%    y <- y[-(1:burn.in)]
%
%
%    # Specify the model
%    skel = mskf.model.skeleton(y, ne=1, nm=2)
%    start = c(`1` = .05, `2` = 3, `3` = .51, `4` = .22, `5` = 1.5, `6` = 2,
%    		`7` = 1, `8` = .5)
%    mod = mskf.model(skel, pap = diag(7:8), paK = 5:6, pac = 1:2, paH = 3:4, start = start,
%    	upper = 50, lower = -50)
%
%    # fit the model
%    fitres <- mskf(mod,control=list(trace=3,maxit=2000),debug=0)


\examples{
### Univariate AR(1) with two regimes

c = array(c(0,3), dim = c(1,2))
H = array(c(0.5, 0.3), dim = c(1,1,2))
G = array(c(1, 1), dim = c(1,1,2))

p = matrix(c(0.7, 0.3, 0.2, 0.8), 2, 2)

# specify model
skel = mskfModelSkeleton(rnorm(100), ne = 1, nm = 2)
mdl = mskfModel(skel, map = p, maH = H, maG = G)

# simulate from model
sim = simulate(mdl, nt=100)
y = sim$y;

## estimate model parameters

# set up model skeleton for simulated data
skel2 = mskfModelSkeleton(y, ne = 1, nm = 2)

# build model and initialize starting values
mdl2 = mskfModel(skel2, pap = diag(7:8), paK = 5:6, pac=1:2, paH=3:4, upper = 50, lower = -50)
startValues(mdl2, 'pac') = c(`1` = .05, `2` = 3)
startValues(mdl2, 'paH') = c(`3` = .51, `4` = .22)
startValues(mdl2, 'paK') = c(`5` = 1.5, `6` = 2)
startValues(mdl2, 'pap') = diag(c(1,.5))

# fit the model
mskf(mdl2)




### Univeriate AR(1) with two regimes and 4 indicators

B = array(c(10), dim = c(1, 1, 2))
W = array(rep(c(1, 0.8, 0.6, 0.4), 2), dim = c(4, 1, 2))
R = array(c(diag(1,4), diag(1.2,4)), dim = c(4,4,2))

# simulate from model
y = matrix(rnorm(100*4),100) # data place holder
x = rep(1, 100)
skel3 = mskfModelSkeleton(y, ne = 1, nm = 2, x=x, ipat=c(W=2,B=2,R=2,H=2,c=2,G=1,K=2,p=2))
mdl3 = mskfModel(skel3, map = p, mac = c, maH = H, maB=B, maG = G, maW = W, maR = R)
sim = simulate(mdl3)
y = sim$y

# set up model skeleton for simulated data
skel4 = mskfModelSkeleton(y, ne = 1, nm = 2, x=x, ipat=c(W=2,B=2,R=2,H=2,c=2,G=1,K=2,p=2))

# build model and initialize starting values
mdl4 = mskfModel(skel4, paW = rep(0:3,2), paB = 4:7, paR = c(diag(8:11),diag(12:15)),
	paH = 16:17, pac = c(0,18),paK = 19:20, pap = diag(21:22))

startValues(mdl4, 'paW')[] = c(0,.8, .6, .4)
startValues(mdl4, 'paB')[] = rep(10, 4)
startValues(mdl4, 'paR')[] = rep(1, 8)
startValues(mdl4, 'paH')[] = c(.5, .2)
startValues(mdl4, 'pac')[] = 3
startValues(mdl4, 'paK')[] = c(1,1)
startValues(mdl4, 'pap')[] = c(.8, 1.4)


}
