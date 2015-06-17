\name{mskf.model.mskf.skeleton}
\alias{mskf.model.mskf.skeleton}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mskf.model.mskf.skeleton(mskf.skeleton, ..., start = NA, lower = NA, upper = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mskf.skeleton}{Object returned by mskf.skeleton} 
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
An object of class \code{mskf.skeleton} with components
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
\code{\link{mskf.model}} for full model specification. } 

\seealso{See also \code{\link{mskf}}, \code{\link{mskf.model}} } 

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(mskf.skeleton, ..., start=NA, lower=NA, upper=NA)
{
    # 
    obj = mskf.skeleton
    inits = list(...);
    nmconst = names(obj$const);
    nmpattr = names(obj$pattern);
    if(!all(chk <- names(inits) \%in\% c(nmconst, nmpattr)))
        warning("arguments ", names(inits)[!chk], " unused");

    # "attach" elements of x to the current environment
    for(nme in names(obj)) assign(nme, obj[[nme]]);

    # default values
    defc = list(
        maW = diag(max(ny,ne))[1:ny, 1:ne],
        maB = diag(max(ny,nx))[1:ny, 1:nx],
        maR = diag(ny),
        mac = rep(0, ne),
        maH = diag(ne),
        maG = diag(ne),
        maK = diag(ne),
        map = diag(1, nm)
    )
    npar = 0;
    defp = list(
        paW = {tmp = matrix(npar+1,ny,ne); tmp[upper.tri(tmp)] = 0; npar<-max(tmp,npar); tmp},
        paB = {tmp = matrix(npar+1,ny,nx); tmp[upper.tri(tmp)] = 0; npar<-max(tmp,npar); tmp},
        paR = {tmp = diag(npar+1,ny); npar <- max(tmp,npar); tmp;},
        pac = {tmp = matrix(npar+1:ne,ne, 1); npar<-max(tmp,npar); tmp},
        paH = {tmp = diag(npar+1:ne,ne); npar<-max(tmp,npar); tmp},
        paG = {tmp = diag(npar+1:ne,ne); npar<-max(tmp,npar); tmp},
        paK = {tmp = diag(npar+1,ne); npar<-max(tmp,npar); tmp},
        pap = {tmp = matrix(npar+1,ne,ne); npar<-max(tmp,npar); tmp}
     )
        
    
    # vectors for starting values, lower and upperbounds
    theta <- lobo <- upbo <- c()
    
    # set values of skeleton matrices
    for(nme in nmconst){
        obj$const[[nme]][] = defc[[nme]];
        if(nme \%in\% names(inits))
            obj$const[[nme]][] = inits[[nme]];
    }
    for(nme in nmpattr){
        if(length(obj$pattern[[nme]][]) \%\% length(defp[[nme]]) !=0)
            stop("this shouldn't happen, you discovered a bug (", nme,")")
        else 
            obj$pattern[[nme]][] = defp[[nme]]
        if(nme \%in\% names(inits))
            if(length(obj$pattern[[nme]][]) \%\% length(inits[[nme]]) !=0)
                stop("length of specification of ", nme, " is not a multiple of this model matrix")
            else
                obj$pattern[[nme]][] = inits[[nme]]
    }

    # create named starting parameters vector theta
    theta.nms = na.omit(unique(unlist(obj$pattern)))
    theta.nms = theta.nms[as.character(theta.nms)!="0"]
    npar = length(theta.nms)
    theta = rep(1/npar, npar) 
    names(theta) = theta.nms
    # check provided starting values and set corresponding elements of theta to values provided
    if(!all(chk <- names(start) \%in\% theta.nms))
        warning("start value(s) named ", names(start)[!chk], " not used")
    theta[names(start)[chk]] = start[chk] 

    # construct lower and upper bounds similarly (treat variance parameters special)
    lobo = rep(-10, npar)
    names(lobo) = theta.nms
    diagR = rep(diag(ny)>0,nm);
    diagK = rep(diag(ne)>0,nm);
    varpars = na.omit(unique(c(obj$pattern$paR[diagR], obj$pattern$paK[diagK])))
    varpars = varpars[as.character(varpars)!="0"]
    lobo[as.character(varpars)] = sqrt(.Machine$double.eps) # small number away from zero
    if(!all(chk <- names(lower) \%in\% theta.nms))
        warning("lower bound(s) named ", names(lower)[!chk], " not used")
    lobo[names(lower)[chk]] = lower[chk] 
    
    upbo = rep(10, npar)
    names(upbo) = theta.nms
    if(!all(chk <- names(upper) \%in\% theta.nms))
        warning("upper bound(s) named ", names(upper)[!chk], " not used")
    lobo[names(upper)[chk]] = upper[chk] 

    
    # initial state and its covariance matrix
    a0 = rep(0, ne)           # LET OP: DIT MOET JE EIGENLIJK VOOR ALLE REGIMES EEN KEER OPGEVEN
    ##a0 = array(a0, c(ne,  1, nm))   # doen dan dus maar? lijkt een probleem te veroorzaken
    P0 = diag(100, ne)        # IDEM
    ##P0 = array(P0, c(ne, ne, nm))   # ook maar doen dan?? lijkt een probleem te veroorzaken

    obj = c(obj, list(theta = theta, lobo = lobo, upbo = upbo, a0 = a0, P0 = P0));
    class(obj) = 'mskf.model'
    obj
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
