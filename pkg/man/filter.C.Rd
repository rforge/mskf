\name{filter.C}
\alias{filter.C}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
filter.C(theta, ny, ne, nx, nt, nm, mdim, MA, PA, a0, P0, y, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{ ~~Describe \code{theta} here~~ }
  \item{ny}{ ~~Describe \code{ny} here~~ }
  \item{ne}{ ~~Describe \code{ne} here~~ }
  \item{nx}{ ~~Describe \code{nx} here~~ }
  \item{nt}{ ~~Describe \code{nt} here~~ }
  \item{nm}{ ~~Describe \code{nm} here~~ }
  \item{mdim}{ ~~Describe \code{mdim} here~~ }
  \item{MA}{ ~~Describe \code{MA} here~~ }
  \item{PA}{ ~~Describe \code{PA} here~~ }
  \item{a0}{ ~~Describe \code{a0} here~~ }
  \item{P0}{ ~~Describe \code{P0} here~~ }
  \item{y}{ ~~Describe \code{y} here~~ }
  \item{x}{ ~~Describe \code{x} here~~ }
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

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(theta,ny,ne,nx,nt,nm,mdim,MA,PA,a0,P0,y,x)
{	# START WITH MAKING THE ARRAYS THAT ARE NEEDED
	auc<-array(,c((nt+1),ne,nm))		# a-updated and collapsed to nm posteriors
	au<-array(,c((nt),ne,nm,nm))		# a-updated (not collapsed: S(t-1)=i, S(t)=j)
	ap<-array(,c((nt),ne,nm,nm))		# a-predicted (by definition not collapsed)
	Puc<-array(,c((nt+1),ne,ne,nm))	# covariance matrix of auc
	Pu<-array(,c(nt,ne,ne,nm,nm))		# covariance matrix of au
	Pp<-array(,c(nt,ne,ne,nm,nm))		# covariance matrix of ap
	v<-array(,c(nt,ny,nm,nm))		# one-step-ahead prediction error for S(t-1)=1, S(t)=j
	F<-array(,c(nt,ny,ny,nm,nm))		# covariance matrix of v
	pr<-array(,c((nt+1),nm))		# updated probabilities
	tpr<-array(,c(nt,nm,nm))		# conditional transition probabilites
	utpr<-array(,c(nt,nm,nm))		# updated conditional transition probabilites
	jd<-array(,c(nt,nm,nm))			# joint densities (product of conditional density and transition probability)
	md<-array(,c(nt,1))			# marginal density (used in optimization)
	# MAKE UPDATED PROABILITIES (I.E., VECTOR WITH PROBABILITIES OF
	# OF BEING IN CERTAIN STATE AT OCCASION 0
	pr[1,]<-matrix(c(1/nm),nm,1)
 	# ENTER INITIAL VALUES OF a0 AND P0 IN auc
	auc[1,,]<-a0
	Puc[1,,,]<-P0
	# PLACE THE ELEMENTS FROM x (= vector with free parameters) 
 	# IN THE APPROPRIATE MODEL MATRICES
	for (m in 1:7)
	{	for (s in 1:nm)
		{	for (r in 1:mdim[m,1])
			{	for (c in 1:mdim[m,2])
				{	p<-PA[m,r,c,s] 
					if (p>0) { MA[m,r,c,s]<-theta[p] }
				}
			}
		}
	}
	# TO OBTAIN THE MATRIX WITH TRANSITION PROBABILITIES, 
	# WE HAVE TO TAKE INTO ACCOUNT THAT THE PORBS HAVE TO SUM TO 1 PER COLUMN
	# HENCE, PER COLUMN THERE ARE nm-1 PROBS TO BE ESTIMATED
	# RATHER THAN ESTIMATING THE PROBS WE WILL ESTIMATE A NUMBER
	# WHICH CAN RUN FROM -INFINITY TO +INFINITY
	tp<-matrix(,nm,nm)
	for (c in 1:(nm-1))
	{	tp[nm,c]<-exp(0)
		for (r in 1:(nm-1))
		{	npar<-npar+1
			tp[r,c]<-exp(theta[npar])
		}
	}	
	tp[1,nm]<-exp(0)
	for (r in 2:nm)
	{	npar<-npar+1
		tp[r,nm]<-exp(theta[npar])
	}
	# THE ZEROS ARE INCLUDED FOR IDENTIFICATION
	# TO COME FROM THE ABOVE MATRIX tp TO THE MATRIX WITH ACTUAL
	# TRANSITION PROBS p, WE DO THE FOLLOWING 
	p<-matrix(,nm,nm)
	for (c in 1:nm)
	{	stp<-sum(tp[,c])
		p[,c]<-tp[,c]/stp
	}
	# SELECT THE MODEL MATRICES AND NAME THEM
	W<-MA[1,1:ny,1:ne,,drop=FALSE]
	B<-MA[2,1:ny,1:nx,,drop=FALSE]
	R<-MA[3,1:ny,1:ny,,drop=FALSE]
	c<-MA[4,1:ne,1,,drop=FALSE]		
	H<-MA[5,1:ne,1:ne,,drop=FALSE]
	G<-MA[6,1:ne,1:ne,,drop=FALSE]
	K<-MA[7,1:ne,1:ne,,drop=FALSE]
	# START THE KALMAN FILTER & REGIME SWITCHING
	L<-0
    L <- Kfilter.timeloop(   
            auc, au, ap,
            Puc, Pu, Pp,
            v, F, p,
            pr, tpr, utpr, jd, md,
            W, B, R, c, H, G, K, y, x)

	# for (t in 1:nt)
	# {	# DETERMINE WHETHER THERE ARE MISSING VALUES
		# if(all(!y[t,]==-999))
		# {	for (j in 1:nm)
			# {	for (i in 1:nm)
				# {	# KALMAN FILTER STEP
					## PREDICTION EQUATIONS
					# ap[t,,i,j]<-c[1,,1,j]+H[1,,,j]\%*\%(auc[t,,i])
					# Pp[t,,,i,j]<-H[1,,,j]\%*\%Puc[t,,,i]\%*\%t(H[1,,,j])+G[1,,,j]\%*\%K[1,,,j]\%*\%t(G[1,,,j])
					## ONE-STEP-AHEAD PREDICTION ERROR
					# v[t,,i,j]<-y[t,]-W[1,,,j]\%*\%ap[t,,i,j]-as.matrix(B[1,,,j])\%*\%x[t,]
					# F[t,,,i,j]<-W[1,,,j]\%*\%Pp[t,,,i,j]\%*\%t(W[1,,,j])+R[1,,,j]
					## UPDATE EQUATIONS
					# au[t,,i,j]<-ap[t,,i,j]+Pp[t,,,i,j]\%*\%t(W[1,,,j])\%*\%solve(F[t,,,i,j])\%*\%v[t,,i,j]
					# Pu[t,,,i,j]<-Pp[t,,,i,j]+Pp[t,,,i,j]\%*\%t(W[1,,,j])\%*\%solve(F[t,,,i,j])\%*\%W[1,,,j]\%*\%Pp[t,,,i,j]
					## END OF KALMAN FILTER STEP
					## BEGIN HAMILTON FILTER STEP 
					## CONDITIONAL TRANSITION PROBABILITY FOR S(t-1)=i to S(t)=j
					# tpr[t,j,i]<-p[j,i]*pr[t,i]
					## CONDITIONAL DENSITY 
				# jd[t,j,i]<-tpr[t,j,i]*((2*pi)^(-ny/2)*(det(as.matrix(F[t,,,i,j])))^(-1/2)*exp((-1/2)*t(v[t,,i,j])\%*\%solve(F[t,,,i,j])\%*\%v[t,,i,j]))
				# }
			# }
			# jm<-sum(jd[t,,])
			# L<-L+log(jm)				# LOG LIKELIHOOD FUNCTION
			# utpr[t,,]<-jd[t,,]/jm
			## END OF THE HAMILTON FILTER STEP
			## BEGIN COLLAPSE STEP to reduce the number of posteriors
			# for (j in 1:nm)
			# {	pr[(t+1),j]<-sum(utpr[t,j,])
				# as<-matrix(c(0),ne,1)
				# for (i in 1:nm)
				# {	as<-as+au[t,,i,j]*(utpr[t,j,i]/pr[(t+1),j])
				# }
				# auc[(t+1),,j]<-as
				# Ps<-matrix(c(0),ne,ne)
				# for (i in 1:nm)
				# {	Ps<-Ps+(Pu[t,,,i,j]+(auc[(t+1),,j]-au[t,,i,j])\%*\%(t(auc[(t+1),,j]-au[t,,i,j])))*(utpr[t,j,i]/pr[(t+1),j])
				# }
				# Puc[(t+1),,,j]<-Ps
			# }
		# }
		# if(any(y[t,]==-999))
		# {	for (j in 1:nm)
			# {	for (i in 1:nm)
				# {	# KALMAN FILTER STEP
					## PREDICTION EQUATIONS
					# ap[t,,i,j]<-c[1,,1,j]+H[1,,,j]\%*\%(auc[t,,i])
					# Pp[t,,,i,j]<-H[1,,,j]\%*\%Puc[t,,,i]\%*\%t(H[1,,,j])+G[1,,,j]\%*\%K[1,,,j]\%*\%t(G[1,,,j])
					## ONE-STEP-AHEAD PREDICTION ERROR CAN NOT BE CALCULATED
					## BUT THE ASSOCIATED COVARIANCE MATRIX CAN
					# F[t,,,i,j]<-W[1,,,j]\%*\%Pp[t,,,i,j]\%*\%t(W[1,,,j])+R[1,,,j]
					## THUS, THE UPDATED STATE IS IDENTICAL TO THE PREDICTED STATE
					## AND THE ASSOCIATED COVARIANCE MATRIX IS CALCULATED AS USUAL
					# au[t,,i,j]<-ap[t,,i,j]
					# Pu[t,,,i,j]<-Pp[t,,,i,j]+Pp[t,,,i,j]\%*\%t(W[1,,,j])\%*\%solve(F[t,,,i,j])\%*\%W[1,,,j]\%*\%Pp[t,,,i,j]
					## END OF KALMAN FILTER STEP
					## BEGIN HAMILTON FILTER STEP 
					## CONDITIONAL TRANSITION PROBABILITY FOR S(t-1)=i to S(t)=j
					# tpr[t,j,i]<-p[j,i]*pr[t,i]
					## CONDITIONAL DENSITY CAN NOT BE CALCULATED 
				# }
			# }
			## CURRENT OCCASION DOES NOT CONTRIBUTE TO THE LIKELIHOOD FUNCTION
			# utpr[t,,]<-tpr[t,,]
			## END OF THE HAMILTON FILTER STEP
			## BEGIN COLLAPSE STEP to reduce the number of posteriors
			# for (j in 1:nm)
			# {	pr[(t+1),j]<-sum(utpr[t,j,])
				# as<-matrix(c(0),ne,1)
				# for (i in 1:nm)
				# {	as<-as+au[t,,i,j]*(utpr[t,j,i]/pr[(t+1),j])
				# }
				# auc[(t+1),,j]<-as
				# Ps<-matrix(c(0),ne,ne)
				# for (i in 1:nm)
				# {	Ps<-Ps+(Pu[t,,,i,j]+(auc[(t+1),,j]-au[t,,i,j])\%*\%(t(auc[(t+1),,j]-au[t,,i,j])))*(utpr[t,j,i]/pr[(t+1),j])
				# }
				# Puc[(t+1),,,j]<-Ps
			# }
		# }
	# }	# END OF TIME LOOP
Lp<-if(is.nan(L)) 1e9 else min(-L,1e9)
attr(Lp,"jd") <- jd
attr(Lp,"all") <- list(jd=jd, Puc=Puc, Pu=Pu, Pp=Pp, v=v, F=F, pr=pr, tpr=tpr, utpr=utpr, au=au, ap=ap, auc=auc, md=md) 
Lp
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
