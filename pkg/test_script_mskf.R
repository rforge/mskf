testfunction <- function(){
simres = replicate(100, {

set.seed(123)
# parameters
phi01 = 0       # constant in regime 1
phi02 = 3       # constant in regime 2
phi11 = 0.5     # AR(1) parameter in regime 1
phi12 = 0.3     # AR(1) parameter in regime 2
rsd1  = 1       # residual standard deviation in regime 1
rsd2  = 1       # residual standard deviation in regime 2

# number of time point and burn in
nt = 1350
burn.in = 100

# Markov transition matrix
p = matrix(c(.7,.3,.2,.8), 2, 2)
lp = eigen(p)$vectors[,1]
lp = lp/sum(lp) # equilibrium regime distribution

true.par = c(phi01, phi02, phi11, phi12, rsd1, rsd2, log(p[1,1]/p[2,1]),
						 log(p[2,2]/p[1,2]))

# Create a discrete Markov series
r = sample(1:2, 1, p = lp);
for (t in 2:(nt+burn.in))
	r[t] = sample(1:2, 1, prob = p[,r[t-1]]);
y = rnorm(1)
for (t in 2:(nt+burn.in))
	y[t]  = if (r[t]==1)  phi01 + phi11 * y[t-1] + rnorm(1,0,rsd1) else
		phi02 + phi12 * y[t-1] + rnorm(1,0,rsd2)
y <- y[-(1:burn.in)]


# Specify the model
skel = mskfModelSkeleton(y, ne=1, nm=2)
start = c(`1` = .05, `2` = 3, `3` = .51, `4` = .22, `5` = 1.5, `6` = 2,
					`7` = 1, `8` = .5)
mod = mskfModel(skel, pap = diag(7:8), paK = 5:6, pac = 1:2, paH = 3:4, start = start,
								upper = 50, lower = -50)
fitres <- mskf(mod,control=list(trace=3,maxit=2000),debug=0)
fitres
})
simres
}
