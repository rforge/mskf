sim.msar <- function(
# Simulate data for a Markov switching AR model

nt = 1000,					        # length of time series
p = matrix(c(.7,.3,.6,.4),2),		# matrix with transition probs (rows should sum to 1)
# Choose paramters for the two regimes
H = array(c(.3, .5), c(1, 1, 2)),	# AR parameters for each regime in each row
m = array(0, c(dim(H)[1],1,dim(H)[3])),			# constants in each regime
G = array(diag(dim(H)[1]), dim(H)),		# innovation covariance matrix for each regime
burn.in = 100 + nt / 10)
{
	nt = nt + burn.in;
	c = m;

    # long run regime probabilities
    a = eigen(p)$ve[,1];
    a = a/sum(a);

    # Sample regimes using transition probs p
    x <- rep(NA,nt)
    x[1] <- sample(1:length(a), 1, prob = a)
    for (t in 2:nt)
    	x[t] <- sample(1:length(a), 1, prob = p[,x[t-1]])

    # Sample MS AR

    resid = array(rnorm(nt * prod(dim(H)[-1])), dim = c(nt, dim(H)[-1]))
    for (i in 1:dim(resid)[3]) {
    	resid[,,i] = resid[,,i] %*% chol(G[,, 1 + (i - 1) %% dim(G)[3]]); # multivariate normal, recycle G
    }
    y = matrix(NA, nt, dim(H)[1]);
    y[1, ] = resid[1, , x[1]];
    for (i in 2:nt) {
    	y[i, ] = c[,, x[i]] + H[,, x[i]] %*% y[i - 1, ] + resid[i,, x[i]];
    }
    list(regime = x[-(1:burn.in)], y = as.ts(y[-(1:burn.in),,drop=FALSE]))
}

