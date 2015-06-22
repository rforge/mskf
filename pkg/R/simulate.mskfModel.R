simulate.mskfModel <- function(object, nsim = 1, seed = NA, nt = x$nt, use.start = FALSE, ...) {
	getStart = function(x, pa) na.omit(start(x)[as.character(pa)]);
	x = object;

	c = mac(x);
	H = maH(x);
	p = map(x);
	G = maG(x);
	W = maW(x);
	B = maB(x);
	R = maR(x);

	if (use.start) {
		c[pac(x) != 0] = getStart(x, pac(x));
		H[paH(x) != 0] = getStart(x, paH(x));
		p[pap(x) != 0] = getStart(x, pap(x));
		G[paG(x) != 0] = getStart(x, paG(x));
		W[paW(x) != 0] = getStart(x, paW(x));
		B[paB(x) != 0] = getStart(x, paB(x));
		R[paR(x) != 0] = getStart(x, paR(x));
	}

	X = x$x
	sqrtR = R = if (is.null(R)) array(diag(x$ny), c(x$ny, x$ny, x$nm)) else R;
	for (i in 1:x$ne) {
		sqrtR[,,i] = chol(R[,,i]);
	}

	if (!is.na(seed)) {
		set.seed(seed);
	}
	r = replicate(nsim, {
		a = sim.msar(nt, p = p, H = H, m = c, G = G);
		y = array(NA, dim = c(nt, x$ny));
		e = matrix(rnorm(nt * x$ny), nt);
		for (i in 1:x$nm) {
			y[a$regime == i, ] = a$y[a$regime == i,, drop=FALSE] %*% t(array(W[,,i], dim(W)[-3])) + e[a$regime ==i,, drop = FALSE] %*% sqrtR[,, i];
			if (!is.null(X) && !is.null(B)) {
				ind = a$regime == i;
				y[ind, ] = y[ind, ] + X[ind,] %*% t(array(B[,,i], dim(B)[-3]));
			}
		}
		list(y = as.ts(y), a = a$y, regime = a$regime);
	});
	attr(r, "seed") = seed;
	drop(r);
}


