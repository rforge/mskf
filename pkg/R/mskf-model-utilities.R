

.PATTERNFUNCTION = "function(x, i, value) { if (missing(i)) { x$pattern$paK[] = if (is.array(value)) value else aperm(array(value, dim=dim(x$pattern$paK)), c(2,1,3));} else { x$pattern$paK[,,i] = if (is.array(value)) value else aperm(array(value, dim=dim(x$pattern$paK)[-3]), c(2,1))}; x }"
.MATRIXFUNCTION = "function(x, i, value) { if (missing(i)) { x$const$maK[] = if (is.array(value)) value else aperm(array(value, dim=dim(x$const$maK)), c(2,1,3));} else { x$const$maK[,,i] = if (is.array(value)) value else aperm(array(value, dim=dim(x$const$maK)[-3]), c(2,1))}; x }"
for (a in c("W", "B", "R", "c", "H", "G", "K", "p")) {
	assign(paste("pa", a,"<-", sep=""), eval(parse(text = gsub("K", a, .PATTERNFUNCTION))));
	assign(paste("pa", a,"", sep=""), eval(parse(text = gsub("K", a, "function(x) x$pattern$paK"))));
	assign(paste("ma", a,"<-", sep=""), eval(parse(text = gsub("K", a, .MATRIXFUNCTION))));
	assign(paste("ma", a,"", sep=""), eval(parse(text = gsub("K", a, "function(x) x$const$maK"))));
}
"pap<-" = function(x, value) { x$pattern$pap[] = if (is.array(value)) value else aperm(array(value, dim=dim(x$pattern$pap)), c(2,1)); x }
"map<-"  = function(x, value) { x$const$map[] = if (is.array(value)) value else aperm(array(value, dim=dim(x$const$map)), c(2,1)); x }
startValues = function(x, i, ...) UseMethod("startValues")
startValues.mskfModel = function(x, i, ...) {
	if (missing(i)) {
		x$theta;
	}
	else if (i %in% names(x$pattern)) {
		arr = x$pattern[[i]];
		nms = as.character(unique(arr[arr != 0]));
		x$theta[nms];
		for (i in nms) arr[arr == i] = x$theta[i];
		arr;
	}
	else x$theta[i];
}
"startValues<-" <- function (x, i, value) UseMethod("startValues<-")
"startValues<-.mskfModel" = function(x, i, value) {
	if (missing(i)) {
		if (!is.null(names(value))) {
			x$theta[names(value)] = value;
		}
		else {
			x$theta[] = value;
		}
	}
	else if (length(i)==1 & i %in% names(x$pattern)) {
		# sapply(names(mdl$pattern), function(nm) {x = get(nm)(mdl); unique(x[x!=0])})
		arr = x$pattern[[i]];
		if (is.null(names(value))) {
			if (!is.array(value) || !all.equal(dim(value), dim(arr))) {
				nms = as.character(unique(arr[arr!=0]))
				stop('value [', value, '] should be a vector of the form\n\tc(', paste('`',nms,'` = ...', sep="", collapse=", "),')');
			}
			else {
				value = c(value);
				names(value) = as.character(c(arr));
				value = value[arr != 0];
			}
		}
		nms = as.character(unique(arr[arr != 0]));
		nms = intersect(nms, names(value));
		if (length(nms) > 0) {
			x$theta[nms] = value[nms];
		}
	}
	else x$theta[i] = value;
	x;
}
lower.mskfModel = function(x, i) {
	if (missing(i)) {
		x$lobo;
	}
	else if (i %in% names(x$pattern)) {
		arr = x$pattern[[i]];
		nms = as.character(unique(arr[arr != 0]));
		x$lobo[nms];
		for (i in nms) arr[arr == i] = x$lobo[i];
		arr;
	}
	else x$lobo[i];
}
"lower" <- function (x, i) UseMethod("lower")
"lower<-" <- function (x, i, value) UseMethod("lower<-")
"lower<-.mskfModel" = function(x, i, value) {
	if (missing(i)) {
		if (!is.null(names(value))) {
			x$lobo[names(value)] = value;
		}
		else {
			x$lobo[] = value;
		}
	}
	else if (length(i)==1 & i %in% names(x$pattern)) {
		# sapply(names(mdl$pattern), function(nm) {x = get(nm)(mdl); unique(x[x!=0])})
		arr = x$pattern[[i]];
		if (is.null(names(value))) {
			if (!is.array(value) || !all.equal(dim(value), dim(arr))) {
				nms = as.character(unique(arr[arr!=0]))
				stop('value [', value, '] should be a vector of the form\n\tc(', paste('`',nms,'` = ...', sep="", collapse=", "),')');
			}
			else {
				value = c(value);
				names(value) = as.character(c(arr));
				value = value[arr != 0];
			}
		}
		nms = as.character(unique(arr[arr != 0]));
		nms = intersect(nms, names(value));
		if (length(nms) > 0) {
			x$lobo[nms] = value[nms];
		}
	}
	else x$lobo[i] = value;
	x;
}
upper.mskfModel = function(x, i) {
	if (missing(i)) {
		x$upbo;
	}
	else if (i %in% names(x$pattern)) {
		arr = x$pattern[[i]];
		nms = as.character(unique(arr[arr != 0]));
		x$upbo[nms];
		for (i in nms) arr[arr == i] = x$upbo[i];
		arr;
	}
	else x$upbo[i];
}
"upper" <- function (x, i) UseMethod("upper")
"upper<-" <- function (x, i, value) UseMethod("upper<-")
"upper<-.mskfModel" = function(x, i, value) {
	if (missing(i)) {
		if (!is.null(names(value))) {
			x$upbo[names(value)] = value;
		}
		else {
			x$upbo[] = value;
		}
	}
	else if (length(i)==1 & i %in% names(x$pattern)) {
		# sapply(names(mdl$pattern), function(nm) {x = get(nm)(mdl); unique(x[x!=0])})
		arr = x$pattern[[i]];
		if (is.null(names(value))) {
			if (!is.array(value) || !all.equal(dim(value), dim(arr))) {
				nms = as.character(unique(arr[arr!=0]))
				stop('value [', value, '] should be a vector of the form\n\tc(', paste('`',nms,'` = ...', sep="", collapse=", "),')');
			}
			else {
				value = c(value);
				names(value) = as.character(c(arr));
				value = value[arr != 0];
			}
		}
		nms = as.character(unique(arr[arr != 0]));
		nms = intersect(nms, names(value));
		if (length(nms) > 0) {
			x$upbo[nms] = value[nms];
		}
	}
	else x$upbo[i] = value;
	x;
}
coef.mskfModel = function(object, ...) object$theta
# "$<-.mskf.model" = function(x, i, value) {stop("Not possible to assign values with $ operator. See help(pac) for assigning values to model"); x}
# "$.mskf.model" = function(x, i, value) {x[[i]]}


addmat <- function(A, B){
	if(!is.matrix(A) || !is.matrix(B)) stop('only matrices allowed')
	if (!all(dim(A) == dim(B))) stop("non-conformable matrices")
	nr = nrow(A)
	nc = ncol(A)
	C = A
	res <- .C("addmat", as.integer(nr), as.integer(nc), as.double(A), as.double(B), as.double(C),
		 NAOK = FALSE, PACKAGE = "mskf")
	matrix(res[[length(res)]], nr, nc)
}


mmatmulABAt  <- function(A, B){
	if(!is.matrix(A) || !is.matrix(B)) stop('only matrices allowed')
	if (dim(A)[2] != dim(B)[1] || dim(B)[1] != dim(B)[2]) stop("non-conformable matrices")
	nr = nrow(A)
	nc = ncol(A)
	C = matrix(0, nr, nr)
	res <- .C("mmatmulABAt", as.integer(nr), as.integer(nc), as.double(A), as.double(B), as.double(C),
						NAOK = FALSE, PACKAGE = "mskf")
	matrix(res[[length(res)]], nr, nr)
}


mmatmulASAt <- function(A, S){
	if(!is.matrix(A) || !is.matrix(S)) stop('only matrices allowed')
	if(ncol(S) != nrow(S) || !all(S == t(S))) stop('`S` should be square symmetric.')
	if(ncol(A) != ncol(S)) stop('non-conformable matrices')
	nr = nrow(A)
	nc = ncol(A)
	C = matrix(0, nr, nr)
	res <- .C("mmatmulASAt", as.integer(nr), as.integer(nc), as.double(A), as.double(S), as.double(C),
						NAOK = FALSE, PACKAGE = "mskf")
	matrix(res[[length(res)]], nr, nr)
}

mmatmulAB <- function(A, B){
	if(!is.matrix(A) || !is.matrix(B)) stop('only matrices allowed')
	if (dim(A)[2] != dim(B)[1]) stop("non-conformable matrices")
	nr = nrow(A)
	nc = ncol(A)
	ncb= ncol(B)
	C = matrix(0, nr, ncb)
	res <- .C("mmatmulAB", as.integer(nr), as.integer(nc), as.double(A), as.double(B),
						as.integer(ncb), as.double(C), NAOK = FALSE, PACKAGE = "mskf")
	matrix(res[[length(res)]], nr, ncb)
}


mmatmulABt <- function(A, B){
	if(!is.matrix(A) || !is.matrix(B)) stop('only matrices allowed')
	if (dim(A)[2] != dim(B)[2]) stop("non-conformable matrices")
	nr = nrow(A)
	nc = ncol(A)
	nrb= nrow(B)
	C = matrix(0, nr, nrb)
	res <- .C("mmatmulABt", as.integer(nr), as.integer(nc), as.double(A), as.double(B),
						as.integer(nrb), as.double(C), NAOK = FALSE, PACKAGE = "mskf")
	matrix(res[[length(res)]], nr, nrb)
}

mmatmulAtB <- function(A, B){
	if(!is.matrix(A) || !is.matrix(B)) stop('only matrices allowed')
	if (dim(A)[1] != dim(B)[1]) stop("non-conformable matrices")
	nr = nrow(A)
	nc = ncol(A)
	ncb= ncol(B)
	C = matrix(0, nc, ncb)
	res <- .C("mmatmulAtB", as.integer(nr), as.integer(nc), as.double(A), as.double(B),
						as.integer(ncb), as.double(C), NAOK = FALSE, PACKAGE = "mskf")
	matrix(res[[length(res)]], nc, ncb)
}

mmataddAB <- function(A, B){
	if(!is.matrix(A) || !is.matrix(B)) stop('only matrices allowed')
	if (!all(dim(A) == dim(B))) stop("non-conformable matrices")
	nr = nrow(A)
	nc = ncol(A)
	C = A
	res <- .C("mmataddAB", as.integer(nr), as.integer(nc), as.double(A), as.double(B), as.double(C),
						NAOK = FALSE, PACKAGE = "mskf")
	matrix(res[[length(res)]], nr, nc)
}

mmataddABt <- function(A, B){
	if(!is.matrix(A) || !is.matrix(B)) stop('only matrices allowed')
	if (!all(dim(A) == dim(t(B)))) stop("non-conformable matrices")
	nr = nrow(A)
	nc = ncol(A)
	C = A
	res <- .C("mmataddABt", as.integer(nr), as.integer(nc), as.double(A), as.double(B), as.double(C),
						NAOK = FALSE, PACKAGE = "mskf")
	matrix(res[[length(res)]], nr, nc)
}


mx_plus_alpha_multAb <- function(x, alpha, A, b){
	if(!is.matrix(A)) stop('only matrices allowed')
	if (dim(A)[2] != length(b) || dim(A)[1] != length(x)) stop("non-conformable arrays")
	nr = nrow(A)
	nc = ncol(A)
	C = x
	res <- .C("mx_plus_alpha_multAb", as.double(x), as.integer(nr), as.double(alpha[1]),
						as.double(A), as.integer(nc), as.double(b), as.double(C),
						NAOK = FALSE, PACKAGE = "mskf")
	matrix(res[[length(res)]], nr)
}

mX_plus_alpha_multAB <- function(X, alpha, A, B){
	if(!is.matrix(A) || !is.matrix(B) || !is.matrix(X))
		stop('only matrices allowed')
	if (dim(X)[1] == dim(A)[1] || dim(A)[2] != dim(B)[1] || dim(X)[2] != dim(B)[2])
		stop("non-conformable arrays")
	nr = nrow(X)
	nc = ncol(X)
	nca = ncol(A)
	C = X
	res <- .C("mX_plus_alpha_multAB", as.double(X), as.integer(nr), as.integer(nc),
						as.double(alpha[1]), as.double(A), as.integer(nca), as.double(B), as.double(C),
						NAOK = FALSE, PACKAGE = "mskf")
	matrix(res[[length(res)]], nr, nc)
}

mX_plus_alpha_multABt <- function(X, alpha, A, B){
	if(!is.matrix(A) || !is.matrix(B) || !is.matrix(X))
		stop('only matrices allowed')
	if (dim(X)[1] != dim(A)[1] || dim(A)[2] != dim(B)[2] || dim(X)[2] != dim(B)[1])
		stop("non-conformable arrays")
	nr = nrow(X)
	nc = ncol(X)
	nca = ncol(A)
	C = X
	res <- .C("mX_plus_alpha_multABt", as.double(X), as.integer(nr), as.integer(nc),
						as.double(alpha[1]), as.double(A), as.integer(nca), as.double(B), as.double(C),
						NAOK = FALSE, PACKAGE = "mskf")
	matrix(res[[length(res)]], nr, nc)
}

