

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
start.mskf.model = function(x, i) {
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
"start<-" <- function (x, ...) UseMethod("start<-")
"start<-.mskf.model" = function(x, i, value, ...) {
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
lower.mskf.model = function(x, i) {
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
"lower" <- function (x, ...) UseMethod("lower")
"lower<-" <- function (x, ...) UseMethod("lower<-")
"lower<-.mskf.model" = function(x, i, value) {
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
upper.mskf.model = function(x, i) {
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
"upper" <- function (x, ...) UseMethod("upper")
"upper<-" <- function (x, ...) UseMethod("upper<-")
"upper<-.mskf.model" = function(x, i, value) {
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


