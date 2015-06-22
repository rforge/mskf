print.mskfModel <-
function (x, ...)
{
		obj = x
    print.mskf.skeleton(obj, ...)
    cat("\n\nMarkov switching state space model specification\n\n")
    cat("Fixed values\n")
    for (i in names(obj$const)) {
    	o = obj$const[[i]]
    	cat(i,": num ", paste("[", paste(dim(o), collapse=", "), "] ", sep=""), sep="")
    	cat(c(o)[1:min(7,length(o))], "\n", ...)
    }
    cat("\nFree (to be estimated) parameters\n")
    for (i in names(obj$pattern)) {
    	o = obj$pattern[[i]]
    	cat(i,": num ", paste("[", paste(dim(o), collapse=", "), "] ", sep=""), sep="")
    	cat(c(o)[1:min(7,length(o))], "\n", ...)
    }
    tmp = rbind(obj$theta, obj$lobo, obj$upbo)
    rownames(tmp) = c("theta", "lower bound", "upper bound")
    print(tmp, ...)
}
