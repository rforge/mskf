print.mskfSkeleton <- function(x, ...)
{
		obj = x
    with(obj, {
    cat("\nMarkov switching state space model skeleton\n\n", ...);
    cat("\nNumber of endogenous variables : ", ny, ...);
    cat("\nNumber of exogenous variables  : ", nx, ...);
    cat("\nNumber of latent variables     : ", ne, ...);
    cat("\nNumber of states               : ", nm, ...);
    cat("\nNumber of time points          : ", nt, ...);
    tmp = lapply(const,dim)
    tmp = sapply(tmp, paste, collapse=" x ")
    tmp = paste(c("\n","",""), names(tmp), " (", tmp, ")", sep="", collapse=", ")
    cat("\n\nMatrices used                : ", tmp, ...);
    tmp = lapply(pattern,dim)
    tmp = sapply(tmp, paste, collapse=" x ")
    tmp = paste(c("\n","",""), names(tmp), " (", tmp, ")", sep="", collapse=", ")
    cat("\n\nMatrices with free parameters used : ", tmp, "\n",...);
    })
    invisible(obj)
}
