mskfModel.mskfSkeleton <- function(object, ..., start=NA, lower=NA, upper=NA)
{
    #
    obj = object
    inits = list(...);
    nmconst = names(obj$const);
    nmpattr = names(obj$pattern);
    if (!all(chk <- names(inits) %in% c(nmconst, nmpattr)))
        warning("arguments ", names(inits)[!chk], " unused");


    # "attach" elements of x to the current environment ## <-- obsolete
    for (nme in names(obj)) assign(nme, obj[[nme]]);
  missingStart = missing(start)
  missingUpper = missing(upper)
  missingLower = missing(lower)
	with(obj, {
    # default values
    defc = list(
        maW = diag(max(ny,ne))[1:ny, 1:ne],
        maB = diag(max(ny,nx))[1:ny, 1:nx],
        maR = diag(ny),
        mac = rep(0, ne),
        maH = diag(0.7, ne),
        maG = diag(ne),
        maK = diag(ne),
        map = diag(1, nm)
    )
    maxn = max(unlist(inits[nmpattr]),0) # max of labels used in inits
    npar = maxn;
    defp = list(
        paW = {tmp = matrix(npar+1,ny,ne); tmp[upper.tri(tmp)] = 0; npar<-max(tmp,npar); tmp},
        paB = {tmp = matrix(npar+1,ny,nx); tmp[upper.tri(tmp)] = 0; npar<-max(tmp,npar); tmp},
        paR = {tmp = diag(npar+1,ny); npar <- max(tmp,npar); tmp;},
        pac = {tmp = matrix(npar+1:ne,ne, 1); npar<-max(tmp,npar); tmp},
        paH = {tmp = diag(npar+1:ne,ne); npar<-max(tmp,npar); tmp},
        paG = {tmp = diag(npar+1:ne,ne); npar<-max(tmp,npar); tmp},
        paK = {tmp = diag(npar+1,ne); npar<-max(tmp,npar); tmp},
        pap = {tmp = matrix(npar+1,nm,nm); npar<-max(tmp,npar); tmp}
     )
     npar = npar - maxn + length(unique(unlist(inits[nmpattr])));


    # vectors for starting values, lower and upperbounds
    theta <- lobo <- upbo <- c()

    # set values of skeleton matrices
    for(nme in nmconst){
        obj$const[[nme]][] = defc[[nme]];
        if(nme %in% names(inits))
            obj$const[[nme]][] = inits[[nme]];
    }
    for(nme in nmpattr){
        len = length(obj$pattern[[nme]][])
        if(len & length(obj$pattern[[nme]][]) %% length(defp[[nme]]) !=0)
            stop("this shouldn't happen, you discovered a bug (", nme,")")
        else
            obj$pattern[[nme]][] = defp[[nme]];
        if(nme %in% names(inits))
            if(length(obj$pattern[[nme]][]) %% length(inits[[nme]]) !=0)
                stop("length of specification of ", nme, " is not a multiple of this model matrix")
            else
                obj$pattern[[nme]][] = inits[[nme]]
    }

    # create named starting parameters vector theta
    theta.nms = na.omit(unique(unlist(obj$pattern)))
    theta.nms = theta.nms[as.character(theta.nms)!="0"]
    npar = length(theta.nms)
    theta = rep(if (missingStart) 1/npar else start, len = npar)
    names(theta) = theta.nms
    # check provided starting values and set corresponding elements of theta to values provided
    if (!missingStart && is.null(names(start))) {
    	warning("start should be a named vector to initialize the estimated parameters");
    	names(start) = 1:length(start); # default names
    }
    if(!all(chk <- names(start) %in% theta.nms))
        warning("start value(s) named ", names(start)[!chk], " not used")
    theta[names(start)[chk]] = start[chk]

    # construct lower and upper bounds similarly (treat variance parameters special)
    lobo = rep(if (missingLower) -10 else lower, len = npar)
    names(lobo) = theta.nms
    diagR = rep(diag(ny)>0,nm);
    diagK = rep(diag(ne)>0,nm);
    varpars = na.omit(unique(c(obj$pattern$paR[diagR], obj$pattern$paK[diagK])))
    varpars = varpars[as.character(varpars)!="0"]
    lobo[as.character(varpars)] = sqrt(.Machine$double.eps) # small number away from zero
    if(!all(chk <- names(lower) %in% theta.nms))
        warning("lower bound(s) named ", names(lower)[!chk], " not used")
    lobo[names(lower)[chk]] = lower[chk]

    upbo = rep(if (missingUpper) 10 else upper, len = npar)
    names(upbo) = theta.nms
    if(!all(chk <- names(upper) %in% theta.nms))
        warning("upper bound(s) named ", names(upper)[!chk], " not used")
    upbo[names(upper)[chk]] = upper[chk]


    # initial state and its covariance matrix
    a0 = rep(0, ne)           # LET OP: DIT MOET JE EIGENLIJK VOOR ALLE REGIMES EEN KEER OPGEVEN
    ##a0 = array(a0, c(ne,  1, nm))   # doen dan dus maar? lijkt een probleem te veroorzaken
    P0 = diag(100, ne)        # IDEM
    ##P0 = array(P0, c(ne, ne, nm))   # ook maar doen dan?? lijkt een probleem te veroorzaken

    obj = c(obj, list(theta = theta, lobo = lobo, upbo = upbo, a0 = a0, P0 = P0));
    class(obj) = 'mskfModel'
    obj
	})
}
