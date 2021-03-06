mskf <- function(model,  method="L-BFGS-B", hessian=TRUE, control=list(trace=3,maxit=500),
    debug=FALSE)
{
		const = model$const
		pattern = model$pattern

  #  for(nme in names(model))   assign(nme, model[nme][[1]])
  #  for(nme in names(const))   assign(nme, const[nme][[1]])
  #  for(nme in names(pattern)) assign(nme, pattern[nme][[1]])
	with(c(model, const, pattern), {
    maxnpar <- length(model$theta)
    maxv <- max(ny,ne,nx,nm)
    MA <- array(0,c(8,maxv,maxv,nm))
    PA <- array(0,c(8,maxv,maxv,nm))
    if (ipat[1]>0) {MA[1,1:ny,1:ne,1:nm]<-maW}
    if (ipat[1]==0) {MA[1,1:ny,1:ne,1:nm]<-array(0,c(ny,ne,nm))}
    if (ipat[1]<2) {PA[1,1:ny,1:ne,1:nm]<-array(0,c(ny,ne,nm))}
    if (ipat[1]==2) {PA[1,1:ny,1:ny,1:nm]<-paW}
    if (ipat[2]>0) {MA[2,1:ny,1:nx,1:nm]<-maB}
    if (ipat[2]==0) {MA[2,1:ny,1,1:nm]<-array(0,c(ny,1,nm))}
    if (ipat[2]<2) {PA[2,1:ny,1,1:nm]<-array(0,c(ny,1,nm))}
    if (ipat[2]==2) {PA[2,1:ny,1:nx,1:nm]<-paB}
    if (ipat[3]>0) {MA[3,1:ny,1:ny,1:nm]<-maR}
    if (ipat[3]==0) {MA[3,1:ny,1:ny,1:nm]<-array(0,c(ny,ny,nm))}
    if (ipat[3]<2) {PA[3,1:ny,1:ny,1:nm]<-array(0,c(ny,ny,nm))}
    if (ipat[3]==2) {PA[3,1:ny,1:ny,1:nm]<-paR}
    if (ipat[4]>0) {MA[4,1:ne,1,1:nm]<-mac}
    if (ipat[4]==0) {MA[4,1:ne,1,1:nm]<-array(0,c(ny,1,nm))}
    if (ipat[4]<2) {PA[4,1:ne,1,1:nm]<-array(0,c(ny,1,nm))}
    if (ipat[4]==2) {PA[4,1:ne,1,1:nm]<-pac}
    if (ipat[5]>0) {MA[5,1:ne,1:ne,1:nm]<-maH}
    if (ipat[5]==0) {MA[5,1:ne,1:ne,1:nm]<-array(0,c(ne,ne,nm))}
    if (ipat[5]<2) {PA[5,1:ne,1:ne,1:nm]<-array(0,c(ne,ne,nm))}
    if (ipat[5]==2) {PA[5,1:ne,1:ne,1:nm]<-paH}
    if (ipat[6]>0) {MA[6,1:ne,1:ne,1:nm]<-maG}
    if (ipat[6]==0) {MA[6,1:ne,1:ne,1:nm]<-array(0,c(ne,ne,nm))}
    if (ipat[6]<2) {PA[6,1:ne,1:ne,1:nm]<-array(0,c(ne,ne,nm))}
    if (ipat[6]==2) {PA[6,1:ne,1:ne,1:nm]<-paG}
    if (ipat[7]>0) {MA[7,1:ne,1:ne,1:nm]<-maK}
    if (ipat[7]==0) {MA[7,1:ne,1:ne,1:nm]<-array(0,c(ne,ne,nm))}
    if (ipat[7]<2) {PA[7,1:ne,1:ne,1:nm]<-array(0,c(ne,ne,nm))}
    if (ipat[7]==2) {PA[7,1:ne,1:ne,1:nm]<-paK}
    # In addition the elements of the transition probabilities matrix can
    # be fixed or we can indicate that they should be estimated freely.
    if (ipat[8]>0) {MA[8,1:nm,1:nm,1]<-map}
    if (ipat[8]==2) {PA[8,1:nm,1:nm,1]<-pap}

    mdim<-matrix(,8,2)
    mdim[1,1]<-ny
    mdim[1,2]<-ne
    mdim[2,1]<-ny
    if (nx==0) {
        nx <-1
    	x  <-matrix(0,nt,1)
    }
    mdim[2,2] <-nx
    mdim[3,1] <-ny
    mdim[3,2] <-ny
    mdim[4,1] <-ne
    mdim[4,2] <-1
    mdim[5,1] <-ne
    mdim[5,2] <-ne
    mdim[6,1] <-ne
    mdim[6,2] <-ne
    mdim[7,1] <-ne
    mdim[7,2] <-ne
    mdim[8,1] <-nm
    mdim[8,2] <-nm

    # Before starting to run the filter, the free parameters have to be
    # renumbered so their numbers actually run from 1 to the total number of
    # free parameters. This only needs to be done once!
    iw <-matrix(0,maxnpar,1)
    npar <-0
    for (m in 1:7)
    {	for (s in 1:nm)
    	{	for (r in 1:mdim[m,1])
    		{	for (c in 1:mdim[m,2])
    			{	pr <-PA[m,r,c,s]
    				if (pr>0) {
    				    if (! any(iw[,1]==pr)) {
    				        npar <-npar+1
    						iw[npar,1] <-pr
    				    }
    				}
    			}
    		}
    	}
    }

    sx <-sort(iw[1:npar,])
    for (m in 1:7)
    {	for (s in 1:nm)
    	{	for (r in 1:mdim[m,1])
    		{	for (c in 1:mdim[m,2])
    			{	pr <-PA[m,r,c,s]
    				if (pr>0)
    				{	for (k in 1:npar)
    					{	if (pr==sx[k]) {PA[m,r,c,s] <-k}
    					}
    				}
    			}
    		}
    	}
    }

    opt<-optim(par=theta, fn=filter.C, method=method,npar=npar,
        ny=ny, ne=ne, nx=nx, nt=nt, nm=nm, mdim=mdim,
        PA=PA, MA=MA, y=y, x=x, a0=a0, P0=P0,
        hessian=hessian, control=control, lower=lobo, upper=upbo, debug=debug)
    opt
	})
}

