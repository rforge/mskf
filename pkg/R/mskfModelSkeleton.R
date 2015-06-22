mskfModelSkeleton <-
function(y, ne=1, nm=2, x = NA, ipat = c(W=1,B=nx,R=0,c=2,H=2,G=1,K=2,p=2))
{
    y = as.matrix(y)
    ny = ncol(y)
    nt = nrow(y)
    nx = if(missing(x)) 0 else NCOL(x)
    x  = as.matrix(x)

    val.ipat.nms = c("W", "B", "R", "c", "H", "G", "K", "p")
    if(!all(names(ipat) %in% val.ipat.nms)||!all(val.ipat.nms %in% names(ipat)))
        stop('ipat should have elements named "W", "B", "R", "c", "H", "G", "K", "p"')

	dims = list(
		W = c(ny, ne, nm),
		B = c(ny, nx, nm),
		R = c(ny, ny, nm),
		c = c(ne,  1, nm),
		H = c(ne, ne, nm),
		G = c(ne, ne, nm),
		K = c(ne, ne, nm),
		p = c(nm, nm)
	)
	const = list();
	pattr = list();
	for(nam in names(ipat)){
		pa = paste("pa",nam,sep="")
		ma = paste("ma",nam,sep="")
		d = dims[[nam]];
		dnam = if(length(d)>2) list( seq_len(d[1]), seq_len(d[2]), paste('regime',1:d[3])) else NULL;
		if(ipat[nam]>=1)
			const[[ma]] = array(NA, d, dnam);
		if(ipat[nam]==2)
			pattr[[pa]] = array(NA, d, dnam);
	}
	mdl <- list(const = const, pattern = pattr, y=y, x=x, ny=ny, nx=nx, nt=nt, ne=ne, nm=nm, ipat=ipat)
	class(mdl) <- "mskfSkeleton"
	mdl
}
