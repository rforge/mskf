`Kfilter.timeloop` <-
function(   auc, au, ap,
            Puc, Pu, Pp,
            v, F, p,
            pr, tpr, utpr, jd, md,
            W, B, R, c, H, G,K,
            y, x, debug=FALSE){
#	auc[is.na(auc)] <- -9999
#	au[is.na(au)] <- -9999
#	ap[is.na(ap)] <- -9999
#	Puc[is.na(Puc)] <- -9999
#	Pu[is.na(Pu)] <- -9999
#	Pp[is.na(Pp)] <- -9999
#	v[is.na(v)] <- -9999
#	F[is.na(F)] <- -9999
#	pr[is.na(pr)] <- -9999
#	tpr[is.na(tpr)] <- -9999
#	utpr[is.na(utpr)] <- -9999
#	jd[is.na(jd)] <- -9999
#	md[is.na(md)] <- -9999
	if(any(c(is.na(W),is.na(B),is.na(R),is.na(c),is.na(H),is.na(G),is.na(K))))
		stop('Arrays W, B, R, c, H, G, K should not contain NAs')

    nt = dim(v)[1];
    ny = dim(v)[2];
    nm = dim(v)[3];
    ne = dim(au)[2];
    nx = dim(x)[2];

    auc  = aperm(auc,  c(2:3,1))
    au   = aperm(au,   c(2:4,1))
    ap   = aperm(ap,   c(2:4,1))
    Puc  = aperm(Puc,  c(2:4,1))
    Pu   = aperm(Pu,   c(2:5,1))
    Pp   = aperm(Pp,   c(2:5,1))
    v    = aperm(v,    c(2:4,1))
    F    = aperm(F,    c(2:5,1))
    pr   = aperm(pr,   c(2:2,1))
    tpr  = aperm(tpr,  c(2:3,1))
    utpr = aperm(utpr, c(2:3,1))
    jd   = aperm(jd,   c(2:3,1))
    md   = aperm(md,   c(2:2,1))
    y    = t(y)
    x    = t(x)
#nl="\n"
#cat(nt,nl,ny,nl,nm,nl,ne,nl,nx,nl,nl)
#cat(auc,nl,au,nl,ap,nl,nl)
#cat(Puc,nl,Pu,nl,Pp,nl,nl)
#cat(v,F,nl,nl);
#cat(pr,tpr,utpr,nl,nl)
#cat(jd,md)
retval = list(empty="failed")
#on.exit(if(is.nan(retval$L)) print(retval));

    retval <-
   .C("kfilter_timeloop",
      as.integer(nt),
      as.integer(nm),
      as.integer(ne),
      as.integer(ny),
      as.integer(nx),
      c=as.double(c),
      H=as.double(H),
      B=as.double(B),
      K=as.double(K),
      G=as.double(G),
      R=as.double(R),
      W=as.double(W),
      v=as.double(v),
      F=as.double(F),
      pr=as.double(pr),
      tpr=as.double(tpr),
      utpr=as.double(utpr),
      jd=as.double(jd),
      md=as.double(md),
      auc=as.double(auc),
      au=as.double(au),
      ap=as.double(ap),
      p=as.double(p),
      Puc=as.double(Puc),
      Pu=as.double(Pu),
      Pp=as.double(Pp),
      y=as.double(y),
      x=as.double(x),
	L=double(1),
      debug=as.integer(debug),
      NAOK = TRUE, PACKAGE = "mskf")

	pf <- parent.frame();
#    auc   <<- aperm(array(retval$auc,     c(ne,nm,nt+1)), c(3,1:2))
#    au    <<- aperm(array(retval$au,     c(ne,nm,nm,nt)), c(4,1:3))
#    ap    <<- aperm(array(retval$ap,     c(ne,nm,nm,nt)), c(4,1:3))
#    Puc   <<- aperm(array(retval$Puc,    c(ne,ne,nm,nt)), c(4,1:3))
#    Pu    <<- aperm(array(retval$Pu,  c(ne,ne,nm,nm,nt)), c(5,1:4))
#    Pp    <<- aperm(array(retval$Pp,  c(ne,ne,nm,nm,nt)), c(5,1:4))
#    v     <<- aperm(array(retval$v,      c(ny,nm,nm,nt)), c(4,1:3))
#    F     <<- aperm(array(retval$F,   c(ne,ne,nm,nm,nt)), c(5,1:4))
#    pr    <<- aperm(array(retval$pr,           c(nm,nt)), c(2,1:1))
#    tpr   <<- aperm(array(retval$tpr,       c(nm,nm,nt)), c(3,1:2))
#    utpr  <<- aperm(array(retval$utpr,      c(nm,nm,nt)), c(3,1:2))
#    jd    <<- aperm(array(retval$jd,        c(nm,nm,nt)), c(3,1:2))
#    md    <<- aperm(array(retval$md,            c(nt,1)), c(1,2:2))
    assign("auc",  aperm(array(retval$auc,     c(ne,nm,nt+1)), c(3,1:2)), envir=pf)
    assign("au",   aperm(array(retval$au,     c(ne,nm,nm,nt)), c(4,1:3)), envir=pf)
    assign("ap",   aperm(array(retval$ap,     c(ne,nm,nm,nt)), c(4,1:3)), envir=pf)
    assign("Puc",  aperm(array(retval$Puc,    c(ne,ne,nm,nt)), c(4,1:3)), envir=pf)
    assign("Pu",   aperm(array(retval$Pu,  c(ne,ne,nm,nm,nt)), c(5,1:4)), envir=pf)
    assign("Pp",   aperm(array(retval$Pp,  c(ne,ne,nm,nm,nt)), c(5,1:4)), envir=pf)
    assign("v",    aperm(array(retval$v,      c(ny,nm,nm,nt)), c(4,1:3)), envir=pf)
	assign("F",    aperm(array(retval$F,   c(ne,ne,nm,nm,nt)), c(5,1:4)), envir=pf)
    assign("pr",   aperm(array(retval$pr,           c(nm,nt)), c(2,1:1)), envir=pf)
    assign("tpr",  aperm(array(retval$tpr,       c(nm,nm,nt)), c(3,1:2)), envir=pf)
    assign("utpr", aperm(array(retval$utpr,      c(nm,nm,nt)), c(3,1:2)), envir=pf)
    assign("jd",   aperm(array(retval$jd,        c(nm,nm,nt)), c(3,1:2)), envir=pf)
    assign("md",   aperm(array(retval$md,            c(nt,1)), c(1,2:2)), envir=pf)

    retval$L;
}

