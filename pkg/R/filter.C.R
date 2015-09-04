`filter.C` <-
function(theta,npar,ny,ne,nx,nt,nm,mdim,MA,PA,a0,P0,y,x, debug=FALSE)
{
  # START WITH MAKING THE ARRAYS THAT ARE NEEDED
  auc  <- array(0,c((nt + 1),ne,nm))     # a-updated and collapsed to nm posteriors
  au   <- array(0,c((nt),ne,nm,nm))      # a-updated (not collapsed: S(t-1)=i, S(t)=j)
  ap   <- array(0,c((nt),ne,nm,nm))      # a-predicted (by definition not collapsed)
  Puc  <- array(0,c((nt + 1),ne,ne,nm))  # covariance matrix of auc
  Pu   <- array(0,c(nt,ne,ne,nm,nm))     # covariance matrix of au
  Pp   <- array(0,c(nt,ne,ne,nm,nm))     # covariance matrix of ap
  v    <- array(0,c(nt,ny,nm,nm))        # one-step-ahead prediction error for S(t-1)=1, S(t)=j
  F    <- array(0,c(nt,ny,ny,nm,nm))     # covariance matrix of v
  pr   <- array(0,c((nt + 1),nm))        # updated probabilities
  tpr  <- array(0,c(nt,nm,nm))           # conditional transition probabilites
  utpr <- array(0,c(nt,nm,nm))           # updated conditional transition probabilites
  jd   <- array(0,c(nt,nm,nm))           # joint densities (product of conditional density and transition probability)
  md   <- array(0,c(nt,1))               # marginal density (used in optimization)

  # MAKE UPDATED PROABILITIES (I.E., VECTOR WITH PROBABILITIES OF
  # OF BEING IN CERTAIN STATE AT OCCASION 0)
  pr[1,] <- matrix(c(1 / nm),nm,1)
  # ENTER INITIAL VALUES OF a0 AND P0 IN auc
  auc[1,,] <- a0
  Puc[1,,,] <- P0
  # PLACE THE ELEMENTS FROM x (= vector with free parameters)
  # IN THE APPROPRIATE MODEL MATRICES
  for (m in 1:7)
  {
    for (s in 1:nm)
    {
      for (r in 1:mdim[m,1])
      {
        for (c in 1:mdim[m,2])
        {
          p <- PA[m,r,c,s]
          if (p > 0) {
            MA[m,r,c,s] <- theta[p]
          }
        }
      }
    }
  }
  # TO OBTAIN THE MATRIX WITH TRANSITION PROBABILITIES,
  # WE HAVE TO TAKE INTO ACCOUNT THAT THE PORBS HAVE TO SUM TO 1 PER COLUMN
  # HENCE, PER COLUMN THERE ARE nm-1 PROBS TO BE ESTIMATED
  # RATHER THAN ESTIMATING THE PROBS WE WILL ESTIMATE A NUMBER
  # WHICH CAN RUN FROM -INFINITY TO +INFINITY
  tp <- matrix(0, nm, nm)
  for (c in 1:(nm - 1))
  {
    tp[nm,c] <- exp(0)
    for (r in 1:(nm - 1))
    {
      npar <- npar + 1
      tp[r,c] <- exp(theta[npar])
    }
  }
  tp[1,nm] <- exp(0)
  for (r in 2:nm)
  {
    npar <- npar + 1
    tp[r,nm] <- exp(theta[npar])
  }
  # THE ZEROS ARE INCLUDED FOR IDENTIFICATION
  # TO COME FROM THE ABOVE MATRIX tp TO THE MATRIX WITH ACTUAL
  # TRANSITION PROBS p, WE DO THE FOLLOWING
  p <- matrix(0,nm,nm)
  for (c in 1:nm)
  {
    stp <- sum(tp[,c])
    p[,c] <- tp[,c] / stp
  }
  # SELECT THE MODEL MATRICES AND NAME THEM
  W <- MA[1,1:ny,1:ne,,drop = FALSE]
  B <- MA[2,1:ny,1:nx,,drop = FALSE]
  R <- MA[3,1:ny,1:ny,,drop = FALSE]
  c <- MA[4,1:ne,1,,drop = FALSE]
  H <- MA[5,1:ne,1:ne,,drop = FALSE]
  G <- MA[6,1:ne,1:ne,,drop = FALSE]
  K <- MA[7,1:ne,1:ne,,drop = FALSE]
  # START THE KALMAN FILTER & REGIME SWITCHING
  L <- 0
  L <- Kfilter.timeloop(auc, au, ap,
                        Puc, Pu, Pp,
                        v, F, p,
                        pr, tpr, utpr, jd, md,
                        W, B, R, c, H, G, K, y, x, debug = debug)
  # for (t in 1:nt)
  # {  # DETERMINE WHETHER THERE ARE MISSING VALUES
    # if(all(!y[t,]==-999))
    # {  for (j in 1:nm)
      # {  for (i in 1:nm)
        # {  # KALMAN FILTER STEP
          ## PREDICTION EQUATIONS
          # ap[t,,i,j]<-c[1,,1,j]+H[1,,,j]%*%(auc[t,,i])
          # Pp[t,,,i,j]<-H[1,,,j]%*%Puc[t,,,i]%*%t(H[1,,,j])+G[1,,,j]%*%K[1,,,j]%*%t(G[1,,,j])
          ## ONE-STEP-AHEAD PREDICTION ERROR
          # v[t,,i,j]<-y[t,]-W[1,,,j]%*%ap[t,,i,j]-as.matrix(B[1,,,j])%*%x[t,]
          # F[t,,,i,j]<-W[1,,,j]%*%Pp[t,,,i,j]%*%t(W[1,,,j])+R[1,,,j]
          ## UPDATE EQUATIONS
          # au[t,,i,j]<-ap[t,,i,j]+Pp[t,,,i,j]%*%t(W[1,,,j])%*%solve(F[t,,,i,j])%*%v[t,,i,j]
          # Pu[t,,,i,j]<-Pp[t,,,i,j]+Pp[t,,,i,j]%*%t(W[1,,,j])%*%solve(F[t,,,i,j])%*%W[1,,,j]%*%Pp[t,,,i,j]
          ## END OF KALMAN FILTER STEP
          ## BEGIN HAMILTON FILTER STEP
          ## CONDITIONAL TRANSITION PROBABILITY FOR S(t-1)=i to S(t)=j
          # tpr[t,j,i]<-p[j,i]*pr[t,i]
          ## CONDITIONAL DENSITY
        # jd[t,j,i]<-tpr[t,j,i]*((2*pi)^(-ny/2)*(det(as.matrix(F[t,,,i,j])))^(-1/2)*exp((-1/2)*t(v[t,,i,j])%*%solve(F[t,,,i,j])%*%v[t,,i,j]))
        # }
      # }
      # jm<-sum(jd[t,,])
      # L<-L+log(jm)        # LOG LIKELIHOOD FUNCTION
      # utpr[t,,]<-jd[t,,]/jm
      ## END OF THE HAMILTON FILTER STEP
      ## BEGIN COLLAPSE STEP to reduce the number of posteriors
      # for (j in 1:nm)
      # {  pr[(t+1),j]<-sum(utpr[t,j,])
        # as<-matrix(c(0),ne,1)
        # for (i in 1:nm)
        # {  as<-as+au[t,,i,j]*(utpr[t,j,i]/pr[(t+1),j])
        # }
        # auc[(t+1),,j]<-as
        # Ps<-matrix(c(0),ne,ne)
        # for (i in 1:nm)
        # {  Ps<-Ps+(Pu[t,,,i,j]+(auc[(t+1),,j]-au[t,,i,j])%*%(t(auc[(t+1),,j]-au[t,,i,j])))*(utpr[t,j,i]/pr[(t+1),j])
        # }
        # Puc[(t+1),,,j]<-Ps
      # }
    # }
    # if(any(y[t,]==-999))
    # {  for (j in 1:nm)
      # {  for (i in 1:nm)
        # {  # KALMAN FILTER STEP
          ## PREDICTION EQUATIONS
          # ap[t,,i,j]<-c[1,,1,j]+H[1,,,j]%*%(auc[t,,i])
          # Pp[t,,,i,j]<-H[1,,,j]%*%Puc[t,,,i]%*%t(H[1,,,j])+G[1,,,j]%*%K[1,,,j]%*%t(G[1,,,j])
          ## ONE-STEP-AHEAD PREDICTION ERROR CAN NOT BE CALCULATED
          ## BUT THE ASSOCIATED COVARIANCE MATRIX CAN
          # F[t,,,i,j]<-W[1,,,j]%*%Pp[t,,,i,j]%*%t(W[1,,,j])+R[1,,,j]
          ## THUS, THE UPDATED STATE IS IDENTICAL TO THE PREDICTED STATE
          ## AND THE ASSOCIATED COVARIANCE MATRIX IS CALCULATED AS USUAL
          # au[t,,i,j]<-ap[t,,i,j]
          # Pu[t,,,i,j]<-Pp[t,,,i,j]+Pp[t,,,i,j]%*%t(W[1,,,j])%*%solve(F[t,,,i,j])%*%W[1,,,j]%*%Pp[t,,,i,j]
          ## END OF KALMAN FILTER STEP
          ## BEGIN HAMILTON FILTER STEP
          ## CONDITIONAL TRANSITION PROBABILITY FOR S(t-1)=i to S(t)=j
          # tpr[t,j,i]<-p[j,i]*pr[t,i]
          ## CONDITIONAL DENSITY CAN NOT BE CALCULATED
        # }
      # }
      ## CURRENT OCCASION DOES NOT CONTRIBUTE TO THE LIKELIHOOD FUNCTION
      # utpr[t,,]<-tpr[t,,]
      ## END OF THE HAMILTON FILTER STEP
      ## BEGIN COLLAPSE STEP to reduce the number of posteriors
      # for (j in 1:nm)
      # {  pr[(t+1),j]<-sum(utpr[t,j,])
        # as<-matrix(c(0),ne,1)
        # for (i in 1:nm)
        # {  as<-as+au[t,,i,j]*(utpr[t,j,i]/pr[(t+1),j])
        # }
        # auc[(t+1),,j]<-as
        # Ps<-matrix(c(0),ne,ne)
        # for (i in 1:nm)
        # {  Ps<-Ps+(Pu[t,,,i,j]+(auc[(t+1),,j]-au[t,,i,j])%*%(t(auc[(t+1),,j]-au[t,,i,j])))*(utpr[t,j,i]/pr[(t+1),j])
        # }
        # Puc[(t+1),,,j]<-Ps
      # }
    # }
  # }  # END OF TIME LOOP
Lp<-if(is.nan(L)) 1e9 else min(-L,1e9)
attr(Lp,"jd") <- jd
attr(Lp,"all") <- list(jd=jd, Puc=Puc, Pu=Pu, Pp=Pp, v=v, F=F, pr=pr, tpr=tpr, utpr=utpr, au=au, ap=ap, auc=auc, md=md)
Lp
}



`filter.old` <-
  function(theta,npar,ny,ne,nx,nt,nm,mdim,MA,PA,a0,P0,y,x, debug=FALSE)
  {
  	Debug = debug
    # START WITH MAKING THE ARRAYS THAT ARE NEEDED
    auc  <- array(0,c((nt + 1),ne,nm))     # a-updated and collapsed to nm posteriors
    au   <- array(0,c((nt),ne,nm,nm))      # a-updated (not collapsed: S(t-1)=i, S(t)=j)
    ap   <- array(0,c((nt),ne,nm,nm))      # a-predicted (by definition not collapsed)
    Puc  <- array(0,c((nt + 1),ne,ne,nm))  # covariance matrix of auc
    Pu   <- array(0,c(nt,ne,ne,nm,nm))     # covariance matrix of au
    Pp   <- array(0,c(nt,ne,ne,nm,nm))     # covariance matrix of ap
    v    <- array(0,c(nt,ny,nm,nm))        # one-step-ahead prediction error for S(t-1)=1, S(t)=j
    F    <- array(0,c(nt,ny,ny,nm,nm))     # covariance matrix of v
    pr   <- array(0,c((nt + 1),nm))        # updated probabilities
    tpr  <- array(0,c(nt,nm,nm))           # conditional transition probabilites
    utpr <- array(0,c(nt,nm,nm))           # updated conditional transition probabilites
    jd   <- array(0,c(nt,nm,nm))           # joint densities (product of conditional density and transition probability)
    md   <- array(0,c(nt,1))               # marginal density (used in optimization)

    # MAKE UPDATED PROABILITIES (I.E., VECTOR WITH PROBABILITIES OF
    # OF BEING IN CERTAIN STATE AT OCCASION 0)
    pr[1,] <- matrix(c(1 / nm),nm,1)
    # ENTER INITIAL VALUES OF a0 AND P0 IN auc
    auc[1,,] <- a0
    Puc[1,,,] <- P0
    # PLACE THE ELEMENTS FROM x (= vector with free parameters)
    # IN THE APPROPRIATE MODEL MATRICES
    for (m in 1:7)
    {
      for (s in 1:nm)
      {
        for (r in 1:mdim[m,1])
        {
          for (c in 1:mdim[m,2])
          {
            p <- PA[m,r,c,s]
            if (p > 0) {
              MA[m,r,c,s] <- theta[p]
            }
          }
        }
      }
    }
    # TO OBTAIN THE MATRIX WITH TRANSITION PROBABILITIES,
    # WE HAVE TO TAKE INTO ACCOUNT THAT THE PORBS HAVE TO SUM TO 1 PER COLUMN
    # HENCE, PER COLUMN THERE ARE nm-1 PROBS TO BE ESTIMATED
    # RATHER THAN ESTIMATING THE PROBS WE WILL ESTIMATE A NUMBER
    # WHICH CAN RUN FROM -INFINITY TO +INFINITY
    tp <- matrix(0, nm, nm)
    for (c in 1:(nm - 1))
    {
      tp[nm,c] <- exp(0)
      for (r in 1:(nm - 1))
      {
        npar <- npar + 1
        tp[r,c] <- exp(theta[npar])
      }
    }
    tp[1,nm] <- exp(0)
    for (r in 2:nm)
    {
      npar <- npar + 1
      tp[r,nm] <- exp(theta[npar])
    }
    # THE ZEROS ARE INCLUDED FOR IDENTIFICATION
    # TO COME FROM THE ABOVE MATRIX tp TO THE MATRIX WITH ACTUAL
    # TRANSITION PROBS p, WE DO THE FOLLOWING
    p <- matrix(0,nm,nm)
    for (c in 1:nm)
    {
      stp <- sum(tp[,c])
      p[,c] <- tp[,c] / stp
    }
    # SELECT THE MODEL MATRICES AND NAME THEM
    W <- MA[1,1:ny,1:ne,,drop = FALSE]
    B <- MA[2,1:ny,1:nx,,drop = FALSE]
    R <- MA[3,1:ny,1:ny,,drop = FALSE]
    c <- MA[4,1:ne,1,,drop = FALSE]
    H <- MA[5,1:ne,1:ne,,drop = FALSE]
    G <- MA[6,1:ne,1:ne,,drop = FALSE]
    K <- MA[7,1:ne,1:ne,,drop = FALSE]
    # START THE KALMAN FILTER & REGIME SWITCHING
    L <- 0
#     L <- Kfilter.timeloop(auc, au, ap,
#                           Puc, Pu, Pp,
#                           v, F, p,
#                           pr, tpr, utpr, jd, md,
#                           W, B, R, c, H, G, K, y, x, debug = debug)

    ### ----------------- debug code ------------------- ###
    if (debug) {
      cat("x[0,1] = ", x[0,1], "\n")
      cat(paste(c("nt","ne","nm","ny","nx"),c(nt,ne,nm,ny,nx),sep=" = ", collapse = "\n"))
      cat("\np =\n"); print(p)
      for (j in 1:nm) {
        cat(sprintf("\n\nRegime %d\n========\n\n\n\ta[t] = c + H a[t-1] + G u[t],     u[t] ~ N(0, K)\n\ty[t] = W a[t] + B x[t] + e[t],    e[t] ~ N(0, R)\n\n", j-1));
        cat("\nH = \n"); print(as.matrix(H[1,,,j]));
        cat("\nK =\n"); print(as.matrix(K[1,,,j]));
        cat("\nG =\n"); print(as.matrix(G[1,,,j]));
        cat("\nc =\n"); print(as.matrix(c[1,,,j]));
        cat("\nB =\n"); print(as.matrix(B[1,,,j]));
        cat("\nR =\n"); print(as.matrix(R[1,,,j]));
        cat("\nW =\n"); print(as.matrix(W[1,,,j]));
      }
      cat("-----------------------------------------------------------------------------------------------------------\n");
      Debug = debug;
    }
    ### ----------------- end debug ------------------- ###

    for (t in 1:nt)
    {
      ### ----------------- debug code ------------------- ###
      if (t >= debug) debug = 0;
      if (debug) {
         cat(sprintf("\nt = %d\n", t-1));
      }
      ### ----------------- end debug ------------------- ###


      # DETERMINE WHETHER THERE ARE MISSING VALUES
      if (all(!is.na(y[t,])))
      {
        for (j in 1:nm)
        {
          for (i in 1:nm)
          {
            # KALMAN FILTER STEP
            # PREDICTION EQUATIONS
            ### ----------------- debug code ------------------- ###
              if(debug && i==1 && j==1) {
                cat("\np = \n"); print(p);
              }
            ### ----------------- end debug ------------------- ###

            ap[t,,i,j] <- c[1,,1,j] + H[1,,,j] %*% (auc[t,,i])
            ### ----------------- debug code ------------------- ###
              if(debug) {
                cat("\nap = c + H %*% auc\n");
                cat("\nc = \n"); print(as.matrix(c[1,,,j]));
                cat("\nH = \n"); print(as.matrix(H[1,,,j]));
                cat("\nauc = \n"); print(as.matrix(auc[t,,j]));
                cat("\nap = \n"); print(as.matrix(ap[t,,i,j]));
              }
            ### ----------------- end debug ------------------- ###

            Pp[t,,,i,j] <-
              H[1,,,j] %*% Puc[t,,,i] %*% t(H[1,,,j]) + G[1,,,j] %*% K[1,,,j] %*% t(G[1,,,j])
            ### ----------------- debug code ------------------- ###
            if(debug) cat(sprintf("\n S(i=%d) --> S(j=%d) at t = %d\t\t***********************************************\n\nPp = H Puc H' + G K G'", i-1, j-1, t-1)); # states are labeled 0, 1, 2, etc
            if(debug){
              cat("\n");
              cat("\nPuc = \n"); print(as.matrix(Puc[t,,,j]));
              cat("\nG = \n"); print(as.matrix(G[1,,,j]));
              cat("\nK = \n"); print(as.matrix(K[1,,,j]));
              cat("\nPp = \n"); print(as.matrix(Pp[t,,,i,j]));
            }
            ### ----------------- end debug ------------------- ###

            # ONE-STEP-AHEAD PREDICTION ERROR
            v[t,,i,j] <- y[t,] - W[1,,,j] %*% ap[t,,i,j] - as.matrix(B[1,,,j]) %*% x[t,]
            F[t,,,i,j] <- W[1,,,j] %*% Pp[t,,,i,j] %*% t(W[1,,,j]) + R[1,,,j]
            ### ----------------- debug code ------------------- ###
            if(debug){
              cat("\nv1 = y - W ap\n");
              cat("\ny = \n");     print(as.matrix(y[t,]))
              cat("\nW = \n");     print(as.matrix(W[1,,,j]));
              cat("\nap = \n");    print(as.matrix(ap[t,,i,j]));
              cat("\nv1 = \n");    print(matrix(y[t,] - W[1,,,j] %*% ap[t,,i,j]));
              cat("\nv = v1 - B x\n");
              cat("\nB = \n");     print(as.matrix(B[1,,,j]));
              cat("\nx = \n");     print(as.matrix(x[t,]))
              cat("\nv = \n");    print(as.matrix(v[t,,i,j]));
              cat("\nF = W Pp W' + R\n");
              cat("\nPp = \n");    print(as.matrix(Pp[t,,,i,j]));
              cat("\nPpW = \n");   print(as.matrix(Pp[t,,,i,j] %*% t(W[1,,,j])));
              cat("\nR = \n");     print(as.matrix(R[1,,,j]));
              cat("\nF = \n");     print(as.matrix(F[t,,,i,j]));
            }
            ### ----------------- end debug ------------------- ###


            # UPDATE EQUATIONS
            au[t,,i,j] <-
              ap[t,,i,j] + Pp[t,,,i,j] %*% t(W[1,,,j]) %*% solve(F[t,,,i,j]) %*% v[t,,i,j]
            Pu[t,,,i,j] <-
              Pp[t,,,i,j] + Pp[t,,,i,j] %*% t(W[1,,,j]) %*% solve(F[t,,,i,j]) %*% W[1,,,j] %*% Pp[t,,,i,j]

            ### ----------------- debug code ------------------- ###
            if(debug){
            	cat("\n");
            	cat("\nFinv =\n");    print(as.matrix(solve(F[t,,,i,j])))
            	cat(sprintf("\n|F| = %12.6g\n", det(F[t,,,i,j])));
              cat("\nau = ap + Pp W Finv v\n");
            	cat("\nPpWFinv = \n");   print(as.matrix(Pp[t,,,i,j] %*% t(W[1,,,j]) %*% solve(F[t,,,i,j])));
            	cat("\nau = \n");   print(as.matrix(au[t, , i, j]));
	            cat("Pu = Pp - Pp W Finv W' Pp\n");
	            cat("\nPu = \n");   print(as.matrix(Pu[t,,,i,j]));
            }
            ### ----------------- end debug ------------------- ###


            # END OF KALMAN FILTER STEP
            # BEGIN HAMILTON FILTER STEP
            # CONDITIONAL TRANSITION PROBABILITY FOR S(t-1)=i to S(t)=j
            tpr[t,j,i] <- p[j,i] * pr[t,i]
            # CONDITIONAL DENSITY
            jd[t,j,i] <-
              tpr[t,j,i] * ((2 * pi) ^ (-ny / 2) * det(as.matrix(F[t,,,i,j]))^(-1/2) * exp((-1/2) * t(v[t,,i,j]) %*% solve(F[t,,,i,j]) %*% v[t,,i,j]))
            ### ----------------- debug code ------------------- ***/
            if (Debug) {
            	cat(sprintf("\n        v = %12.6g\nv' Finv v = %12.6g\n", v[t,1,i,j], t(v[t,,i,j]) %*% solve(F[t,,,i,j]) %*% v[t,,i,j]));
            	cat(sprintf("\njd[%d,%d] = tpr[%d,%d] * (2 * pi)^(ny/2) * exp(-0.5 * v' Finv v) / sqrt(det(F))", j, i, j, i));
            	cat(sprintf("\n\t             p = %12.6g\n\t            pr = %12.6g\n\t           tpr = %12.6g\n\t(2*pi)^(-ny/2) = %12.6g\n\t    vtFinvv[0] = %12.6g\n\t          detF = %12.6g\n\t       jd[%d,%d] = %12.6g\n", p[j, i], pr[t, i], tpr[t,j,i], (2*pi)^(-ny/2.0), t(v[t,,i,j]) %*% solve(F[t,,,i,j]) %*% v[t,,i,j], det(as.matrix(F[t,,,i,j])), j, i, jd[t, j, i]));
            }
            ### ----------------- end debug ------------------- ***/
          }
        }
        jm <- sum(jd[t,,])
        L <- L + log(jm)        # LOG LIKELIHOOD FUNCTION
        utpr[t,,] <- jd[t,,] / jm

        ### ----------------- debug code ------------------- ***/
        if (debug) {
          cat(sprintf("\n\n              +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\t       missing = %12s\n\tlog-Likelihood = %12.6g\n\tjm             = %12.6g\n\n", if (any(is.na(y[t,]))) "TRUE" else "FALSE", L, jm));
        }
        ### ----------------- end debug ------------------- ***/


        # END OF THE HAMILTON FILTER STEP
        # BEGIN COLLAPSE STEP to reduce the number of posteriors
        for (j in 1:nm)
        {
          pr[(t + 1),j] <- sum(utpr[t,j,])
          as <- matrix(c(0),ne,1)
          for (i in 1:nm)
          {
            as <- as + au[t,,i,j] * (utpr[t,j,i] / pr[(t + 1),j])
          }
          auc[(t + 1),,j] <- as
          Ps <- matrix(c(0),ne,ne)
          for (i in 1:nm)
          {
            Ps <-
              Ps + (Pu[t,,,i,j] + (auc[(t + 1),,j] - au[t,,i,j]) %*% (t(auc[(t + 1),,j] -
                                                                          au[t,,i,j]))) * (utpr[t,j,i] / pr[(t + 1),j])
          }
          Puc[(t + 1),,,j] <- Ps
        }
      }

    ###
    ### If there are missing values we have a separate time loop... (for no reason)
    ###

      if (any(is.na(y[t,])))
      {
        ### ----------------- debug code ------------------- ###
        if(debug)
          cat(sprintf("missing value(s) at t = %d\n", t));
        ### ----------------- end  debug ------------------- ###

        for (j in 1:nm)
        {
          for (i in 1:nm)
          {
            # KALMAN FILTER STEP
            # PREDICTION EQUATIONS
            ap[t,,i,j] <- c[1,,1,j] + H[1,,,j] %*% (auc[t,,i])
            Pp[t,,,i,j] <-
              H[1,,,j] %*% Puc[t,,,i] %*% t(H[1,,,j]) + G[1,,,j] %*% K[1,,,j] %*% t(G[1,,,j])
            # ONE-STEP-AHEAD PREDICTION ERROR CAN NOT BE CALCULATED
            # BUT THE ASSOCIATED COVARIANCE MATRIX CAN
            F[t,,,i,j] <- W[1,,,j] %*% Pp[t,,,i,j] %*% t(W[1,,,j]) + R[1,,,j]
            # THUS, THE UPDATED STATE IS IDENTICAL TO THE PREDICTED STATE
            # AND THE ASSOCIATED COVARIANCE MATRIX IS CALCULATED AS USUAL
            au[t,,i,j] <- ap[t,,i,j]
            Pu[t,,,i,j] <-
              Pp[t,,,i,j] + Pp[t,,,i,j] %*% t(W[1,,,j]) %*% solve(F[t,,,i,j]) %*% W[1,,,j] %*% Pp[t,,,i,j]
            # END OF KALMAN FILTER STEP
            # BEGIN HAMILTON FILTER STEP
            # CONDITIONAL TRANSITION PROBABILITY FOR S(t-1)=i to S(t)=j
            tpr[t,j,i] <- p[j,i] * pr[t,i]
            # CONDITIONAL DENSITY CAN NOT BE CALCULATED
          }
        }
        # CURRENT OCCASION DOES NOT CONTRIBUTE TO THE LIKELIHOOD FUNCTION
        utpr[t,,] <- tpr[t,,]
        # END OF THE HAMILTON FILTER STEP
        # BEGIN COLLAPSE STEP to reduce the number of posteriors
        for (j in 1:nm)
        {
          pr[(t + 1),j] <- sum(utpr[t,j,])
          as <- matrix(c(0),ne,1)
          for (i in 1:nm)
          {
            as <- as + au[t,,i,j] * (utpr[t,j,i] / pr[(t + 1),j])
          }
          auc[(t + 1),,j] <- as
          Ps <- matrix(c(0),ne,ne)
          for (i in 1:nm)
          {
            Ps <-
              Ps + (Pu[t,,,i,j] + (auc[(t + 1),,j] - au[t,,i,j]) %*% (t(auc[(t + 1),,j] -
                                                                          au[t,,i,j]))) * (utpr[t,j,i] / pr[(t + 1),j])
          }
          Puc[(t + 1),,,j] <- Ps
        }
      }
    }  # END OF TIME LOOP
    Lp<-if(is.nan(L)) 1e9 else min(-L,1e9)
    attr(Lp,"jd") <- jd
    attr(Lp,"all") <- list(jd=jd, Puc=Puc, Pu=Pu, Pp=Pp, v=v, F=F, pr=pr, tpr=tpr, utpr=utpr, au=au, ap=ap, auc=auc, md=md)
    Lp
  }
