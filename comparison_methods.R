
########################################################
########################################################
# Proj-Huber extended from Wong and Lee (2017)
########################################################
########################################################
R_Lee <-function(x2, r, lambda= 0.2 ){
  
  Psi  = function(M,c){
    M.2 = 2*M
    ind1 =( M>c)
    M.2[ind1]= 2*c
    ind2 = (-M>c)
    M.2[ind2]=-2*c
    
    return(M.2)
  }
  
  X=x2
  c= 1.345*sqrt(var(as.vector(x2)) )
  eps = 10^-5
  maxiter = 200
  cond=TRUE
  k=0
  r.huber = r
  nrob = softImpute(x2,rank=r.huber,lambda=lambda)
  nrob.recon = nrob$u%*%diag(nrob$d)%*%t(nrob$v)
  Y.1 = nrob.recon
  err=c()
  while(cond){
    k=k+1
    Y.0=Y.1
    E = X - Y.0 
    Z = Y.0 + Psi(E, c=c)/2
    tmp = softImpute(Z,rank=r.huber,lambda=lambda)
    Y.1 = tmp$u%*%diag(tmp$d)%*%t(tmp$v)
    
    cond = !(sum((Y.0 - Y.1)^2)/sum((Y.0^2)) < eps | k > maxiter)  
    
    err=c(err,sum((Y.0 - Y.1)^2)/sum((Y.0^2)))
  }
  
  output<-list()
  output$U = tmp$u
  output$V = tmp$v
  output$Dsq = tmp$d
  output$rec = Y.1
  
  return(output)
  
}

########################################################
########################################################
# Robust FPCA by Boente and Salibian Barrera (2021) 
########################################################
########################################################

cov_boente <- function(x,gr, alpha, bw.mu, bw.cov, cv = FALSE, ncores = 1, seed = 123) {
  X <- list(x = x$Ly, pp = x$Lt)
  
  start_time <- Sys.time()
  
  # Start cluster
  if (isTRUE(cv)) {
    
    hs.mu <- seq(0.02, 0.3, length.out = 5)   # candidate of bw of mu
    hs.cov <- seq(0.02, 0.3, length.out = 5)   # candidate of bw of cov
    k.cv <- 5
    rho.param <- 1e-3
    k <- 3
    s <- k
    ncov <- length(gr)
    
    n_cores <- ncores
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    
    # run CV to find smoothing parameters for mean and covariance function
    aa <- cv.mu.par(X, alpha=alpha, hs=hs.mu, seed=seed, k.cv=k.cv)
    bw.mu <- aa$h[ which.min(aa$tmse) ]
    mh <- mu.hat3.lin(X=X, h=bw.mu)
    
    # covariance
    bb <- cov.fun.cv.par(X=X, muh=mh, ncov=ncov, k.cv=k.cv, hs=hs.cov,
                         alpha=alpha, seed=seed, k=k, s=s, reg.rho=rho.param)[1:2]
    bw.cov <- bb$h[ which.min(bb$tmspe) ]
    ma <- matrixx(X, mh)
    
    stopCluster(cl)
  } else {
    mh <- mu.hat3.lin(X=X, h=bw.mu)
    ma <- matrixx(X, mh)
  }
  
  # Compute the estimated cov function
  cov.fun2 <- cov.fun.hat2(X=X, h=bw.cov, mh=mh, ma=ma, ncov=length(gr), trace=FALSE)
  
  end_time <- Sys.time()
  print(paste0("boente stage : ", 
               round(difftime(end_time, 
                              start_time, 
                              units = "secs"), 3),
               " secs"))
  
  df <- data.frame(t = unlist(X$pp),
                   mu = unlist(mh))
  df <- unique(df)
  idx <- sort(df$t, index.return = T)$ix
  # mu <- df$mu[idx]
  mu <- ConvertSupport(fromGrid = df$t[idx], 
                       toGrid = gr,
                       mu = df$mu[idx])
  
  # noise variance
  noise_var <- eigen(cov.fun2$G)$values[1] / (1e3 - 1)
  
  
  return(list(mu = mu,
              cov = cov.fun2$G,
              noise_var = noise_var))
}


#################################################################
#################################################################
# partial FPCA by Kraus (2015) - code borrowed by the author
################################################################
#################################################################

### (1) mean function, covariance operator and its eigendecomposition


mean.missfd = function(x)
{
  colMeans(x,na.rm=TRUE)
}

var.missfd = function(x, make.pos.semidef=TRUE)
{
  R = var(x,use="pairwise.complete.obs")
  if (make.pos.semidef) R = make.var.pos.semidef(R)
  R
}

cov.missfd = var.missfd

make.var.pos.semidef = function(R)
{
  eig.R = eigen(R,symmetric=TRUE)
  w = which(eig.R$values>0)
  R = eig.R$vectors[,w]%*%(eig.R$values[w]*t(eig.R$vectors[,w]))
  R
}

eigen.missfd = function(R,d.t,...)
{
  # d.t is the distance between the points of the equidistant grid
  out = eigen(R,symmetric=TRUE,...)
  if (missing(d.t)) d.t = 1/ncol(R)
  out$values = out$values*d.t
  out$vectors = out$vectors/sqrt(d.t)
  out
}

### (2) prediction of scores of a partially observed function


pred.score.missfd = function(x1,phi,x,R,mu,n,C,ncomp,alpha,d.t,
                             gcv.df.factor=1,gcv.plot=FALSE,
                             gcv.print=FALSE)
{
  # x1 is one partially observed curve
  # this function predicts the Fourier scores of x1
  # with respect to phi (i.e., inner products of x1 and phi)
  # phi is a matrix containing functions in columns
  # if missing, phi defaults to the first four principal components
  # x is a functional data set (functions in rows)
  # R is the covar oper, mu the mean fun, n the sample size of x
  # C, ncomp are the covar oper and sample size for the set of
  # complete functions
  # R, mu, n, C, ncomp are optional, may be precomputed
  # beforehand, e.g., to save time when this function is to be
  # called repeatedly (for different incomplete x1)
  # alpha is the ridge regularisation parameter, determined by gcv
  # from complete curves if missing
  # d.t is the distance between the points of the equidistant grid
  miss = is.na(x1)
  obs = !miss
  k.O = sum(obs)
  k.M = sum(miss)
  k = k.O + k.M
  if (missing(d.t)) d.t = 1/k
  if (missing(phi)) {
    R = var.missfd(x)
    phi = eigen.missfd(R)$vectors[,1:4]
  }
  phi = as.matrix(phi)
  nscores = ncol(phi)
  if (missing(mu)) {
    mu = mean.missfd(x)
    n = nrow(x)
  }
  if (missing(R)) {
    R = var.missfd(x)
    n = nrow(x)
  }
  if (missing(C)) {
    comp = apply(!is.na(x),1,all)
    ncomp = sum(comp)
    C = var.missfd(x[comp,]) # covariance operator of complete curves
  }
  
  if (k.M==0) { # complete observation, no prediction needed
    b = rep(NA,nscores)
    for (j in 1:nscores) {
      b[j] = sum((x1[obs]-mu[obs])*phi[obs,j])*d.t
    }
  } else {
    b = rep(NA,nscores)
    b.details = matrix(NA,7,nscores)
    rownames(b.details) = c("score.obs","score.miss","se","relerr","alpha","df.alpha","varprop.alpha")
    
    # precompute quantities that don't depend on alpha
    # to speed up the optimisation in gcv
    cMM = colSums(phi[miss,,drop=FALSE]*(C[miss,miss,drop=FALSE]%*%phi[miss,,drop=FALSE]))*(d.t^2)
    eig.ROO = eigen.missfd(R[obs,obs,drop=FALSE],d.t=d.t)
    lambda.ROO = eig.ROO$values
    psi.ROO = eig.ROO$vectors
    # rO = crossprod(R[obs,miss,drop=FALSE],phi[miss,,drop=FALSE])*d.t # columns = righ-hand sides of the inverse problem
    psirO = crossprod(psi.ROO,(R[obs,miss,drop=FALSE]%*%phi[miss,,drop=FALSE])*d.t)*d.t # projection of the right-hand sides on psi.ROO
    psirO.psicO = psirO*crossprod(psi.ROO,(C[obs,miss,drop=FALSE]%*%phi[miss,,drop=FALSE])*d.t)*d.t
    psirOrOpsi.psiCOOpsi = array(0,c(k.O,k.O,nscores))
    for (j in 1:nscores) {
      psirOrOpsi.psiCOOpsi[,,j] = tcrossprod(psirO[,j])*crossprod(psi.ROO,C[obs,obs,drop=FALSE]%*%psi.ROO)*d.t^2
    }
    
    if (missing(alpha)) { # gcv selection of alpha for each score
      alpha = rep(NA,nscores)
      if (gcv.plot) {
        a = ceiling(sqrt(nscores)) # number of columns in the plot
        op = par(mfrow=c(ceiling(nscores/a),a))
      }
      for (j in 1:nscores) {
        # gcv function for scores and for functions is the same
        # just use it with the right input parameters
        alpha[j] = alpha.gcv.pred.missfd(trace.CMM=cMM[j],
                                         diag.psiROMCMOpsi=psirO.psicO[,j],
                                         psiROMMOpsi.psiCOOpsi=psirOrOpsi.psiCOOpsi[,,j],
                                         lambda.ROO=lambda.ROO,n=n,
                                         ncomp=ncomp,df.factor=gcv.df.factor,
                                         gcv.plot=gcv.plot,gcv.print=gcv.print)
        if (gcv.plot) title(main=j)
      }
      if (gcv.plot) par(op)
    } else {
      if (length(alpha)==1) alpha = rep(alpha,nscores)
    }
    
    rMM = colSums(phi[miss,,drop=FALSE]*(R[miss,miss,drop=FALSE]%*%phi[miss,,drop=FALSE]))*(d.t^2)
    for (j in 1:nscores) {
      b.details["score.obs",j] = sum((x1[obs]-mu[obs])*phi[obs,j])*d.t
      b.details["score.miss",j] = sum(psirO[,j]*crossprod(psi.ROO,x1[obs]-mu[obs])/(lambda.ROO+alpha[j]))*d.t
      b[j] = b.details["score.obs",j] + b.details["score.miss",j]
      b.details["se",j] = sqrt( rMM[j] - sum(psirO[,j]^2*lambda.ROO/(lambda.ROO+alpha[j])^2) )
      b.details["relerr",j] = b.details["se",j]/sqrt(sum(phi[,j]*(R%*%phi[,j]))*(d.t)^2)
      b.details["alpha",j] = alpha[j]
      b.details["df.alpha",j] = sum(lambda.ROO/(lambda.ROO+alpha[j]))
      b.details["varprop.alpha",j] = sum(lambda.ROO^3/(lambda.ROO+alpha[j])^2)/sum(lambda.ROO)
    }
    attr(b,"details") = b.details
  }
  
  b  
}

### (3) prediction of the missing part of an incomplete function

pred.missfd = function(x1,x,R,mu,n,C,ncomp,alpha,d.t,
                       gcv.df.factor=1,gcv.plot=FALSE,gcv.print=FALSE)
{
  # x1 is one partially observed curve
  # this function predicts the missing part of x1
  # x is a functional data set
  # R is the covar oper, mu the mean fun, n the sample size of x
  # C, ncomp covar oper and sample size for the set of
  # complete functions
  # R, mu, n, C, ncomp are optional, may be precomputed
  # beforehand, e.g., to save time when this function is to be
  # called repeatedly (for different incomplete x1)
  # alpha is the ridge regularisation parameter, determined by gcv
  # from complete curves if missing
  # d.t is the distance between the points of the equidistant grid
  miss = is.na(x1)
  obs = !miss
  k.O = sum(obs)
  k.M = sum(miss)
  k = k.O + k.M
  if (missing(d.t)) d.t = 1/k
  if (missing(mu)) {
    mu = mean.missfd(x)
    n = nrow(x)
  }
  if (missing(R)) {
    R = var.missfd(x)
    n = nrow(x)
  }
  if (missing(C)) {
    comp = apply(!is.na(x),1,all)
    ncomp = sum(comp)
    C = var.missfd(x[comp,]) # covariance operator of complete curves
  }
  
  x1.pred = rep(NA,k) # this will contain predicted x1 on miss and NAs on obs
  x1.pred.covar = matrix(NA,k,k) # covar oper of the predictive distribution
  
  if (k.M>0) {
    # precompute quantities that don't depend on alpha
    # to speed up the optimisation in gcv
    trace.CMM = sum(diag(C)[miss])*d.t
    eig.ROO = eigen.missfd(R[obs,obs,drop=FALSE],d.t=d.t)
    lambda.ROO = eig.ROO$values
    psi.ROO = eig.ROO$vectors
    psiROM = crossprod(psi.ROO,R[obs,miss,drop=FALSE])*d.t
    diag.psiROMCMOpsi = rowSums(psiROM*crossprod(psi.ROO,C[obs,miss,drop=FALSE])*d.t)*d.t
    psiROMMOpsi.psiCOOpsi = tcrossprod(psiROM)*d.t*crossprod(psi.ROO,C[obs,obs,drop=FALSE]%*%psi.ROO)*d.t^2
    
    # if alpha not provided on input, use gcv
    if (missing(alpha)) {
      alpha = alpha.gcv.pred.missfd(trace.CMM=trace.CMM,
                                    diag.psiROMCMOpsi=diag.psiROMCMOpsi,
                                    psiROMMOpsi.psiCOOpsi=psiROMMOpsi.psiCOOpsi,
                                    lambda.ROO=lambda.ROO,n=n,
                                    ncomp=ncomp,df.factor=gcv.df.factor,
                                    gcv.plot=gcv.plot,gcv.print=gcv.print)
    }
    
    # compute the regularised prediction and predictive covariance operator
    x1.pred[miss] = crossprod(psiROM,crossprod(psi.ROO,x1[obs]-mu[obs])/(lambda.ROO+alpha))*d.t + mu[miss]
    x1.pred.covar[miss,miss] = R[miss,miss] - crossprod(psiROM*(lambda.ROO/(lambda.ROO+alpha)^2),psiROM)
    
    attributes(x1.pred) = list(covar=x1.pred.covar,se=sqrt(diag(x1.pred.covar)),
                               relerr=sqrt(sum(diag(x1.pred.covar)[miss])/sum(diag(R))),
                               alpha=alpha,df.alpha=sum(lambda.ROO/(lambda.ROO+alpha)),
                               varprop.alpha=sum(lambda.ROO^3/(lambda.ROO+alpha)^2)/sum(lambda.ROO))
  } else { # complete curve, no prediction
    x1.pred = x1
  }
  
  x1.pred
}

gcv.pred.missfd = function(log.alpha, trace.CMM, diag.psiROMCMOpsi,
                           psiROMMOpsi.psiCOOpsi, lambda.ROO, ncomp,
                           df.factor=1, prn=FALSE)
{
  # gcv on complete curves
  # this version is fully based on eigendecomposition
  alpha = exp(log.alpha)
  gof.term = trace.CMM - 2*sum(diag.psiROMCMOpsi/(lambda.ROO+alpha)) + sum(psiROMMOpsi.psiCOOpsi*tcrossprod(1/(lambda.ROO+alpha)))
  df.alpha = sum(lambda.ROO/(lambda.ROO+alpha))
  out = log(gof.term) - 2*log(1-df.factor*df.alpha/ncomp)
  if (prn) print(c(alpha=alpha,gof.term=gof.term,df.alpha=df.alpha,log.gcv=out))
  out
}

alpha.gcv.pred.missfd = function(trace.CMM, diag.psiROMCMOpsi,
                                 psiROMMOpsi.psiCOOpsi, lambda.ROO,
                                 n, ncomp, df.factor=1, gcv.plot=FALSE,
                                 gcv.print=FALSE)
{
  # scale (standardise) everything
  std = lambda.ROO[1]
  trace.CMM = trace.CMM/std
  lambda.ROO = lambda.ROO/std
  diag.psiROMCMOpsi = diag.psiROMCMOpsi/(std^2)
  psiROMMOpsi.psiCOOpsi = psiROMMOpsi.psiCOOpsi/(std^3)
  # find minimum alpha so that df <= df.max
  nn = sum(lambda.ROO>.Machine$double.eps^.5) # rank of ROO
  if (nn>=length(lambda.ROO)) { # ROO full rank, alpha small OK
    log.alpha.min = log(.Machine$double.eps^.5)
  } else {
    df.max = min(nn/4) # maximum df that we will allow in optim
    # now need to find alpha corresponding to df.max
    log.alpha.min = log(.Machine$double.eps^.5)
    m1 = log.alpha.min
    m2 = log(n/df.max)
    if (df.alpha.eq(log.alpha=log.alpha.min,deg.fr=df.max,lambda=lambda.ROO)>0) {
      # find log alpha between m1, m2 such that df equals df.max
      # df at m2 certainly < df.max
      # if df at m1 > df.max, there exists alpha such that df equals df.max
      # (opposite signs at end points), we find it using uniroot
      # (if same signs, we keep m1)
      log.alpha.min = ( uniroot(df.alpha.eq,c(m1,m2),deg.fr=df.max,lambda=lambda.ROO)$root )
    }
  }
  # minimise gcv
  log.alpha = optim(max(log(mean(lambda.ROO)),log.alpha.min),
                    gcv.pred.missfd,NULL,
                    trace.CMM=trace.CMM,
                    diag.psiROMCMOpsi=diag.psiROMCMOpsi,
                    psiROMMOpsi.psiCOOpsi=psiROMMOpsi.psiCOOpsi,
                    lambda.ROO=lambda.ROO,ncomp=ncomp,
                    df.factor=df.factor,method="L-BFGS-B",
                    # control=list(trace=0),
                    lower=log.alpha.min)$par
  # plot gcv on a grid of alpha values
  if (gcv.plot) {
    log.alpha.max = log.alpha+3
    log.alpha.grid = seq(log.alpha.min,log.alpha.max,len=40)
    log.gcv.grid = double(length(log.alpha.grid))
    for (i in 1:length(log.gcv.grid)) {
      log.gcv.grid[i] = gcv.pred.missfd(log.alpha=log.alpha.grid[i],
                                        trace.CMM=trace.CMM,
                                        diag.psiROMCMOpsi=diag.psiROMCMOpsi,
                                        psiROMMOpsi.psiCOOpsi=psiROMMOpsi.psiCOOpsi,
                                        lambda.ROO=lambda.ROO,ncomp=ncomp,
                                        df.factor=df.factor,
                                        prn=gcv.print) #+log(std)
    }
    plot(log.alpha.grid+log(std),log.gcv.grid+log(std),type="l",xlab="log(alpha)",ylab="log(gcv(alpha))")
    abline(v=log.alpha.min+log(std),lty=3)
    abline(v=log.alpha+log(std),col=2)
  }
  exp(log.alpha)*std
}

df.alpha.eq = function(log.alpha,deg.fr,lambda)
{
  sum(lambda/(lambda+exp(log.alpha))) - deg.fr
}

### (4) prediction bands for the missing part of an incomplete function

qsupgp2 = function(p,R,h=pmax(sqrt(diag(R)),.2*sqrt(max(diag(R)))),nsim=2500)
{
  # simulate the p-th quantile of sup(|X|/h) where X is a mean zero
  # Gaussian process with covariance function R
  eig.R = eigen(R,symmetric=TRUE)
  npos = sum(eig.R$values>0)
  y = (eig.R$vectors[,1:npos]*rep(sqrt(eig.R$values[1:npos]),each=nrow(R))) %*% matrix(rnorm(npos*nsim),npos,nsim) # simulated curves are in columns
  quantile(apply(abs(y)/h,2,max), probs=p)
}

gpband = function(mu,R,h=pmax(sqrt(diag(R)),.2*sqrt(max(diag(R),na.rm=TRUE))),coverage=.95)
{
  # for a Gaussian process with mean mu and covariance R
  # this function computes a band of the form mu +- u*h that
  # contains a trajectory with probability coverage,
  # h is the boundary function for the band (h=1 for const width,
  # default h reflects the sd but is bounded away from 0)
  # the band is constructed only in regions where mu is not NA
  w = !is.na(mu)
  if (length(h)==1) h = rep(h,length(mu))
  u = qsupgp2(p=coverage,R=R[w,w,drop=FALSE],h=h[w])
  out = cbind(lower=mu-u*h,upper=mu+u*h)
}

#' @export
gtshat.ls <-function(X, t0, s0, h, mh, matx, eps=1e-6){
  n <- length(X$x)
  err <- 1+eps
  j <- 0
  M <- matx$m
  MT <- matx$mt
  w2 <- k.epan((MT[,1]-t0)/h)
  w3 <- k.epan((MT[,2]-s0)/h)
  we <- w2*w3
  M <- M[ we > 0, ]
  # MM <- matrix(NA, nrow(M), 3)
  # MM[, 1] <- M[, 1]
  # MM[, 2] <- M[, 1] * (MT[ we > 0, 1]-t0)
  # MM[, 3] <- M[, 1] * (MT[ we > 0, 2]-s0)
  we <- we[ we > 0]
  if( length(we)==0 ) return(NA)
  # B <- as.numeric( coef( lm( M[, 2] ~ MM - 1, weights = we) )[1] )
  B <- as.numeric( coef( lm(M[ ,2] ~ M[, 1] - 1, weights = we) ) )
  # # get rid of points on the diagonal
  # #   xs <- M[,1]
  # #   ys <- M[,2]
  # #   tmp <- (xs != ys)
  # #   B <- median( M[tmp,2] / M[tmp,1 ])
  # #   sigma.hat <- mad( M[tmp,2] - B * M[tmp,1])
  # B <- median( M[,2] / M[,1 ])
  # sigma.hat <- mad( M[,2] - B * M[,1])
  # # B <- med.w( x=M[,2] / M[,1 ], w=we)
  # # sigma.hat <- mad.w( x=M[,2] - B * M[,1], w=we)
  # while ( ( (j <- j+1) < 1000 ) && (err > eps) ){
  #   w1 <- psi((M[,2]-B*M[,1])/sigma.hat,cc)/(M[,2]-B*M[,1])
  #   w1[ is.na(w1) ] <- 1/sigma.hat # missing values are "psi(0)/0", that lim is 1/sigma.hat
  #   w <- w1 * we
  #   B2 <- sum(w*M[,1]*M[,2])/sum(w*(M[,1]^2))
  #   err <- abs(B2/B - 1)
  #   if( is.na(err) || is.nan(err) ) return(NA)
  #   #   print('gtshat')
  #   #   print(B)
  #   #   print(summary(w1))
  #   #   print(summary(w))
  #   #   print(sigma.hat)
  #   # }
  #   B <- B2
  # }
  return(B)
}


#' @export
gtthat.ls <- function(X, t0, h, muhat, max.it=300, eps=1e-10) {
  # find the scale, full iterations
  # muhat should be a list as returned by mu.hat3
  i <- 0
  err <- 1 + eps
  t <- unlist(X$pp)
  x <- unlist(X$x)
  if(missing(muhat))
    muhat <- mu.hat3.lin.ls(X=X, h=h, ep=eps)
  muhat <- unlist(muhat)
  kerns <- k.epan((t-t0)/h)
  sc <- sqrt( sum(kerns * (x - muhat)^2 ) / (sum(kerns) - 1) )
  # if(missing(initial.sc)) sc <- mad(x - muhat) else sc <- initial.sc
  # while( ( (i <- i+1) < max.it ) && (err > eps) ) {
  #   kerns <- k.epan((t-t0)/h)
  #   sc2 <- sqrt(sc^2 * sum(kerns * rho((x-muhat)/sc,cc)) / (b * sum(kerns)))
  #   err <- abs(sc2/sc - 1)
  #   if( is.na(err) || is.nan(err) ) return(NA)
  #   #   print(summary(kerns))
  #   #   print(sc2)
  #   #   print(sc)
  #   # }
  #   sc <- sc2
  # }
  return(list(gtt=sc, muhat=muhat))
}

#' @export
uhat.ls <- function(X, t0, h=0.1, ep=1e-6, max.it=100){
  x <- unlist(X$x)
  t <- unlist(X$pp)
  # s <- localMAD(X,t0,h)
  # oldu <- (u <- median(x)) + 100*ep
  # it <- 0
  kw <- k.epan((t-t0)/h)
  return( sum( kw*x ) / sum(kw) )
  # while( ((it <- it + 1) < max.it ) && ( (abs(oldu) - u) > ep ) ){
  #   w <- psi((x-u)/s,cc)/(x-u)
  #   w[ is.nan(w) ] <- 1
  #   w <- w*kw
  #   oldu <- u
  #   u <- sum(w*x)/sum(w)
  # }
  # return(u)
}

#' @export
uhat.lin.ls <- function(X, t0, h=0.1, ep=1e-6, max.it=100){
  x <- unlist(X$x)
  t <- unlist(X$pp)
  # s <- localMAD(X,t0,h)
  # oldu <- (u <- median(x)) + 100*ep
  it <- 0
  tt <- cbind(rep(1, length(t)), t-t0)
  wk <- k.epan((t-t0)/h)
  return( coef( lm(x ~ I(t - t0), weights=wk) )[1] )
  #
  #   beta <- rlm(x=tt[wk>0,], y=x[wk>0])$coef
  # while( ((it <- it + 1) < max.it ) ){
  #   re <- as.vector(x - tt %*% beta)/s
  #   w <- ( psi(re,cc)/re )
  #   w[ is.nan(w) ] <- 1
  #   w <- w * wk
  #   beta.n <- solve( t( tt * w ) %*% tt, t(tt * w) %*% x )
  #   if( any( is.na(beta.n) || is.nan(beta.n) ) ) return(NA)
  #   if( sum( (beta.n - beta)^2 ) < ep ) it = max.it
  #   beta <- beta.n
  # }
  # return(beta.n[1])
}

#' @export
mu.hat2.ls <- function(X, h=0.1, ep=1e-6) {
  tt <- unlist(X$pp)
  nt <- length(tt)
  us <- rep(0, nt)
  for(i in 1:nt)
    us[i] <- uhat.ls(X=X, t0=tt[i], h=h, ep=ep)
  return(us)
}

#' @export
mu.hat2.lin.ls <- function(X, h=0.1, ep=1e-6) {
  tt <- unlist(X$pp)
  nt <- length(tt)
  us <- rep(0, nt)
  for(i in 1:nt)
    us[i] <- uhat.lin.ls(X=X, t0=tt[i], h=h, ep=ep)
  return(us)
  
}

#' @export
mu.hat3.ls <- function(X, h=0.1, ep=1e-6) {
  us=relist(mu.hat2.ls(X=X, h=h, ep=ep),X$x)
  return(us)
  #return(list(x=X$x,pp=X$pp,u=us))
}

#' @export
mu.hat3.lin.ls <- function(X, h=0.1, ep=1e-6) {
  us <- relist(mu.hat2.lin.ls(X=X, h=h, ep=ep), X$x)
  return(us)
  #return(list(x=X$x,pp=X$pp,u=us))
}

#' @export
cov.fun.hat2.ls <- function(X, h, mh, ma, ncov=50, trace=FALSE) {
  # this function uses the diagonal
  if(trace) print("Computing cov function")
  if(missing(ma)) ma <- matrixx(X=X, mh=mh)
  mii <- min( ti <- unlist(X$pp) )
  maa <- max( ti )
  tt <- seq(mii, maa, length=ncov)
  ss <- seq(mii, maa, length=ncov)
  pps <- as.matrix(expand.grid(tt, ss))
  np <- nrow(pps)
  betahat <- rep(0, np)
  sigmahat <- rep(0, ncov)
  for(j in 1:ncov) sigmahat[j] <- gtthat.ls(X=X, t0=tt[j], h=h, muhat=mh)$gtt
  for(j in 1:np) {
    t0 <- pps[j, 1]
    s0 <- pps[j, 2]
    betahat[j] <- gtshat.ls(X=X, t0=t0, s0=s0, h=h, mh=mh, matx=ma, eps=1e-6)
    # gamma[s0, t0] / gamma[t0, t0]
  }
  G <- betahat <- matrix(betahat, ncov, ncov)
  for(i in 1:ncov)
    for(j in 1:ncov)
      G[i,j] <- betahat[i,j] * sigmahat[i]^2
  G <- ( G + t(G) ) / 2
  if(trace) print("Done computing cov function")
  return(list(G=G, grid=pps))
}

#' @export
cv.mu.par.ls <- function(X, k.cv=5, k = k.cv, hs=exp(seq(-4, 0, by=.6)), seed=123) {
  # parallel computing version
  # Cross validation for the mean function
  lh <- length(hs)
  n <- length(X$x)
  if(exists(".Random.seed", where=.GlobalEnv)) old.seed <- .Random.seed
  set.seed(seed)
  fl <- sample( (1:n) %% k.cv + 1)
  tmses <- rep(NA, lh)
  Xtmp <- vector('list', 2)
  names(Xtmp) <- c('x', 'pp')
  tmp.par <- foreach(h=hs, .combine=c, .inorder=FALSE, .packages='MASS',
                     .export=c('uhat.lin.ls', 'k.epan')) %dopar% {
                       Xh <- relist(NA, X$x) # predictions go here
                       for(i in 1:k.cv) {
                         this <- (1:n)[ fl != i ]
                         Xtmp$x <- X$x[ this ] # training obs
                         Xtmp$pp <- X$pp[ this ] # training times
                         ps <- unlist(X$pp[ -this ]) # times where we need to predict
                         xhats <- rep(NA, length(ps))
                         for(l in 1:length(ps)) {
                           tmp2 <- try(uhat.lin.ls(X=Xtmp, t0=ps[l], h=h, ep=1e-6, max.it=100), silent=TRUE)
                           if( class(tmp2) != 'try-error' )
                             xhats[l] <- tmp2
                         }
                         Xh[-this] <- relist(xhats, X$x[ -this ] ) # fill predictions
                       }
                       # tmp <- sapply(mapply('-', X$x, Xh), function(a) a^2) # squared residuals, list-wise
                       tmp <- mapply('-', X$x, Xh)  # squared residuals, list-wise
                       # # tmses <- tm(unlist(tmp), alpha=alpha) # 40% trimming, use 60% smallest resids
                       # tmses <- RobStatTM::mscale(unlist(tmp)) #, delta=.3, tuning.chi=2.560841)
                       tmp2 <- unlist(tmp)
                       if(any(is.na(tmp2))) {
                         tmses <- NA } else {
                           # tmp2 <- tmp2[ !is.na(tmp2) ]
                           # if(length(tmp2) > 0) {
                           tmses <- sqrt(mean(tmp2^2)) #RobStatTM::mscale(tmp2) #, delta=.3, tuning.chi=2.560841)
                         }
                       
                     }
  if(exists('old.seed')) assign('.Random.seed', old.seed, envir=.GlobalEnv)
  return(list(tmse=tmp.par, h=hs))
}


#' @export
cov.fun.cv.par.ls <- function(X, muh, ncov=50, k.cv=5, hs=exp(seq(-4, 0, by=.6)),
                              seed=123, k=2, s=20, reg.rho=1e-5) {
  # parallel computing version
  # CV for the covariance function
  # muh = estimated mean function
  mii <- min( ti <- unlist(X$pp) )
  maa <- max( ti )
  tt <- seq(mii, maa, length=ncov)
  lh <- length(hs)
  n <- length(X$x)
  if(exists(".Random.seed", where=.GlobalEnv)) old.seed <- .Random.seed
  set.seed(seed)
  fl <- sample( (1:n) %% k.cv + 1)
  tmspe <- rep(NA, lh)
  Xtest <- Xtrain <- vector('list', 2)
  names(Xtrain) <- c('x', 'pp')
  names(Xtest) <- c('x', 'pp')
  tmp.par <- foreach(h=hs, .combine=c, .inorder=FALSE, .packages=c('MASS', 'mgcv'),
                     .export=c('cov.fun.hat2.ls', 'matrixx', 'subsets',
                               'gtthat.ls', 'gtshat.ls', 'mu.hat3.lin.ls',
                               'k.epan', 'pred.cv', 'L2.norma.mesh', 'L2.dot.product.mesh',
                               'integral')) %dopar% {
                                 Xhat <- relist(NA, X$x)
                                 for(i in 1:k.cv) {
                                   this <- (1:n)[ fl != i ] # training set
                                   Xtrain$x <- X$x[ this ]
                                   Xtrain$pp <- X$pp[ this ]
                                   Xtest$x <- X$x[ -this ]
                                   Xtest$pp <- X$pp[ -this ]
                                   ma <- matrixx(Xtrain, muh[ this ])
                                   cov.fun <- try(cov.fun.hat2.ls(X=Xtrain, h=h, mh=muh[ this ],
                                                                  ma=ma, ncov=50, trace=FALSE)) # $G
                                   if( class(cov.fun) != 'try-error') {
                                     if(!any(is.na(cov.fun$G))) {
                                       uu <- as.vector(cov.fun$G) #
                                       ttx <- cov.fun$grid #
                                       cov.fun <- matrix(fitted(gam(uu ~ s(ttx[,1], ttx[,2]))), length(tt), length(tt)) #
                                       cov.fun <- ( cov.fun + t(cov.fun) ) / 2
                                       tmp <- try( pred.cv(X=Xtrain, muh=muh[ this ], X.pred=Xtest,
                                                           muh.pred=muh[ -this ], cov.fun=cov.fun, tt=tt,
                                                           k=k, s=s, rho=reg.rho) )
                                       if( class(tmp) != 'try-error') Xhat[ -this ] <- tmp
                                     }
                                   }
                                 }
                                 # tmp <- sapply(mapply('-', X$x, Xhat), function(a) a^2) # squared residuals, list-wise
                                 tmp <- mapply('-', X$x, Xhat) # squared residuals, list-wise
                                 # tmspe <- tm(unlist(tmp), alpha=alpha) # 40% trimming, use 60% smallest resids
                                 # tmspe <- RobStatTM::mscale(unlist(tmp)) #, delta=.3, tuning.chi=2.560841)
                                 tmp2 <- unlist(tmp)
                                 if(any(is.na(tmp2))) {
                                   tmspe <- NA } else {
                                     # tmp2 <- tmp2[ !is.na(tmp2) ]
                                     # if(length(tmp2) > 0) {
                                     tmspe <- sqrt(mean(tmp2^2)) #RobStatTM::mscale(tmp2) #, delta=.3, tuning.chi=2.560841)
                                   }
                                 # } else {
                                 #   tmspe <- NA
                                 # }
                               }
  if(exists('old.seed')) assign('.Random.seed', old.seed, envir=.GlobalEnv)
  return(list(tmspe=tmp.par, h=hs, ncov=ncov, k=k, s=s, rho=reg.rho))
}

# Elliptical-FPCA main function
#
#
#' @export
lsfpca <- function(X, ncpus=4, opt.h.mu, opt.h.cov, hs.mu=seq(10, 25, by=1), hs.cov=hs.mu,
                   rho.param=1e-5, k = 3, s = k, trace=FALSE, seed=123, k.cv=5, ncov=50,
                   max.kappa=1e3) {
  
  # X is a list of 2 named elements "x", "pp"
  # which contain the lists of observations and times for each item
  # X$x[[i]] and X$pp[[i]] are the data for the i-th individual
  
  # Start cluster
  if( missing(opt.h.mu) || missing(opt.h.cov) ) {
    cl <- makeCluster(ncpus) # stopCluster(cl)
    registerDoParallel(cl)
  }
  # run CV to find smoothing parameters for
  # mean and covariance function
  if(missing(opt.h.mu)) {
    aa <- cv.mu.par.ls(X, hs=hs.mu, seed=seed, k.cv=k.cv)
    opt.h.mu <- aa$h[ which.min(aa$tmse) ]
  }
  mh <- mu.hat3.lin.ls(X=X, h=opt.h.mu)
  if(missing(opt.h.cov)) {
    bb <- cov.fun.cv.par.ls(X=X, muh=mh, ncov=ncov, k.cv=k.cv, hs=hs.cov, seed=seed,
                            k=k, s=s, reg.rho=rho.param)[1:2]
    opt.h.cov <- bb$h[ which.min(bb$tmspe) ]
  }
  if( exists('cl', inherits=FALSE) ) {
    stopCluster(cl)
  }
  ma <- matrixx(X, mh)
  # Compute the estimated cov function
  cov.fun2 <- cov.fun.hat2.ls(X=X, h=opt.h.cov, mh=mh, ma=ma, ncov=ncov, trace=FALSE)
  # smooth it
  yy <- as.vector(cov.fun2$G)
  xx <- cov.fun2$grid
  tmp <- fitted(mgcv::gam(yy ~ s(xx[,1], xx[,2]), family='gaussian'))
  cov.fun2$G <- matrix(tmp, length(unique(xx[,1])), length(unique(xx[,1])))
  cov.fun2$G <- ( cov.fun2$G + t(cov.fun2$G) ) / 2
  ours <- list(mh=mh, ma=ma, cov.fun2 = cov.fun2)
  
  # predicted scores, fitted values
  Xpred.fixed <- Xpred <- pred(X=X, muh=ours$mh, cov.fun=ours$cov.fun2$G,
                               tt=unique(ours$cov.fun2$grid[,1]),
                               ss=unique(ours$cov.fun2$grid[,1]), k=k, s=s, rho=rho.param)
  
  
  la1 <- eigen(ours$cov.fun2$G)$values[1]
  
  rho.param <- la1/(max.kappa-1)
  Xpred <- pred(X=X, muh=ours$mh, cov.fun=ours$cov.fun2$G,
                tt=unique(ours$cov.fun2$grid[,1]),
                ss=unique(ours$cov.fun2$grid[,1]), k=k, s=s, rho=rho.param)
  return(list(cov.fun = ours$cov.fun$G, muh=ours$mh, tt=unique(ours$cov.fun2$grid[,1]),
              ss=unique(ours$cov.fun2$grid[,2]), ma=ma, xis=Xpred$xis, pred=Xpred$pred,
              xis.fixed = Xpred.fixed$xis, pred.fixed = Xpred.fixed$pred,
              opt.h.mu=opt.h.mu, opt.h.cov=opt.h.cov, rho.param=rho.param))
}
