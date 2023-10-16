######################################################
# Huber loss functions
######################################################
psi<-function(xx,tau=.5){
  cut= 1.345 
  result =matrix(nrow=nrow(xx), ncol=ncol(xx))
  for(i in 1:nrow(xx)){
    for(j in 1:ncol(xx)){
      x<-xx[i,j]
      if(x<(-cut)){u=(tau-1)/2}
      if(x>=(-cut) & x<0 ){ u=(1-tau)*x}
      if(x>=0 & x<cut){u=tau*x}
      if(x>=(cut)) {u=tau/2}
      result[i,j]<-u
    }}
  return(result)
}   

######################################################
# Composite loss functions (Lim and Oh, 2016)
######################################################
com_psi<-function(x, opt_w=NULL,b, tau_range, plot=FALSE){
  
  
  if(is.null(opt_w)){
    
    
    opt_w<-matrix(nrow=length(tau_range), ncol=nrow(x))
    for(i in 1: nrow(x)){
      a<-vector()
      for(k in 1:length(tau_range)){	a[k]<- dens(quantile(x[i,],tau_range[k]), x[i,])  }
      opt_w[,i] =(a)
      if(max(a)>1){ opt_w[,i]=a/sum(a)
      }
    }
    
    if(plot==TRUE){
      plot(tau_range ,apply(opt_w,1, mean) ) 
    }
  }
  
  A <- matrix(0, nrow=nrow(x), ncol=ncol(x))
  for(j in 1:nrow(x)){
    for(i in 1:length(tau_range)){				
      z=opt_w[i,j]*psi(as.matrix(x[j,]-b[i,j]),tau_range[i])
      A[j,]= A[j,] + z 
    }
  }
  return(A)
}

dens<-function(x,y){
  index=density(y)$x 
  fit=density(y)$y[which.min(abs(index-x))]
  return(fit)
}

norm_eigenfunction<-function(est.pc,t ){
  
  est.pc <- apply(est.pc, 2, function(x) {
    x <- x/sqrt(trapzRcpp(t, x^2))
    return(x)
  })
  
  return(est.pc)
  
}

sign_eigenfunction<-function(x, true){
  
  
  rmse1=  mean( ( true - x )^2  )
  rmse2= mean( ( true +x )^2  )
  
  if(rmse2 <rmse1){ x = - x }
  
  return(x)
  
}


######################################################
# help functions for the proposed method
######################################################

####################################################
# code related to the softImpute-ALS algorithm (Hastie et al., 2015)
clean.warm.start=function(a){
  if(is.null(a))return(NULL)
  d=a$d
  if(is.null(d))return(NULL)
  if(any(d>0)){
    if(length(d)==1){
      a$u=matrix(a$u,ncol=1)
      a$v=matrix(a$v,ncol=1)
    }
    a
  }
  else NULL
}

###################
softImpute.x.matrix=function(x, J,lambda,type, thresh,maxit,trace.it,warm.start,final.svd){
  switch(type,
         "als"=simpute.als(x, J,thresh, lambda,maxit,trace.it,warm.start,final.svd),
         "svd"=simpute.svd(x, J,thresh, lambda,maxit,trace.it,warm.start,final.svd)
  )
}

softImpute.x.Incomplete=function(x, J,lambda,type, thresh,maxit,trace.it,warm.start,final.svd){
  switch(type,
         "als"=Ssimpute.als(x, J,thresh, lambda,maxit,trace.it,warm.start,final.svd),
         "svd"=Ssimpute.svd(x, J,thresh, lambda,maxit,trace.it,warm.start,final.svd)
  )
}

setGeneric("softImpute.x",softImpute.x.matrix)

#####################
Frob=function(Uold,Dsqold,Vold,U,Dsq,V){
  denom=sum(Dsqold^2)
  utu=Dsq* (t(U)%*%Uold)
  vtv=Dsqold* (t(Vold)%*%V)
  uvprod= sum(diag(utu%*%vtv))
  num=denom+sum(Dsq^2) -2*uvprod
  num/max(denom,1e-9)
}


######################
simpute.als<-function (x, J = 2, thresh = 1e-05,lambda=0,maxit=100,trace.it=TRUE, warm.start=NULL, final.svd=TRUE) 
{### x is a matrix, possibly with NAs
  n <- dim(x)
  m <- n[2]
  n <- n[1]
  this.call=match.call()
  a=names(attributes(x))
  binames=c("biScale:row","biScale:column")
  if(all(match(binames,a,FALSE))){
    biats=attributes(x)[binames]
  } else biats=NULL
  
  xnas <- is.na(x)
  nz=m*n-sum(xnas)
  xfill <- x
  if(!is.null(warm.start)){
    #must have u,d and v components
    if(!all(match(c("u","d","v"),names(warm.start),0)>0))stop("warm.start does not have components u, d and v")
    warm=TRUE
    D=warm.start$d
    JD=sum(D>0)
    if(JD >= J){
      U=warm.start$u[,seq(J),drop=FALSE]
      V=warm.start$v[,seq(J),drop=FALSE]
      Dsq=D[seq(J)]
    }
    else{
      Dsq=c(D,rep(D[JD],J-JD))
      Ja=J-JD
      U=warm.start$u
      Ua=matrix(rnorm(n*Ja),n,Ja)
      Ua=Ua-U%*% (t(U)%*%Ua)
      Ua=svd(Ua)$u
      U=cbind(U,Ua)
      V=cbind(warm.start$v,matrix(0,m,Ja))
    }
    xfill[xnas]=(U%*%(Dsq*t(V)))[xnas]  
  }
  else
  {
    V=matrix(0,m,J)
    U=matrix(rnorm(n*J),n,J)
    U=svd(U)$u
    Dsq=rep(1,J)# we call it Dsq because A=UD and B=VD and AB'=U Dsq V^T
    xfill[xnas]=0
  }
  ratio <- 1
  iter <- 0
  while ((ratio > thresh)&(iter<maxit)) {
    iter <- iter + 1
    U.old=U
    V.old=V
    Dsq.old=Dsq
    
    ## U step
    B=t(U)%*%xfill
    if(lambda>0)B=B*(Dsq/(Dsq+lambda))
    Bsvd=svd(t(B))
    V=Bsvd$u
    Dsq=(Bsvd$d)
    U=U%*%Bsvd$v
    xhat=U %*%(Dsq*t(V))
    xfill[xnas]=xhat[xnas]
    ###The next line we could have done later; this is to match with sparse version
    if(trace.it) obj=(.5*sum( (xfill-xhat)[!xnas]^2)+lambda*sum(Dsq))/nz
    ## V step
    A=t(xfill%*%V)
    if(lambda>0)A=A*(Dsq/(Dsq+lambda))
    Asvd=svd(t(A))
    U=Asvd$u
    Dsq=Asvd$d
    V=V %*% Asvd$v
    xhat=U %*%(Dsq*t(V))
    xfill[xnas]=xhat[xnas]
    ratio=Frob(U.old,Dsq.old,V.old,U,Dsq,V)
    if(trace.it) cat(iter, ":", "obj",format(round(obj,5)),"ratio", ratio, "\n")
  }
  if(iter==maxit)warning(paste("Convergence not achieved by",maxit,"iterations"))
  
  if(lambda>0&final.svd){
    U=xfill%*%V
    sU=svd(U)
    U=sU$u
    Dsq=sU$d
    V=V%*%sU$v
    Dsq=pmax(Dsq-lambda,0)
    if(trace.it){
      xhat=U %*%(Dsq*t(V))
      obj=(.5*sum( (xfill-xhat)[!xnas]^2)+lambda*sum(Dsq))/nz
      cat("final SVD:", "obj",format(round(obj,5)),"\n")
    }
    
  }
  J=min(sum(Dsq>0)+1,J)
  out=list(u=U[,seq(J)],d=Dsq[seq(J)],v=V[,seq(J)])
  attributes(out)=c(attributes(out),list(lambda=lambda,call=this.call),biats)
  out
  
}

######################
simpute.als.yj.L1<-function (x, J = 2, thresh = 1e-05,lambda=0,maxit=100,trace.it=TRUE, warm.start=NULL, 
                             final.svd=TRUE, alpha=0.5) {
  n <- dim(x)
  m <- n[2]
  n <- n[1]
  this.call=match.call()
  a=names(attributes(x))
  binames=c("biScale:row","biScale:column")
  if(all(match(binames,a,FALSE))){
    biats=attributes(x)[binames]
  } else biats=NULL
  
  xnas <- is.na(x)
  nz=m*n-sum(xnas)
  xfill <- x
  if(!is.null(warm.start)){
    #must have u,d and v components
    if(!all(match(c("u","d","v"),names(warm.start),0)>0))stop("warm.start does not have components u, d and v")
    warm=TRUE
    D=warm.start$d
    JD=sum(D>0)
    if(JD >= J){
      U=warm.start$u[,seq(J),drop=FALSE]
      V=warm.start$v[,seq(J),drop=FALSE]
      Dsq=D[seq(J)]
    }else{
      Dsq=c(D,rep(D[JD],J-JD))
      Ja=J-JD
      U=warm.start$u
      Ua=matrix(rnorm(n*Ja),n,Ja)
      Ua=Ua-U%*% (t(U)%*%Ua)
      Ua=svd(Ua)$u
      U=cbind(U,Ua)
      V=cbind(warm.start$v,matrix(0,m,Ja))
    }
    xfill[xnas]=(U%*%(Dsq*t(V)))[xnas]  
  }else{
    V=matrix(0,m,J)
    U=matrix(rnorm(n*J),n,J)
    U=svd(U)$u
    Dsq=rep(1,J)# we call it Dsq because A=UD and B=VD and AB'=U Dsq V^T
    xfill[xnas]=0
  }
  ratio <- 1
  iter <- 0
  while ((ratio > thresh)&(iter<maxit)) {
    iter <- iter + 1
    U.old=U
    V.old=V
    Dsq.old=Dsq
    
    
    X.tmp=UD(U,sqrt(Dsq),n)
    
    B.tmp = c()
    for(j in 1:ncol(x)){
      B.tmp = cbind(B.tmp, as.matrix(coef(glmnet(X.tmp, xfill[,j], lambda=lambda/sd(xfill[,j]), family="gaussian", intercept = F, alpha=alpha, 
                                                 thresh = 1e-06, standardize = FALSE))))
    }
    
    B = diag(sqrt(Dsq),J,J)%*%B.tmp[-1,]
    
    Bsvd=svd(t(B))
    V=Bsvd$u
    Dsq=(Bsvd$d)
    U=U%*%Bsvd$v
    xhat=U %*%(Dsq*t(V))
    xfill[xnas]=xhat[xnas]
    
    if(trace.it) obj=(.5*sum( (xfill-xhat)[!xnas]^2)+lambda*sum(Dsq))/nz
    
    A=t(xfill%*%V)
    if(lambda>0)A=A*(Dsq/(Dsq+lambda))
    
    Asvd=svd(t(A))
    U=Asvd$u
    Dsq=Asvd$d
    V=V %*% Asvd$v
    xhat=U %*%(Dsq*t(V))
    xfill[xnas]=xhat[xnas]
    ratio=Frob(U.old,Dsq.old,V.old,U,Dsq,V)
    if(trace.it) cat(iter, ":", "obj",format(round(obj,5)),"ratio", ratio, "\n")
  }
  if(iter==maxit)warning(paste("Convergence not achieved by",maxit,"iterations"))
  
  if(lambda>0&final.svd){
    U=xfill%*%V
    sU=svd(U)
    U=sU$u
    Dsq=sU$d
    V=V%*%sU$v
    Dsq=pmax(Dsq-lambda,0)
    if(trace.it){
      xhat=U %*%(Dsq*t(V))
      obj=(.5*sum( (xfill-xhat)[!xnas]^2)+lambda*sum(Dsq))/nz
      cat("final SVD:", "obj",format(round(obj,5)),"\n")
    }
    
  }
  out=list(u=U[,seq(J)],d=Dsq[seq(J)],v=V[,seq(J)])
  attributes(out)=c(attributes(out),list(lambda=lambda,call=this.call),biats)
  out
  
}

######################
svd.als.yj.L1=function(x, rank.max=2, lambda=0,thresh = 1e-05, maxit=100,trace.it=FALSE,warm.start=NULL,final.svd=TRUE, alpha=0.5){
  if(rank.max>(rmax<-min(dim(x)))){
    rank.max=rmax
    warning(paste("rank.max should not exceed min(dim(x)); changed to ",rmax))
  }
  ismiss=is.na(x)
  if(any(ismiss))stop("NAs in x; use softImpute instead")
  this.call=match.call()
  out=simpute.als.yj.L1(x,J=rank.max,thresh,lambda,maxit,trace.it,warm.start,final.svd)
  attr(out,"call")=this.call
  attr(out,"lambda")=lambda
  out
}


setGeneric("svd.als",svd.als)

######################################

UD=function(U,D,n=nrow(U)){
  U*outer(rep(1,n),D,"*")
}

##################################
suvC <-
  function(u,v,irow,pcol){
    dd=dim(u)
    nnrow=as.integer(dd[1])
    nncol=as.integer(nrow(v))
    nrank=dd[2]
    storage.mode(u)="double"
    storage.mode(v)="double"
    storage.mode(irow)="integer"
    storage.mode(pcol)="integer"
    nomega=as.integer(length(irow))
    .Fortran("suvC",
             nnrow,nncol,nrank,u,v,irow,pcol,nomega,
             r=double(nomega),
             PACKAGE="softImpute"
    )$r
  }

##################################
suv <-
  function(u,v,irow,jcol){
    ###  computes (u%*%t(v))[i,j]= sum(u[i,]*v[j,]) for all pairs in irow,jcol
    dd=dim(u)
    nnrow=as.integer(dd[1])
    nncol=as.integer(nrow(v))
    nrank=dd[2]
    storage.mode(u)="double"
    storage.mode(v)="double"
    storage.mode(irow)="integer"
    storage.mode(jcol)="integer"
    nomega=as.integer(length(jcol))
    .Fortran("suv",
             nnrow,nncol,nrank,u,v,irow,jcol,nomega,
             r=double(nomega),
             PACKAGE="softImpute"
    )$r
  }


##################################
dens<-function(x,y){
  index=density(y)$x 
  fit=density(y)$y[which.min(abs(index-x))]
  return(fit)
}

##################################
norm_eigenfunction<-function(est.pc,t ){
  
  est.pc <- apply(est.pc, 2, function(x) {
    x <- x/sqrt(trapzRcpp(t, x^2))
    return(x)
  })
  
  return(est.pc)
  
}

##################################
sign_eigenfunction<-function(x, true){
  
  
  rmse1=  mean( ( true - x )^2  )
  rmse2= mean( ( true +x )^2  )
  
  if(rmse2 <rmse1){ x = - x }
  
  return(x)
  
}
