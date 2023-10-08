##################################################
# Data generation 
#      
######################################################


data_generation <-function(t=seq(0.01,1,by = 0.01 ), n=100, r.true=4 , setting=1, sigma.u = 1, 
                           m.prob=0.2, dd=0.7, f=0.2, symm=0.2, model=1, sigma.con=6){
  
  
  if(model==1){
    phi.mat=c()
    for (i in 1:r.true){
      phi.mat = cbind(phi.mat, sqrt(2)*cos(2*i*pi*t))
    }

    eval.lambda = 5*(c(1:r.true))^{-2}

    
    true.dat=c()   # data without noise (smooth curves without measurement errors)

    for(i in 1:n){
      coef = rmvnorm(1,rep(0,length(eval.lambda)), diag(eval.lambda))
      coef.mat = matrix(rep(coef, each=length(t)), nrow=length(t), ncol=r.true, byrow = FALSE)
      true.dat = rbind(true.dat, rowSums(coef.mat * phi.mat)) # each row is each curve
    }
    
  } else if (model==2){
    
    dist.mat<-matrix(nrow=length(t),ncol=length(t))
    for(i in 1:length(t)){
      for(j in 1:length(t)){
        dist.mat[i,j] = abs(t[i]-t[j])
      }
    }
    
    matern.cor=matrix(Matern(as.vector(dist.mat), range=0.2, smoothness=1),nrow=nrow(dist.mat))
    tmp = svd(matern.cor)
    
    phi.mat = 10*tmp$v[,c(1:r.true)]
    eval.lambda = tmp$d[c(1:r.true)]/10
    
    
    true.dat=c()  
    
    for(i in 1:n){
      coef = rmvnorm(1,rep(0,length(eval.lambda)), diag(eval.lambda))
      coef.mat = matrix(rep(coef, each=length(t)), nrow=length(t), ncol=r.true, byrow = FALSE)
      true.dat = rbind(true.dat, rowSums(coef.mat * phi.mat)) # each row is each curve
    }
    
    
  }
  
  if(setting==1){
    # (setting 1) gaussian error
    sim.dat = true.dat +sigma.u* matrix(rnorm(length(true.dat),0,sigma.u),ncol=ncol(true.dat), nrow=nrow(true.dat))
  }
  #  (setting 2):contaminated gaussian error
  if(setting==2){
    con.id=sample(1:(length(true.dat)), size=length(true.dat)*0.3, replace=F)
    tmp.vec <- as.vector(true.dat)
    sim.vec <- tmp.vec + rnorm(length(tmp.vec),0,sigma.u)
    sim.vec[con.id] <- sim.vec[con.id] + rnorm(length(con.id),0,sigma.con)
    sim.dat = matrix(sim.vec, nrow=nrow(true.dat), ncol=ncol(true.dat))
  }
  #  (setting 3):contaminated skewed error
  if(setting==3){ 
    con.id1=sample(1:(length(true.dat)), size=length(true.dat)*0.15, replace=F)
    con.id2=sample(1:(length(true.dat)), size=length(true.dat)*0.15, replace=F)
    tmp.vec <- as.vector(true.dat)
    sim.vec <- tmp.vec + rnorm(length(tmp.vec),0,sigma.u)
    sim.vec[con.id1] <- sim.vec[con.id1] + rnorm(length(con.id1),0,sigma.con)
    sim.vec[con.id2] <- sim.vec[con.id2] + rALD(length(con.id2),mu=0,sigma=1,p= symm)
    sim.dat = matrix(sim.vec, nrow=nrow(true.dat), ncol=ncol(true.dat))
  }
  
  

  ### Make the data partially observed 
  delta=matrix(1,nr=nrow(sim.dat), nc=ncol(sim.dat))
  
  if(m.prob > 0){  

    for(jj in 1:nrow(delta)){

      if(rbinom(n=1, size=1, prob=m.prob)==1){
        
        u1=sqrt(runif(1)); u2=runif(1)
        Ci=dd*u1; Ei = f*u2
        if((Ci-Ei)<0.01){
          l.b=1
        }else if((Ci-Ei)>0.99){
          l.b=length(t)
        } else {l.b=which(round(t,2)==round(Ci-Ei,2)) }
        
        if((Ci+Ei)>= 0.99){
          u.b=length(t)
        }else if((Ci+Ei)<0.01){
          u.b=1
        } else {u.b=which(round(t,2)==round(Ci+Ei,2)) }
        
        if(l.b==u.b & l.b==99){
          delta[jj,] = 1
        }else if(l.b==u.b & u.b==1){
          delta[jj,] = 1
        }else delta[jj,c(l.b:u.b)] = 0

        
      }
    }
  } 
  par.dat0=sim.dat*delta  # par.dat0: missing is replaced by 0
  par.dat = par.dat0
  par.dat[par.dat==0]=NA       # par.dat: missing is prelaced by NA
  
  
  output<-list()
  output$true.dat= true.dat
  output$sim.dat= sim.dat;  output$par.dat=par.dat;  output$par.dat0 = par.dat0 ; 
  output$delta=delta
  output$phi.mat = phi.mat
  return(output)
}	


##################################################
# Proposed method
#      
######################################################



proposed_yj_V2<-function( x, r=15,  lambda.seq = seq(0,2, by = 0.5), maxiter = 200, maxiter.cv=200,
                       cv=FALSE, adhoc.cv=FALSE,  op.lambda, thresh = 5e-05, K.cv=5,sparsity.thres=0.8){
  
  
  tau_range<-seq(0.05, 0.95, by=0.05)
  b_tau =matrix(nrow=length(tau_range), ncol=(nrow(x)),0)
  
 
  alpha=0.5  # for elastic net
  
  if(adhoc.cv==TRUE){	
    lambda.sparsity <- 0; j<-0
    
    tmp.V.list=list(); tmp.U.list=list(); tmp.Dsq.list=list()

    while (lambda.sparsity < sparsity.thres & j < length(lambda.seq)){
     
      j<- j+1
      print(paste('lambda',j))  
      lambda = lambda.seq[j]
      x.train <- x
      
      r.closs=r
      iter=0; ratio=1; cond=TRUE
      
      while(cond){
        
        if(iter==0){
          nrob = softImpute(x.train,rank.max=r.closs,lambda=lambda)
          U = nrob$u
          V = nrob$v
          Dsq = nrob$d
          Y.1 = U%*%diag(Dsq)%*%t(V)
          
          iiter=1
          if(iiter<10){
            iiter=iiter+1
            U.old=U
            V.old=V
            Dsq.old=Dsq
            
            Y.0=Y.1
            E = x - Y.0 
            Z = Y.0 + com_psi( x=E, b= b_tau, tau_range =tau_range, plot=FALSE)/2  
            tmp.als = svd.als.yj.L1(Z,rank.max =r,lambda=(lambda-0.1), final.svd = FALSE, maxit = maxiter, alpha=alpha)
            
            U = tmp.als$u
            V = tmp.als$v
            Dsq = tmp.als$d
            
            Y.1 = U%*%diag(Dsq)%*%t(V)
          }
          
        }
        
        iter=iter+1

        tryCatch({
          
          U.old=U
          V.old=V
          Dsq.old=Dsq
          
          Y.0=Y.1
          E = x.train - Y.0 
          Z = Y.0 + com_psi( x=E, b= b_tau, tau_range =tau_range)/2 
          tmp.als = svd.als.yj.L1(Z,rank.max =r,lambda=lambda, final.svd = FALSE, maxit = maxiter, alpha=alpha)
          
          U = tmp.als$u
          V = tmp.als$v
          Dsq = tmp.als$d
          
          Y.1 = U%*%diag(Dsq)%*%t(V)
          
          ratio=Frob(U.old,Dsq.old,V.old,U,Dsq,V)
          cond = !(ratio < thresh| iter > maxiter.cv)
        
          
          
          
        }, error = function(e) { 
          print("error")
          print(e)
          skip_sim <<- TRUE
        })
   
      }
      
 
      lambda.sparsity  <- sum(V==0)/length(V)
      tmp.V.list[[j]] <- V; tmp.U.list[[j]] <- U; tmp.Dsq.list[[j]] <- Dsq
      
      print(lambda.sparsity)
    }
    

    op.lambda  = lambda

    
    print(paste('cv selected=', op.lambda) )
  }
  
  ###########################################

  
  output <-list()
  
  output$Y.1= tmp.U.list[[j]] %*% diag(tmp.Dsq.list[[j]])%*%t(tmp.V.list[[j]])
  output$U = tmp.U.list[[j]]
  output$V = tmp.V.list[[j]]
  output$Dsq = tmp.Dsq.list[[j]]
  output$op.lambda = op.lambda
  
  return(output)
  
  
}


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








dens<-function(x1,y){
  index=density(y)$x 
  fit=density(y)$y[which.min(abs(index-x1))]
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

