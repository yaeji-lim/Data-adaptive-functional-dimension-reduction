##################################################
# Data generation 
#      
######################################################


data_generation <-function(t=seq(0.01,1,by = 0.01 ), n=100, K.true=4 , setting=1, sigma.e = 1, 
                           m.prob=0.2, dd=0.7, f=0.2, symm=0.2, model=1, sigma.con=6){
  
  
  if(model==1){
    phi.mat=c()
    for (i in 1:K.true){
      phi.mat = cbind(phi.mat, sqrt(2)*cos(2*i*pi*t))
    }
    
    eval.lambda = 5*(c(1:K.true))^{-2}
    
    
    true.dat=c()   # data without noise (smooth curves without measurement errors)
    
    for(i in 1:n){
      coef = rmvnorm(1,rep(0,length(eval.lambda)), diag(eval.lambda))
      coef.mat = matrix(rep(coef, each=length(t)), nrow=length(t), ncol=K.true, byrow = FALSE)
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
    
    phi.mat = 10*tmp$v[,c(1:K.true)]
    eval.lambda = tmp$d[c(1:K.true)]/10
    
    
    true.dat=c()  
    
    for(i in 1:n){
      coef = rmvnorm(1,rep(0,length(eval.lambda)), diag(eval.lambda))
      coef.mat = matrix(rep(coef, each=length(t)), nrow=length(t), ncol=K.true, byrow = FALSE)
      true.dat = rbind(true.dat, rowSums(coef.mat * phi.mat)) # each row is each curve
    }
    
    
  }
  
  if(setting==1){
    # (setting 1): without contamination
    sim.dat = true.dat +sigma.e* matrix(rnorm(length(true.dat),0,sigma.e),ncol=ncol(true.dat), nrow=nrow(true.dat))
  }
  #  (setting 2): contaminated by additive gaussian measurement errors
  if(setting==2){
    con.id=sample(1:(length(true.dat)), size=length(true.dat)*0.3, replace=F)
    tmp.vec <- as.vector(true.dat)
    sim.vec <- tmp.vec + rnorm(length(tmp.vec),0,sigma.e)
    sim.vec[con.id] <- sim.vec[con.id] + rnorm(length(con.id),0,sigma.con)
    sim.dat = matrix(sim.vec, nrow=nrow(true.dat), ncol=ncol(true.dat))
  }
  #  (setting 3): contaminated by additive asymmetric measurement errors
  if(setting==3){ 
    con.id1=sample(1:(length(true.dat)), size=length(true.dat)*0.15, replace=F)
    con.id2=sample(1:(length(true.dat)), size=length(true.dat)*0.15, replace=F)
    tmp.vec <- as.vector(true.dat)
    sim.vec <- tmp.vec + rnorm(length(tmp.vec),0,sigma.e)
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
  
  par.dat0=sim.dat*delta  # par.dat0: missing represented by 0
  par.dat = par.dat0
  par.dat[par.dat==0]=NA  # par.dat: missing represented by NA
  
  
  output<-list()
  output$true.dat= true.dat
  output$sim.dat= sim.dat;  output$par.dat=par.dat;  output$par.dat0 = par.dat0 ; 
  output$delta=delta
  output$phi.mat = phi.mat
  return(output)
}	