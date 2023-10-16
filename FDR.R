

##################################################
# Proposed method
###################################################

# need source files "others_function.R"


proposed_FDR<-function( x, r=15,  lambda.seq = seq(0,2, by = 0.5), maxiter = 200, maxiter.cv=200,
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


######################################################################
# non-robust version of proposed method using squared loss (L2) errors
#####################################################################3


proposed_FDR_L2<-function( x, r=15,  lambda.seq = seq(0,2, by = 0.5), maxiter = 200, maxiter.cv=200,
                          cv=FALSE, adhoc.cv=FALSE,  op.lambda, thresh = 5e-05, K.cv=5,sparsity.thres=0.8){
  
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
            Z = Y.0 + E 
            tmp.als = svd.als.yj.L1(Z,rank.max =r,lambda=(lambda-0.1), final.svd = FALSE, maxit = maxiter, alpha=alpha)
            
            U = tmp.als$u
            V = tmp.als$v
            Dsq = tmp.als$d
            
            Y.1 = U%*%diag(Dsq)%*%t(V)
          }
          
        }
        
        iter=iter+1
        #print(iter)
        tryCatch({
          
          U.old=U
          V.old=V
          Dsq.old=Dsq
          
          Y.0=Y.1
          E = x.train - Y.0 
          Z = Y.0 + E 
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
        
        if(iter %% 100 == 0) (c(iter, ratio))  
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



