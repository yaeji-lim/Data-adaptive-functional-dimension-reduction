
#devtools::install_github('msalibian/sparseFPCA', ref = "master")
#devtools::install_github("statKim/robfpca")

library(sparseFPCA)
library(evd)
library(fGarch)
library(fMultivar)
library(elasticnet)
library(adegenet)
library(mvtnorm)
library(fdapace)
library(rospca)
library(mvnfast) 
library(ald)
library(truncnorm)
library(skewt)
library(fields)
library(spam)
library(elasticnet)
library(rrcov)
library(maps)
library(glmnet)
library(waveslim)
library(ggplot2)
library(LatticeKrig) 
library(pracma)
library(fda)
library(far)
library(softImpute)
library(rrcovNA)
library(ExtDist)
library(glmnet)
library(irlba)
library(misty)
library(robfpca)

source("simData_gen.R") # data generation
source("comparison_methods.R") # comparison methods
source("help_functions.R")  # necessary files to run FDR and simulation
source("FDR.R") # proposed method

#################################################
# simulation settings described in Section 3.1
#################################################
model.id=1   # model = 1 or 2
setting=1    # setting = 1, 2 or 3
m.prob=0.1   # m.prob = 0, 0.1, or 0.2
K.true=5     # K.true = 5 or 8
n=100        # sample size
t=seq(0.01,1,by = 0.01 )  # grids
sim.itr=c(1:100) 
n.sim= length(sim.itr)

rec_error_full <-eigen_error_angle <- eigen_error_angle_full <- matrix(nrow=n.sim, ncol=5)

for(sim in sim.itr){
  
  
  print(paste(sim,'th iteration'))
  
  set.seed(343*(sim*1)+setting+5*K.true)
  output =data_generation(t=t, n=n, K.true= K.true , setting=setting , sigma.e=1, m.prob=m.prob,symm=0.8, model=model.id,
                          sigma.con=5)
  ########################
  ### Data visualization
  ########################
  par(mfrow=c(1,3))
  sam.id <- sample(1:n, size=3)
  for (i in sam.id){
    plot(output$par.dat[i,], type='l', main='', xlab='t', ylab='')
    lines(output$true.dat[i,],col=2)
  }
  
  true.dat =output$true.dat ;sim.dat= output$sim.dat;  par.dat=output$par.dat; par.dat0= output$par.dat0  ; phi.mat =output$phi.mat  
  phi.mat=output$phi.mat	# true PC functions
  delta =output$delta     # missing indicators (0 for missing; 1 for observed)
  
  ###############################################################
  # 1. Proposed
  ###################################################################
  
  p=51  # number of basis functions/ fourier requires the odd number/ any large number can be used
  fr_obj = create.fourier.basis(rangeval = c(0,1), nbasis = (p)) 
  fbasisevals = eval.basis(t , fr_obj) 
  fr.gs = orthonormalization(fbasisevals,basis=FALSE ,norm=TRUE) # orthonormalize the basis
  
  x.fr = matrix(nrow = n, ncol = p)
  for (i in 1:n){
    for (j in 1:p){
      x.fr[i, j] = sum(par.dat0[i, ]*fr.gs[,j]) 
    }
  }
  
  x=x.fr; gs=fr.gs
  r.prop = 15# when K.true=5 we set r.prop=15 
             # when K.true=8 we set r.porp=20
  
  # sparsity.thres can be set 0.5 to 0.8
  # lambda.seq can be properly selected achieving good convergence or reasonable sparsity rate
  # may take some time to converge..
  fit_FDR = proposed_FDR(x=x, r= r.prop, cv=FALSE, adhoc.cv=TRUE, op.lambda = 0.2, thresh = 5e-05, 
                          lambda.seq = seq(0.35, 0.75, length=10), maxiter=200,maxiter.cv=200, K.cv=5, 
                          sparsity.thres=0.85)
  
  
  # estimated PC functions
  est.pc<- matrix(0,nrow=length(t), ncol= r.prop)
  for(ii in 1: r.prop){
    pcf <- rep(0, length.out=length(t))
    for(j in 1:p){
      pcf<- pcf+ fit_FDR$V[j,ii]*gs[,j]
    }
    
    est.pc[,ii]=pcf
  }
  
  # reconstruction using chosen PC functions
  n.pc=K.true
  closs.pca_rec=(fit_FDR$U[,1:n.pc]%*% diag(fit_FDR$Dsq[1:n.pc])) %*%t(est.pc[,1:n.pc])
  
  
  ################################################################################
  ## 2.  non-robust FDR: the non-robust version of our method using squared errors
  ###############################################################################
  
  fit_FDR_L2 = proposed_FDR_L2(x=x, r= r.prop, cv=FALSE, adhoc.cv=TRUE, op.lambda = 0.2, thresh = 5e-05, 
                             lambda.seq = seq(0.5,1, length=10), maxiter=200,maxiter.cv=200, K.cv=5, 
                             sparsity.thres=0.85) 
  
  
  # estimated PC functions
  est.pc_L2<- matrix(0,nrow=length(t), ncol= r.prop)
  for(ii in 1: r.prop){
    pcf <- rep(0, length.out=length(t))
    for(j in 1:p){
      pcf<- pcf+ fit_FDR_L2$V[j,ii]*gs[,j]
    }
    
    est.pc_L2[,ii]=pcf
  }
  
  # reconstruction using chosen PC functions
  n.pc=K.true
  closs.pca_rec_L2=(fit_FDR_L2$U[,1:n.pc]%*% diag(fit_FDR_L2$Dsq[1:n.pc])) %*%t(est.pc_L2[,1:n.pc])
  
  
  ######################################################################
  # 3. Robust FPCA by Boente and Salibian-Barrera (2021)
  #####################################################################
  r.pca= K.true
  x2= par.dat
  x.obs=na.indicator(x2)
  gt=t(matrix(t, nrow=length(t), ncol=n))
  gt[which(x.obs==0)]=NA
  if(length(which(x.obs==0))!=0){
    Bonte_data <- list(Ly = apply(x2, 1, function(y){ y[!is.na(y)] }),
                       Lt = apply(gt, 1, function(y){ y[!is.na(y)] }))}
  if(length(which(x.obs==0))==0){
    Bonte_data <- list( Ly = lapply(1:n, function(i) { x2[i, ] }), 
                        Lt = lapply(1:n, function(i) { gt[i, ] }) )  }
  
  bwmu_boente <- bwcov_boente <- 0.02;
  alpha=0.2
  cov.boente.obj <- cov_boente(Bonte_data,t, alpha=alpha, bw.mu = bwmu_boente, bw.cov = bwcov_boente,
                               seed = seed)
  mu.boente <- cov.boente.obj$mu
  cov.boente <- cov.boente.obj$cov
  boente.noise.est <- cov.boente.obj$noise_var
  pca.boente.obj <- funPCA(Bonte_data$Lt, Bonte_data$Ly, 
                           mu.boente, cov.boente, sig2 = boente.noise.est, 
                           t,  K = r.pca)
  boente.eigen=pca.boente.obj$eig.fun
  Boente_rec=predict(pca.boente.obj, K = r.pca)
  
  
  #######################################################
  ## 4.  Proj-Huber extended from Raymond and Lee (2017)
  ######################################################
  
  x2= x.fr
  RL_res = R_Lee(x2, r=r.prop, lambda= 20  )
  RL_pca_temp = RL_res$V
  
  RL_pca2 <- matrix(0,nrow=length(t), ncol= r.prop)
  for(ii in 1: r.pca){
    pcf <- rep(0, length.out=length(t))
    for(j in 1:p){
      pcf<- pcf+ RL_pca_temp[j,ii]*gs[,j]
    }
    
    RL_pca2[,ii]=pcf
  }
  
  
  RL_res_rec2=( RL_res$U[,1:n.pc]%*% diag( RL_res$Dsq[1:n.pc])) %*%t( RL_pca2[,1:n.pc])
  
  
  ####################################
  ## 5.  Partial FPCA by Kraus (2015)
  ###################################
  x2= par.dat   
  mu.kraus <- mean.missfd(x2)
  cov.kraus <- var.missfd(x2)
  eig.kraus	<- eigen.missfd(cov.kraus)$vectors
  K_kraus <- r.pca
  kraus.pca <- eig.kraus[, 1:K_kraus]
  
  
  cand <- which(
    (apply(x2, 1, function(x){ sum(is.na(x)) }) > 0) 
  )
  
  
  if(length(cand)!=0){ 
    rec_Kraus=matrix(nrow=nrow(x2), ncol=ncol(x2))
    for (i in 1:length(cand)) {
      ind <- cand[i]
      pred_comp <-  pred.missfd(x2[ind, ], x2)
      NA_ind <- which(is.na(x2[ind, ]))   # index of missing periods
      rec_Kraus[ind,NA_ind ] <-  pred_comp[NA_ind]
    }
  }else{
    
  }
  
  
  Kraus.score=matrix(nrow=nrow(x2), ncol=r.pca)
  for(i in 1:nrow(x2)){	
    tryCatch({
      Kraus.score[i,]=pred.score.missfd(x2[i,], kraus.pca,x2 ,gcv.df.factor=2)
    }, error=function(e){})
  }
  
  
  
  rec_Kraus_full= mu.kraus + Kraus.score %*% t(kraus.pca )
  
  
  #####################################
  # PC function estimation performance
  ######################################
  
  # normalize estimated functions
  est.pc =norm_eigenfunction(est.pc, t)
  est.pc_L2 =norm_eigenfunction(est.pc_L2, t)
  boente.eigen =norm_eigenfunction(boente.eigen, t)
  RL_pca2 =norm_eigenfunction(RL_pca2, t)
  kraus.pca=norm_eigenfunction(kraus.pca, t)
  
  # match the sign
  for(l in c(1: K.true)){    
    est.pc[,l]=sign_eigenfunction((est.pc[,l]),(phi.mat[,l])   )
    est.pc_L2[,l]=sign_eigenfunction((est.pc_L2[,l]),(phi.mat[,l])   )
    boente.eigen[,l] = sign_eigenfunction((boente.eigen[,l]),(phi.mat[,l])   )
    RL_pca2[,l] = sign_eigenfunction((RL_pca2[,l]),(phi.mat[,l])   )
    kraus.pca[,l] = sign_eigenfunction((kraus.pca[,l]),(phi.mat[,l])   )
  } 
  
  ###################################
  # Visualize estimated PC functions 
  ###################################

  par(mfrow=c(2,2))
  for(l in c(1:4)){
    plot(t,(phi.mat[,l]), type="l", ylab="", xlab="", ylim=c(-2,2))
    lines(t,(est.pc[,l]),col="blue", lwd = 2)
    lines(t,(est.pc_L2[,l]),col="red", lwd = 2)
    lines(t,(boente.eigen[,l]),col="purple", lwd = 2)
    lines(t,(RL_pca2[,l]),col="green", lwd = 2)
    lines(t,(kraus.pca[,l]),col="pink", lwd = 2)
    legend("topright", col=c("black","blue","red","purple","green",'pink'), lwd=2, legend=c("true",'proposed FDR','non-robust FDR','Robust FPCA','Proj-Huber','Kraus'),cex=1)
  }
  
  
  ####################################
  #  Visualize reconstruction performance 
  ####################################
  
  par(mfrow=c(2,2))
  tmp.id = sample(1:n, 4)
  for(l in tmp.id){
    plot(t,sim.dat[l,], type="p", ylab="", xlab="")
    lines(t,true.dat[l,],col="black", lwd=2)
    lines(t, closs.pca_rec[l,],col="blue", lwd = 2)
    lines(t, closs.pca_rec_L2[l,],col="red", lwd = 2)
    lines(t, Boente_rec[l,],col="purple", lwd = 2)
    lines(t, RL_res_rec2[l,],col="green", lwd = 2)
    lines(t, rec_Kraus_full[l,],col="pink", lwd = 2)
    
  }  
  
  ###################################################
  # calculate recovery errors (reconstruction MISE)
  ###################################################
  rec_error_full[which(sim.itr==sim),]=
    c(
      mean((true.dat - closs.pca_rec)^2) ,# MSE between true (before contamination) and reconstructed curves
      mean((true.dat - closs.pca_rec_L2)^2),
      mean((true.dat - Boente_rec)^2),
      mean((true.dat - RL_res_rec2)^2),
      mean((true.dat - rec_Kraus_full)^2)
    )
  
  ###################################################
  # calculate PC estimation errors (eigen-angle 90%)
  ###################################################
  eva.pca=3 # first three explain 90% of total variations
  eigen_error_angle[which(sim.itr==sim), ] <-c(
    mean(
      sapply(1:eva.pca, function(i){
        subspace(est.pc[, i], phi.mat[, i])
      })
    ),
    mean(
      sapply(1:eva.pca, function(i){
        subspace(est.pc_L2[, i], phi.mat[, i])
      })
    ),
    mean(
      sapply(1:eva.pca, function(i){
        subspace(boente.eigen[, i], phi.mat[, i])
      })
    ),
    mean(
      sapply(1:eva.pca, function(i){
        subspace(RL_pca2[, i], phi.mat[, i])
      })
    ),
    mean(
      sapply(1:eva.pca, function(i){
        subspace(kraus.pca[, i], phi.mat[, i])
      })
    )
  )  

  ###################################################
  # calculate PC estimation errors (eigen-angle for all)
  ###################################################
  
  eigen_error_angle_full[which(sim.itr==sim), ] <-c(
    mean(
      sapply(1:K.true, function(i){
        subspace(est.pc[, i], phi.mat[, i])
      })
    ),
    mean(
      sapply(1:K.true, function(i){
        subspace(est.pc_L2[, i], phi.mat[, i])
      })
    ),
    mean(
      sapply(1:K.true, function(i){
        subspace(boente.eigen[, i], phi.mat[, i])
      })
    ),
    mean(
      sapply(1:K.true, function(i){
        subspace(RL_pca2[, i], phi.mat[, i])
      })
    ),
    mean(
      sapply(1:K.true, function(i){
        subspace(kraus.pca[, i], phi.mat[, i])
      })
    )
  )  
  
  print(which(sim.itr==sim))
  
}


t1=cbind( round( apply( eigen_error_angle  ,2, mean, na.rm=T) ,3),
          round( apply( eigen_error_angle_full  ,2, mean, na.rm=T) ,3),
          round( apply( rec_error_full  ,2, mean, na.rm=T) ,3)
          )

t2=cbind( round( apply( eigen_error_angle  ,2, sd, na.rm=T) ,3),
          round( apply( eigen_error_angle_full  ,2, sd, na.rm=T) ,3),
          round( apply( rec_error_full  ,2, sd, na.rm=T) ,3)
          )


temp2=NA
for(i in 1:nrow(t1)){
  temp1=NA
  for(j in 1:ncol(t1)){
    temp1=c( temp1 , paste(t1[i,j], '(', t2[i,j] , ') & ' , sep='') )
  }
  
  temp2= rbind(temp2, temp1)
  
}

result=temp2[-1,-1]

rownames(result)=c('Proposed FDR', 'non-robust FDR','Robust FPCA','Huber-Proj','Partial FPCA')
colnames(result)=c('Proposed FDR', 'non-robust FDR','Robust FPCA','Huber-Proj','Partial FPCA')



