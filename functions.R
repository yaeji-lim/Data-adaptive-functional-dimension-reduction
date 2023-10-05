
data_generation <-function(t=seq(0.01,1,by = 0.01 ), n=100, r.true=4 , setting=1, sigma.u = 1, 
                           m.prob=0.2, dd=0.7, f=0.2, symm=0.2, model=1, sigma.con=6){
  
  
  ##################################################
  # Yeonjoo: true rank 5, 10 정도 두가지 경우?
  #       true rank 크던 작던 잘 작동하는걸 보여주고싶음
  ######################################################
  
  if(model==1){
    phi.mat=c()
    for (i in 1:r.true){
      phi.mat = cbind(phi.mat, sqrt(2)*cos(2*i*pi*t))
    }
    # 082222 - YJ added: make the scale 1
    
    
    eval.lambda = 5*(c(1:r.true))^{-2}
    #eval.lambda = 10*(c(1:r.true))^{-2}
    #cumsum(eval.lambda)/sum(eval.lambda)
    
    
    true.dat=c()   # data without noise (smooth curves without measurement errors)
    # variance for measurement error
    
    
    for(i in 1:n){
      #set.seed(23+89*i)
      coef = rmvnorm(1,rep(0,length(eval.lambda)), diag(eval.lambda))
      coef.mat = matrix(rep(coef, each=length(t)), nrow=length(t), ncol=r.true, byrow = FALSE)
      true.dat = rbind(true.dat, rowSums(coef.mat * phi.mat)) # each row is each curve
    }
    
  } else if (model==2){
    m = length (t)
    
    phi.mat = get_delaigle_eigen(t, model = 2)
    cov_sim <- get_delaigle_cov(t, model = 2)
    true.dat <- mvtnorm::rmvnorm(n, rep(0, m), cov_sim)
    
    
  } else if (model==3){
    
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
    # (setting 1) gaussian measurement error over all grids
    sim.dat = true.dat +sigma.u* matrix(rnorm(length(true.dat),0,sigma.u),ncol=ncol(true.dat), nrow=nrow(true.dat))
  }
  #  (setting 2):t(3) measurement error over all grids
  if(setting==2){
    sim.dat = true.dat + (sigma.u/sqrt(3))*matrix(rt(length(true.dat),df=3),ncol=ncol(true.dat), nrow=nrow(true.dat))
  }
  
  if(setting==3){
    # (setting 3): asymmetric Laplace -- skewed
    noise= matrix(rALD(length(true.dat),mu=0,sigma=1,p= symm),ncol=ncol(true.dat), nrow=nrow(true.dat))
    #sim.dat = true.dat + sigma.u* noise	
    sim.dat = true.dat + sigma.u* noise /sqrt( var(as.vector(noise))) # stablize the variance
  }
  
  if(setting==4){
    con.id=sample(1:(length(true.dat)), size=length(true.dat)*0.3, replace=F)
    tmp.vec <- as.vector(true.dat)
    sim.vec <- tmp.vec + rnorm(length(tmp.vec),0,sigma.u)
    sim.vec[con.id] <- sim.vec[con.id] + rnorm(length(con.id),0,sigma.con)
    sim.dat = matrix(sim.vec, nrow=nrow(true.dat), ncol=ncol(true.dat))
  }
  
  if(setting==5){
    con.id1=sample(1:(length(true.dat)), size=length(true.dat)*0.15, replace=F)
    con.id2=sample(1:(length(true.dat)), size=length(true.dat)*0.15, replace=F)
    tmp.vec <- as.vector(true.dat)
    sim.vec <- tmp.vec + rnorm(length(tmp.vec),0,sigma.u)
    sim.vec[con.id1] <- sim.vec[con.id1] + rnorm(length(con.id1),0,sigma.con)
    sim.vec[con.id2] <- sim.vec[con.id2] + rALD(length(con.id2),mu=0,sigma=1,p= symm)
    sim.dat = matrix(sim.vec, nrow=nrow(true.dat), ncol=ncol(true.dat))
  }
  
  

  
  delta=matrix(1,nr=nrow(sim.dat), nc=ncol(sim.dat))
  
  if(m.prob > 0){  
    #dd=1.4; f=0.2 # setting from David Kraus paper
    
    for(jj in 1:nrow(delta)){
      
      # For example, for 20% of observations have partial missing
      if(rbinom(n=1, size=1, prob=m.prob)==1){
        
        u1=sqrt(runif(1)); u2=runif(1)
        Ci=dd*u1; Ei = f*u2
        
        #l.b = Ci-Ei; u.b = Ci+Ei  # lower and upper bound
        
        if((Ci-Ei)<0.01){
          l.b=1
        }else if((Ci-Ei)>0.99){
          l.b=length(t)
        } else {l.b=which(round(t,2)==round(Ci-Ei,2)) }
        #else {l.b=which(round(t,2)==round(Ci-Ei,2)*100) }
        
        
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
        # }  
        
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
  #print(k)
  
  output<-list()
  output$U = tmp$u
  output$V = tmp$v
  output$Dsq = tmp$d
  output$rec = Y.1
  
  return(output)
  
}



proposed_yj_V2<-function( x, r=15,  lambda.seq = seq(0,2, by = 0.5), maxiter = 200, maxiter.cv=200,
                       cv=FALSE, adhoc.cv=FALSE,  op.lambda, thresh = 5e-05, K.cv=5,sparsity.thres=0.8){
  
  
  tau_range<-seq(0.05, 0.95, by=0.05)
  b_tau =matrix(nrow=length(tau_range), ncol=(nrow(x)),0)
  
  # 현재 true.r = 5 일때 r=15 setting은 좋아보임
  
  
  ##############################
  # for CV
  
  # Yeonjoo: 미팅때 얘기했듯.. 모든 simulation set에 대해 다 CV 돌리면 너무 시간이 오래 걸리므로
  #         주어진 세팅에서 10-20개정도 돌려보고 대충 비슷한 lamba 고르면 그걸로 fix 하고 돌려되 될듯요 
  
  # lambda search 할때 어떤 lambda에서 에러가 나서 코드가 안돌아가는경우가 있는데. 
  # 그건 lambda가 너무 커서 다 0으로 보내버렸을경우라 그럴땐 lambda candidate을 좀 작게 잡아주면 되요
  
  
  alpha=0.5  # for elastic net
  
  if(adhoc.cv==TRUE){	
    lambda.sparsity <- 0; j<-0
    
    tmp.V.list=list(); tmp.U.list=list(); tmp.Dsq.list=list()
    #for(j in 1:length(lambda.seq)){
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
            Z = Y.0 + com_psi( x=E, b= b_tau, tau_range =tau_range, plot=FALSE)/2  #Psi(E, c=c)/2
            #tmp.als = svd.als(Z,rank.max =r.closs,lambda=lambda, final.svd = TRUE)
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
          Z = Y.0 + com_psi( x=E, b= b_tau, tau_range =tau_range)/2  #Psi(E, c=c)/2
          #tmp.als = svd.als(Z,rank.max =r.closs,lambda=lambda, final.svd = TRUE)
          tmp.als = svd.als.yj.L1(Z,rank.max =r,lambda=lambda, final.svd = FALSE, maxit = maxiter, alpha=alpha)
          
          U = tmp.als$u
          V = tmp.als$v
          Dsq = tmp.als$d
          
          Y.1 = U%*%diag(Dsq)%*%t(V)
          
          ratio=Frob(U.old,Dsq.old,V.old,U,Dsq,V)
          cond = !(ratio < thresh| iter > maxiter.cv)
          #cond = !(sum((Y.0 - Y.1)^2)/sum((Y.0^2)) < thresh | k > maxiter)  
          
          #err=c(err,sum((Y.0 - Y.1)^2)/sum((Y.0^2)))
          
          
          
        }, error = function(e) { 
          print("error")
          print(e)
          skip_sim <<- TRUE
        })
        
        #if(iter %% 100 == 0) (c(iter, ratio))  
      }
      
      
      #lambda.sparsity <- c(lambda.sparsity, sum(V==0)/length(V))
      lambda.sparsity  <- sum(V==0)/length(V)
      tmp.V.list[[j]] <- V; tmp.U.list[[j]] <- U; tmp.Dsq.list[[j]] <- Dsq
      
      print(lambda.sparsity)
    }
    
    
    #op.candidate = min(which(lambda.sparsity>0.73))
    op.lambda  = lambda
    #ind = min(which(lambda.sparsity>0.73))
    #op.lambda = lambda.seq[op.candidate]
    #if (lambda.sparsity[(op.candidate)]> 0.85){
    #  op.lambda = (lambda.seq[op.candidate-1]+ lambda.seq[op.candidate])/2
    #}else {
    #  op.lambda = lambda.seq[op.candidate]
    #}
    
    print(paste('cv selected=', op.lambda) )
  }
  
  ###########################################

  
  output <-list()
  
  output$Y.1= tmp.U.list[[j]] %*% diag(tmp.Dsq.list[[j]])%*%t(tmp.V.list[[j]])
  output$U = tmp.U.list[[j]]
  output$V = tmp.V.list[[j]]
  output$Dsq = tmp.Dsq.list[[j]]
  output$op.lambda = op.lambda
  #output$cv.rmse = mean(cv.rmse) 
  #output$cv.mad = mean(cv.mad)
  #output$sparsity = sparsity
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









surface.Krig<-
  function (object, grid.list = NULL, extrap = FALSE, graphics.reset = NULL, 
            xlab = NULL, ylab = NULL, main = NULL, zlab = NULL, zlim = NULL, 
            levels = NULL, type = "C", nx = 80, ny = 80, ...) 
  {
    out.p <- predictSurface(object, grid.list = grid.list, extrap = extrap, 
                            nx = nx, ny = ny, drop.Z = TRUE)
    if (!is.null(ylab)) 
      out.p$ylab <- ylab
    if (!is.null(xlab)) 
      out.p$xlab <- xlab
    if (!is.null(zlab)) 
      out.p$zlab <- zlab
    if (!is.null(main)) 
      out.p$main <- main
    plot.surface_yj(out.p, type = type, graphics.reset = graphics.reset, 
                    levels = levels, zlim = zlim, ...)
    invisible()
  }





plot.surface_yj<-
  function (x, main = NULL, type = "C", zlab = NULL, xlab = NULL, 
            ylab = NULL, levels = NULL, zlim = NULL, graphics.reset = NULL, 
            labcex = 0.6, add.legend = TRUE, ...) 
  {
    obj <- x
    old.par <- par(no.readonly = TRUE)
    if (is.na(match(type, c("b", "c", "C", "I", "p")))) {
      stop("plot type does not match b, C, I, or p.")
    }
    if (is.null(zlim)) {
      zlim = range(obj$z, na.rm = TRUE)
    }
    if (is.null(graphics.reset) & (type == "b")) {
      graphics.reset <- TRUE
    }
    else {
      graphics.reset <- FALSE
    }
    if (graphics.reset) {
      on.exit(par(old.par))
    }
    if (is.null(xlab)) {
      if (is.null(obj$xlab)) 
        xlab <- "X"
      else xlab <- obj$xlab
    }
    if (is.null(ylab)) {
      if (is.null(obj$ylab)) 
        ylab <- "Y"
      else ylab <- obj$ylab
    }
    if (is.null(zlab)) {
      if (is.null(obj$zlab)) 
        zlab <- "Z"
      else zlab <- obj$zlab
    }
    if (is.null(main)) 
      if (!is.null(obj$main)) 
        main <- obj$main
    if (type == "b") 
      set.panel(1, 2, TRUE)
    if (type == "p" | type == "b") {
      if (type == "b") {
        add.legend <- FALSE
        old.mar <- par()$mar
        par(mar = c(0, 5, 0, 0))
      }
      drape.plot(obj, xlab = xlab, ylab = ylab, zlab = zlab, 
                 zlim = zlim, add.legend = add.legend, ...)
      if (!is.null(main)) 
        title(main)
    }
    if (type == "I") {
      image.plot(obj$x, obj$y, obj$z, xlab = xlab, ylab = ylab, 
                 zlim = zlim, ...)
      if ((!is.null(main)) & type != "b") 
        title(main)
    }
    if (type == "c") {
      if (is.null(levels)) 
        levels <- pretty(obj$z[!is.na(obj$z)], 5)
      contour(obj$x, obj$y, obj$z, levels = levels, labcex = labcex, 
              lwd = .5, ...)
      if ((!is.null(main)) & type != "b") 
        title(main)
    }
    if (type == "b" | type == "C") {
      if (type == "b") {
        par(mar = old.mar)
      }
      image.plot(obj$x, obj$y, obj$z, xlab = xlab, ylab = ylab, 
                 graphics.reset = graphics.reset, zlim = zlim, ...)
      if (is.null(levels)) 
        levels <- pretty(obj$z[!is.na(obj$z)], 5)
      contour(obj$x, obj$y, obj$z, add = TRUE, levels = levels, 
              labcex = labcex, col = "dark grey", lwd = 1)
      if ((!is.null(main)) & type != "b") 
        title(main)
    }
    invisible()
  }









grid_finder<-function(lon, lat, full.lon, full.lat,plot=TRUE){
  
  
  library(fields)
  log.grid <- seq.default(0,360,, full.lon)
  lat.grid <- seq.default(-90,90,, full.lat)
  xy.grid <- expand.grid(log.grid,lat.grid)
  
  
  start = which.min(abs(xy.grid[,1]-lon[1] ) +   abs(xy.grid[,2]-lat[1] )    )
  mid1 = which.min(abs(xy.grid[,1]-lon[2] ) +   abs(xy.grid[,2]-lat[1] )    )
  mid2 = which.min(abs(xy.grid[,1]-lon[1] ) +   abs(xy.grid[,2]-lat[2] )    )
  final=  which.min(abs(xy.grid[,1]-lon[2] ) +   abs(xy.grid[,2]-lat[2] )    )
  
  
  st=ASIA= c(start:mid1 )
  for(m in 1:( (mid2-start)/full.lon)){
    ASIA<-  cbind( ASIA, st+ full.lon*m   )}
  
  
  if( plot==TRUE){
    A<-matrix(5,nrow=1, ncol= full.lon* full.lat)
    A[1,ASIA]<- 1
    image.fit <- as.image(A,x=xy.grid, grid=list(x=log.grid, y=lat.grid))
    image.plot(image.fit, main='Region', zlim=c(1:2))
    map("world2", ylim=c(-90,90), xlim = c(0,360), add = TRUE)
    
  }
  
  print( paste('number of lon : '  ,length(which(xy.grid[ASIA,]==xy.grid[ASIA[1],1])) ) )
  print(paste('number of lat : '  , length(which(xy.grid[ASIA,]==xy.grid[ASIA,2][1]))))
  return (as.vector(ASIA))
  
}




norm_eigenfunction<-function(est.pc,t ){
  
  est.pc <- apply(est.pc, 2, function(x) {
    x <- x/sqrt(trapzRcpp(t, x^2))
    return(x)
  })
  
  return(est.pc)
  
}


dens<-function(x1,y){
  index=density(y)$x 
  fit=density(y)$y[which.min(abs(index-x1))]
  return(fit)
}


sign_eigenfunction<-function(x, true){
  
  
  rmse1=  mean( ( true - x )^2  )
  rmse2= mean( ( true +x )^2  )
  
  if(rmse2 <rmse1){ x = - x }
  
  return(x)
  
}

