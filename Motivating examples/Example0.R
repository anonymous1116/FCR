  #### High-dimensional Cost-constrained Regression via Non-convex Optimization #####
  #### Regression Examples
  #### 2/15/2021
  source("./CCOfunctions.R")
  library(MASS)
  library(adagio)
  library(glmnet)
  
  maxiter = 1000
  nsim=100
  p0 = 8
  tol = 1e-5
  
  nums = c(50)
  true.total.cost = c()
  # three methods: Ours, HCR, LASSO 
  #Four measure: L2 estimation, Prediction, FDR, FNR, time
  Example.result.all = array(NA,dim=c(length(nums),nsim,3,5)) 
  beta.result.all = array(NA, dim=c(length(nums), nsim,4, 1000))
  
  n = 2000
  p = 1000
  n.test=10000
  set.seed(1000)
  budget=nums[1]
  
  costs = rep(20, p0/2)
  costs = c(costs,rep(10, p0/2))
  costs=c(costs,sample(1:2,p-p0,replace=T))
  
  #beta=c(rnorm(p0/4,2,0.1))
  #beta=c(beta, rnorm(p0/4,2,0.1))
  #beta=c(beta, rnorm(p0/4,2,0.1))
  #beta=c(beta, rnorm(p0/4,2,0.1))
  
  beta = rep(2, p0)
  beta=c(beta, rep(0, p-length(beta)))
  
  true.total.cost = c(true.total.cost, sum(costs[1:p0]) )
  
  
  SIGMA11 = matrix(0,3*p0,3*p0)
  diag(SIGMA11[1:p0,(p0+1):(2*p0)])=0.8
  diag(SIGMA11[1:p0,(2*p0+1):(3*p0)])=0.8
  
  diag(SIGMA11[(2*p0+1):(3*p0),1:p0])=0.8
  diag(SIGMA11[(2*p0+1):(3*p0),1:p0])=0.8
  diag(SIGMA11)=1
  
  SIGMA12 = matrix(0, 3*p0, p - 3*p0)
  SIGMA21 = t(SIGMA12)
  SIGMA22 = matrix(0, p-3*p0, p-3*p0);diag(SIGMA22)<-1
  
  SIGMA = rbind(cbind(SIGMA11, SIGMA12),
                cbind(SIGMA21, SIGMA22))
  
  
  ## Search the paramter of interest
  beta.interest=beta
  all.subset=as.matrix(expand.grid(rep(list(0:1),3*p0))) 
  feasible.subset.index=which(all.subset%*%costs[1:(3*p0)]<=budget)
  true.obj.values=rep(NA,length(feasible.subset.index))
  true.obj.values[1]=as.vector(t(beta)%*%SIGMA%*%beta)
  for(i in 2:length(feasible.subset.index)){
    a.set=as.vector(which(all.subset[feasible.subset.index[i],]==1))
    a.set.c=setdiff(1:(3*p0),a.set)
    true.obj.values[i]=as.vector(t(beta[a.set.c])%*%(SIGMA[a.set.c,a.set.c]-
                                                       SIGMA[a.set.c,a.set]%*%solve(SIGMA[a.set,a.set])%*%SIGMA[a.set,a.set.c])%*%beta[a.set.c])
  }
  opt.index=which(true.obj.values==min(true.obj.values))[1]
  opt.subset=as.vector(which(all.subset[feasible.subset.index[opt.index],]==1))
  opt.subset.c=setdiff(1:(3*p0),opt.subset)
  beta.interest[opt.subset.c]=0
  beta.interest[opt.subset]=beta[opt.subset]+as.vector(solve(SIGMA[opt.subset,opt.subset])%*%
                                                         SIGMA[opt.subset,opt.subset.c]%*%beta[opt.subset.c])
  rm(all.subset)
  rm(feasible.subset.index)
  
  for (j in 1:(length(nums)) ){
    
    for(sim in 1:nsim)
    { 
      time_result <-c()
      set.seed(1000*sim)
      X = mvrnorm(n,rep(0,p),SIGMA)
      
      # For checking 
      #apply(X,2,f<-function(tmp){ sum(tmp^2 )})
      Y=X%*%beta+0.5*rnorm(n)
      
      X.test = mvrnorm(n.test,rep(0,p),SIGMA)
      #X.test = scale(X.test)/sqrt(n.test-1)
      Y.test=X.test%*%beta+0.5*rnorm(n.test)
      
      ## FCR
      ptm=proc.time()[3]
      FCR.result=FCR(X,Y,costs,budget,maxiter,tol, rep(0,p))
      time_result=c(time_result, proc.time()[3] - ptm)
      
      beta.FCR=FCR.result$beta
      L2.error.FCR=sqrt(sum((beta.interest-beta.FCR)^2))
      Prediction.FCR=mean((Y.test-X.test%*%beta.FCR)^2)
      FPR.FCR=sum(beta.interest==0&beta.FCR!=0)/sum(beta.interest==0)
      FNR.FCR=sum(beta.interest!=0&beta.FCR==0)/sum(beta.interest!=0)
      
      ## HCR
      ptm=proc.time()[3]
      HCR.result=HCR(X,Y,costs,budget,maxiter,tol, rep(0,p))
      time_result=c(time_result, proc.time()[3] - ptm)
      
      beta.HCR=HCR.result$beta
      L2.error.HCR=sqrt(sum((beta.interest-beta.HCR)^2))
      Prediction.HCR=mean((Y.test-X.test%*%beta.HCR)^2)
      FPR.HCR=sum(beta.interest==0&beta.HCR!=0)/sum(beta.interest==0)
      FNR.HCR=sum(beta.interest!=0&beta.HCR==0)/sum(beta.interest!=0)
      
      ## Lasso method 2
      ptm=proc.time()[3]
      beta.lasso2=Lasso.search2(X,Y,costs,budget)
      time_result=c(time_result, proc.time()[3] - ptm)
      
      L2.error.lasso2=sqrt(sum((beta.interest-beta.lasso2)^2))
      Prediction.lasso2=mean((Y.test-X.test%*%beta.lasso2)^2)
      FPR.lasso2=sum(beta.interest==0&beta.lasso2!=0)/sum(beta.interest==0)
      FNR.lasso2=sum(beta.interest!=0&beta.lasso2==0)/sum(beta.interest!=0)
      
      
      Example.result=rbind(c(L2.error.FCR,Prediction.FCR,FPR.FCR,FNR.FCR),
                           c(L2.error.HCR,Prediction.HCR,FPR.HCR,FNR.HCR),
                           c(L2.error.lasso2,Prediction.lasso2,FPR.lasso2,FNR.lasso2)
                          )
      betas = rbind(beta.interest, beta.FCR, beta.HCR, beta.lasso2)
      
      beta.result.all[j,sim,,] = betas
      Example.result = cbind(Example.result, time_result)
      Example.result.all[j,sim,,]=Example.result
      print(c("n", nums[j], "nsim", sim))
    }  
  }
  
  tmp = list(true.total.cost, nums, Example.result.all, beta.result.all, X[1:(3*p0),1:(3*p0)])
  save(tmp,file="Example0_tmp.RData")
  
  