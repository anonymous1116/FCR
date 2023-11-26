#### High-dimensional Cost-constrained Regression via Non-convex Optimization #####
#### Regression Examples
#### 2/15/2021
source("../CCOfunctions.R")
library(MASS)
library(adagio)
library(glmnet)

maxiter = 1000
nsim=100
p0 = 16
tol = 1e-5
nums = c(500, 1000, 3000) 
budgets = c(80)
rhos = c(0.2, 0.4, 0.6, 0.8)

true.total.cost = c()
# three methods: Ours, HCR, Lasso
#Four measure: L2 estimation, Prediction, FDR, FNR, time
Example.result.all = array(NA,dim=c(length(nums),length(rhos),nsim,3,5)) 

for (j in 1:(length(nums)) ){
  for (k in 1:(length(rhos)) ){
    n = nums[j]
    p= 1000 
    n.test=10000
    set.seed(1000)
    costs=sample(1:10,p,replace=T)
    #budget=12
    budget=budgets[1]
    beta=c(rnorm(p0/4,4,0.5))
    beta=c(beta, rnorm(p0/4,3,0.5))
    beta=c(beta, rnorm(p0/4,2,0.5))
    beta=c(beta, rnorm(p0/4,1,0.5))
    beta=c(beta, rep(0, p-length(beta)))
    
    
    SIGMA11 = matrix(0,p0,p0)
    
    rho = rhos[k]
    for(i1 in 1:p0){for (i2 in 1:p0){SIGMA11[i1,i2]=rho^abs(i1-i2)}}
    SIGMA12 = matrix(0, p0, p -p0)
    SIGMA21 = t(SIGMA12)
    SIGMA22 = matrix(0, p-p0, p-p0);diag(SIGMA22)<-1
    
    SIGMA = rbind(cbind(SIGMA11, SIGMA12),
                  cbind(SIGMA21, SIGMA22))
    
    ## Search the paramter of interest
    beta.interest=beta
    all.subset=as.matrix(expand.grid(rep(list(0:1),p0))) 
    feasible.subset.index=which(all.subset%*%costs[1:p0]<=budget)
    true.obj.values=rep(NA,length(feasible.subset.index))
    true.obj.values[1]=as.vector(t(beta)%*%SIGMA%*%beta)
    for(i in 2:length(feasible.subset.index))
    {
      a.set=as.vector(which(all.subset[feasible.subset.index[i],]==1))
      a.set.c=setdiff(1:p0,a.set)
      true.obj.values[i]=as.vector(t(beta[a.set.c])%*%(SIGMA[a.set.c,a.set.c]-
                                                         SIGMA[a.set.c,a.set]%*%solve(SIGMA[a.set,a.set])%*%SIGMA[a.set,a.set.c])%*%beta[a.set.c])
    }
    opt.index=which(true.obj.values==min(true.obj.values))[1]
    opt.subset=as.vector(which(all.subset[feasible.subset.index[opt.index],]==1))
    opt.subset.c=setdiff(1:p0,opt.subset)
    beta.interest[opt.subset.c]=0
    beta.interest[opt.subset]=beta[opt.subset]+as.vector(solve(SIGMA[opt.subset,opt.subset])%*%
                                                           SIGMA[opt.subset,opt.subset.c]%*%beta[opt.subset.c])
    
    
    true.total.cost = c(true.total.cost, sum(costs[beta.interest !=0]))
    
    for(sim in 1:nsim)
    { 
      time_result <-c()
      set.seed(1000*sim)
      X=mvrnorm(n,rep(0,p),SIGMA)
      #X = scale(X)/sqrt(n-1) #sqrt(n)-normalized column
      #X = scale(X) #sqrt(n)-normalized column
      
      # For checking 
      #apply(X,2,f<-function(tmp){ sum(tmp^2 )})
      Y=X%*%beta+0.5*rnorm(n)
      
      X.test = mvrnorm(n.test,rep(0,p),SIGMA)
      #X.test = scale(X.test)/sqrt(n.test-1)
      #X.test = scale(X.test)
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
      
      ## Lasso
      ptm=proc.time()[3]
      Lasso.result=Lasso.search2(X,Y,costs,budget)
      time_result=c(time_result, proc.time()[3] - ptm)
      
      beta.Lasso=Lasso.result
      L2.error.Lasso=sqrt(sum((beta.interest-beta.Lasso)^2))
      Prediction.Lasso=mean((Y.test-X.test%*%beta.Lasso)^2)
      FPR.Lasso=sum(beta.interest==0&beta.Lasso!=0)/sum(beta.interest==0)
      FNR.Lasso=sum(beta.interest!=0&beta.Lasso==0)/sum(beta.interest!=0)
      
      
      ## Report Results
      Example.result=rbind(c(L2.error.FCR,Prediction.FCR,FPR.FCR,FNR.FCR),
                            c(L2.error.HCR,Prediction.HCR,FPR.HCR,FNR.HCR),
                            c(L2.error.Lasso,Prediction.Lasso,FPR.Lasso,FNR.Lasso)
      )
      
      Example.result = cbind(Example.result, time_result)
      Example.result.all[j,k,sim,,]=Example.result
      print(c("n", nums[j], "rhos", rhos[k], "nsim", sim))
    }  
  }
}
tmp = list(true.total.cost, nums, rhos, Example.result.all)
#save(tmp,file="simulation2.RData")
#save(tmp,file="simulation2_tmp.RData")

#Example.result.all
#save.image("Example_2_15_2021.Rdata")
#obj_list2 = obj_list
