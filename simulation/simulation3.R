#### High-dimensional Cost-constrained Regression via Non-convex Optimization #####
#### Regression Examples
#### 2/15/2021
source("../CCOfunctions.R")
library(MASS)
library(adagio)
library(glmnet)

maxiter = 1000
nsim=100
p0 = 24
tol = 1e-5
nums = c(500, 1000, 3000) 
budgets = c(25, 50, 75, 100)
true.total.cost = c()
# three methods: Ours, HCR, Lasso
#Four measure: Pred Error, time
Example.result.all = array(NA,dim=c(length(nums),length(budgets),nsim,3,4)) 

for (j in 1:(length(nums)) ){
  for (k in 1:(length(budgets)) ){
    n = nums[j]
    p= 1000 
    n.test=10000
    set.seed(1000)
    costs=sample(1:50,p,replace=T)
    #budget=12
    budget=budgets[k]
    beta=c(rnorm(p0/4,4,0.5))
    beta=c(beta, rnorm(p0/4,3,0.5))
    beta=c(beta, rnorm(p0/4,2,0.5))
    beta=c(beta, rnorm(p0/4,1,0.5))
    beta=c(beta, rep(0, p-length(beta)))
    
    SIGMA11 = matrix(0,p0,p0)
    
    rho = .3
    for(i1 in 1:p0){for (i2 in 1:p0){SIGMA11[i1,i2]=rho^abs(i1-i2)}}
    SIGMA12 = matrix(0, p0, p -p0)
    SIGMA21 = t(SIGMA12)
    SIGMA22 = matrix(0, p-p0, p-p0);diag(SIGMA22)<-1
    
    SIGMA = rbind(cbind(SIGMA11, SIGMA12),
                  cbind(SIGMA21, SIGMA22))
    
    true.total.cost = c(true.total.cost, sum(costs[1:p0]) )
    
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
      Prediction.FCR=mean((Y.test-X.test%*%beta.FCR)^2)
      FPR.FCR=sum(beta==0&beta.FCR!=0)/sum(beta==0)
      FNR.FCR=sum(beta!=0&beta.FCR==0)/sum(beta!=0)
      
      
      ## HCR 
      ptm=proc.time()[3]
      HCR.result=HCR(X,Y,costs,budget,maxiter,tol, rep(0,p))
      time_result=c(time_result, proc.time()[3] - ptm)
      
      beta.HCR=HCR.result$beta
      Prediction.HCR=mean((Y.test-X.test%*%beta.HCR)^2)
      FPR.HCR=sum(beta==0&beta.HCR!=0)/sum(beta==0)
      FNR.HCR=sum(beta!=0&beta.HCR==0)/sum(beta!=0)
      
      ## Lasso
      ptm=proc.time()[3]
      Lasso.result=Lasso.search2(X,Y,costs,budget)
      time_result=c(time_result, proc.time()[3] - ptm)
      
      beta.Lasso=Lasso.result
      Prediction.Lasso=mean((Y.test-X.test%*%beta.Lasso)^2)
      FPR.Lasso=sum(beta==0&beta.Lasso!=0)/sum(beta==0)
      FNR.Lasso=sum(beta!=0&beta.Lasso==0)/sum(beta!=0)
      
      
      ## Report Results
      Example.result=rbind(c(Prediction.FCR,FPR.FCR,FNR.FCR),
                            c(Prediction.HCR,FPR.HCR,FNR.HCR),
                            c(Prediction.Lasso,FPR.Lasso,FNR.Lasso))
      
      Example.result = cbind(Example.result, time_result)
      Example.result.all[j,k,sim,,]=Example.result
      print(c("n", nums[j], "budget", budgets[k], "nsim", sim))
    }  
  }
}
tmp = list(true.total.cost, nums, budgets, Example.result.all)
#save(tmp,file="simulation3.RData")
#save(tmp,file="simulation3_tmp.RData")

#Example.result.all
#save.image("Example_2_15_2021.Rdata")
#obj_list2 = obj_list
