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
#Four measure: L2 estimation, Prediction, FDR, FNR, time
Example1.result.all = array(NA,dim=c(length(nums),length(budgets),nsim,3,5)) 

for (j in 1:(length(nums)) ){
  for (k in 1:(length(budgets)) ){
  n = nums[j]
  p= 1000 
  n.test=10000
  set.seed(1000)
  costs=sample(1:10,p,replace=T)
  #budget=12
  budget=budgets[k]
  beta=c(rnorm(p0/4,4,0.5))
  beta=c(beta, rnorm(p0/4,3,0.5))
  beta=c(beta, rnorm(p0/4,2,0.5))
  beta=c(beta, rnorm(p0/4,1,0.5))
  beta=c(beta, rep(0, p-length(beta)))
  
  SIGMA=diag(p)
  beta.interest=CE.projection(beta,costs,budget)
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
    Example1.result=rbind(c(L2.error.FCR,Prediction.FCR,FPR.FCR,FNR.FCR),
                          c(L2.error.HCR,Prediction.HCR,FPR.HCR,FNR.HCR),
                          c(L2.error.Lasso,Prediction.Lasso,FPR.Lasso,FNR.Lasso)
    )
    
    Example1.result = cbind(Example1.result, time_result)
    Example1.result.all[j,k,sim,,]=Example1.result
    print(c("n", nums[j], "budget", budgets[k], "nsim", sim))
  }  
}
}
tmp = list(true.total.cost, nums, budgets, Example1.result.all)
#save(tmp,file="simulation1.RData")
#save(tmp,file="simulation1_tmp.RData")

#Example1.result.all
#save.image("Example1_2_15_2021.Rdata")
#obj_list2 = obj_list
