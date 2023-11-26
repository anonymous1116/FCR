#### High-dimensional Cost-constrained Regression via Non-convex Optimization #####
#### Regression Examples
#### 2/15/2021
source("../CCOfunctions.R")
library(MASS)
library(adagio)
library(glmnet)

maxiter = 1000
nsim=100
p0 = 32
tol = 1e-5

nums = c(500, 1000, 3000)
budgets = c(25, 50, 75, 100)
true.total.cost = c()
# three methods: Ours, HCR
#Four measure: L2 estimation, Prediction, FDR, FNR, time
results = array(NA,dim=c(length(nums),length(budgets),nsim,3,3)) 

for (j in 1:(length(nums)) ){
  n = nums[j]
  p=1000
  n.test=10000
  set.seed(1000)
  
  beta<-c()
  for(pnum in 1:(p0/2)){ beta <- c(beta, rnorm(2,mean = 2, sd = .5), rep(0,2) ) }
  
  beta <- c(beta, rep(0, p-length(beta)))
  
  group_size = 4
  group = rep(1:(p/group_size), each = group_size)
  group.costs= sample(1:10, size = length(unique(group)), T)
  
  #costs = c(); costs[1] = group.costs[1]; tmp = group.costs[1]
  #l=1
  #for (i in 2:length(group)) {
  #  if (group[i] > group[i-1]) {
  #    l = l+1
  #  }
  #  costs[i] = group.costs[l]
  #}
  
  true.total.cost = sum(group.costs[1:(p0/2)])
  cat("total cost:",  true.total.cost, ", p0:", sum(beta>0), "\n")
  for (k in 1:(length(budgets)) ){
  #budget=12
  budget=budgets[k]
  SIGMA=diag(p)
  
  for(sim in 1:nsim)
  { 
    time_result <-c()
    set.seed(1000*sim)
    X=mvrnorm(n,rep(0,p),SIGMA)
    
    # For checking 
    #apply(X,2,f<-function(tmp){ sum(tmp^2 )})
    Y=X%*%beta+0.5*rnorm(n)
    
    X.test = mvrnorm(n.test,rep(0,p),SIGMA)
    Y.test=X.test%*%beta+0.5*rnorm(n.test)
    
    
    #nlambda = 100
    #lambda.all=exp(seq(-5,-1,length=nlambda))
    
    ## FCR
    ptm=proc.time()[3]
    FCR.result=FCR.group(X,Y,group=group,group.costs,budget,maxiter,tol,rep(0,p))
    time.FCR = proc.time()[3] - ptm
    logl.FCR = FCR.result$obj
    
    pred.FCR = mean((Y.test - X.test%*%FCR.result$beta)^2)
    
    ## HCR method
    ptm=proc.time()[3]
    HCR.result = CEVS.group(X,Y,group,group.costs,budget,maxiter,tol,lambda = 0,alpha=1)
    time.HCR = proc.time()[3] - ptm
    logl.HCR = HCR.result$obj.values[HCR.result$niter]
    pred.HCR = mean((Y.test - X.test%*%HCR.result$beta)^2)
    
    ## lasso
    ptm=proc.time()[3]
    beta.lasso=GroupLasso.search2(X,Y,group,group.costs,budget)
    time.lasso =proc.time()[3] - ptm
    
    logl.lasso= mean((Y - X%*%beta.lasso)^2)/2
    pred.lasso = mean((Y.test - X.test%*%beta.lasso)^2)
    
    results[j,k,sim,,] = rbind(c(pred.FCR, logl.FCR, time.FCR), 
                               c(pred.HCR, logl.HCR, time.HCR),
                               c(pred.lasso, logl.lasso, time.lasso))
    
    print(c("n", nums[j], "budget", budgets[k], "nsim", sim))
  }  
  }
}

tmp = list(true.total.cost, nums, budgets, results)
save(tmp,file="simulation4.RData")
#load("simulation6.RData")
#save(tmp,file = "simulation6_tmp.RData")
results = tmp[[4]]
dim(results)
apply(results[1,1,,,], c(2,3), mean) 
apply(results[2,1,,,], c(2,3), mean) 
apply(results[3,1,,,], c(2,3), mean) 

apply(results[1,2,,,], c(2,3), mean) 
apply(results[2,2,,,], c(2,3), mean) 
apply(results[3,2,,,], c(2,3), mean) 

apply(results[1,3,,,], c(2,3), mean) 
apply(results[2,3,,,], c(2,3), mean) 
apply(results[3,3,,,], c(2,3), mean) 

apply(results[1,4,,,], c(2,3), mean) 
apply(results[2,4,,,], c(2,3), mean) 
apply(results[3,4,,,], c(2,3), mean) 

