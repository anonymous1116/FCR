source("./CCOfunctions.R")
library(MASS)
library(adagio)
library(glmnet)

maxiter = 1000
nsim=50
p0 = 24
tol = 1e-5

nums = c(500, 1000, 3000)
budgets = c(25,50,75,100)
true.total.cost = c()
# Five methods: FCR, HCR, Lasso 
#Four measure: Prediction,FNR, FPR, time
Example.result.all = array(NA,dim=c(length(nums),length(budgets),nsim,3,4))
obj_list<-c()

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
    
    SIGMA11 = matrix(0,p0,p0)
    
    rho = .5
    for(i1 in 1:p0){for (i2 in 1:p0){SIGMA11[i1,i2]=rho^abs(i1-i2)}}
    SIGMA12 = matrix(0, p0, p -p0)
    SIGMA21 = t(SIGMA12)
    SIGMA22 = matrix(0, p-p0, p-p0);diag(SIGMA22)<-1
    
    SIGMA = rbind(cbind(SIGMA11, SIGMA12),
                  cbind(SIGMA21, SIGMA22))
    
    true.total.cost = c(true.total.cost, sum(costs[1:p0]) )
    for(sim in 1:nsim){ 
      time_result <-c()
      set.seed(1000*sim)
      
      X=mvrnorm(n,rep(0,p),SIGMA)
      #X = scale(X)#/sqrt(n-1)
      Y=2*rbinom(n,1,exp(X%*%beta)/(1+exp(X%*%beta)))-1
      X.test=mvrnorm(n.test,rep(0,p),SIGMA)
      #X.test = scale(X.test)#/sqrt(n.test-1)
      Y.test=2*rbinom(n.test,1,exp(X.test%*%beta)/(1+exp(X.test%*%beta)))-1
      
      ## FCR
      ptm=proc.time()[3]
      FCR.result=FCR.logistic(X,Y,costs,budget,maxiter,tol,rep(0,p))
      time_result<- c(time_result, proc.time()[3] - ptm)
      beta.FCR = FCR.result$beta
      Prediction.FCR=mean(Y.test!=sign(X.test%*%FCR.result$beta))
      FPR.FCR=sum(beta==0&beta.FCR!=0)/sum(beta==0)
      FNR.FCR=sum(beta!=0&beta.FCR==0)/sum(beta!=0)
      
      ## HCR method
      ptm=proc.time()[3]
      HCR.result=CEVS.logistic(X,Y,costs,budget,maxiter,tol,rep(0,p))
      time_result<- c(time_result, proc.time()[3] - ptm)
      beta.HCR = HCR.result$beta
      Prediction.HCR=mean(Y.test!=sign(X.test%*%HCR.result$beta))
      FPR.HCR=sum(beta==0&beta.HCR!=0)/sum(beta==0)
      FNR.HCR=sum(beta!=0&beta.HCR==0)/sum(beta!=0)
      
      ## Lasso method
      ptm=proc.time()[3]
      beta.lasso=Lasso.logistic.search(X,Y,costs,budget)
      time_result<- c(time_result, proc.time()[3] - ptm)
      Prediction.lasso=mean(Y.test!=sign(X.test%*%beta.lasso))
      FPR.lasso=sum(beta==0&beta.lasso!=0)/sum(beta==0)
      FNR.lasso=sum(beta!=0&beta.lasso==0)/sum(beta!=0)
      
      ## Report Results
      Example.result=rbind(c(Prediction.FCR,FPR.FCR,FNR.FCR),
                            c(Prediction.HCR,FPR.HCR,FNR.HCR),
                            c(Prediction.lasso,FPR.lasso,FNR.lasso))
      
      Example.result = cbind(Example.result, time_result)
      Example.result.all[j,k,sim,,]=Example.result
      print(c("n", nums[j], "budget", budgets[k], "nsim", sim))
    }  
  }
}
tmp = list(true.total.cost, nums, budgets, Example.result.all)
#save(tmp,file="simulation5.RData")
