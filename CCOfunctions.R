HCR<-function(X,Y,costs,budget,maxiter,epsilon,beta.initial)
{
  ######### This function performs cost effective variable 
  ######### selection by iterative dynamic programming
  # X:  nxp design matrix
  # Y:  nx1 response variable
  # costs: px1 vector. The j-th element is the cost to collect the j-th predictor
  # budget: the budget to build the model
  # maxiter: the maximum number of interations for the iterative dynamic programming
  # epsilon: the algorithm will terminate if the change of the objective function is smaller than epsilon
  n=nrow(X)
  p=ncol(X)
  if (n>p){L=max(eigen(t(X)%*%X/n)$values) + 0.1}
  else {L=max(eigen(X%*%t(X)/n)$values) + 0.1}
  
  beta.m=beta.initial #rep(0,p)
  obj.value=rep(NA,maxiter)
  iter.index=0
  for(m in 1:maxiter)
  {
    iter.index=iter.index+1
    temp=sum((Y-X%*%beta.m)^2)/(2*n) 
    a.vector=as.vector(beta.m+t(X)%*%(Y-X%*%beta.m)/(n*L))
    result=knapsack(costs,(a.vector^2),budget)
    z.vector=rep(0,p)
    z.vector[result$indices]=1
    beta.m=a.vector*z.vector
    obj.value[m]=sum((Y-X%*%beta.m)^2)/(2*n)  
    if(abs(temp-obj.value[m])<=epsilon){break}
  }
  if(iter.index==maxiter)
  {warning(paste("Algorithm does not converge in ", maxiter, " iterations!", sep=""))}
  
  list(obj.values=obj.value[1:iter.index], beta=beta.m, niter=iter.index)
}

FCR<-function(X,Y,costs,budget,maxiter,epsilon,beta.initial)
{
  ######### This function performs cost effective variable 
  ######### selection by iterative dynamic programming
  # X:  nxp design matrix
  # Y:  nx1 response variable
  # costs: px1 vector. The j-th element is the cost to collect the j-th predictor
  # budget: the budget to build the model
  # maxiter: the maximum number of interations for the iterative dynamic programming
  # epsilon: the algorithm will terminate if the change of the objective function is smaller than epsilon
  n=nrow(X)
  p=ncol(X)
  
  beta.m=beta.initial #rep(0,p)
  obj.value=rep(NA,maxiter)
  iter.index=0
  obj.min=c(sum((Y-X%*%beta.m)^2)/(2*n))
  
  W = 1/diag(t(X) %*% X)
  d.m = diag(W) %*% t(X) %*% (Y - X %*% beta.m)
  
  
  for(m in 1:maxiter)
  {
    iter.index=iter.index+1
    temp=sum((Y-X%*%beta.m)^2)/(2*n) 
    a.vector=as.vector(beta.m+ d.m)
    
    result=knapsack(costs,(a.vector^2),budget)
    
    beta.m <-rep(0,p)
    d.m <- rep(0,p)
    
    X.A.m = X[,result$indices]
    X.I.m = X[,-result$indices]
    
    beta.A.m = solve(t(X.A.m) %*% X.A.m) %*% t(X.A.m) %*% Y
    beta.m[result$indices] = beta.A.m
    
    W.I.m = 1/diag(t(X.I.m) %*% X.I.m)
    d.I.m = diag(W.I.m) %*% t(X.I.m) %*% (Y - X.A.m %*% beta.A.m)
    d.m[-result$indices] = d.I.m
    
    obj.value[m]=sum((Y-X%*%beta.m)^2)/(2*n)
    
    obj.min = min(obj.min,obj.value[m])
    if (obj.min == obj.value[m]){beta.min = beta.m; obj.min = obj.min}
    if(abs(temp-obj.value[m])<=epsilon){break}
    if( (m>2) & (min(abs(obj.value[1:(m-1)]-obj.value[m]))<=1e-15) ){break}
  }
  if(iter.index==maxiter)
  {warning(paste("Algorithm does not converge in ", maxiter, " iterations! when n=", n, sep=""))}
  
  list(obj.values=obj.value[1:iter.index], beta=beta.min, niter=iter.index, obj = obj.min)
}


Lasso.search<-function(X,Y,costs,budget)
{
  lasso.fit=glmnet(X,Y,family="gaussian",alpha=1,nlambda=1000,intercept=F)
  beta.all=t(as.matrix(lasso.fit$beta))
  avail.index=which((beta.all!=0)%*%costs<=budget)
  
  cv.lasso=cv.glmnet(X,Y,family="gaussian",alpha=1,type.measure="mse",
                     lambda=lasso.fit$lambda,intercept=F)
  opt.index=which.min(cv.lasso$cvm[avail.index])
  beta.est=as.vector(beta.all[opt.index,])
  beta.est
}


Lasso.search2<-function(X,Y,costs,budget)
{
  lasso.fit=glmnet(X,Y,family="gaussian",alpha=1,nlambda=1000,intercept=F)
  beta.all=t(as.matrix(lasso.fit$beta))
  avail.index=which((beta.all!=0)%*%costs<=budget)
  
  cv.lasso=cv.glmnet(X,Y,family="gaussian",alpha=1,type.measure="mse",
                     lambda=lasso.fit$lambda,intercept=F)
  opt.index=which.min(cv.lasso$cvm[avail.index])
  beta.est=as.vector(beta.all[opt.index,])
  a.set =(beta.est!=0) 
  tmp = X[,a.set]
  
  if (sum(a.set) == 1){
    beta.est[a.set] = 1/sum(tmp*tmp) * sum(tmp*Y) 
  }
  else if (sum(a.set) == 0) {beta.est = rep(0,p)}
  else{
    tmp = as.matrix(tmp)
    beta.est.A = solve(t(tmp)%*% tmp) 
    beta.est.A = beta.est.A %*% t(tmp) %*% Y
    beta.est[a.set] = beta.est.A
  }
  return (beta.est)
}

Lasso.search3<-function(X,Y,costs,budget)
{
  cv.lasso=cv.glmnet(X,Y,family="gaussian",alpha=1,type.measure="mse",nlambda= 100,intercept=F)
  lasso.fit=glmnet(X,Y,family="gaussian",alpha=1,lambda = cv.lasso$lambda.1se,intercept=F)
  
  beta.all = as.matrix(coef(lasso.fit,s="lambda.1se"))[-1]
  p0 = sum(beta.all != 0)
  ind = (beta.all != 0)
  all.subset=as.matrix(expand.grid(rep(list(0:1),p0)))
  feasible.subset.index=which(all.subset%*%costs[1:p0]<=budget)
  obj.values=rep(NA,length(feasible.subset.index))
  
  X_reduce = X[,ind]
  for(i in 2:length(feasible.subset.index)){
    a.set=as.vector(which(all.subset[feasible.subset.index[i],]==1))
    
    if (length(a.set) == 1){
      tmp = X_reduce[,a.set]
      beta.est = 1/sum(tmp*tmp) * sum(tmp*Y) 
      obj.values[i] = sum((Y-tmp*beta.est)^2)/(2*n)  
    }else{
      tmp = X_reduce[,a.set]
      tmp = as.matrix(tmp)
      beta.est = solve(t(tmp)%*% tmp) 
      beta.est = beta.est %*% t(tmp) %*% Y
      obj.values[i] = sum((Y-tmp%*%beta.est)^2)/(2*n) 
    }
  }
  
  ind_min = feasible.subset.index[which.min(obj.values)]
  obj.values = min(obj.values, na.rm = TRUE)
  a.set = all.subset[ind_min,]
  ind[ind==TRUE]<-a.set
  tmp = X[,ind == 1]
  if (sum(a.set) == 1){
    beta.est = 1/sum(tmp*tmp) * sum(tmp*Y) 
  }else{
    tmp = as.matrix(tmp)
    beta.est = solve(t(tmp)%*% tmp) 
    beta.est = beta.est %*% t(tmp) %*% Y
  }
  beta.m<-rep(0, dim(X)[2])
  beta.m[ind==1] = beta.est
  return(beta.m)
}

CE.projection<-function(beta0,costs,budget)
{
  result=knapsack(costs,(beta0^2),budget)
  z.vector=rep(0,length(beta0))
  z.vector[result$indices]=1
  beta0*z.vector
}

##########################################################################################
########################## Logistic Regression Functions #################################
##########################################################################################

FCR.logistic<-function(X,Y,costs,budget,maxiter,epsilon,beta.initial){
  ######### This function performs cost effective variable 
  ######### selection by SDAR
  # X:  nxp design matrix
  # Y:  nx1 response variable
  # costs: px1 vector. The j-th element is the cost to collect the j-th predictor
  # budget: the budget to build the model
  # maxiter: the maximum number of interations for the iterative dynamic programming
  # epsilon: the algorithm will terminate if the change of the objective function is smaller than epsilon
  n=nrow(X)
  p=ncol(X)
  
  beta.m=beta.initial #rep(0,p)
  obj.min=mean(log(1+exp((-Y)*(X%*%beta.m))))
  
  obj.value=rep(NA,maxiter)
  iter.index=0
  W = 1/diag(t(X) %*% X) * 4
  d.m = diag(W) %*% t(X) %*% (Y*(exp((-Y)*(X%*%beta.m))/(1+exp((-Y)*(X%*%beta.m)))))
  
  for(m in 1:maxiter){
    iter.index=iter.index+1
    temp=mean(log(1+exp((-Y)*(X%*%beta.m))))
    
    a.vector=as.vector(beta.m+ d.m)
    a.vector[is.na(a.vector)] = 0
    result=knapsack(costs,(a.vector^2),budget)
    
    beta.m <-rep(0,p)
    d.m <- rep(0,p)
    
    X.A.m = X[,result$indices]
    X.I.m = X[,-result$indices]
    beta.A.m = as.vector(glmnet(X.A.m,Y,family="binomial",alpha=1,lambda=0,intercept=F)$beta)
    
    
    beta.m[result$indices]<-beta.A.m
    W.I.m = as.vector(1/diag(t(X.I.m) %*% X.I.m))
    if (length(W.I.m) == 1){
      d.I.m = t(X.I.m) * W.I.m
      } else {
      d.I.m = diag(c(W.I.m)) %*% t(X.I.m) 
      }
    
    d.I.m = d.I.m %*% (Y*exp(-Y*(X.A.m%*%beta.A.m)/(1+exp(-Y*(X.A.m%*%beta.A.m))) ) ) * 4
    d.m[-result$indices] = d.I.m
    obj.value[m]=mean(log(1+exp((-Y)*(X%*%beta.m))))
    
    obj.min = min(obj.min,obj.value[m])
    if (obj.min == obj.value[m]){beta.min = beta.m; obj.min = obj.min}
    if (obj.value[m] == Inf){break}
    if(abs(temp-obj.value[m])<=epsilon| obj.value){break}
    if( (m>2) & (min(abs(obj.value[1:(m-1)]-obj.value[m]))<=1e-15) ){break}
  }
  if(iter.index==maxiter)
  {warning(paste("Algorithm does not converge in ", maxiter, " iterations! when n=", n, sep=""))}
  
  list(obj.values=obj.value[1:iter.index], beta=beta.min, niter=iter.index, obj = obj.min)
}


CEVS.logistic<-function(X,Y,costs,budget,maxiter,epsilon,beta.initial)
{
  ######### This function performs cost effective variable 
  ######### selection by iterative dynamic programming
  # X:  nxp design matrix
  # Y:  nx1 response variable
  # costs: px1 vector. The j-th element is the cost to collect the j-th predictor
  # budget: the budget to build the model
  # maxiter: the maximum number of interations for the iterative dynamic programming
  # epsilon: the algorithm will terminate if the change of the objective function is smaller than epsilon
  n=nrow(X)
  p=ncol(X)
  L=0.25*max(eigen(t(X)%*%X/n)$values)+0.1
  
  beta.m=beta.initial #rep(0,p)
  obj.value=rep(NA,maxiter)
  iter.index=0
  for(m in 1:maxiter)
  {
    iter.index=iter.index+1
    temp=mean(log(1+exp((-Y)*(X%*%beta.m))))
    a.vector=as.vector(beta.m+t(X)%*%(Y*(exp((-Y)*(X%*%beta.m))/(1+exp((-Y)*(X%*%beta.m)))))/(n*L))
    result=knapsack(costs,(a.vector^2),budget)
    z.vector=rep(0,p)
    z.vector[result$indices]=1
    beta.m=a.vector*z.vector
    obj.value[m]=mean(log(1+exp((-Y)*(X%*%beta.m))))
    if(abs(temp-obj.value[m])<=epsilon){break}
  }
  if(iter.index==maxiter)
  {warning(paste("Algorithm does not converge in ", maxiter, " iterations!", sep=""))}
  
  list(obj.values=obj.value[1:iter.index], beta=beta.m, niter=iter.index)
}

Lasso.logistic.search<-function(X,Y,costs,budget)
{
  lasso.fit=glmnet(X,Y,family="binomial",alpha=1,nlambda=1000,intercept=F)
  beta.all=t(as.matrix(lasso.fit$beta))
  cv.lasso=cv.glmnet(X,Y,family="binomial",alpha=1,type.measure="deviance",
                     lambda=lasso.fit$lambda,intercept=F)
  avail.index=which((beta.all!=0)%*%costs<=budget)
  opt.index=which.min(cv.lasso$cvm[avail.index])
  beta.est=as.vector(beta.all[opt.index,])
  beta.est
}



######################################################################################
################################ "Group" functions ###################################
######################################################################################

gknapsack<-function(weight,profit,cap)
{
  # This function uses dynamic programming to solve the 0-1 knapsack problem
  # Maximize sum(profit*beta) such that sum(weight*beta) <= cap, 
  # where beta is a vector with beta[i] == 0 or 1.
  p=length(profit)
  beta=rep(NA,p)
  nn=sum(weight>cap)
  if(nn==0)
  {
    result=knapsack(weight,profit,cap)
    z.vector=rep(0,p)
    z.vector[result$indices]=1
    beta=z.vector
  }
  if(nn>0&nn<p)
  {
    index=which(weight>cap)
    beta[index]=0
    result=knapsack(weight[-index],profit[-index],cap)
    z.vector=rep(0,p-nn)
    z.vector[result$indices]=1
    beta[-index]=z.vector
  }
  if(nn==p){beta=rep(0,p)}
  beta
}

Elastic.group.search<-function(X,Y,group,group.costs,budget,alpha)
{
  ######### This function searches the best elastic net solution that satisfies the budget 
  ######### constraint
  # X:  nxp design matrix
  # Y:  nx1 response variable
  # costs: px1 vector. The j-th element is the cost to collect the j-th predictor
  # budget: the budget to build the model
  # alpha: the parameter used in the elastic net
  # alpha=1 for Lasso and alpha=0 for the ridge regression
  n.group=length(group.costs)
  elastic.fit=glmnet(X,Y,family="gaussian",alpha=alpha,nlambda=100,intercept=F)
  beta.all=as.matrix(elastic.fit$beta)
  nfolds=5
  cvfold.out=cvFolds(nrow(X),K=nfolds,type="interleaved")
  cv.elastic=cv.glmnet(X,Y,family="gaussian",alpha=alpha,type.measure="mse",
                       lambda=elastic.fit$lambda,intercept=F,foldid=cvfold.out$which)
  model.cost=rep(NA,ncol(beta.all))
  for(i in 1:ncol(beta.all)){
    temp=beta.all[,i]
    group.indicator=rep(NA,n.group)
    for(j in 1:n.group){
      group.indicator[j]=as.integer(sum(temp[group==j]^2)!=0)
    }
    model.cost[i]=sum(group.indicator*group.costs)
  }
  avail.index=which(model.cost<=budget)
  opt.index=which.min(cv.elastic$cvm[avail.index])
  beta.est=as.vector(beta.all[,opt.index])
  beta.est
}

CEVS.group<-function(X,Y,group,group.costs,budget,maxiter,epsilon,lambda,alpha)
{
  ######### This function performs cost effective variable selection by iterative DP
  ######### (We consider the group cost rather than the separae cost for each variable)
  # X:  nxp design matrix
  # Y:  nx1 response variable
  # group:  px1 vector (e.g., 1,1,2,3,4,... where 1 is the 1st group, 2 is the 2nd group...)
  # group.costs: gx1 vector. The j-th element is the cost to collect the j-th group
  # budget: the budget to build the model
  # maxiter: the maximum number of interations for the iterative dynamic programming
  # epsilon: the algorithm will terminate if the change of the objective function is 
  # smaller than epsilon
  # lambda and alpha: parameters in the following elastic net penalty: 
  #                lambda*(alpha*|beta|+0.5*(1-alpha)*beta^2)
  # alpha=1 for Lasso and alpha=0 for the ridge regression
  n=nrow(X)
  p=ncol(X)
  n.group=length(group.costs)
  if (n>p){L=max(eigen(t(X)%*%X/n)$values) + 0.1}
  else {L=max(eigen(X%*%t(X)/n)$values) + 0.1}
  beta.m=rep(0,p)
  obj.value=rep(NA,maxiter)
  logl.value=rep(NA,maxiter)
  
  iter.index=0
  for(m in 1:maxiter)
  {
    iter.index=iter.index+1
    temp=sum((Y-X%*%beta.m)^2)/(2*n)+lambda*alpha*sum(abs(beta.m))
    +lambda*((1-alpha)/2)*sum(beta.m^2) 
    a.vector=as.vector(beta.m+t(X)%*%(Y-X%*%beta.m)/(n*L))
    a.vector.soft=(L/(L+lambda*(1-alpha)))*sign(a.vector-alpha*lambda/L)*
      pmax(abs(a.vector)-alpha*lambda/L,0)
    temp1=(a.vector^2)*(L^2)/(2*(L+lambda*(1-alpha)))
    temp2=-abs(a.vector)*alpha*lambda*L/(L+lambda*(1-alpha))
    temp3=(alpha^2)*(lambda^2)/(2*(L+lambda*(1-alpha)))
    obj.coefficients.seperate=(temp1+temp2+temp3)*(1+sign(abs(a.vector)-alpha*lambda/L))/2
    obj.coefficients=rep(NA,n.group)
    for(i in 1:n.group){obj.coefficients[i]=sum(obj.coefficients.seperate[group==i])}
    z.vector=gknapsack(group.costs,obj.coefficients,budget)
    beta.m=rep(NA,p)
    for(i in 1:n.group){beta.m[group==i]=a.vector.soft[group==i]*z.vector[i]}
    obj.value[m]=sum((Y-X%*%beta.m)^2)/(2*n)+lambda*alpha*sum(abs(beta.m))
    +lambda*((1-alpha)/2)*sum(beta.m^2)   
    logl.value[m] = sum((Y-X%*%beta.m)^2)/(2*n)
    if(abs(temp-obj.value[m])<=epsilon){break}
  }
  if(iter.index==maxiter)
  {warning(paste("Algorithm does not converge in ", maxiter, " iterations!", sep=""))}
  
  list(obj.values=obj.value[1:iter.index], logl.value = logl.value[1:iter.index], beta=beta.m, niter=iter.index)
}


FCR.group<-function(X,Y,group,group.costs,budget,maxiter,epsilon,beta.initial){
  ######### This function performs cost effective variable selection by iterative DP
  ######### (We consider the group cost rather than the separae cost for each variable)
  # X:  nxp design matrix
  # Y:  nx1 response variable
  # group:  px1 ordered vector (e.g., 1,1,2,3,4,... where 1 is the 1st group, 2 is the 2nd group...)
  # group.costs: gx1 vector. The j-th element is the cost to collect the j-th group
  # budget: the budget to build the model
  # maxiter: the maximum number of interations for the iterative dynamic programming
  # epsilon: the algorithm will terminate if the change of the objective function is 
  # smaller than epsilon
  n=nrow(X)
  p=ncol(X)
  n.group=length(unique(group))
  
  
  beta.m=beta.initial #rep(0,p)
  obj.value=rep(NA,maxiter)
  obj.min=sum((Y-X%*%beta.m)^2)/(2*n) 
  
  iter.index=0
  W = 1/diag(t(X) %*% X)
  d.m = diag(W) %*% t(X) %*% (Y - X %*% beta.m)
  
  for(m in 1:maxiter){
    iter.index = iter.index+1
    temp=sum((Y-X%*%beta.m)^2)/(2*n)
    
    tmp = beta.m + d.m
    
    l = 1
    a.vector = rep(0,n.group); a.vector[1] = tmp[1]
    for (j in 2:(length(group)) ){
      if (group[j] > group[j-1]) {
        l = l+1
      }
      a.vector[l] = a.vector[l] + tmp[j]
    }
    result = knapsack(group.costs,(a.vector^2),budget)
    ind.group = result$indices
    
    ind = rep(FALSE, p)
    for (j in 1:(length(ind.group)) ){
      ind = ind + (group == ind.group[j])
    }
    ind = (ind==1)
    
    beta.m <-rep(0,p); d.m <- rep(0,p)
    
    X.A.m = X[,ind]
    X.I.m = X[,-ind]
    
    beta.A.m = solve(t(X.A.m) %*% X.A.m) %*% t(X.A.m) %*% Y
    beta.m[ind] = beta.A.m
    
    W.I.m = 1/diag(t(X.I.m) %*% X.I.m)
    d.I.m = diag(W.I.m) %*% t(X.I.m) %*% (Y - X.A.m %*% beta.A.m)
    d.m[-ind] = d.I.m
    
    obj.value[m]=sum((Y-X%*%beta.m)^2)/(2*n)
    
    obj.min = min(obj.min,obj.value[m])
    if (obj.min == obj.value[m]){beta.min = beta.m; obj.min = obj.min}
    if(abs(temp-obj.value[m])<=epsilon){break}
    if( (m>2) & (min(abs(obj.value[1:(m-1)]-obj.value[m]))<=1e-10) ){break}
  }
  if(iter.index==maxiter)
  {warning(paste("Algorithm does not converge in ", maxiter, " iterations! when n=", n, sep=""))}
  
  list(obj.values=obj.value[1:iter.index], beta=beta.min, niter=iter.index, obj = obj.min)
}


FCR.logistic.group<-function(X,Y,group,group.costs,budget,maxiter,epsilon,lambda,beta.initial){
  # X:  nxp design matrix
  # Y:  nx1 response variable
  # group:  px1 ordered vector (e.g., 1,1,2,3,4,... where 1 is the 1st group, 2 is the 2nd group...)
  # group.costs: gx1 vector. The j-th element is the cost to collect the j-th group
  # budget: the budget to build the model
  # maxiter: the maximum number of interations for the iterative dynamic programming
  # epsilon: the algorithm will terminate if the change of the objective function is 
  # smaller than epsilon
  n=nrow(X)
  p=ncol(X)
  n.group=length(unique(group))
  
  beta.m=beta.initial #rep(0,p)
  obj.value=rep(NA,maxiter)
  obj.min=mean(log(1+exp((-Y)*(X%*%beta.m))))
  obj.min = Inf
  iter.index=0
  W = 1/diag(t(X) %*% X) * 4
  d.m = diag(W) %*% t(X) %*% (Y*(exp((-Y)*(X%*%beta.m))/(1+exp((-Y)*(X%*%beta.m)))))
  
  for(m in 1:maxiter){
    iter.index = iter.index+1
    temp=mean(log(1+exp((-Y)*(X%*%beta.m))))
    
    tmp = beta.m + d.m
    
    l = 1
    a.vector = rep(0,n.group); a.vector[1] = tmp[1]
    for (j in 2:(length(group)) ){
      if (group[j] > group[j-1]) {
        l = l+1
      }
      a.vector[l] = a.vector[l] + tmp[j]
    }
    a.vector[is.na(a.vector)] = 0
    result = knapsack(group.costs,(a.vector^2),budget)
    ind.group = result$indices
    
    ind = rep(FALSE, p)
    for (j in 1:(length(ind.group)) ){
      ind = ind + (group == ind.group[j])
    }
    ind = (ind==1)
    
    beta.m <-rep(0,p); d.m <- rep(0,p)
    
    X.A.m = X[,ind]
    X.I.m = X[,-ind]
    
    beta.A.m = as.vector(glmnet(X.A.m,Y,family="binomial",alpha=1,lambda=lambda,intercept=F)$beta)
    beta.m[ind] = beta.A.m
    
    W.I.m = as.vector(1/diag(t(X.I.m) %*% X.I.m))
    if (length(W.I.m) == 1){
      d.I.m = t(X.I.m) * W.I.m
    } else {
      d.I.m = diag(c(W.I.m)) %*% t(X.I.m) 
      }
    d.I.m = d.I.m %*% (Y*exp(-Y*(X.A.m%*%beta.A.m)/(1+exp(-Y*(X.A.m%*%beta.A.m))) ) ) * 4
    d.m[-ind] = d.I.m
    obj.value[m]=mean(log(1+exp((-Y)*(X%*%beta.m))))
    
    obj.min = min(obj.min,obj.value[m])
    if (obj.min == obj.value[m]){beta.min = beta.m; obj.min = obj.min}
    if(abs(temp-obj.value[m])<=epsilon){break}
    if( (m>2) & (min(abs(obj.value[1:(m-1)]-obj.value[m]))<=1e-10) ){break}
  }
  if(iter.index==maxiter)
  {warning(paste("Algorithm does not converge in ", maxiter, " iterations! when n=", n, sep=""))}
  
  list(obj.values=obj.value[1:iter.index], beta=beta.min, niter=iter.index, obj = obj.min)
}

FCR.logistic.group.cv <-function(nfolds,X,Y,group,group.costs,budget,maxiter,epsilon,beta.initial){
  lasso.fit=glmnet(X,Y,family="binomial",alpha=1,nlambda=100,intercept=F)
  lambda.all=lasso.fit$lambda
  nlambda=length(lambda.all)
  nX=nrow(X)
  
  ## 5-fold CV to choose the best lambda
  cvfold=cvFolds(nX,K=nfolds,type="interleaved")
  valid.error=rep(NA,nlambda)
  for(i in 1:nlambda){
    Y.est=rep(NA,length(y))
    cv.score = c()
    for(nn in 1:nfolds){
      X.train=X[cvfold$subsets[cvfold$which!=nn],]
      X.valid=X[cvfold$subsets[cvfold$which==nn],]
      Y.train=Y[cvfold$subsets[cvfold$which!=nn]]
      Y.valid=Y[cvfold$subsets[cvfold$which==nn]]
      
      FCR.group.result= FCR.logistic.group(X.train,Y.train,group,group.costs,budget,1000,epsilon,
                                           lambda.all[i],rep(0,p))
      cv.score = c(cv.score, mean(log(1+exp((-Y.valid)*(X.valid%*%FCR.group.result$beta)))))
    }
    valid.error[i]=mean(cv.score)
  }
  
  opt.lambda=lambda.all[which.min(valid.error)[1]]
  return(opt.lambda)
}



HCR.logistic.group<-function(X,Y,group,group.costs,budget,maxiter,epsilon, lambda,alpha, beta.initial)
{
  ######### This function performs cost effective variable selection by iterative DP
  ######### (We consider the group cost rather than the separae cost for each variable)
  # X:  nxp design matrix
  # Y:  nx1 response variable
  # group:  px1 vector (e.g., 1,1,2,3,4,... where 1 is the 1st group, 2 is the 2nd group...)
  # group.costs: gx1 vector. The j-th element is the cost to collect the j-th group
  # budget: the budget to build the model
  # maxiter: the maximum number of interations for the iterative dynamic programming
  # epsilon: the algorithm will terminate if the change of the objective function is 
  # smaller than epsilon
  # lambda and alpha: parameters in the following elastic net penalty: 
  #                lambda*(alpha*|beta|+0.5*(1-alpha)*beta^2)
  # alpha=1 for Lasso and alpha=0 for the ridge regression
  n=nrow(X)
  p=ncol(X)
  n.group=length(group.costs)
  if (n>p){L=max(eigen(t(X)%*%X/n)$values) + 0.1}
  else {L=max(eigen(X%*%t(X)/n)$values) + 0.1}
  beta.m=beta.initial
  obj.value=rep(NA,maxiter)
  logl.value=rep(NA,maxiter)
  
  iter.index=0
  for(m in 1:maxiter)
  {
    iter.index=iter.index+1
    temp=mean(log(1+exp((-Y)*(X%*%beta.m))))+lambda*alpha*sum(abs(beta.m))
    +lambda*((1-alpha)/2)*sum(beta.m^2) 
    a.vector = as.vector(beta.m+t(X)%*%(Y*(exp((-Y)*(X%*%beta.m))/(1+exp((-Y)*(X%*%beta.m)))))/(n*L))
    a.vector.soft=(L/(L+lambda*(1-alpha)))*sign(a.vector-alpha*lambda/L)*
      pmax(abs(a.vector)-alpha*lambda/L,0)
    temp1=(a.vector^2)*(L^2)/(2*(L+lambda*(1-alpha)))
    temp2=-abs(a.vector)*alpha*lambda*L/(L+lambda*(1-alpha))
    temp3=(alpha^2)*(lambda^2)/(2*(L+lambda*(1-alpha)))
    obj.coefficients.seperate=(temp1+temp2+temp3)*(1+sign(abs(a.vector)-alpha*lambda/L))/2
    obj.coefficients=rep(NA,n.group)
    for(i in 1:n.group){obj.coefficients[i]=sum(obj.coefficients.seperate[group==i])}
    z.vector=gknapsack(group.costs,obj.coefficients,budget)
    beta.m=rep(NA,p)
    for(i in 1:n.group){beta.m[group==i]=a.vector.soft[group==i]*z.vector[i]}
    obj.value[m]=mean(log(1+exp((-Y)*(X%*%beta.m))))+lambda*alpha*sum(abs(beta.m))
    +lambda*((1-alpha)/2)*sum(beta.m^2)   
    logl.value[m] = mean(log(1+exp((-Y)*(X%*%beta.m))))
    if(abs(temp-obj.value[m])<=epsilon){break}
  }
  if(iter.index==maxiter)
  {warning(paste("Algorithm does not converge in ", maxiter, " iterations!", sep=""))}
  
  list(obj.values=obj.value[1:iter.index], logl.value = logl.value[1:iter.index], beta=beta.m, niter=iter.index)
}


Elastic.group.logistic.search<-function(X,Y,group,group.costs,budget,alpha)
{
  ######### This function searches the best elastic net solution that satisfies the budget 
  ######### constraint
  # X:  nxp design matrix
  # Y:  nx1 response variable
  # costs: px1 vector. The j-th element is the cost to collect the j-th predictor
  # budget: the budget to build the model
  # alpha: the parameter used in the elastic net
  # alpha=1 for Lasso and alpha=0 for the ridge regression
  n.group=length(group.costs)
  elastic.fit=glmnet(X,Y,family="binomial",alpha=1,intercept=F)
  beta.all=as.matrix(elastic.fit$beta)
  nfolds=5
  cvfold.out=cvFolds(nrow(X),K=nfolds,type="interleaved")
  cv.elastic=cv.glmnet(X,Y,family="binomial",alpha=alpha,type.measure="deviance",
                       lambda=elastic.fit$lambda,intercept=F,foldid=cvfold.out$which)
  model.cost=rep(NA,ncol(beta.all))
  for(i in 1:ncol(beta.all)){
    temp=beta.all[,i]
    group.indicator=rep(NA,n.group)
    for(j in 1:n.group){
      group.indicator[j]=as.integer(sum(temp[group==j]^2)!=0)
    }
    model.cost[i]=sum(group.indicator*group.costs)
  }
  avail.index=which(model.cost<=budget)
  opt.index=which.min(cv.elastic$cvm[avail.index])
  beta.est=as.vector(beta.all[,opt.index])
  beta.est
}


GroupLasso.logistic.search<-function(X,Y,group,group.costs,budget)
{
  n.group=length(group.costs)
  grouplasso.fit=gglasso(X,Y,group,loss="logit",nlambda=100,intercept=F)
  beta.all=as.matrix(grouplasso.fit$beta)
  nfolds=5
  cvfold.out=cvFolds(nrow(X),K=nfolds,type="interleaved")
  cv.grouplasso=cv.gglasso(X,Y,group,lambda=grouplasso.fit$lambda,loss="logit",
                         pred.loss="misclass",foldid=cvfold.out$which, intercept = F)
  model.cost=rep(NA,ncol(beta.all))
  for(i in 1:ncol(beta.all)){
    temp=beta.all[,i]
    group.indicator=rep(NA,n.group)
    for(j in 1:n.group){
      group.indicator[j]=as.integer(sum(temp[group==j]^2)!=0)
    }
    model.cost[i]=sum(group.indicator*group.costs)
  }
  avail.index=which(model.cost<=budget)
  opt.index=which.min(cv.grouplasso$cvm[avail.index])
  beta.est=as.vector(beta.all[,opt.index])
  beta.est
}

