source("../CCOfunctions.R")
library(MASS)
library(adagio)
library(glmnet)
library(data.table)
library(ggplot2)

y <- fread("./realdata/arthritis_y.csv", header = T)[,2]
X <- fread("./realdata/arthritis_X.csv", header = T)
X <- X[,2:length(X)]
colnames(X)<-as.character(as.numeric(colnames(X))+1)

dim(X) # 49509 * 48 
costs<-fread("./realdata/arthritis_c.csv")
costs<-as.vector(t(costs)[,2])
costs <- costs[-1]
#cost<-cost[-c(22,23,25)]

y<-as.matrix(y)
X<-as.matrix(X)


ind = c(30,35, 39,43,44)
#ind_row = (X[,30] == 1 | X[,35] ==1| X[,39] == 1| X[,43] == 1| X[,44] ==1)
#X = X[!ind_row,]
X = X[,-ind]
#y = y[!ind_row]
costs = costs[-ind]
summary(X)
y = 2*y -1

#X = scale(X)
group = c(1, # gender
          2, # age at the time of screening
          rep(3, 6), # age at the time of screening
          rep(4, 5), # race / ethnicity
          seq(5, 17, by =1), # 14- 26
          rep(18, 3),   # 27-29 Doctor told overweight (risk factor) 
          19, # sleep 31
          rep(20, 3), # 32-34
          21, # 36 smoking
          rep(22, 3), # 37-38, 40 # Blood relatives with arthritis
          rep(23, 2), # 41-42 # joint pain/aching/stiffness in past year
          rep(24, 3), # 45 - 47 # symptoms began only because of injury
          25) # how long experiencing pain
l = 1
n.group = length(unique(group))
group.costs = rep(0,n.group); group.costs[1] = costs[1]
for (j in 2:(length(group)) ){
  if (group[j] > group[j-1]) {
    l = l+1
  }
  group.costs[l] = costs[j]
}

sum(group.costs)


nsim = 100
maxiter = 1000
tol = 1e-5
budget_seq = c(20, 40, 60, 80, 100)

# 2 methods: Ours,Ours2, HCR
# 3 measures: prediction, loglikelihood, time
results_arthritis = array(NA,dim=c(length(budget_seq),nsim,2,4))
maxiter = 1000
obj_results = array(NA, dim= c(length(budget_seq), nsim, 2, maxiter+1))
p = dim(X)[2]
lambda = 0
set.seed(2726)
for (j in 1:length(budget_seq) ){
  budget = budget_seq[j]
  for (sim in 1:nsim){
    ind = sample(1:dim(X)[1], size = floor(dim(X)[1] * 0.8))
    train_y = y[ind]
    train_X = X[ind,]
    
    test_y = y[-ind]
    test_X = X[-ind,]
    obj.init<-mean(log(1+exp((-train_y)*(train_X%*%rep(0,dim(X)[2])))))
    
    ptm=proc.time()[3]
    FCR.result = FCR.logistic.group(train_X,train_y,group=group,group.costs,budget,maxiter=1000,tol,lambda= lambda, rep(0,p))
    time.FCR = proc.time()[3] - ptm
    logl.FCR = FCR.result$obj
    pred.FCR = mean(test_y != sign(test_X%*%FCR.result$beta))
    obj.FCR = c(obj.init,FCR.result$obj.values)
    obj_results[j,sim,1,1:length(obj.FCR)]<-obj.FCR
    a<-roc(test_y, test_X%*%FCR.result$beta, quiet = TRUE)
    auc.FCR = as.numeric(a$auc)
    
    
    ptm=proc.time()[3]
    HCR.result = HCR.logistic.group(train_X,train_y,group,group.costs,budget,maxiter=1000,tol,lambda,0,rep(0,p))
    time.HCR = proc.time()[3] - ptm
    logl.HCR = HCR.result$logl.value[HCR.result$niter]
    pred.HCR = mean(test_y != sign(test_X%*%HCR.result$beta))
    obj.HCR = c(obj.init,HCR.result$obj.values)
    obj_results[j,sim,2,1:length(obj.HCR)]<-obj.HCR
    a<-roc(test_y, test_X%*%HCR.result$beta, quiet = TRUE)
    auc.HCR = as.numeric(a$auc)
    
    
    
    results_arthritis[j,sim,,] = rbind(c(pred.FCR, logl.FCR, auc.FCR, time.FCR), 
                                           c(pred.HCR, logl.HCR, auc.HCR, time.HCR))
    print(c("budget", budget_seq[j], "nsim: ", sim))
  }  
}

results = results_arthritis
tmp = list(group.costs, budget_seq, results, obj_results)
save(tmp,file="results_arthritis2.RData")

measures = c("pred", "logl", "AUC", "time")
names = c("FCR",  "HCR")
#load("results_arthritis.RData")
#results = tmp[[3]]
#budget_seq = tmp[[2]]
df_n = prod(dim(results)[1:3])
dim_df = dim(results)
plot_list<-list()
for (j in 1:length(measures)){
  n_clus = rep(as.character(budget_seq), dim_df[2] * dim_df[3])
  methods = rep(names, each = dim_df[1] * dim_df[2])
  value = as.vector(results[,,,j])
  
  df = data.frame(n_clus = n_clus, methods = methods, value = value)
  # I reorder the groups order : I change the order of the factor data$names
  #df$n_clus <- factor(df$n_clus , levels=c("100", "300", "500", "700", "900"))
  df$n_clus <- factor(df$n_clus , levels=as.character(budget_seq))
  
  plot_list[[j]]<-ggplot(df, aes(x = n_clus, y = value, color = methods)) +  # ggplot function
    geom_boxplot() + labs(x = "C") + labs(title = measures[j])
  
}

library(grid)
grid.arrange(grobs=plot_list,ncol=2, top=textGrob("arthritis"))

