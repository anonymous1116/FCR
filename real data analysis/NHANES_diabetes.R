setwd("/home/hyun18/R/CCO")
source("./CCOfunctions.R")
library(MASS)
library(adagio)
library(glmnet)
library(gglasso)
library(cvTools)
library(data.table)
library(ggplot2)
library(pROC)
library(gridExtra)
library(grid)
library(lattice)
library(corrplot)

#diabetes
y<-fread("./realdata/diabetes_y.csv",header = T)[,2]
y[y==2] =1
X<-fread("./realdata/diabetes_X.csv",header = T)
X<-X[,2:length(X)]
#X<-X[,-c(7,22,23,25) ] # All components of each column has the same values
dim(X) 
costs<-fread("./realdata/diabetes_c.csv")
costs<-as.vector(t(costs)[,2])
costs <- costs[-1]
#cost<-cost[-c(7,22,23,25)]
total.cost=sum(costs);total.cost

y<-as.matrix(y)
X<-as.matrix(X)
y = 2* y - 1

group = c(1, # gender
          2, # age at the time of screening
          rep(3, 6), # age at the time of screening
          rep(4, 5), # race / ethnicity
          seq(5, 29, by =1),
          rep(30, 4), # smoking
          31,
          32,
          33)

category = c(paste0(rep("Demo_",13),seq(1,13), "(2)"), 
             paste0(c("Demo_"), 15, "(4)"),  
             paste0(c("Demo_"), 14, "(2)"),  
             paste0(rep("Exam_", 13),seq(1,13), "(5)"), 
             paste0(rep("Lab_", 4),seq(1,4), "(9)") , 
             paste0(rep("Ques_", 13),seq(1,13), "(4)")
)
colnames(X)<-category

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
sum(group.costs == 2)/ length(group.costs)
sum(group.costs == 4)/ length(group.costs)
sum(group.costs == 5)/ length(group.costs)
sum(group.costs == 9)/ length(group.costs)

nsim = 100
maxiter = 1000
tol = 1e-5
budget_seq = c(30, 50, 70, 90, 110)

p = dim(X)[2]
results_diabetes = array(NA,dim=c(length(budget_seq),nsim,3,4)) 
obj_results = array(NA, dim= c(length(budget_seq), nsim, 2, maxiter+1))
set.seed(2726)
lambda= 0
alpha = 1
for (sim in 1:nsim){
  for (j in 1:length(budget_seq) ){
    budget = budget_seq[j]
    ind = sample(1:dim(X)[1], size = floor(dim(X)[1] * 0.8))
    train_y = y[ind]
    train_X = X[ind,]
    
    test_y = y[-ind]
    test_X = X[-ind,]
    obj.init<-log(2)
    
    ptm=proc.time()[3]
    FCR.result = FCR.logistic.group(train_X,train_y,group=group,group.costs,budget,maxiter=1000,tol,lambda= lambda, rep(0,p))
    time.FCR = proc.time()[3] - ptm
    logl.FCR = FCR.result$obj
    pred.FCR = mean(test_y != sign(test_X%*%FCR.result$beta))
    obj.FCR = c(obj.init,FCR.result$obj.values)
    obj_results[j,sim,1,1:length(obj.FCR)]<-obj.FCR
    a<-roc(test_y, as.vector(test_X%*%FCR.result$beta), quiet = TRUE)
    auc.FCR = as.numeric(a$auc)
    
    ptm=proc.time()[3]
    HCR.result = HCR.logistic.group(train_X,train_y,group,group.costs,budget,maxiter=1000,tol,lambda,0,rep(0,p))
    time.HCR = proc.time()[3] - ptm
    logl.HCR = HCR.result$logl.value[HCR.result$niter]
    pred.HCR = mean(test_y != sign(test_X%*%HCR.result$beta))
    obj.HCR = c(obj.init,HCR.result$obj.values)
    obj_results[j,sim,2,1:length(obj.HCR)]<-obj.HCR
    a<-roc(test_y, as.vector(test_X%*%HCR.result$beta), quiet = TRUE)
    auc.HCR = as.numeric(a$auc)
    
    ptm=proc.time()[3]
    Lasso.result = Elastic.group.logistic.search(train_X, train_y, group, group.costs,budget, alpha)
    time.Lasso = proc.time()[3] - ptm
    logl.Lasso = mean(log(1+exp((-train_y)*(train_X%*%Lasso.result))))
    pred.Lasso = mean(test_y != sign(test_X%*%Lasso.result))
    a<-roc(test_y, as.vector(test_X%*%Lasso.result), quiet = TRUE)
    auc.Lasso = as.numeric(a$auc)
  
    results_diabetes[j,sim,,] = rbind(c(pred.FCR, logl.FCR, auc.FCR, time.FCR), 
                                      c(pred.HCR, logl.HCR, auc.HCR, time.HCR),
                                      c(pred.Lasso, logl.Lasso, auc.Lasso, time.Lasso)
    )
    print(c("budget", budget_seq[j], "nsim: ", sim))
  }  
}


results = results_diabetes
tmp = list(group.costs, budget_seq, results, obj_results)
save(tmp,file="revision_results/results_diabetes.RData")
load("revision_results/results_diabetes.RData")

results = results[,,1:3,]

apply(results[1,,,], c(2,3), mean) 
apply(results[2,,,], c(2,3), mean) 
apply(results[3,,,], c(2,3), mean) 

apply(results[1,,,], c(2,3), mean) 
apply(results[2,,,], c(2,3), mean) 
apply(results[3,,,], c(2,3), mean) 


measures = c("pred", "logl", "AUC", "time")
names = c("FCR",  "HCR", "Lasso")
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

library(gridExtra)
grid.arrange(grobs=plot_list,ncol=2, top=textGrob("Diabetes"))


{
  budgets = seq(10, 150, by =5)
  betas = matrix(NA, nrow = dim(X)[2], ncol = length(budgets))
  for (i in 1:length(budgets)){
    FCR.result=FCR.logistic.group(X,y,group=group,group.costs,budgets[i],maxiter=1000,epsilon=1e-5,lambda,rep(0,p))
    betas[,i] = FCR.result$beta  
  }
  
  
  colnames(betas) <- budgets
  rownames(betas) <- category
  betas
  
  betas_grid <- expand.grid(X=as.character(budgets), Y=category)
  betas_grid$Z = as.numeric(abs(as.vector(t(betas)))>0)
  
  tmp = sort(category)
  tmp = c(tmp[c(1,8:15)], tmp[c(2:7)], tmp[c(33,38:45,34:37)], tmp[c(16,21:28,17:20)], tmp[c(29:32)])
  
  ggplot(betas_grid, aes(X, Y, fill= Z)) + geom_tile() +
    scale_fill_gradient(low="yellow", high="red") + 
    theme(legend.position="none")+ theme(
      plot.title = element_text(color="black", size=14, face="bold"),
      axis.title.x = element_text(color="black", size=14, face="bold"),
      axis.title.y = element_text(color="black", size=14, face="bold")
    )+ scale_x_discrete(name="C",limits=as.character(budgets)) + 
    scale_y_discrete(name = "Variables", limits = tmp)
  
  sort(glm((y+1)/2 ~ X-1, family = "binomial")$coefficients, decreasing = T)
  summary(lm(y~X-1))
}
