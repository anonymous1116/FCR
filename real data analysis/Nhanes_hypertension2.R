source("../CCOfunctions.R")
library(MASS)
library(adagio)
library(glmnet)
library(data.table)

library(pROC)

#Hypertension
y <- fread("./realdata/hypertension_y.csv", header = T)[,2]
X <- fread("./realdata/hypertension_X.csv", header = T)
X <- X[,2:length(X)]
colnames(X)<-as.character(as.numeric(colnames(X))+1)
dim(X)
costs<-fread("./realdata/hypertension_c.csv")
costs<-as.vector(t(costs)[,2])
costs <- costs[-1]

y<-as.matrix(y)
X<-as.matrix(X)
#result<-glm(y~X, family = binomial)

y = 2*y -1

ind = c(22,23,25) # Same value in this variables
ind_row = (X[,28] == 1 | X[,29]==1)
X = X[!ind_row,]
X = X[,-c(28,29)]
y = y[!ind_row]

X = X[,-ind]
costs = costs[-ind]
total.cost = sum(costs);total.cost


group = c(1, # gender -cate
          2, # age at the time of screening  -conti
          rep(3, 5), # race / ethnicity -cate
          seq(4, 18, by =1), 
          rep(19, 2), # smoking -cate
          20, 
          21)



category = c(paste0(rep("Demo_", 7),seq(1,7), "(2)"), 
             paste0(c("Demo"), 9, "(4)"), 
             paste0(c("Demo_"), 8, "(2)"),  
             paste0(c("Diet"), "(4)"), 
             paste0(rep("Exam_", 5),seq(1,5), "(5)"), 
             paste0(rep("Lab_", 4),seq(1,4), "(9)") , 
             paste0(rep("Ques_", 7),seq(1,7), "(4)")
)
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
budget_seq = c(10, 20, 30, 40, 50)
#budget_seq = c(20, 40, 60, 80, 100)

p = dim(X)[2]
results_hypertension2 = array(NA,dim=c(length(budget_seq),nsim,2,4)) 
obj_results = array(NA, dim= c(length(budget_seq), nsim, 2, maxiter+1))
set.seed(2726)
lambda=0
for (sim in 1:nsim){
  for (j in 1:length(budget_seq) ){
    budget = budget_seq[j]
    ind = sample(1:dim(X)[1], size = floor(dim(X)[1] * 0.8))
    train_y = y[ind]
    train_X = X[ind,]
    
    test_y = y[-ind]
    test_X = X[-ind,]
    obj.init<-mean((train_y)^2)
    
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
    
    
    results_hypertension2[j,sim,,] = rbind(c(pred.FCR, logl.FCR, auc.FCR, time.FCR), 
                                       c(pred.HCR, logl.HCR, auc.HCR, time.HCR))
    print(c("budget", budget_seq[j], "nsim: ", sim))
  }  
}


results = results_hypertension2
tmp = list(group.costs, budget_seq, results, obj_results)
save(tmp,file="results_hypertension2.RData")

apply(results[1,,,], c(2,3), mean) 
apply(results[2,,,], c(2,3), mean) 
apply(results[3,,,], c(2,3), mean) 

apply(results[1,,,], c(2,3), mean) 
apply(results[2,,,], c(2,3), mean) 
apply(results[3,,,], c(2,3), mean) 


#measures = c("pred", "logl",  "time")
measures = c("pred", "logl", "AUC", "time")

df_n = prod(dim(results)[1:3])
dim_df = dim(results)
plot_list<-list()
for (j in 1:length(measures)){
  n_clus = rep(as.character(budget_seq), dim_df[2] * dim_df[3])
  names = c("FCR", "HCR")
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
grid.arrange(grobs=plot_list,ncol=2, top=textGrob("hypertension"))

