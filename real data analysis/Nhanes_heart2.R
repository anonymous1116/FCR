source("../CCOfunctions.R")
library(MASS)
library(adagio)
library(glmnet)
library(data.table)

#diabetes
y<-fread("./realdata/heart_y.csv",header = T)[,2]
X<-fread("./realdata/heart_X.csv",header = T)
X<-X[,2:length(X)]
dim(X) 
costs<-fread("./realdata/heart_c.csv")
costs<-as.vector(t(costs)[,2])
costs <- costs[-1]
total.cost = sum(costs);total.cost

y<-as.matrix(y)
y = 2* y - 1

X<-as.matrix(X)
colnames(X)<-as.character(as.numeric(colnames(X))+1)

ind = c(seq(49,181,by=1),
        217, 221, 222, 234, 239, 242, 246, 247, 250, 254, 258, 262, 266, 270, 274, 278, 282) # Same value in this variables
#ind_row = (X[,28] == 1 | X[,29]==1)
#X = X[!ind_row,]
X = X[,-ind]
#y = y[!ind_row]

#X = X[,-ind]
costs = costs[-ind]
total.cost = sum(costs);total.cost


group = c(rep(1,3), #1-3
          rep(2,9), #4-12
          rep(3,2), #13-14
          rep(4,2), #15-16
          rep(5,1), # 17
          rep(6,5), # 18- 22
          rep(7,2), #23-24
          rep(8,4), #25-28
          rep(9,7), #29-35
          rep(10,8), #36-43
          rep(11,1), #44
          rep(12,1), # 45
          rep(13,3), #46-48
          rep(14,7), # 182-188
          rep(15,2), # 189-190
          rep(16,1), # 191
          rep(17,7), #192 - 198
          rep(18,8), # 199-206
          rep(19,3), # 207-209
          rep(20,4), # 210-213
          rep(21,4), # 214-218       -217
          rep(22,2), # 219- 222.     -221,222
          rep(23,1), # 223
          rep(24,7), # 224-230      
          rep(25,4), # 231-235      -234
          rep(26,3), # 236-239      -239
          rep(27,3), # 240-243      -242
          rep(28,2), # 244-247.     -246,247
          rep(29,3), # 248-251      -250
          rep(30,3),  # 252-255     -254
          rep(31,3), #256-259       -258
          rep(32,3), # 260-263      -262
          rep(33,3), # 264-267      -266
          rep(34,3), # 268-271      -270
          rep(35,3), # 272-275      -274
          rep(36,3), # 276-279      -278
          rep(37,3), # 280-283      -282
          rep(38,5), # 284-288
          rep(39,15), #289-303
          rep(40,5) #304-308
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
budget_seq = c(20, 40, 60, 80, 100)

# 2 methods: Ours,Ours2, HCR
# 3 measures: prediction, loglikelihood, time
results_heart = array(NA,dim=c(length(budget_seq),nsim,2,4))
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
    
    
    
    results_heart[j,sim,,] = rbind(c(pred.FCR, logl.FCR, auc.FCR, time.FCR), 
                                       c(pred.HCR, logl.HCR, auc.HCR, time.HCR))
    print(c("budget", budget_seq[j], "nsim: ", sim))
  }  
}

results = results_heart
tmp = list(group.costs, budget_seq, results, obj_results)
save(tmp,file="results_heart2.RData")

measures = c("pred", "logl", "AUC", "time")
names = c("FCR",  "HCR")
#load("results_heart.RData")
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
grid.arrange(grobs=plot_list,ncol=2, top=textGrob("heart"))

