library(ggplot2)
setwd("/home/hyun18/R/CCO/revision_results")

plot_list_diabetes<-list()
plot_list_hypertension<-list()
plot_list_arthritis<-list()
plot_list_heart<-list()

# diabetes
{ 
  a<-load("results_diabetes.RData")
  group.costs<- tmp[[1]]
  budget_seq = tmp[[2]]
  results= tmp[[3]][,,1:3,]
  sum(group.costs)
  measures = c("Pred Error", "-logl", "AUC",  "time(s)")
  names = c("FCR", "HCR", "Lasso")
  
  df_n = prod(dim(results)[1:3])
  dim_df = dim(results)
  
  print(apply(results[1,,,], c(2,3), mean))
  print(apply(results[2,,,], c(2,3), mean))
  print(apply(results[3,,,], c(2,3), mean))
  print(apply(results[4,,,], c(2,3), mean))
  print(apply(results[5,,,], c(2,3), mean))
  print(apply(results, c(3,4), mean))
  
  for (j in 1:length(measures)){
    n_clus = rep(as.character(budget_seq), dim_df[2] * dim_df[3])
    methods = rep(names, each = dim_df[1] * dim_df[2])
    value = as.vector(results[,,,j])
    
    df = data.frame(n_clus = n_clus, methods = methods, value = value)
    # I reorder the groups order : I change the order of the factor data$names
    #df$n_clus <- factor(df$n_clus , levels=c("100", "300", "500", "700", "900"))
    df$n_clus <- factor(df$n_clus , levels=as.character(budget_seq))
    
    p<- ggplot(df, aes(x = n_clus, y = value, color = methods)) +  # ggplot function
      geom_boxplot() + labs(x = "C" ) + 
      labs(y = measures[j]) + theme_bw()  + theme(plot.margin = margin(1,1,1,1)) +
      scale_color_brewer(palette="Dark2")
    if (j ==3 | j==4){
      p<- p+  ggtitle("Diabetes")
    }
    plot_list_diabetes[[j]]<- p
  }
}

# hypertension
{ 
  a<-load("results_hypertension.RData")
  group.costs<- tmp[[1]]
  budget_seq = tmp[[2]]
  results= tmp[[3]][,,1:3,]
  sum(group.costs)
  measures = c("Pred Error", "-logl", "AUC",  "time(s)")
  names = c("FCR", "HCR", "Lasso")
  df_n = prod(dim(results)[1:3])
  dim_df = dim(results)
  
  print(apply(results[1,,,], c(2,3), mean))
  print(apply(results[2,,,], c(2,3), mean))
  print(apply(results[3,,,], c(2,3), mean))
  print(apply(results[4,,,], c(2,3), mean))
  print(apply(results[5,,,], c(2,3), mean))
  print(apply(results, c(3,4), mean))
  
  for (j in 1:length(measures)){
    n_clus = rep(as.character(budget_seq), dim_df[2] * dim_df[3])
    methods = rep(names, each = dim_df[1] * dim_df[2])
    value = as.vector(results[,,,j])
    
    df = data.frame(n_clus = n_clus, methods = methods, value = value)
    # I reorder the groups order : I change the order of the factor data$names
    #df$n_clus <- factor(df$n_clus , levels=c("100", "300", "500", "700", "900"))
    df$n_clus <- factor(df$n_clus , levels=as.character(budget_seq))
    
    p<- ggplot(df, aes(x = n_clus, y = value, color = methods)) +  # ggplot function
      geom_boxplot() + labs(x = "C" ) + 
      labs(y = measures[j]) + theme_bw()  +
      #scale_fill_manual(values=c("blue","green"))
      scale_color_brewer(palette="Dark2")
    if (j ==3| j==4){
      p<- p+  ggtitle("Hypertension")
    }
    plot_list_hypertension[[j]]<- p
  }
}

# arthritis
{ 
  a<-load("results_arthritis.RData")
  group.costs<- tmp[[1]]
  budget_seq = tmp[[2]]
  results= tmp[[3]][,,1:3,]
  sum(group.costs)
  measures = c("Pred Error", "-logl", "AUC",  "time(s)")
  names = c("FCR", "HCR", "Lasso")
  
  df_n = prod(dim(results)[1:3])
  dim_df = dim(results)
  
  for (j in 1:length(measures)){
    n_clus = rep(as.character(budget_seq), dim_df[2] * dim_df[3])
    methods = rep(names, each = dim_df[1] * dim_df[2])
    value = as.vector(results[,,,j])
    
    df = data.frame(n_clus = n_clus, methods = methods, value = value)
    # I reorder the groups order : I change the order of the factor data$names
    #df$n_clus <- factor(df$n_clus , levels=c("100", "300", "500", "700", "900"))
    df$n_clus <- factor(df$n_clus , levels=as.character(budget_seq))
    
    p<- ggplot(df, aes(x = n_clus, y = value, color = methods)) +  # ggplot function
      geom_boxplot() + labs(x = "C" ) + 
      labs(y = measures[j]) + theme_bw()  +
      #scale_fill_manual(values=c("blue","green"))
      scale_color_brewer(palette="Dark2")
    if (j ==3| j==4){
      p<- p+  ggtitle("Arthritis")
    }
    plot_list_arthritis[[j]]<- p
  }
}

# heart
{ 
  a<-load("results_heart.RData")
  group.costs<- tmp[[1]]
  budget_seq = tmp[[2]]
  results= tmp[[3]][,,1:3,]
  sum(group.costs)
  measures = c("Pred Error", "-logl", "AUC",  "time(s)")
  names = c("FCR", "HCR", "Lasso")
  
  df_n = prod(dim(results)[1:3])
  dim_df = dim(results)
  
  print(apply(results[1,,,], c(2,3), mean))
  print(apply(results[2,,,], c(2,3), mean))
  print(apply(results[3,,,], c(2,3), mean))
  print(apply(results[4,,,], c(2,3), mean))
  print(apply(results[5,,,], c(2,3), mean))
  print(apply(results, c(3,4), mean))
  
  for (j in 1:length(measures)){
    n_clus = rep(as.character(budget_seq), dim_df[2] * dim_df[3])
    methods = rep(names, each = dim_df[1] * dim_df[2])
    value = as.vector(results[,,,j])
    
    df = data.frame(n_clus = n_clus, methods = methods, value = value)
    # I reorder the groups order : I change the order of the factor data$names
    #df$n_clus <- factor(df$n_clus , levels=c("100", "300", "500", "700", "900"))
    df$n_clus <- factor(df$n_clus , levels=as.character(budget_seq))
    
    p<- ggplot(df, aes(x = n_clus, y = value, color = methods)) +  # ggplot function
      geom_boxplot() + labs(x = "C" ) + 
      labs(y = measures[j]) + theme_bw()  +
      #scale_fill_manual(values=c("blue","green"))
      scale_color_brewer(palette="Dark2")
    if (j ==3| j==4){
      p<- p+  ggtitle("Heart")
    }
    plot_list_heart[[j]]<- p
  }
}


library(grid)
library(patchwork)

combined <- plot_list_diabetes[[3]] + plot_list_diabetes[[1]] + 
  plot_list_diabetes[[2]] +
  plot_list_hypertension[[3]] + plot_list_hypertension[[1]] + 
  plot_list_hypertension[[2]]  +
  plot_list_arthritis[[3]] + plot_list_arthritis[[1]] + 
  plot_list_arthritis[[2]]  +
  plot_list_heart[[3]] + plot_list_heart[[1]] + 
  plot_list_heart[[2]]  & theme(legend.position = "right")

combined_time <- plot_list_diabetes[[4]] + plot_list_hypertension[[4]] + 
  plot_list_arthritis[[4]] +  plot_list_heart[[4]]

pdf("RA.pdf", width = 12, height = 9)
print(combined+ plot_layout(ncol = 3, widths = c(4, 4, 4), heights = c(2,2,2,2),  guides = "collect"))
dev.off()

pdf("RA2.pdf",width = 12, height = 3)
print(combined_time + plot_layout(guides = "collect", ncol = 4)& theme(legend.position = "right"))
dev.off()

# diabetes obj trajectory
red_col = c("coral", "coral1", "coral2", "coral3", "coral4")
green_col = c("chartreuse", "chartreuse1", "chartreuse2", "chartreuse3", "chartreuse4")
{
  pdf("RA2_2.pdf",width = 10, height = 3)
  a<-load("results_diabetes.RData")
  a<-tmp[[4]]
  budgets<-tmp[[2]]
  a[,,,1]<-log(2)
  
  legend_FCR = paste0("FCR C=", budgets)
  legend_HCR = paste0("HCR C=", budgets)
  obj_trajectory<-a[1,1,2,]
  obj = obj_trajectory[!is.na(obj_trajectory)] 
  par(mar = c(4, 4, 2, 2))
  
  plot(log(1:length(obj), 10), obj, xlim = c(0,3) , 
       xlab= expression(log[10]~(iteration) ),
       ylim= c(0.35,log(2)), ylab = "-logl",
       type = "l", col = "red", lty = 3, main = "Loss trajectory")
  
  for(i in 1:(dim(a)[1])){
    for(j in 1:(dim(a)[2])){
      obj_trajectory<-a[i,j,2,]
      obj = obj_trajectory[!is.na(obj_trajectory)] 
      points(log(1:length(obj), 10), obj, type = "l", col = red_col[i], lty = 3, lwd = 2)
      obj_trajectory<-a[i,j,1,]
      obj = obj_trajectory[!is.na(obj_trajectory)] 
      points(log(1:length(obj), 10), obj, type = "l", col = green_col[i], lty = 2, lwd = 2)
    }
  }
  legend("topright", col = c(green_col, red_col), lty = c(rep(2,5),rep(3,5)) ,
         cex = .6, legend = c(legend_FCR, legend_HCR), lwd = rep(2,10))
  dev.off() 
}

# hypertension, arthritis, and hearts
{
  pdf("RA2_appendx.pdf",width = 10, height = 9)
  layout(matrix(c(1,2,3), ncol = 1, nrow = 3),
         heights=c(1,1,1))
  {
    a<-load("results_hypertension.RData")
    a<-tmp[[4]]
    budgets<-tmp[[2]]
    a[,,,1]<-log(2)
    
    legend_FCR = paste0("FCR C=", budgets)
    legend_HCR = paste0("HCR C=", budgets)
    obj_trajectory<-a[1,1,2,]
    obj = obj_trajectory[!is.na(obj_trajectory)] 
    par(mar = c(4, 4, 2, 2))
    
    plot(log(1:length(obj), 10), obj, xlim = c(0,3) , 
         xlab= expression(log[10]~(iteration) ),
         ylim= c(0.2,log(2)), ylab = "-logl",
         type = "l", col = red_col[1], lty = 3, main = "Hypertension Loss trajectory")
    
    for(i in 1:(dim(a)[1])){
      for(j in 1:(dim(a)[2])){
        obj_trajectory<-a[i,j,2,]
        obj = obj_trajectory[!is.na(obj_trajectory)] 
        points(log(1:length(obj), 10), obj, type = "l", col = red_col[i], lty = 3, lwd = 2)
        obj_trajectory<-a[i,j,1,]
        obj = obj_trajectory[!is.na(obj_trajectory)] 
        points(log(1:length(obj), 10), obj, type = "l", col = green_col[i], lty = 2, lwd = 2)
      }
    }
    legend("topright", col = c(green_col, red_col), lty = c(rep(2,5),rep(3,5)) ,
           cex = .6, legend = c(legend_FCR, legend_HCR), lwd = rep(2,10))
  }
  {
    a<-load("results_arthritis.RData")
    a<-tmp[[4]]
    budgets<-tmp[[2]]
    a[,,,1]<-log(2)
    
    legend_FCR = paste0("FCR C=", budgets)
    legend_HCR = paste0("HCR C=", budgets)
    obj_trajectory<-a[1,1,2,]
    obj = obj_trajectory[!is.na(obj_trajectory)] 
    par(mar = c(4, 4, 2, 2))
    
    plot(log(1:length(obj), 10), obj, xlim = c(0,3) , 
         xlab= expression(log[10]~(iteration) ),
         ylim= c(0.35,log(2)), ylab = "-logl",
         type = "l", col = red_col[1], lty = 3, main = "Arthritis Loss trajectory")
    
    for(i in 1:(dim(a)[1])){
      for(j in 1:(dim(a)[2])){
        obj_trajectory<-a[i,j,2,]
        obj = obj_trajectory[!is.na(obj_trajectory)] 
        points(log(1:length(obj), 10), obj, type = "l", col = red_col[i], lty = 3, lwd = 2)
        obj_trajectory<-a[i,j,1,]
        obj = obj_trajectory[!is.na(obj_trajectory)] 
        points(log(1:length(obj), 10), obj, type = "l", col = green_col[i], lty = 2, lwd = 2)
      }
    }
    legend("topright", col = c(green_col, red_col), lty = c(rep(2,5),rep(3,5)) ,
           cex = .6, legend = c(legend_FCR, legend_HCR), lwd = rep(2,10))
  }
  {
    a<-load("results_heart.RData")
    a<-tmp[[4]]
    budgets<-tmp[[2]]
    a[,,,1]<-log(2)
    
    legend_FCR = paste0("FCR C=", budgets)
    legend_HCR = paste0("HCR C=", budgets)
    obj_trajectory<-a[1,1,2,]
    obj = obj_trajectory[!is.na(obj_trajectory)] 
    par(mar = c(4, 4, 2, 2))
    
    plot(log(1:length(obj), 10), obj, xlim = c(0,3) , 
         xlab= expression(log[10]~(iteration) ),
         ylim= c(0.2,log(2)), ylab = "-logl",
         type = "l", col = red_col[1], lty = 3, main = "Heart Loss trajectory")
    
    for(i in 1:(dim(a)[1])){
      for(j in 1:(dim(a)[2])){
        obj_trajectory<-a[i,j,2,]
        obj = obj_trajectory[!is.na(obj_trajectory)] 
        points(log(1:length(obj), 10), obj, type = "l", col = red_col[i], lty = 3, lwd = 2)
        obj_trajectory<-a[i,j,1,]
        obj = obj_trajectory[!is.na(obj_trajectory)] 
        points(log(1:length(obj), 10), obj, type = "l", col = green_col[i], lty = 2, lwd = 2)
      }
    }
    legend("topright", col = c(green_col, red_col), lty = c(rep(2,5),rep(3,5)) ,
           cex = .6, legend = c(legend_FCR, legend_HCR), lwd = rep(2,10))
  }
  
  
  dev.off() 
  
}
