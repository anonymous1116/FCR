#install.packages("ggplot2")                                  
#install.packages("gridExtra")

library(ggplot2)    
library(forcats)
library(dplyr)
library(gridExtra)
#load("Example0_tmp.RData")
setwd("/Users/hyeongjin/Library/CloudStorage/OneDrive-purdue.edu/Research/CostConstrainedOptimization/R_codes/CCO")

load("simulation5.RData")
#load("COO2_tmp.RData")



total.cost = tmp[[1]]
total.cost
nums = tmp[[2]]
rhos = tmp[[3]]
results = tmp[[4]]
dim(results)

# n = 200

#mar1 = c(2, 3, 2, 1)
#par(mar = mar1)
{
  pdf(file = "./sim5.pdf", width = 12, height = 7.5) 
  #par(mfrow = c(2,3))
  layout(matrix(c(1,2,3,4,5,6, 7, 8, 9),3, 3, byrow = TRUE),
         widths=c(1,1,1), heights=c(1,1, 1))
  par(mar = c(2, 4.5, 2, 1.5))
  
  for(num in 1:length(nums)){
    results_tmp = results[num,,,,]
    mean_mat = apply(results_tmp[,1:2,,1], c(1,3), mean) #1: prediction
    upper_mat = apply(results_tmp[,1:2,,1], c(1,3), function(a) quantile(a, probs = .975)) #1: prediction
    lower_mat = apply(results_tmp[,1:2,,1], c(1,3), function(a) quantile(a, probs = .025)) #1: prediction
    
    plot(x=rhos,y =  mean_mat[,1], type="l", axes = F, ylim = c(min(lower_mat),max(upper_mat)), 
         xlim = c(min(rhos)-.1, max(rhos)+.1), lwd = 2, ylab ="Pred Error",
         main = paste0("n=", nums[num]))
    # x axis
    xaxes=c(0, rhos, rhos[length(rhos)]+20)
    axis(1, at= xaxes, labels=xaxes, col.axis="black", las=1)
    # y axis
    axis(2, col.axis="black", las=1)
    for(b in 1:length(rhos)){
      lines(x = rep(rhos[b],2), y = c(upper_mat[b,1],lower_mat[b,1]), lwd = 2)
    }
    lines(x=rhos, y = mean_mat[,2], lty = 2, col = 2, lwd = 2)
    for(b in 1:length(rhos)){
      lines(x = rep(rhos[b],2), y = c(upper_mat[b,2],lower_mat[b,2]), lty =2, col = "red", lwd = 2)
    }
    
    lines(x=rhos, y = mean_mat[,3], lty = 4, col = "blue", lwd = 2)
    for(b in 1:length(rhos)){
      lines(x = rep(rhos[b],2), y = c(upper_mat[b,3],lower_mat[b,3]), lty=4, col = "blue", lwd = 2)
    }
    if(num == 1){
      legend("bottomleft", legend=c("FCR", "HCR", "Lasso"), lty=c(1,2,4), 
             lwd = c(1,1,1), cex = .7, col = c("black", "red", "blue"))
      #legend("bottomleft", legend=c("FCR", "HCR"), lty=c(1,2), 
      #       lwd = c(1,1), cex = .7, col = c("black", "red"))
      
    }
  }
  par(mar = c(2, 4.5, 2, 1.5))
  for(num in 1:length(nums)){
    results_tmp = results[num,,,,]
    mean_mat = apply(results_tmp[,,,2], c(1,3), mean) #2: FPR
    upper_mat = apply(results_tmp[,,,2], c(1,3), function(a) quantile(a, probs = .975)) #4: FPR
    lower_mat = apply(results_tmp[,,,2], c(1,3), function(a) quantile(a, probs = .025)) #4: FPR
    
    plot(x=rhos,y =  mean_mat[,1], type="l", axes = F, ylim = c(min(lower_mat),max(upper_mat)), 
         xlim = c(min(rhos)-.1, max(rhos)+.1), lwd = 2, ylab ="AUC",
         main = paste0("n=", nums[num]))
    # x axis
    xaxes=c(0, rhos, rhos[length(rhos)]+.1)
    axis(1, at= xaxes, labels=xaxes, col.axis="black", las=1)
    # y axis
    axis(2, col.axis="black", las=1)
    for(b in 1:length(rhos)){
      lines(x = rep(rhos[b],2), y = c(upper_mat[b,1],lower_mat[b,1]), lwd = 2)
    }
    lines(x=rhos, y = mean_mat[,2], lty = 2, col = 2, lwd = 2)
    for(b in 1:length(rhos)){
      lines(x = rep(rhos[b],2), y = c(upper_mat[b,2],lower_mat[b,2]), lty =2, col = "red", lwd = 2)
    }
    
    lines(x=rhos, y = mean_mat[,3], lty = 4, col = "blue", lwd = 2)
    for(b in 1:length(rhos)){
      lines(x = rep(rhos[b],2), y = c(upper_mat[b,3],lower_mat[b,3]), lty=4, col = "blue", lwd = 2)
    }
    if(num == 1){
      legend("topleft", legend=c("FCR", "HCR", "Lasso"), lty=c(1,2,4), 
             lwd = c(1,1,1), cex = .7, col = c("black", "red", "blue"))
    }
  }
  
  par(mar = c(2, 4.5, 2, 1.5))
  for(num in 1:length(nums)){
    results_tmp = results[num,,,,]
    mean_mat = apply(log(results_tmp[,,,5]), c(1,3), mean) #2: FPR
    upper_mat = apply(log(results_tmp[,,,5]), c(1,3), function(a) quantile(a, probs = .975)) #4: FPR
    lower_mat = apply(log(results_tmp[,,,5]), c(1,3), function(a) quantile(a, probs = .025)) #4: FPR
    
    plot(x=rhos,y =  mean_mat[,1], type="l", axes = F, ylim = c(min(lower_mat),max(upper_mat)), 
         xlim = c(min(rhos)-.1, max(rhos)+.1), lwd = 2, ylab ="log(time)",
         main = paste0("n=", nums[num]))
    # x axis
    xaxes=c(0, rhos, rhos[length(rhos)]+.1)
    axis(1, at= xaxes, labels=xaxes, col.axis="black", las=1)
    # y axis
    axis(2, col.axis="black", las=1)
    for(b in 1:length(rhos)){
      lines(x = rep(rhos[b],2), y = c(upper_mat[b,1],lower_mat[b,1]), lwd = 2)
    }
    lines(x=rhos, y = mean_mat[,2], lty = 2, col = 2, lwd = 2)
    for(b in 1:length(rhos)){
      lines(x = rep(rhos[b],2), y = c(upper_mat[b,2],lower_mat[b,2]), lty =2, col = "red", lwd = 2)
    }
    
    lines(x=rhos, y = mean_mat[,3], lty = 4, col = "blue", lwd = 2)
    for(b in 1:length(rhos)){
      lines(x = rep(rhos[b],2), y = c(upper_mat[b,3],lower_mat[b,3]), lty=4, col = "blue", lwd = 2)
    }
    if(num == 1){
      legend("topleft", legend=c("FCR", "HCR", "Lasso"), lty=c(1,2,4), 
             lwd = c(1,1,1), cex = .7, col = c("black", "red", "blue"))
    }
  }
  
  dev.off()
}
