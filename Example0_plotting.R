#install.packages("ggplot2")                                  
#install.packages("gridExtra")

library(ggplot2)    
library(forcats)
library(dplyr)
library(grid)
load("Example0_tmp.RData")
total.cost = tmp[[1]]
total.cost
nums = tmp[[2]]
results = tmp[[3]]
beta.result.all = tmp[[4]]
X_p0 = tmp[[5]]

measures = c("L2", "Pred", "FPR", "FNR", "time")
names = c("FCR", "HCR", "LASSO")

#measures = c("Pred", "time")

apply(results[1,,,],c(2,3),mean)
apply(results[1,,,],c(2,3),median)
apply(results[1,,,],c(2,3),sd)


df_n = prod(dim(results)[1:3])

dim_df = dim(results)
plot_list<-list()
for (j in 1:length(measures)){
  #for (j in 1:2){
  n_clus = rep(as.character(nums), dim_df[2] * dim_df[3])
  methods = rep(names, each = dim_df[1] * dim_df[2])
  value = as.vector(results[,,,j])
  
  df = data.frame(n_clus = n_clus, methods = methods, value = value)
  #df$n_clus <- factor(df$n_clus , levels=c("100", "300", "500", "700", "900"))
  df$n_clus <- factor(df$n_clus , levels=as.character(nums))
  
  plot_list[[j]]<-ggplot(df, aes(x = n_clus, y = value, color = methods)) +  # ggplot function
    geom_boxplot() + labs(x = "n") + labs(title = measures[j])
  
}
library(ggplot2)
library(gridExtra)
library(patchwork)
grid.arrange(grobs=plot_list,ncol=2, top=textGrob("Example0"))


library(corrplot)
{
  # Step 1: Call the pdf command to start the plot
pdf(file = "./intro.pdf", width = 12, height = 5) 
layout(matrix(c(1,2,4,1,3, 5), 2, 3, byrow = TRUE),
       widths=c(1,0.8,0.8), heights=c(1,1))
mgp1 = c(1.3,0.3,0.5)
mar1 = c(2.2, 3, 2.5, 2.5)
mar2 = c(2.5, 3, 2.2, 2.5)
#mar3 = c(2.5,1,2.5,1)
mar3 = c(2.5,1,2.5,1)
corrplot(cor(X_p0))
par(mar = mar1)

plot(1:24, beta.result.all[1,1,1,1:24], col = "black", cex = 1, lwd = 4, type ="b", 
     ylim = c(-.2,2.3), xlim = c(0,25), ylab = expression(hat(beta)), 
        xlab = paste0(expression(beta)),  mgp=mgp1)
y1 = apply(beta.result.all[1,,4,1:24], 2, function(a) quantile(a, p=.95))
y2 = apply(beta.result.all[1,,4,1:24], 2, function(a) quantile(a, p=.05))
y3 = apply(beta.result.all[1,,4,1:24], 2, function(a) quantile(a, p=.5))

points(1:24, y1, col = "blue", cex = 1, lwd = 1, type = "l", lty = 3)
points(1:24, y2, col = "blue", cex = 1, lwd = 1, type = "l", lty = 3)
# Fill area between lines
polygon(c(1:24, rev(1:24)), c(y2, rev(y1)),
        col = "blue",density = 30, angle = 45)

points(1:24, y3, col = "blue", cex = 1.5, lwd = 1, type = "b", pch = 5)
legend("topright", legend = c("LASSO"), col = c("blue"), pch = 5, cex = 1,
       box.lty = 0)

par(mar = mar2)
plot(1:24, beta.result.all[1,1,1,1:24], col = "black", cex = 1, lwd = 4, type ="b", 
     ylim = c(-.2,2.3), xlim = c(0,25), ylab = expression(hat(beta)), xlab = paste0(expression(beta)),
     mgp=mgp1)
y1 = apply(beta.result.all[1,,2,1:24], 2, function(a) quantile(a, p=.95))
y2 = apply(beta.result.all[1,,2,1:24], 2, function(a) quantile(a, p=.05))
y3 = apply(beta.result.all[1,,2,1:24], 2, function(a) quantile(a, p=.5))

points(1:24, y1, col = "red", cex = 1, lwd = 1, type = "l", lty = 2)
points(1:24, y2, col = "red", cex = 1, lwd = 1, type = "l", lty = 2)
# Fill area between lines
polygon(c(1:24, rev(1:24)), c(y2, rev(y1)),
        col = "red",density = 30, angle = 45)
points(1:24, y3, col = "red", cex = 1.5, lwd = 1, type = "b", pch = 2)
legend("topright", legend = c("FCR"), col = c("red"), pch = 2, cex=1,
       box.lty=0)

par(mar = mar1)

plot(1:24, beta.result.all[1,1,1,1:24], col = "black", cex = 1, lwd = 4, type ="b", 
     ylim = c(-.2,2.3), xlim = c(0,25), ylab = expression(hat(beta)), xlab = paste0(expression(beta)),
     mgp=mgp1)
y1 = apply(beta.result.all[1,,3,1:24], 2, function(a) quantile(a, p=.95))
y2 = apply(beta.result.all[1,,3,1:24], 2, function(a) quantile(a, p=.05))
y3 = apply(beta.result.all[1,,3,1:24], 2, function(a) quantile(a, p=.5))


points(1:24, y1, col = "green", cex = 1, lwd = 1, type = "l", lty = 3)
points(1:24, y2, col = "green", cex = 1, lwd = 1, type = "l", lty = 3)
# Fill area between lines
polygon(c(1:24, rev(1:24)), c(y2, rev(y1)),
        col = "green",density = 30, angle = 45)
points(1:24, y3, col = "green", cex = 1.5, lwd = 1, type = "b", pch = 3)
legend("topright", legend = c("HCR"), col = c("green"), pch = 3, cex=1,
       box.lty=0)

dev.off() 
  }

#combined <- plot_list[[1]] + plot_list[[2]] + plot_list[[3]] & theme(legend.position = "bottom")
#combined + plot_layout(guides = "collect", ncol = 3)
