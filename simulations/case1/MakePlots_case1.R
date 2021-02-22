###########

rm(list=ls())


library(lattice)
library(ggplot2)

load("sim1res.rda")

esc <- expand.grid(vars=c(50, 100), prop.miss=c(0.1,0.2,0.3), sigma.noise=c(0.125, 0.250))

taula <- rbind(
  data.frame("method"="MGCCA","dist"=sapply(1:length(res), function(i) mean(res[[i]][1,])), esc),
  data.frame("method"="IMPUTE","dist"=sapply(1:length(res), function(i) mean(res[[i]][2,])), esc),
  data.frame("method"="COMPLETE","dist"=sapply(1:length(res), function(i) mean(res[[i]][3,])), esc)
)

png("plot1.png", 480*10/7, 480)
ggplot(data = taula, aes(x = factor(prop.miss), y = dist, colour = method)) +       
  geom_line(aes(group = method)) + geom_point() +
  facet_grid(rows= vars(vars), cols=vars(sigma.noise)) +
  xlab("proportions of missing individuals") +
  ylab("mean square-distance")
dev.off()

######

nsim <- ncol(res[[1]])
esc2 <- data.frame(esc=rep(1:nrow(esc), nsim), rep=rep(1:nsim, each=nrow(esc)), do.call("rbind", replicate(nsim, esc, simplify = FALSE)))
esc2 <- esc2[order(esc2$esc, esc2$rep),]

taula <- rbind(
  data.frame("method"="MGCCA","dist"=unlist(lapply(1:length(res), function(i) res[[i]][1,])), esc2),
  data.frame("method"="IMPUTE","dist"=unlist(lapply(1:length(res), function(i) res[[i]][2,])), esc2),
  data.frame("method"="COMPLETE","dist"=unlist(lapply(1:length(res), function(i) res[[i]][3,])), esc2)
)

png("plot2.png", 480*10/7, 480)
ggplot(taula, aes(x = factor(prop.miss), y = dist))+ 
  geom_boxplot(aes(fill = method), alpha = 0.5) + 
  facet_grid(rows= vars(vars), cols=vars(sigma.noise)) +
  xlab("proportions of missing individuals") +
  ylab("mean square-distance")
dev.off()

