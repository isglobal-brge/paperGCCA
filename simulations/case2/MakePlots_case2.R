rm(list=ls())


library(lattice)

load("sim2res.rda")

esc <- expand.grid(difer=c(0, 0.1, 0.2, 0.3, 0.4, 0.5), prop.miss=c(0.1,0.3))


pow <- sapply(1:length(res), function(i) rowMeans(res[[i]]<0.05))
taula <- cbind(t(pow), esc)
colnames(taula)[1:5] <- c("true","full","impute","complete","mgcca")

               
ccc <- c("red","blue","darkgreen")
lll <- c(1,1,1)
list.xyplot<-list(title="",
  space="top", columns=3,
  text=list(c("MGCCA","IMPUTE","COMPLETE")),
            lines=list(col=ccc, lty=lll, lwd=3),
            cex.title=1, cex=1.1)
               
#png("plot1.png")
print(xyplot(mgcca + impute + complete ~ difer | prop.miss, data=taula,
             panel=function(...) {panel.xyplot(...); },
             xlab="prop missing",ylab="dist",type=c("p","l"),pch=19,cex=0.7,main="Power",col=ccc,lty=lll,lwd=1.7,
             strip=strip.custom(strip.names=TRUE,strip.levels=TRUE, sep=": ", var.name=c("prop miss")),
             par.strip.text = list(cex = 0.85), key=list.xyplot, ylim=c(0,max(taula[,c("mgcca","impute","complete")])),
             scales=list(x=list(tick.number=4, cex=0.8, rot=45), y= list(tick.number=10, cex=0.8))))
#dev.off()



#####

rm(list=ls())


library(ggplot2)

load("sim2res.rda")

esc <- expand.grid(difer=c(0, 0.25, 0.5), prop.miss=c(0.1,0.2,0.3))

alpha <- 0.01
taula <- rbind(
  data.frame("method"="MGCCA","pow"=sapply(1:length(res), function(i) mean(res[[i]][5,]<alpha)), esc),
  data.frame("method"="IMPUTE","pow"=sapply(1:length(res), function(i) mean(res[[i]][3,]<alpha)), esc),
  data.frame("method"="COMPLETE","pow"=sapply(1:length(res), function(i) mean(res[[i]][4,]<alpha)), esc)
)

taula$method <- factor(taula$method, levels=c("MGCCA","IMPUTE","COMPLETE"))

png("plot2.png", 0.5*21*480/7, 0.5*480)
ggplot(taula, aes(x = difer, y = pow, colour = method), size=0.5)+ 
  geom_line(aes(group = method)) + geom_point() + 
  facet_grid(cols=vars(prop.miss)) +
  ylab("Power") + 
  xlab("Difference") + 
  xlim(0, 0.5)
dev.off()

               