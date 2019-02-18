rm(list=ls())

############################################
############################################
## Probs Cancer Death in 4ER- / 10
############################################
############################################


load(file="../Models/INTCLUST_AllProbs.RData")
library(survival)

dr <- mapply(function(x,y) x[,-1] + y[,-1], x=psd, y=psdc)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psdo)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psld)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psldc)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psldo)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psc)

r0 <- do.call("cbind", dr)
rownames(r0) <- psc[[1]][,1]

load(file="../Models/INTCLUST_AllProbs_FROM5.RData")
library(survival)

dr <- mapply(function(x,y) x[,-1] + y[,-1], x=psd, y=psdc)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psdo)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psld)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psldc)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psldo)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psc)

r5 <- do.call("cbind", dr)
rownames(r5) <- psc[[1]][,1]

load(file="../Models/INTCLUST_AllProbs_FROM10.RData")
library(survival)

dr <- mapply(function(x,y) x[,-1] + y[,-1], x=psd, y=psdc)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psdo)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psld)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psldc)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psldo)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psc)


r10 <- do.call("cbind", dr)
rownames(r10) <- psc[[1]][,1]


Clinical <- read.table("../../TableS6.txt", header=T, sep="\t", stringsAsFactors=F)
Clinical$iC10 <- factor(Clinical$iC10, levels=c(1:3, "4ER+", "4ER-", 5:10))
Clinical <- Clinical[which(Clinical$METABRIC.ID %in% colnames(r0)),]
Clinical <- Clinical[match(colnames(r0), Clinical$METABRIC.ID),]
mean(Clinical$METABRIC.ID==colnames(r0))


ICM.r0.INTCLUST.mean <- t(apply(r0, 1, function(y) tapply(y, Clinical$iC10, mean)))
ICM.r0.INTCLUST.SE <- t(apply(r0, 1, function(y) tapply(y, Clinical$iC10, function(z) sd(z)/sqrt(length(z)))))
lr.INTCLUST.size <- t(apply(r0, 1, function(y) tapply(y, Clinical$iC10, length)))
ss <- lr.INTCLUST.size[1,]
ss <- paste("(n=", ss, ")", sep="")
ICM.r5.INTCLUST.mean <- t(apply(r5, 1, function(y) tapply(y, Clinical$iC10, mean)))
ICM.r5.INTCLUST.SE <- t(apply(r5, 1, function(y) tapply(y, Clinical$iC10, function(z) sd(z)/sqrt(length(z)))))
ICM.r5.INTCLUST.mean <- ICM.r5.INTCLUST.mean[c("6", "10", "15", "20"),]
ICM.r5.INTCLUST.SE <- ICM.r5.INTCLUST.SE[c("6", "10", "15", "20"),]
ICM.r10.INTCLUST.mean <- t(apply(r10, 1, function(y) tapply(y, Clinical$iC10, mean)))
ICM.r10.INTCLUST.SE <- t(apply(r10, 1, function(y) tapply(y, Clinical$iC10, function(z) sd(z)/sqrt(length(z)))))

load(file="../Models/FourGroups_AllProbs.RData")
library(survival)

dr <- mapply(function(x,y) x[,-1] + y[,-1], x=psd, y=psdc)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psdo)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psld)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psldc)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psldo)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psc)

r0 <- do.call("cbind", dr)
rownames(r0) <- psc[[1]][,1]

load(file="../Models/FourGroups_AllProbs_FROM5.RData")
library(survival)

dr <- mapply(function(x,y) x[,-1] + y[,-1], x=psd, y=psdc)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psdo)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psld)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psldc)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psldo)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psc)



r5 <- do.call("cbind", dr)
rownames(r5) <- psc[[1]][,1]

load(file="../Models/FourGroups_AllProbs_FROM10.RData")
library(survival)

dr <- mapply(function(x,y) x[,-1] + y[,-1], x=psd, y=psdc)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psdo)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psld)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psldc)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psldo)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psc)

r10 <- do.call("cbind", dr)
rownames(r10) <- psc[[1]][,1]


Clinical <- read.table("../../TableS6.txt", header=T, sep="\t", stringsAsFactors=F)
Clinical$iC10 <- factor(Clinical$iC10, levels=c(1:3, "4ER+", "4ER-", 5:10))
Clinical <- Clinical[which(Clinical$METABRIC.ID %in% colnames(r0)),]
Clinical <- Clinical[match(colnames(r0), Clinical$METABRIC.ID),]
mean(Clinical$METABRIC.ID==colnames(r0))


FOURGROUPS.r0.INTCLUST.mean <- t(apply(r0, 1, function(y) tapply(y, Clinical$iC10, mean)))
FOURGROUPS.r0.INTCLUST.SE <- t(apply(r0, 1, function(y) tapply(y, Clinical$iC10, function(z) sd(z)/sqrt(length(z)))))
FOURGROUPS.r5.INTCLUST.mean <- t(apply(r5, 1, function(y) tapply(y, Clinical$iC10, mean)))
FOURGROUPS.r5.INTCLUST.SE <- t(apply(r5, 1, function(y) tapply(y, Clinical$iC10, function(z) sd(z)/sqrt(length(z)))))
FOURGROUPS.r10.INTCLUST.mean <- t(apply(r10, 1, function(y) tapply(y, Clinical$iC10, mean)))
FOURGROUPS.r10.INTCLUST.SE <- t(apply(r10, 1, function(y) tapply(y, Clinical$iC10, function(z) sd(z)/sqrt(length(z)))))



coliClust <- c('#FF5500', '#00EE76', '#CD3278','#00C5CD', '#B5D0D2', '#8B0000',
               '#FFFF40', '#0000CD', '#FFAA00', '#EE82EE', '#7D26CD')
res <- data.frame(X=c(ICM.r5.INTCLUST.mean, FOURGROUPS.r5.INTCLUST.mean), Year=c(rep(rownames(ICM.r5.INTCLUST.mean), 11),
                                                                                  rep(rownames(FOURGROUPS.r5.INTCLUST.mean), 11)),
                  Model=rep(c("IntClust Model", "IHC Model"), c(11*nrow(ICM.r5.INTCLUST.mean), 11*nrow(FOURGROUPS.r5.INTCLUST.mean))),
                  IntClust=c(rep(colnames(ICM.r5.INTCLUST.mean), rep(nrow(FOURGROUPS.r5.INTCLUST.mean), 11))),
                  li=c(ICM.r5.INTCLUST.mean, FOURGROUPS.r5.INTCLUST.mean) - 1.96*c(ICM.r5.INTCLUST.SE, FOURGROUPS.r5.INTCLUST.SE),
                  ui=c(ICM.r5.INTCLUST.mean, FOURGROUPS.r5.INTCLUST.mean) + 1.96*c(ICM.r5.INTCLUST.SE, FOURGROUPS.r5.INTCLUST.SE))

res$IntClust <- factor(res$IntClust, levels=c(1:3, "4ER+", "4ER-", 5:10))
names(coliClust) <- levels(res$IntClust)
res$pch <- 18
res$pch[which(res$Model=="IntClust Model")] <- 19
res$colr <- "black"
res$colr[which(res$Model=="IHC Model")] <- "olivedrab"
res$Year <- as.numeric(as.character(res$Year))

res <- subset(res, IntClust %in% c("10", "4ER-"))
ss <- ss[c(11,5)]
coliClust <- coliClust[c(5, 11)]
library(lattice)
pp <- xyplot(X ~ Year| IntClust, groups=Model, data=res,
             layout=c(2, 1, 1),
             ylab=list(label="Probability of distant relapse/cancer death", cex=1.2),
             xlab=list(label="Years after surgery", cex=1.2),
             ylim=c(0, 0.4), xlim=c(4, 22),
             scales=list(x=list(at=c(2, 5, 10, 15, 20),
                             cex=1.1), y=list(cex=1.1)),
             par.settings = list(superpose.symbol = list(pch=c(18, 19),
                                     col=c("olivedrab", "black")),
                 superpose.line = list(
                     col = rep(c("olivedrab", "black"), 4),
                     lwd = c(1.5, 1.5, 1.5),
                     lty = c(1, 2))
                                 ),
             panel=function(x,y,col, pch,...,subscripts) {
                 panel.fill(col=adjustcolor(coliClust[panel.number()],
                                alpha=0.4))
                 panel.xyplot(x,y, lwd=1.5, col=res$col[subscripts], pch=res$pch[subscripts])
                 panel.superpose(x,y,subscripts,...)
                 panel.text(12.5, 0.35, labels=ss[panel.number()], outer=T, cex=1.2)
                 panel.arrows(x, y, x, res$ui[subscripts], length = 0.05, lwd=0.5,
                      angle = 90, col=res$col[subscripts])
                 panel.arrows(x, y, x, res$li[subscripts], length = 0.05, lwd=0.5,
                      angle = 90, col=res$col[subscripts])
             },panel.groups=function(x, y, ..., subscripts) {
                 panel.lines(approx(x,y)$x, approx(x,y)$y, lty=c(1, 2),
                             col=res$col[subscripts])
                 }, auto.key=T)
print(pp)



