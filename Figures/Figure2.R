rm(list=ls())
## First panel:

pdf("Panel2a.pdf", width=12, height=6)
library(Hmisc)
library(survival)
load(file="../Models/INTCLUST_AllProbs.RData")
r <- mapply(function(x,y) x[,-1] + y[,-1], x=psd, y=psdc)
r <- mapply(function(x,y) x + y[,-1], x=r, y=psdo)
r <- mapply(function(x,y) x + y[,-1], x=r, y=psld)
r <- mapply(function(x,y) x + y[,-1], x=r, y=psldc)
r <- mapply(function(x,y) x + y[,-1], x=r, y=psldo)
r <- mapply(function(x,y) x + y[,-1], x=r, y=psl)
r <- mapply(function(x,y) x + y[,-1], x=r, y=pslc)
r <- mapply(function(x,y) x + y[,-1], x=r, y=pslo)
r <- mapply(function(x,y) x + y[,-1], x=r, y=psc)

r0 <- do.call("cbind", r)
rownames(r0) <- psc[[1]][,1]

load(file="../Models/INTCLUST_AllProbs_FROM5.RData")
library(survival)

r <- mapply(function(x,y) x[,-1] + y[,-1], x=psd, y=psdc)
r <- mapply(function(x,y) x + y[,-1], x=r, y=psdo)
r <- mapply(function(x,y) x + y[,-1], x=r, y=psld)
r <- mapply(function(x,y) x + y[,-1], x=r, y=psldc)
r <- mapply(function(x,y) x + y[,-1], x=r, y=psldo)
r <- mapply(function(x,y) x + y[,-1], x=r, y=psl)
r <- mapply(function(x,y) x + y[,-1], x=r, y=pslc)
r <- mapply(function(x,y) x + y[,-1], x=r, y=pslo)
r <- mapply(function(x,y) x + y[,-1], x=r, y=psc)

r5 <- do.call("cbind", r)
rownames(r5) <- psc[[1]][,1]
r5 <- r5[c("6", "10", "15", "20"),]

load(file="../Models/INTCLUST_AllProbs_FROM10.RData")
library(survival)

r <- mapply(function(x,y) x[,-1] + y[,-1], x=psd, y=psdc)
r <- mapply(function(x,y) x + y[,-1], x=r, y=psdo)
r <- mapply(function(x,y) x + y[,-1], x=r, y=psld)
r <- mapply(function(x,y) x + y[,-1], x=r, y=psldc)
r <- mapply(function(x,y) x + y[,-1], x=r, y=psldo)
r <- mapply(function(x,y) x + y[,-1], x=r, y=psl)
r <- mapply(function(x,y) x + y[,-1], x=r, y=pslc)
r <- mapply(function(x,y) x + y[,-1], x=r, y=pslo)
r <- mapply(function(x,y) x + y[,-1], x=r, y=psc)

r10 <- do.call("cbind", r)
rownames(r10) <- psc[[1]][,1]


Clinical <- read.table("../../TableS6.txt", header=T, sep="\t", stringsAsFactors=F)
Clinical$iC10 <- factor(Clinical$iC10, levels=c(1:3, "4ER+", "4ER-", 5:10))
Clinical <- Clinical[which(Clinical$METABRIC.ID %in% colnames(r0)),]
Clinical <- Clinical[match(colnames(r0), Clinical$METABRIC.ID),]
mean(Clinical$METABRIC.ID==colnames(r0))



r0.INTCLUST.mean <- t(apply(r0, 1, function(y) tapply(y, Clinical$iC10, mean)))
r0.INTCLUST.SE <- t(apply(r0, 1, function(y) tapply(y, Clinical$iC10, function(z) sd(z)/sqrt(length(z)))))
r5.INTCLUST.mean <- t(apply(r5, 1, function(y) tapply(y, Clinical$iC10, mean)))
r5.INTCLUST.SE <- t(apply(r5, 1, function(y) tapply(y, Clinical$iC10, function(z) sd(z)/sqrt(length(z)))))
r10.INTCLUST.mean <- t(apply(r10, 1, function(y) tapply(y, Clinical$iC10, mean)))
r10.INTCLUST.SE <- t(apply(r10, 1, function(y) tapply(y, Clinical$iC10, function(z) sd(z)/sqrt(length(z)))))


coliClust <- c('#FF5500', '#00EE76', '#CD3278','#00C5CD', '#B5D0D2', '#8B0000',
               '#FFFF40', '#0000CD', '#FFAA00', '#EE82EE', '#7D26CD')
ss <- paste0("n=", sapply(r, function(x) ncol(x)))
res <- data.frame(X=c(r0.INTCLUST.mean, r5.INTCLUST.mean, r10.INTCLUST.mean), Year=c(rep(rownames(r0.INTCLUST.mean), 11),
                                                                                  rep(rownames(r5.INTCLUST.mean), 11),
                                                                                  rep(rownames(r10.INTCLUST.mean), 11)),
                  Predictions=rep(c("After Surgery", "5 years disease-free",
                      "10 years disease-free"), c(11*nrow(r0.INTCLUST.mean), 11*nrow(r5.INTCLUST.mean), 11*nrow(r10.INTCLUST.mean))),
                  IntClust=c(rep(colnames(r0.INTCLUST.mean), rep(nrow(r0.INTCLUST.mean), 11)),
                      rep(colnames(r5.INTCLUST.mean), rep(nrow(r5.INTCLUST.mean), 11)),
                      rep(colnames(r10.INTCLUST.mean), rep(nrow(r10.INTCLUST.mean), 11))),
                  li=c(r0.INTCLUST.mean, r5.INTCLUST.mean, r10.INTCLUST.mean) - 1.96*c(r0.INTCLUST.SE, r5.INTCLUST.SE, r10.INTCLUST.SE),
                  ui=c(r0.INTCLUST.mean, r5.INTCLUST.mean, r10.INTCLUST.mean) + 1.96*c(r0.INTCLUST.SE, r5.INTCLUST.SE, r10.INTCLUST.SE))

res$Predictions <- factor(res$Predictions, levels=c("After Surgery",
                                               "5 years disease-free",
                                               "10 years disease-free"))

res$IntClust <- factor(res$IntClust, levels=c(1:3, "4ER+", "4ER-", 5:10))
names(ss) <- levels(res$IntClust)
names(coliClust) <- levels(res$IntClust)
res$pch <- 17
res$pch[which(res$Predictions=="5 years disease-free")] <- 18
res$pch[which(res$Predictions=="10 years disease-free")] <- 19
res$colr <- "black"
res$colr[which(res$Predictions=="5 years disease-free")] <- "darkred"
res$colr[which(res$Predictions=="10 years disease-free")] <- "olivedrab"
res$IntClust <- factor(res$IntClust, levels=c(3, 7, 8, "4ER+", 10, "4ER-", 1,6,9,2,5))
coliClust <- coliClust[levels(res$IntClust)]
ss <- ss[levels(res$IntClust)]

res$Year <- as.numeric(as.character(res$Year))
pp <- xyplot(X ~ Year| IntClust, groups=Predictions, data=res,
             layout=c(11, 1, 1), ylab="Probability of Relapse", xlab="Years after surgery", ylim=c(0, 0.8),
             scales=list(x=list(at=c(2, 5, 10, 15, 20))),
             par.settings = list(superpose.symbol = list(pch=c(17, 18, 19),
                                     col=c("black", "darkred", "olivedrab")),
                 superpose.line = list(
                     col = rep(c("black", "darkred", "olivedrab"), 4),
                     lwd = c(1.5, 1.5, 1.5),
                     lty = c(1, 1, 1))
                                 ),
             panel=function(x,y,col, pch,...,subscripts) {
                 panel.fill(col=adjustcolor(coliClust[panel.number()],
                                alpha=0.4))
                 panel.xyplot(x,y, lwd=1.5, col=res$col[subscripts], pch=res$pch[subscripts])
                 panel.superpose(x,y,subscripts,...)
                 panel.arrows(x, y, x, res$ui[subscripts], length = 0.05, lwd=0.5,
                      angle = 90, col=res$col[subscripts])
                 panel.arrows(x, y, x, res$li[subscripts], length = 0.05, lwd=0.5,
                      angle = 90, col=res$col[subscripts])
                 panel.text(10, 0.7, labels=ss[panel.number()], outer=T, cex=0.75)
             },panel.groups=function(x, y, ..., subscripts) {
                 panel.lines(approx(x,y)$x, approx(x,y)$y, lty=2,
                             col=res$col[subscripts])
                 }, auto.key=T)
print(pp)

dev.off()




## Second panel:


rm(list=ls())

## Get CN from CBioportal, for example
CN <- read.table("CURTIS_data_CNA.txt", header=T, sep="\t")
Clinical <- read.table("../../TableS6.txt", header=T, sep="\t")
genelist <- list()
genelist[[1]] <- c("HSF5", "PPM1E", "PRR11", "DHX40", "TUBD1", "RPS6KB1", "CA4", "C17orf64", "BCAS3", "TBX2", "BRIP1", "TBC1D3P2")
genelist[[2]] <- c("FGF3", "CCND1", "CTTN", "FLJ42102", "CLPB", "P2RY2", "UCP2", "CHRDL2",
"MAP6", "EMSY", "OMP", "PAK1", "RSF1", "NARS2")
genelist[[3]] <- c("ZNF703", "EIF4EBP1", "LETM2", "STAR", "FGFR1")
genelist[[4]] <- c("FBXO32", "TMWM65", "SQLE", "LINC00861", "PCAT1", "MYC", "LINC00977", "MIR5194", "ADCY8")

pdf("Panel2b.pdf", width=12, height=4)
layout(rbind(rep(5, 4), c(1, 2, 3, 4)), height=c(0.2, 0.8))
par(mar=c(7, 1, 3, 3))
groups <- c("1", "2", "6", "9")
genes <- c(4, 5, 2, 3)
nams <- c("17q23Amp", "11q13Amp", "8p12Amp", "8q24Amp")
coliClust <- c('#FF5500', '#00EE76', '#CD3278','#00C5CD', '#B5D0D2', '#8B0000',
               '#FFFF40', '#0000CD', '#FFAA00', '#EE82EE', '#7D26CD')
names(coliClust) <- c(1:3, "4ER+", "4ER-", 5:10)
for (i in c(1, 3, 4, 2)) {
    X <- CN[which(CN[,1] %in% genelist[[i]]),]
    found <- genelist[[i]][which(genelist[[i]] %in% X[,1])]
    X <- X[match(found, X[,1]),]
    rownames(X) <- X[,1]
    X <- X[,-c(1:2)]
    X <- t(X)
    X <- data.frame(METABRIC.ID=sub(".", "-", rownames(X), fixed=T), X)
    X[,-1] <- apply(X[,-1], 2, function(x) factor(x, levels=seq(-2, 2, 1), labels=c("HOMD", "HETD", "NEUT", "GAIN", "AMP")))
    X[,-1] <- apply(X[,-1], 2, function(x) as.character(x))
    X <- data.frame(rep(X[,1], ncol(X)-1), Gene=rep(colnames(X)[-1], rep(1980, ncol(X)-1)), do.call("c", X[,-1]))
    colnames(X)[c(1,3)] <- c("METABRIC.ID", "CN")
    X$CN <- factor(X$CN)
    X$Gene <- factor(X$Gene, levels=found)
    levels(X$CN) <- list("YES"=c("GAIN", "AMP"), "NO"=c("HOMD", "HETD", "NEUT"))
    X <- merge(X, Clinical[,c('METABRIC.ID', 'iC10')])
    X$iC10 <- factor(X$iC10, levels=c(1:3, "4ER+", "4ER-", 5:10))
    tmp <- t(prop.table(table(X$iC10, X$Gene, X$CN)[groups[i],,], 1))
    tt <- barplot(tmp, las=2, beside=F, col=c(coliClust[groups[i]], adjustcolor(coliClust[groups[i]], 0.1)), axes=F, cex.names=1.3, font=3)
    mtext(paste("IC", groups[i], sep=""), side = 3, line = 1.3)
    mtext(nams[i], side = 3, line = 0, cex=1.1)
    if (i %in% c(3, 4, 2)) axis(2, las=2)
    if (i %in% c(1, 3, 4)) axis(4, labels=F)
    if (i==1) mtext("*", side=1, line=0.75, at=tt[c(6)], cex=1.5)
    if (i==2) mtext("*", side=1, line=0.75, at=tt[c(1, 2, 10, 11, 12)], cex=1.5)
    if (i==3) mtext("*", side=1, line=0.75, at=tt[c(1, 5)], cex=1.5)
    if (i==4) mtext("*", side=1, line=0.75, at=tt[c(5)], cex=1.5)
}

par(mar=c(2, 1, 2, 3))
Clinical$iC10 <- factor(Clinical$iC10, levels=c(1:3, "4ER+", "4ER-", 5:10))
res <- prop.table(table(Clinical$iC10))
res <- matrix(res, ncol=1)
perms <- c(3, 8, 9, 4, 11, 5, 1, 7, 10, 2, 6)
res <- res[perms,,drop=F]
barplot(res, col=coliClust[perms], horiz=T, beside=F, axes=F)
start <- c(0, cumsum(res))[1:11]
end <- cumsum(res)[1:11]

cols <- c("grey", "grey", "black", "black", "grey", "black", "grey", "black", "black", "black", "grey")

for (i in 1:11) {
    mtext(paste(round(res[i]*100, 0), "%", sep=""), side=1, at=start[i] + (end[i]-start[i])/2, line=1)
    text(x=start[i] + (end[i]-start[i])/2, y=0.7, labels=paste("IC", levels(Clinical$iC10)[perms[i]], sep=""), col=cols[i])
}
dev.off()


## Third panel:
pdf("Panel2c.pdf", width=8, height=8)
res <- prop.table(table(Clinical$iC10, Clinical$ER.Status), 2)
res <- res[perms,2]
res <- res[c(7, 8, 9, 10)]
barplot(res, beside=T,
        ylab="ER+ proportion",
        names.arg=c("IC1", "IC6", "IC9", "IC2"), cex.names=1.8, cex.axis=1.8, cex.lab=1.8, col=coliClust[c(1, 7, 10, 2)])
dev.off()
