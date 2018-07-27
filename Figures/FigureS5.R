rm(list=ls())

library(Hmisc)
load(file="../Models/IntClust_AllProbs.RData")
library(survival)

lr <- mapply(function(x,y) x[,-1] + y[,-1], x=plc, y=pld)
lr <- mapply(function(x,y) x + y[,-1], x=lr, y=pldc)

dr <- mapply(function(x,y) x[,-1] + y[,-1], x=psd, y=psdc)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psdo)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psld)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psldc)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psldo)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psc)



times <-  ps[[1]][,1]

lr.INTCLUST.mean <- sapply(lr, function(x) apply(x, 1, mean))
lr.INTCLUST.SE <- sapply(lr, function(x) apply(x, 1, function(y) sd(y)/sqrt(length(y))))
dr.INTCLUST.mean <- sapply(dr, function(x) apply(x, 1, mean))
dr.INTCLUST.SE <- sapply(dr, function(x) apply(x, 1, function(y) sd(y)/sqrt(length(y))))

all.lr <- do.call("cbind", lr)
all.dr <- do.call("cbind", dr)
Clinical <- read.table(file="~/Documents/Projects/Metastasis/Metastasis/Paper/Draft/MetabricClinicalMolecularDataset.txt", header=T, sep="\t", quote="", comment.char="", stringsAsFactors=FALSE)
Clinical$ER.Status <- factor(Clinical$ER.Status,
                             levels=c("neg", "pos"), labels=c("ER-", "ER+"))
Clinical$ER.Status[which(is.na(Clinical$ER.Status) & Clinical$ER.Expr=="+")] <- "ER+"
Clinical$ER.Status[which(is.na(Clinical$ER.Status) & Clinical$ER.Expr=="-")] <- "ER-"
Clinical$Group <- NA
Clinical$Group[which(Clinical$ER.Status=="ER-" & Clinical$Her2.Expr=="-")] <- "ER-/HER2-"
Clinical$Group[which(Clinical$ER.Status=="ER-" & Clinical$Her2.Expr=="+")] <- "ER-/HER2+"
Clinical$Group[which(Clinical$ER.Status=="ER+" & Clinical$Her2.Expr=="+")] <- "ER+/HER2+"
Clinical$Group[which(Clinical$ER.Status=="ER+" & Clinical$Her2.Expr=="-")] <- "ER+/HER2-"
Clinical$Group <- factor(Clinical$Group)
Clinical <- Clinical[which(Clinical$Group=="ER+/HER2-"),]
all.lr <- all.lr[,which(colnames(all.lr) %in% Clinical$METABRIC.ID)]
all.dr <- all.dr[,which(colnames(all.dr) %in% Clinical$METABRIC.ID)]
Clinical.LR <- Clinical[which(Clinical$METABRIC.ID %in% colnames(all.lr)),]
Clinical.LR <- Clinical.LR[match(colnames(all.lr), Clinical.LR$METABRIC.ID),]
Clinical.DR <- Clinical[which(Clinical$METABRIC.ID %in% colnames(all.dr)),]
Clinical.DR <- Clinical.DR[match(colnames(all.dr), Clinical.DR$METABRIC.ID),]




coliClust <- c('#FF5500', '#00EE76', '#CD3278','#00C5CD', '#B5D0D2', '#8B0000',
               '#FFFF40', '#0000CD', '#FFAA00', '#EE82EE', '#7D26CD')

ss <- paste0("n=", sapply(lr, function(x) ncol(x)-1))
TOT <- ncol(lr.INTCLUST.mean)
res <- data.frame(X=c(lr.INTCLUST.mean, dr.INTCLUST.mean), Year=rep(times, TOT*2), Relapse=rep(c("Loco-regional Relapse",
                                                                                      "Distant Relapse"), c(TOT*5, TOT*5)),
                  IntClust=rep(rep(colnames(lr.INTCLUST.mean), rep(5, TOT)), 2), li=c(lr.INTCLUST.mean, dr.INTCLUST.mean) - 1.96*c(lr.INTCLUST.SE, dr.INTCLUST.SE), ui=c(lr.INTCLUST.mean, dr.INTCLUST.mean) + 1.96*c(lr.INTCLUST.SE, dr.INTCLUST.SE))
res$IntClust <- factor(res$IntClust, levels=c(1:3, "4ER+", "4ER-", 5:10),
                       labels=c(1:3, "4ER+", "4ER-", 5:10))
res$pch <- 19
res$colr <- coliClust[as.numeric(res$IntClust)]

res <- subset(res, IntClust %in% c(3, 1, 2, 6, 9))
res$IntClust <- factor(res$IntClust, levels=c(3, 1, 6, 9, 2))
ICModel <- res

## Top Plot

pdf("FigureS5a.pdf", width=12, height=7)
library(Hmisc)
library(survival)


par(mfrow=c(1, 2))
sub.res <- subset(ICModel, Relapse=="Distant Relapse")
Y <- sub.res[order(sub.res$Year, sub.res$IntClust),]
Y$Dist <- 1:nrow(Y)
Y$Dist <- Y$Dist + rep(c(0, 2, 4, 6, 8), rep(5,5))
plot(X ~ Dist, data=Y, ylim=c(0,1), axes=F, col=Y$col, pch=Y$pch,
     ylab="Probability of having distant relapse", xlab="Years after surgery", main="a",
     cex.axis=1.5, cex.lab=1.5)
usr <- par('usr')
rect(Y$Dist[seq(from=1, by=5, length=5)]-0.85, usr[3],
     Y$Dist[seq(from=5, by=5, length=5)]+0.85, usr[4],
     col=adjustcolor('grey', alpha=0.4), border=NA, cex=1.5)
points(li ~ Dist, data=Y, pch="-", col=Y$col, cex=2)
points(ui ~ Dist, data=Y, pch="-", col=Y$col, cex=2)
segments(Y$Dist, Y$li, Y$Dist, Y$ui, col=Y$col, cex=2)
axis(2)
axis(1, at=unique(Y$Dist)[seq(from=3, by=5, length=5)] - 0.5, labels=c(2, 5, 10, 15, 20))
box()

legend("topleft", bty="n", lwd=2, pch=19, col=coliClust[c(1, 2, 3, 7, 10)],
       legend=paste("IC", c(1, 2, 3, 6, 9), sep=""), cex=1.2)


library(Hmisc)
load(file="../Models/IntClust_AllProbs.RData")
library(survival)

lr <- mapply(function(x,y) x[,-1] + y[,-1], x=plc, y=pld)
lr <- mapply(function(x,y) x + y[,-1], x=lr, y=pldc)

dr <- mapply(function(x,y) x[,-1] + y[,-1], x=psd, y=psdc)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psdo)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psld)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psldc)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psldo)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psc)

all.lr <- do.call("cbind", lr)
all.dr <- do.call("cbind", dr)


times <-  ps[[1]][,1]


Clinical <- read.table(file="../../TableS6.txt", header=T, sep="\t", quote="", comment.char="", stringsAsFactors=FALSE)
Clinical$ER.Status <- factor(Clinical$ER.Status,
                             levels=c("neg", "pos"), labels=c("ER-", "ER+"))
Clinical$ER.Status[which(is.na(Clinical$ER.Status) & Clinical$ER.Expr=="+")] <- "ER+"
Clinical$ER.Status[which(is.na(Clinical$ER.Status) & Clinical$ER.Expr=="-")] <- "ER-"
Clinical$Group <- NA
Clinical$Group[which(Clinical$ER.Status=="ER-" & Clinical$Her2.Expr=="-")] <- "ER-/HER2-"
Clinical$Group[which(Clinical$ER.Status=="ER-" & Clinical$Her2.Expr=="+")] <- "ER-/HER2+"
Clinical$Group[which(Clinical$ER.Status=="ER+" & Clinical$Her2.Expr=="+")] <- "ER+/HER2+"
Clinical$Group[which(Clinical$ER.Status=="ER+" & Clinical$Her2.Expr=="-")] <- "ER+/HER2-"
Clinical$Group <- factor(Clinical$Group)
Clinical <- Clinical[which(Clinical$Group=="ER+/HER2-"),]
all.lr <- all.lr[,which(colnames(all.lr) %in% Clinical$METABRIC.ID)]
all.dr <- all.dr[,which(colnames(all.dr) %in% Clinical$METABRIC.ID)]
Clinical.LR <- Clinical[which(Clinical$METABRIC.ID %in% colnames(all.lr)),]
Clinical.LR <- Clinical.LR[match(colnames(all.lr), Clinical.LR$METABRIC.ID),]
Clinical.DR <- Clinical[which(Clinical$METABRIC.ID %in% colnames(all.dr)),]
Clinical.DR <- Clinical.DR[match(colnames(all.dr), Clinical.DR$METABRIC.ID),]


IC.lr <- all.lr[,match(Clinical.LR$METABRIC.ID, colnames(all.lr))]
IC.dr <- all.dr[,match(Clinical.DR$METABRIC.ID, colnames(all.dr))]
IC.lr <- IC.lr[,which(Clinical.LR$iC10 %in% c(1, 2, 6, 9))]
IC.dr <- IC.dr[,which(Clinical.DR$iC10 %in% c(1, 2, 6, 9))]
Clinical.LR <- subset(Clinical.LR, iC10 %in% c(1, 2, 6, 9))
Clinical.DR <- subset(Clinical.DR, iC10 %in% c(1, 2, 6, 9))
## m1 <- manova(t(IC.lr) ~ Clinical.LR$iC10)
## m3 <- manova(t(IC.dr) ~ Clinical.DR$iC10)

library(multcomp)
res <- matrix(NA, nrow=6, ncol=5)
for (i in 1:5) {
    tmp <- data.frame(y=asin(IC.lr[i,]), iC10=Clinical.LR$iC10)
    m <- lm(y ~ iC10, data=tmp)
    ct <- print(summary(glht(m, linfct = mcp(iC10 = "Tukey"))))
    res[,i] <- ct$test$pvalues
}
ICM.LR <- data.frame(Contr=names(coef(ct)),res)
library(multcomp)
res <- matrix(NA, nrow=6, ncol=5)
for (i in 1:5) {
    tmp <- data.frame(y=asin(IC.dr[i,]), iC10=Clinical.DR$iC10)
    m <- lm(y ~ iC10, data=tmp)
    ct <- print(summary(glht(m, linfct = mcp(iC10 = "Tukey"))))
    res[,i] <- ct$test$pvalues
}
ICM.DR <- data.frame(Contr=names(coef(ct)),res)



sub.res <- subset(ICModel, Relapse=="Loco-regional Relapse")
Y <- sub.res[order(sub.res$Year, sub.res$IntClust),]
Y$Dist <- 1:nrow(Y)
Y$Dist <- Y$Dist + rep(c(0, 2, 4, 6, 8), rep(5,5))
plot(X ~ Dist, data=Y, ylim=c(0,1), axes=F, col=Y$col, pch=Y$pch, ylab="Probability of distant relapse/cancer death", xlab="Years after loco-regional relapse", main="b", cex.axis=1.5, cex.lab=1.5)
usr <- par('usr')
rect(Y$Dist[seq(from=1, by=5, length=5)]-0.85, usr[3],
     Y$Dist[seq(from=5, by=5, length=5)]+0.85, usr[4],
     col=adjustcolor('grey', alpha=0.4), border=NA, cex=2)
points(li ~ Dist, data=Y, pch="-", col=Y$col, cex=2)
points(ui ~ Dist, data=Y, pch="-", col=Y$col, cex=2)
segments(Y$Dist, Y$li, Y$Dist, Y$ui, col=Y$col, cex=2)
axis(2)
axis(1, at=unique(Y$Dist)[seq(from=3, by=5, length=5)] - 0.5, labels=c(2, 5, 10, 15, 20))
box()
dev.off()

##########################################################
##########################################################
##########################################################
# Bottom panel

pdf("FigureS5b.pdf", width=12, height=5)
library(survival)

Clinical <- read.table(file="../../TableS6.txt", header=T, sep="\t", quote="", comment.char="", stringsAsFactors=FALSE)
Clinical$NaturalDeath <- 1 * (Clinical$Last.Followup.Status %in% c("d", "d-o.c."))

## We remove Samples with no follow-up time

ids <- which(Clinical$T==0)
if (length(ids)>0) Clinical <- Clinical[-ids,]
## We remove Samples with no follow-up time or death known

## We remove samples with stage 4
Clinical <- Clinical[which(Clinical$Stage!=4),]

## We remove benign, DCIS or PHYL
bad.hist <- which(Clinical$Histological.Type %in% c("BENIGN", "PHYL"))
if (length(bad.hist)>0) Clinical <- Clinical[-bad.hist,]

Clinical <- Clinical[which(!is.na(Clinical$T)),]
Clinical <- Clinical[which(!is.na(Clinical$Death)),]

## We move a bit relapses from diagnosis
Clinical$TLR[which(Clinical$TLR==0)] <- 0.1
Clinical$TDR[which(Clinical$TDR==0)] <- 0.1
Clinical$T[which(Clinical$T==Clinical$TDR & Clinical$DR==1)] <- Clinical$T[which(Clinical$T==Clinical$TDR & Clinical$DR==1)] + 0.1
Clinical$T[which(Clinical$T==Clinical$TLR & Clinical$LR==1)] <- Clinical$T[which(Clinical$T==Clinical$TLR & Clinical$LR==1)] + 0.1
Clinical$TDR[which(Clinical$TLR==Clinical$TDR & Clinical$LR==1 & Clinical$DR==1)] <- Clinical$TDR[which(Clinical$TLR==Clinical$TDR & Clinical$LR==1 & Clinical$DR==1)] + 0.1

Clinical$TLR[which(Clinical$LR==0)] <- Clinical$T[which(Clinical$LR==0)]
Clinical$TDR[which(Clinical$DR==0)] <- Clinical$T[which(Clinical$DR==0)]

## If local relapse occured after distant, we dont consider it
Clinical$LR[which(Clinical$TLR>Clinical$TDR)] <- 0
Clinical$TLR[which(Clinical$TLR>Clinical$TDR)] <- Clinical$T[which(Clinical$TLR>Clinical$TDR)]

Clinical$T <- Clinical$T/365.25
Clinical$TLR <- Clinical$TLR/365.25
Clinical$TDR <- Clinical$TDR/365.25

Clinical$LN <- Clinical$Lymph.Nodes.Positive
Clinical$LN[which(Clinical$LN>=10)] <- 10
Clinical$AGE <- Clinical$Age.At.Diagnosis
Clinical$GRADE <- as.numeric(as.character(Clinical$Grade))
Clinical$SIZE <- as.numeric(as.character(Clinical$Size))
Clinical$HT <- 1 * (Clinical$HT!="null")
Clinical$HT <- factor(Clinical$HT, levels=c(0, 1), labels=c("NO", "YES"))
Clinical$CT <- 1 * (Clinical$CT!="null")
Clinical$CT <- factor(Clinical$CT, levels=c(0, 1), labels=c("NO", "YES"))
Clinical$RT[which(Clinical$RT %in% c("null", "NONE RECORDED IN LANTIS", "Nnne"))] <- 0
Clinical$RT[which(Clinical$RT!=0)] <- 1
Clinical$RT <- factor(Clinical$RT, levels=c(0, 1), labels=c("NO", "YES"))
Clinical$BS <- factor(Clinical$Breast.Surgery, levels=c("BREAST CONSERVING", "MASTECTOMY"), labels=c("BC", "M"))
Clinical$iC10 <- factor(Clinical$iC10, levels=c(1:3, "4ER+", "4ER-", 5:10))
Clinical$ClinicalClass <- NA
Clinical$ClinicalClass[which(Clinical$ER.Expr=="-" & Clinical$Her2.Expr=="-")] <- "ER-/HER2-"
Clinical$ClinicalClass[which(Clinical$ER.Expr=="+" & Clinical$Her2.Expr=="-")] <- "ER+/HER2-"
Clinical$ClinicalClass[which(Clinical$ER.Expr=="+" & Clinical$Her2.Expr=="+")] <- "ER+/HER2+"
Clinical$ClinicalClass[which(Clinical$ER.Expr=="-" & Clinical$Her2.Expr=="+")] <- "ER-/HER2+"
Clinical$ClinicalClass <- factor(Clinical$ClinicalClass)
TClinical <- subset(Clinical, LR==1)
TClinical$T2 <- TClinical$T - TClinical$TLR
TClinical$T2[which(TClinical$DR==1)] <- I(TClinical$TDR-TClinical$TLR)[which(TClinical$DR==1)]
TClinical$I2 <- 1 *(TClinical$DR==1 | TClinical$DeathBreast==1)

m <- survfit(Surv(T2, I2) ~ iC10, data=TClinical)
event <- rep(NA, nrow(TClinical))
event[which(TClinical$Death==0 & TClinical$DR==0)] <- "RELAPSE-FREE"
event[which(TClinical$I2==1)] <- "RELAPSE"
event[which(TClinical$I2==0 & TClinical$Death==1)] <- "DEATH/OTHER"
library(cmprsk)
event <-  factor(event, levels=c("RELAPSE-FREE", "RELAPSE", "DEATH/OTHER"))
m1 <- cuminc(ftime=TClinical$T2, fstatus=as.character(event),
             group=TClinical$iC10, cencode="RELAPSE-FREE")
m2 <- survfit(Surv(TClinical$T2, as.numeric(event),type="mstate" ) ~ TClinical$iC10)
## same result
res <- sapply(m1[12:22], function(x) {
           id <- which.min((x$est-0.5)<0)
           x$time[id]
       })
max.obs <- sapply(m1[12:22], function(x) max(x$est))
res[which(max.obs<0.5)] <- NA
ui <- lapply(m1[12:22], function(x) pmax(0, x$est - 1.96*sqrt(x$var)))
max.obs <- sapply(ui, function(x) max(x))
li <- lapply(m1[12:22], function(x) pmin(1, x$est + 1.96*sqrt(x$var)))
li <- lapply(li, function(x) {
           id <- which.min((x-0.5)<0)
       })
li <- mapply(function(x,y) y$time[x], x=li, y=m1[12:22])
ui <- lapply(ui, function(x) {
           id <- which.min((x-0.5)<0)
       })
ui <- mapply(function(x,y) y$time[x], x=ui, y=m1[12:22])
ui[which(max.obs<0.5)] <- NA


par(mfrow=c(1,4))
load(file="../Models/IntClust_AllProbs.RData")


tmp <- mapply(function(x, y) x[,-1]+y[,-1], x=pld, y=plc)
tmp <- mapply(function(x, y) x+y[,-1], x=tmp, y=pldc)
INTCLUST.mean <- sapply(tmp, function(x) apply(x, 1, mean))
INTCLUST.sd <- sapply(tmp, function(x) apply(x, 1, function(z) sd(z)/sqrt(length(z))))

coliClust <- c('#FF5500', '#00EE76', '#CD3278','#00C5CD', '#B5D0D2', '#8B0000',
               '#FFFF40', '#0000CD', '#FFAA00', '#EE82EE', '#7D26CD')
times <-  seq(from=0, to=20, by=0.25)


plot(times, rep(0, length(times)), type="n",
     ylim=c(0,1), xlim=c(0, 10), xlab="Years after loco-regional relapse",
     ylab="", main="")
mtext("Probability of distant relapse\n or cancer-related death",side=2, line=2,
      outer=F, cex=0.6)

for (id in 1:11) {
    lines(times, INTCLUST.mean[,id],
          lwd=2, col=coliClust[id])
    polygon(c(times, rev(times), times[1]),
            c(INTCLUST.mean[,id]- 1.96*INTCLUST.sd[,id],
              rev(INTCLUST.mean[,id]+ 1.96*INTCLUST.sd[,id]),
              INTCLUST.mean[,id][1]- 1.96*INTCLUST.sd[,id][1]),
            col=adjustcolor(coliClust[id], alpha=0.2), border=NA)
    }


dd <- order(res)

library(rms)
ui[which(is.na(ui))] <- Inf
correct.ids <- which(ui>10 & is.finite(res))
ui[correct.ids] <- res[correct.ids]
ui[which(is.na(res))] <- NA
li[which(is.na(res))] <- NA

errbar(1:11, res[dd], ui[dd],
       pmax(0, li[dd]), add=F, col=coliClust[dd],
       errbar.col=coliClust[dd], axes=F, ylab="Years after loco regional relapse",
       cex=2, xlab="", ylim=c(0, 10), lwd=0.5)

if (length(correct.ids>0)) {
        for (i in correct.ids) {
            arrows(which(dd==i), res[i], which(dd==i), 10,
                   col=coliClust[i], length=0.05, lwd=0.5)
        }
    }
axis(2)
axis(1, at=1:11, sub(" RELAPSE", "", names(res))[dd], cex.axis=0.9, las=2)
box()
text(c(10, 11), c(10, 10), c("*", "*"), cex=0.9)


times <- seq(from=0, to=20, by=0.25)
INTCLUST.mean <- sapply(pdc, function(x) apply(x[,-1], 1, mean))
INTCLUST.sd <- sapply(pdc, function(x) apply(x[,-1], 1, function(z) sd(z)/sqrt(length(z))))
plot(times, rep(0, length(times)), type="n",
     ylim=c(0,1), xlim=c(0, 10), xlab="Years after distant relapse",main="",
     ylab="")
mtext("Probability of cancer-related\n death",side=2, line=2,
      outer=F, cex=0.6)

for (id in 1:11) {
    lines(times, INTCLUST.mean[,id],
          lwd=2, col=coliClust[id])
    polygon(c(times, rev(times), times[1]),
            c(INTCLUST.mean[,id]- 1.96*INTCLUST.sd[,id],
              rev(INTCLUST.mean[,id]+ 1.96*INTCLUST.sd[,id]),
              INTCLUST.mean[,id][1]- 1.96*INTCLUST.sd[,id][1]),
            col=adjustcolor(coliClust[id], alpha=0.2), border=NA)
    }

TClinical <- subset(Clinical, DR==1)
TClinical$T2 <- TClinical$T - TClinical$TDR
TClinical$I2 <- 1 *(TClinical$DeathBreast==1)

m <- survfit(Surv(T2, I2) ~ iC10, data=TClinical)
event <- rep(NA, nrow(TClinical))
event[which(TClinical$Death==0)] <- "RELAPSE-FREE"
event[which(TClinical$I2==1)] <- "RELAPSE"
event[which(TClinical$I2==0 & TClinical$Death==1)] <- "DEATH/OTHER"
library(cmprsk)
event <-  factor(event, levels=c("RELAPSE-FREE", "RELAPSE", "DEATH/OTHER"))
m1 <- cuminc(ftime=TClinical$T2, fstatus=as.character(event),
             group=TClinical$iC10, cencode="RELAPSE-FREE")
m2 <- survfit(Surv(TClinical$T2, as.numeric(event),type="mstate" ) ~ TClinical$iC10)
## same result
res <- sapply(m1[12:22], function(x) {
           id <- which.min((x$est-0.5)<0)
           x$time[id]
       })
max.obs <- sapply(m1[12:22], function(x) max(x$est))
res[which(max.obs<0.5)] <- NA
ui <- lapply(m1[12:22], function(x) pmax(0, x$est - 1.96*sqrt(x$var)))
max.obs <- sapply(ui, function(x) max(x))
li <- lapply(m1[12:22], function(x) pmin(1, x$est + 1.96*sqrt(x$var)))
li <- lapply(li, function(x) {
           id <- which.min((x-0.5)<0)
       })
li <- mapply(function(x,y) y$time[x], x=li, y=m1[12:22])
ui <- lapply(ui, function(x) {
           id <- which.min((x-0.5)<0)
       })
ui <- mapply(function(x,y) y$time[x], x=ui, y=m1[12:22])
ui[which(max.obs<0.5)] <- NA

dd <- order(res)
## plot(res$iClust[dd], res$Years[dd], pch=21, col=coliClust[dd])
library(rms)
ui[which(is.na(ui))] <- Inf
correct.ids <- which(ui>10 & is.finite(res))
ui[correct.ids] <- res[correct.ids]
ui[which(is.na(res))] <- NA
li[which(is.na(res))] <- NA

errbar(1:11, res[dd], ui[dd],
       pmax(0, li[dd]), add=F, col=coliClust[dd],
       errbar.col=coliClust[dd], axes=F, ylab="Years after distant relapse",
       cex=2, xlab="", ylim=c(0, 10), lwd=0.5)

if (length(correct.ids>0)) {
        for (i in correct.ids) {
            arrows(which(dd==i), res[i], which(dd==i), 10,
                   col=coliClust[i], length=0.05, lwd=0.5)
        }
    }
axis(2)
axis(1, at=1:11, sub(" RELAPSE", "", names(res))[dd], cex.axis=0.9, las=2)
box()
dev.off()

