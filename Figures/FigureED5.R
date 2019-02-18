## First part
rm(list=ls())

library(Hmisc)
load(file="../Models/IntClust_AllProbs.RData")
library(survival)

lr <- mapply(function(x,y) x[,-1] + y[,-1], x=plc, y=pld)
lr <- mapply(function(x,y) x + y[,-1], x=lr, y=pldc)
lr <- mapply(function(x,y) x + y[,-1], x=lr, y=pldo)

dr <- mapply(function(x,y) x[,-1] + y[,-1], x=psd, y=psdc)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psdo)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psld)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psldc)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psldo)
dr <- mapply(function(x,y) x + y[,-1], x=dr, y=psc)



times.lr <-  pl[[1]][,1]
ids.lr <- which(times.lr %in% c(2, 5, 10, 15, 20))
times <-  ps[[1]][,1]
ids.s <- which(times %in% c(2, 5, 10, 15, 20))

lr.INTCLUST.mean <- sapply(lr, function(x) apply(x, 1, mean))[ids.lr,]
lr.INTCLUST.SE <- sapply(lr, function(x) apply(x, 1, function(y) sd(y)/sqrt(length(y))))[ids.lr,]
dr.INTCLUST.mean <- sapply(dr, function(x) apply(x, 1, mean))[ids.s,]
dr.INTCLUST.SE <- sapply(dr, function(x) apply(x, 1, function(y) sd(y)/sqrt(length(y))))[ids.s,]



all.lr <- do.call("cbind", lr)
all.dr <- do.call("cbind", dr)
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




coliClust <- c('#FF5500', '#00EE76', '#CD3278','#00C5CD', '#B5D0D2', '#8B0000',
               '#FFFF40', '#0000CD', '#FFAA00', '#EE82EE', '#7D26CD')

ss.lr <- paste0("n=", sapply(lr, function(x) ncol(x)))
ss.dr <- paste0("n=", sapply(dr, function(x) ncol(x)))
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

library(Hmisc)
library(survival)


par(mfrow=c(1, 2))
sub.res <- subset(ICModel, Relapse=="Distant Relapse")
Y <- sub.res[order(sub.res$Year, sub.res$IntClust),]
Y$Dist <- 1:nrow(Y)
Y$Dist <- Y$Dist + rep(c(0, 2, 4, 6, 8), rep(5,5))
plot(X ~ Dist, data=Y, ylim=c(0,1), axes=F, col=Y$col, pch=Y$pch,
     ylab="Probability of distant relapse/cancer death", xlab="Years after surgery", main="",
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



sub.res <- subset(ICModel, Relapse=="Loco-regional Relapse")
Y <- sub.res[order(sub.res$Year, sub.res$IntClust),]
Y$Dist <- 1:nrow(Y)
Y$Dist <- Y$Dist + rep(c(0, 2, 4, 6, 8), rep(5,5))
plot(X ~ Dist, data=Y, ylim=c(0,1), axes=F, col=Y$col, pch=Y$pch, ylab="Probability of distant relapse/cancer death", xlab="Years after loco-regional relapse", main="", cex.axis=1.5, cex.lab=1.5)
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

library(survival)

Clinical <- read.table(file="../../TableS6.txt", header=T, sep="\t", quote="", comment.char="", stringsAsFactors=FALSE)
Clinical$NaturalDeath <- 1 * (Clinical$Last.Followup.Status %in% c("d", "d-o.c."))

## We remove Samples with no follow-up time

ids <- which(Clinical$T==0)
if (length(ids)>0) Clinical <- Clinical[-ids,]
## We remove Samples with no follow-up time or death known

## We remove samples with stage 4
Clinical <- Clinical[-which(Clinical$Stage==4),]

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
     ylab="", main="", cex.axis=1.1, cex.lab=1.3)
mtext("Probability of DR/cancer death",side=2, line=2,
      outer=F, cex=0.8)

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
       errbar.col=coliClust[dd], axes=F, ylab="",
       cex=2, xlab="", ylim=c(0, 10), lwd=0.5, cex.lab=1.1)

if (length(correct.ids>0)) {
        for (i in correct.ids) {
            arrows(which(dd==i), res[i], which(dd==i), 10,
                   col=coliClust[i], length=0.05, lwd=0.5)
        }
    }
axis(2, cex.axis=1.1, cex.lab=1.3)
mtext("Years after loco regional relapse", side=2, line=2, cex=0.8)
axis(1, at=1:11, sub(" RELAPSE", "", names(res))[dd], cex.axis=1.1, las=2, cex.lab=1.3)
box()
text(c(10, 11), c(10, 10), c("*", "*"), cex=0.9)


times <- seq(from=0, to=20, by=0.25)
INTCLUST.mean <- sapply(pdc, function(x) apply(x[,-1], 1, mean))
INTCLUST.sd <- sapply(pdc, function(x) apply(x[,-1], 1, function(z) sd(z)/sqrt(length(z))))
plot(times, rep(0, length(times)), type="n",
     ylim=c(0,1), xlim=c(0, 10), xlab="Years after distant relapse",main="",
     ylab="", cex.axis=1.1, cex.lab=1.2)
mtext("Probability of cancer death",side=2, line=2,
      outer=F, cex=0.8)

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
       errbar.col=coliClust[dd], axes=F, ylab="",
       cex=2, xlab="", ylim=c(0, 10), lwd=0.5)

if (length(correct.ids>0)) {
        for (i in correct.ids) {
            arrows(which(dd==i), res[i], which(dd==i), 10,
                   col=coliClust[i], length=0.05, lwd=0.5)
        }
    }
axis(2, cex.axis=1.1, cex.lab=1.3)
mtext("Years after distant relapse",side=2, line=2,
      outer=F, cex=0.8)
axis(1, at=1:11, sub(" RELAPSE", "", names(res))[dd], cex.axis=1.1, cex.lab=1.3, las=2)
box()


## Second part

rm(list=ls())
## Top panel
library(Hmisc)
library(Hmisc)
load(file="../Models/FOURGROUPS_AllProbs.RData")
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

r0 <- do.call("cbind", r)
rownames(r0) <- psc[[1]][,1]

load(file="../Models/FourGroups_AllProbs_FROM5.RData")
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

load(file="../Models/FourGroups_AllProbs_FROM10.RData")
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
Clinical <- Clinical[which(Clinical$METABRIC.ID %in% colnames(r0)),]
Clinical <- Clinical[match(colnames(r0), Clinical$METABRIC.ID),]
mean(Clinical$METABRIC.ID==colnames(r0))



r0.INTCLUST.mean <- t(apply(r0, 1, function(y) tapply(y, Clinical$Group, mean)))
r0.INTCLUST.SE <- t(apply(r0, 1, function(y) tapply(y, Clinical$Group, function(z) sd(z)/sqrt(length(z)))))

r5.INTCLUST.mean <- t(apply(r5, 1, function(y) tapply(y, Clinical$Group, mean)))
r5.INTCLUST.SE <- t(apply(r5, 1, function(y) tapply(y, Clinical$Group, function(z) sd(z)/sqrt(length(z)))))

r10.INTCLUST.mean <- t(apply(r10, 1, function(y) tapply(y, Clinical$Group, mean)))
r10.INTCLUST.SE <- t(apply(r10, 1, function(y) tapply(y, Clinical$Group, function(z) sd(z)/sqrt(length(z)))))



coliClust <- c("lightskyblue", "lightsalmon", "plum2", "wheat4")
names(coliClust) <- c("ER+/HER2+", "ER+/HER2-", "ER-/HER2+", "ER-/HER2-")

ss <- paste0("n=", sapply(r, function(x) ncol(x)))
res <- data.frame(X=c(r0.INTCLUST.mean, r5.INTCLUST.mean, r10.INTCLUST.mean), Year=c(rep(rownames(r0.INTCLUST.mean), 4),
                                                                                  rep(rownames(r5.INTCLUST.mean), 4),
                                                                                  rep(rownames(r10.INTCLUST.mean), 4)),
                  Predictions=rep(c("After Surgery", "5 years disease-free",
                      "10 years disease-free"), c(4*nrow(r0.INTCLUST.mean), 4*nrow(r5.INTCLUST.mean), 4*nrow(r10.INTCLUST.mean))),
                  IntClust=c(rep(colnames(r0.INTCLUST.mean), rep(nrow(r0.INTCLUST.mean), 4)),
                      rep(colnames(r5.INTCLUST.mean), rep(nrow(r5.INTCLUST.mean), 4)),
                      rep(colnames(r10.INTCLUST.mean), rep(nrow(r10.INTCLUST.mean), 4))),
                  li=c(r0.INTCLUST.mean, r5.INTCLUST.mean, r10.INTCLUST.mean) - 1.96*c(r0.INTCLUST.SE, r5.INTCLUST.SE, r10.INTCLUST.SE),
                  ui=c(r0.INTCLUST.mean, r5.INTCLUST.mean, r10.INTCLUST.mean) + 1.96*c(r0.INTCLUST.SE, r5.INTCLUST.SE, r10.INTCLUST.SE))

res$Predictions <- factor(res$Predictions, levels=c("After Surgery",
                                               "5 years disease-free",
                                               "10 years disease-free"))

names(ss) <- levels(res$IntClust)

res$pch <- 17
res$pch[which(res$Predictions=="5 years disease-free")] <- 18
res$pch[which(res$Predictions=="10 years disease-free")] <- 19
res$colr <- "black"
res$colr[which(res$Predictions=="5 years disease-free")] <- "darkred"
res$colr[which(res$Predictions=="10 years disease-free")] <- "olivedrab"


res$Year <- as.numeric(as.character(res$Year))
res$IntClust <- factor(res$IntClust, levels=c("ER+/HER2-", "ER-/HER2-", "ER+/HER2+", "ER-/HER2+"))
ss <- ss[levels(res$IntClust)]
coliClust <- coliClust[levels(res$IntClust)]

pp <- xyplot(X ~ Year| IntClust, groups=Predictions, data=res,
             layout=c(4, 1, 1), ylab="Probability of Relapse", xlab="Years after surgery", ylim=c(0, 0.8),
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



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## Bottom panel
rm(list=ls())
library(survival)

Clinical <- read.table(file="../../TableS6.txt", header=T, sep="\t", quote="", comment.char="", stringsAsFactors=FALSE)
Clinical$NaturalDeath <- 1 * (Clinical$Last.Followup.Status %in% c("d", "d-o.c."))

## We remove Samples with no follow-up time

ids <- which(Clinical$T==0)
if (length(ids)>0) Clinical <- Clinical[-ids,]
## We remove Samples with no follow-up time or death known

## We remove samples with stage 4
Clinical <- Clinical[-which(Clinical$Stage==4),]

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
             group=TClinical$ClinicalClass, cencode="RELAPSE-FREE")
m2 <- survfit(Surv(TClinical$T2, as.numeric(event),type="mstate" ) ~ TClinical$ClinicalClass)
## same result
res <- sapply(m1[5:8], function(x) {
           id <- which.min((x$est-0.5)<0)
           x$time[id]
       })
max.obs <- sapply(m1[5:8], function(x) max(x$est))
res[which(max.obs<0.5)] <- NA
ui <- lapply(m1[5:8], function(x) pmax(0, x$est - 1.96*sqrt(x$var)))
max.obs <- sapply(ui, function(x) max(x))
li <- lapply(m1[5:8], function(x) pmin(1, x$est + 1.96*sqrt(x$var)))
li <- lapply(li, function(x) {
           id <- which.min((x-0.5)<0)
       })
li <- mapply(function(x,y) y$time[x], x=li, y=m1[5:8])
ui <- lapply(ui, function(x) {
           id <- which.min((x-0.5)<0)
       })
ui <- mapply(function(x,y) y$time[x], x=ui, y=m1[5:8])
ui[which(max.obs<0.5)] <- NA


load(file="../Models/FourGroups_AllProbs.RData")

par(mfrow=c(1,4))


INTCLUST.mean <- mapply(function(i,j) i + j[,-1], i=mapply(function(x, y) x[,-1]+y[,-1], x=pld, y=plc), j=pldc)
INTCLUST.mean <- sapply(INTCLUST.mean, function(x) apply(x, 1, mean))
INTCLUST.sd <- mapply(function(i,j) i + j[,-1], i=mapply(function(x, y) x[,-1]+y[,-1], x=pld, y=plc), j=pldc)
INTCLUST.sd <- sapply(INTCLUST.sd, function(x) apply(x, 1, function(z) sd(z)/sqrt(length(z))))

coliClust <- c("lightskyblue", "lightsalmon", "plum2", "wheat4")
names(coliClust) <- c("ER+/HER2+", "ER+/HER2-", "ER-/HER2+", "ER-/HER2-")
times <-  pl[[1]][,1]

coliClust <- coliClust[colnames(INTCLUST.mean)]
plot(times, rep(0, length(times)), type="n",
     ylim=c(0,1), xlim=c(0, 10), xlab="Years after loco-regional relapse", main="",
     ylab="Probability of distant relapse or cancer-related death")
for (id in 1:4) {
    lines(times, INTCLUST.mean[,id],
          lwd=2, col=coliClust[id])
    polygon(c(times, rev(times), times[1]),
            c(INTCLUST.mean[,id]- 1.96*INTCLUST.sd[,id],
              rev(INTCLUST.mean[,id]+ 1.96*INTCLUST.sd[,id]),
              INTCLUST.mean[,id][1]- 1.96*INTCLUST.sd[,id][1]),
            col=adjustcolor(coliClust[id], alpha=0.2), border=NA)
    }




coliClust <- c("lightskyblue", "lightsalmon", "plum2", "wheat4")
names(coliClust) <- c("ER+/HER2+", "ER+/HER2-", "ER-/HER2+", "ER-/HER2-")

dd <- order(res)
library(rms)
ui[which(is.na(ui))] <- Inf
correct.ids <- which(ui>10 & is.finite(res))
ui[correct.ids] <- res[correct.ids]
ui[which(is.na(res))] <- NA
li[which(is.na(res))] <- NA

errbar(1:4, res[dd], ui[dd],
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
axis(1, at=1:4, sub(" RELAPSE", "", names(res))[dd], cex.axis=0.9, las=2)
box()


coliClust <- c("lightskyblue", "lightsalmon", "plum2", "wheat4")
names(coliClust) <- c("ER+/HER2+", "ER+/HER2-", "ER-/HER2+", "ER-/HER2-")

coliClust <- coliClust[colnames(INTCLUST.mean)]
times <- pdc[[1]][,1]
INTCLUST.mean <- sapply(pdc, function(x) apply(x[,-1], 1, mean))
INTCLUST.sd <- sapply(pdc, function(x) apply(x[,-1], 1, function(z) sd(z)/sqrt(length(z))))
plot(times, rep(0, length(times)), type="n",
     ylim=c(0,1), xlim=c(0, 10), xlab="Years after distant relapse",
     ylab="Probability of cancer-related death", main="", )

for (id in 1:4) {
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
             group=TClinical$ClinicalClass, cencode="RELAPSE-FREE")
m2 <- survfit(Surv(TClinical$T2, as.numeric(event),type="mstate" ) ~ TClinical$ClinicalClass)
## same result
res <- sapply(m1[5:8], function(x) {
           id <- which.min((x$est-0.5)<0)
           x$time[id]
       })
max.obs <- sapply(m1[5:8], function(x) max(x$est))
res[which(max.obs<0.5)] <- NA
ui <- lapply(m1[5:8], function(x) pmax(0, x$est - 1.96*sqrt(x$var)))
max.obs <- sapply(ui, function(x) max(x))
li <- lapply(m1[5:8], function(x) pmin(1, x$est + 1.96*sqrt(x$var)))
li <- lapply(li, function(x) {
           id <- which.min((x-0.5)<0)
       })
li <- mapply(function(x,y) y$time[x], x=li, y=m1[5:8])
ui <- lapply(ui, function(x) {
           id <- which.min((x-0.5)<0)
       })
ui <- mapply(function(x,y) y$time[x], x=ui, y=m1[5:8])
ui[which(max.obs<0.5)] <- NA

coliClust <- c("lightskyblue", "lightsalmon", "plum2", "wheat4")
names(coliClust) <- c("ER+/HER2+", "ER+/HER2-", "ER-/HER2+", "ER-/HER2-")
dd <- order(res)
library(rms)
ui[which(is.na(ui))] <- Inf
correct.ids <- which(ui>10 & is.finite(res))
ui[correct.ids] <- res[correct.ids]
ui[which(is.na(res))] <- NA
li[which(is.na(res))] <- NA

errbar(1:4, res[dd], ui[dd],
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
axis(1, at=1:4, sub(" RELAPSE", "", names(res))[dd], cex.axis=0.9, las=2)
box()


## Third part

## Top panel

rm(list=ls())
## First panel:
library(Hmisc)
load(file="../Models/Pam50_AllProbs.RData")
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

r0 <- do.call("cbind", r)
rownames(r0) <- psc[[1]][,1]

load(file="../Models/Pam50_AllProbs_FROM5.RData")
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

load(file="../Models/Pam50_AllProbs_FROM10.RData")
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
Clinical <- Clinical[which(Clinical$METABRIC.ID %in% colnames(r0)),]
Clinical <- Clinical[match(colnames(r0), Clinical$METABRIC.ID),]
mean(Clinical$METABRIC.ID==colnames(r0))



r0.INTCLUST.mean <- t(apply(r0, 1, function(y) tapply(y, Clinical$Pam50, mean)))
r0.INTCLUST.SE <- t(apply(r0, 1, function(y) tapply(y, Clinical$Pam50, function(z) sd(z)/sqrt(length(z)))))

r5.INTCLUST.mean <- t(apply(r5, 1, function(y) tapply(y, Clinical$Pam50, mean)))
r5.INTCLUST.SE <- t(apply(r5, 1, function(y) tapply(y, Clinical$Pam50, function(z) sd(z)/sqrt(length(z)))))

r10.INTCLUST.mean <- t(apply(r10, 1, function(y) tapply(y, Clinical$Pam50, mean)))
r10.INTCLUST.SE <- t(apply(r10, 1, function(y) tapply(y, Clinical$Pam50, function(z) sd(z)/sqrt(length(z)))))

coliClust <- c("#E41A1C", "#FB9A99", "#1F78B4", "#A6CEE3", "#66A61E")
names(coliClust) <- colnames(r0.INTCLUST.mean)

ss <- paste0("n=", sapply(r, function(x) ncol(x)))
res <- data.frame(X=c(r0.INTCLUST.mean, r5.INTCLUST.mean, r10.INTCLUST.mean), Year=c(rep(rownames(r0.INTCLUST.mean), 5),
                                                                                  rep(rownames(r5.INTCLUST.mean), 5),
                                                                                  rep(rownames(r10.INTCLUST.mean), 5)),
                  Predictions=rep(c("After Surgery", "5 years disease-free",
                      "10 years disease-free"), c(5*nrow(r0.INTCLUST.mean), 5*nrow(r5.INTCLUST.mean), 5*nrow(r10.INTCLUST.mean))),
                  IntClust=c(rep(colnames(r0.INTCLUST.mean), rep(nrow(r0.INTCLUST.mean), 5)),
                      rep(colnames(r5.INTCLUST.mean), rep(nrow(r5.INTCLUST.mean), 5)),
                      rep(colnames(r10.INTCLUST.mean), rep(nrow(r10.INTCLUST.mean), 5))),
                  li=c(r0.INTCLUST.mean, r5.INTCLUST.mean, r10.INTCLUST.mean) - 1.96*c(r0.INTCLUST.SE, r5.INTCLUST.SE, r10.INTCLUST.SE),
                  ui=c(r0.INTCLUST.mean, r5.INTCLUST.mean, r10.INTCLUST.mean) + 1.96*c(r0.INTCLUST.SE, r5.INTCLUST.SE, r10.INTCLUST.SE))

res$Predictions <- factor(res$Predictions, levels=c("After Surgery",
                                               "5 years disease-free",
                                               "10 years disease-free"))

names(ss) <- levels(res$IntClust)

res$pch <- 17
res$pch[which(res$Predictions=="5 years disease-free")] <- 18
res$pch[which(res$Predictions=="10 years disease-free")] <- 19
res$colr <- "black"
res$colr[which(res$Predictions=="5 years disease-free")] <- "darkred"
res$colr[which(res$Predictions=="10 years disease-free")] <- "olivedrab"


res$Year <- as.numeric(as.character(res$Year))
res$IntClust <- factor(res$IntClust, levels=c("LumA", "Normal", "Basal", "LumB", "Her2"))
ss <- ss[levels(res$IntClust)]
coliClust <- coliClust[levels(res$IntClust)]


pp <- xyplot(X ~ Year| IntClust, groups=Predictions, data=res,
             layout=c(5, 1, 1), ylab="Probability of Relapse", xlab="Years after surgery", ylim=c(0, 0.8),
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
####################################################
## Bottom panel


library(survival)

## par(mfrow=c(2,2))

Clinical <- read.table(file="../../TableS6.txt", header=T, sep="\t", quote="", comment.char="", stringsAsFactors=FALSE)
Clinical$NaturalDeath <- 1 * (Clinical$Last.Followup.Status %in% c("d", "d-o.c."))

## We remove Samples with no follow-up time

ids <- which(Clinical$T==0)
if (length(ids)>0) Clinical <- Clinical[-ids,]
## We remove Samples with no follow-up time or death known

## We remove samples with stage 4
Clinical <- Clinical[-which(Clinical$Stage==4),]

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
Clinical$Pam50Subtype <- factor(Clinical$Pam50Subtype, levels=c("Basal", "Her2", "LumA", "LumB", "Normal"))
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
             group=TClinical$Pam50, cencode="RELAPSE-FREE")
m2 <- survfit(Surv(TClinical$T2, as.numeric(event),type="mstate" ) ~ TClinical$Pam50)
## same result
res <- sapply(m1[6:10], function(x) {
           id <- which.min((x$est-0.5)<0)
           x$time[id]
       })
max.obs <- sapply(m1[6:10], function(x) max(x$est))
res[which(max.obs<0.5)] <- NA
ui <- lapply(m1[6:10], function(x) pmax(0, x$est - 1.96*sqrt(x$var)))
max.obs <- sapply(ui, function(x) max(x))
li <- lapply(m1[6:10], function(x) pmin(1, x$est + 1.96*sqrt(x$var)))
li <- lapply(li, function(x) {
           id <- which.min((x-0.5)<0)
       })
li <- mapply(function(x,y) y$time[x], x=li, y=m1[6:10])
ui <- lapply(ui, function(x) {
           id <- which.min((x-0.5)<0)
       })
ui <- mapply(function(x,y) y$time[x], x=ui, y=m1[6:10])
ui[which(max.obs<0.5)] <- NA


load(file="../Models/Pam50_AllProbs.RData")

par(mfrow=c(1,4))



INTCLUST.mean <- mapply(function(i,j) i + j[,-1], i=mapply(function(x, y) x[,-1]+y[,-1], x=pld, y=plc), j=pldc)
INTCLUST.mean <- sapply(INTCLUST.mean, function(x) apply(x, 1, mean))
INTCLUST.sd <- mapply(function(i,j) i + j[,-1], i=mapply(function(x, y) x[,-1]+y[,-1], x=pld, y=plc), j=pldc)
INTCLUST.sd <- sapply(INTCLUST.sd, function(x) apply(x, 1, function(z) sd(z)/sqrt(length(z))))

coliClust <- c("#E41A1C", "#FB9A99", "#1F78B4", "#A6CEE3", "#66A61E")
times <-  pl[[1]][,1]


plot(times, rep(0, length(times)), type="n", main="",
     ylim=c(0,1), xlim=c(0, 10), xlab="Years after loco-regional relapse",
     ylab="Probability of DR/cancer death")
for (id in 1:5) {
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

errbar(1:5, res[dd], ui[dd],
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
axis(1, at=1:5, sub(" RELAPSE", "", names(res))[dd], cex.axis=0.9, las=2)
text(c(5), c(10), c("*"), cex=0.9)
box()


times <- pdc[[1]][,1]
INTCLUST.mean <- sapply(pdc, function(x) apply(x[,-1], 1, mean))
INTCLUST.sd <- sapply(pdc, function(x) apply(x[,-1], 1, function(z) sd(z)/sqrt(length(z))))
plot(times, rep(0, length(times)), type="n", main="",
     ylim=c(0,1), xlim=c(0, 10), xlab="Years after distant relapse",
     ylab="Probability of cancer death")

for (id in 1:5) {
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


m <- survfit(Surv(T2, I2) ~ Pam50Subtype, data=TClinical)
event <- rep(NA, nrow(TClinical))
event[which(TClinical$Death==0)] <- "RELAPSE-FREE"
event[which(TClinical$I2==1)] <- "RELAPSE"
event[which(TClinical$I2==0 & TClinical$Death==1)] <- "DEATH/OTHER"
library(cmprsk)
event <-  factor(event, levels=c("RELAPSE-FREE", "RELAPSE", "DEATH/OTHER"))
m1 <- cuminc(ftime=TClinical$T2, fstatus=as.character(event),
             group=TClinical$Pam50Subtype, cencode="RELAPSE-FREE")
m2 <- survfit(Surv(TClinical$T2, as.numeric(event),type="mstate" ) ~ TClinical$Pam50Subtype)
## same result
res <- sapply(m1[6:10], function(x) {
           id <- which.min((x$est-0.5)<0)
           x$time[id]
       })
max.obs <- sapply(m1[6:10], function(x) max(x$est))
res[which(max.obs<0.5)] <- NA
ui <- lapply(m1[6:10], function(x) pmax(0, x$est - 1.96*sqrt(x$var)))
max.obs <- sapply(ui, function(x) max(x))
li <- lapply(m1[6:10], function(x) pmin(1, x$est + 1.96*sqrt(x$var)))
li <- lapply(li, function(x) {
           id <- which.min((x-0.5)<0)
       })
li <- mapply(function(x,y) y$time[x], x=li, y=m1[6:10])
ui <- lapply(ui, function(x) {
           id <- which.min((x-0.5)<0)
       })
ui <- mapply(function(x,y) y$time[x], x=ui, y=m1[6:10])
ui[which(max.obs<0.5)] <- NA

dd <- order(res)

library(rms)
ui[which(is.na(ui))] <- Inf
correct.ids <- which(ui>10 & is.finite(res))
ui[correct.ids] <- res[correct.ids]
ui[which(is.na(res))] <- NA
li[which(is.na(res))] <- NA

errbar(1:5, res[dd], ui[dd],
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
axis(1, at=1:5, sub(" RELAPSE", "", names(res))[dd], cex.axis=0.9, las=2)
box()



