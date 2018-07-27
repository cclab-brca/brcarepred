rm(list=ls())
## Top panel
library(Hmisc)
pdf("FigureS6a.pdf", width=9, height=4)
load(file="../Models/FourGroups_AllProbs.RData")
library(survival)

lr <- mapply(function(x,y) x[,-1] + y[,-1], x=psl, y=pslc)
lr <- mapply(function(x,y) x + y[,-1], x=lr, y=psld)
lr <- mapply(function(x,y) x + y[,-1], x=lr, y=psldc)
lr <- mapply(function(x,y) x + y[,-1], x=lr, y=psldo)
lr <- mapply(function(x,y) x + y[,-1], x=lr, y=pslo)

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

coliClust <- c("lightskyblue", "lightsalmon", "plum2", "wheat4")
names(coliClust) <- c("ER+/HER2+", "ER+/HER2-", "ER-/HER2+", "ER-/HER2-")
ss <- paste0("n=", sapply(lr, function(x) ncol(x)))
names(ss) <- names(lr)
n <- ncol(lr.INTCLUST.mean)
res <- data.frame(X=c(lr.INTCLUST.mean, dr.INTCLUST.mean), Year=rep(times, n*2), Relapse=rep(c("Loco-regional Relapse",
                                                                                      "Distant Relapse"), c(n*5, n*5)),
                  IntClust=rep(rep(colnames(lr.INTCLUST.mean), rep(5, n)), 2), li=c(lr.INTCLUST.mean, dr.INTCLUST.mean) - 1.96*c(lr.INTCLUST.SE, dr.INTCLUST.SE), ui=c(lr.INTCLUST.mean, dr.INTCLUST.mean) + 1.96*c(lr.INTCLUST.SE, dr.INTCLUST.SE))
res$pch <- 19
res$pch[which(res$Relapse=="Distant Relapse")] <- 17
res$colr <- "olivedrab"
res$colr[which(res$Relapse=="Distant Relapse")] <- "black"
res$IntClust <- factor(res$IntClust, levels=c("ER+/HER2-", "ER-/HER2-", "ER+/HER2+", "ER-/HER2+"))
ss <- ss[levels(res$IntClust)]
coliClust <- coliClust[levels(res$IntClust)]
pp <- xyplot(X ~ Year| IntClust, groups=Relapse, data=res,
             layout=c(n, 1, 1), ylab="Probability of relapse", xlab="Years after surgery", ylim=c(0, 0.65),
             par.settings = list(superpose.symbol = list(pch=c(17, 19), col=c("black", "olivedrab")),
                 superpose.line = list(
                     col = rep(c("black", "olivedrab"), 4),
                     lwd = c(1.5, 1.5),
                     lty = c(1, 1))
                                 ),
             panel=function(x,y,col, pch,...,subscripts) {
                 panel.fill(col=adjustcolor(coliClust[panel.number()],
                                alpha=0.4))
                 panel.xyplot(x,y, lwd=1.5, col=res$col[subscripts], pch=res$pch[subscripts])
                 panel.arrows(x, y, x, res$ui[subscripts], length = 0.05, lwd=0.5,
                      angle = 90, col=res$col[subscripts])
                 panel.arrows(x, y, x, res$li[subscripts], length = 0.05, lwd=0.5,
                      angle = 90, col=res$col[subscripts])
                 panel.text(10, 0.6, labels=ss[panel.number()], outer=T, cex=0.75)
             }, auto.key=T)
print(pp)
dev.off()



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## Bottom panel
rm(list=ls())
pdf("FigureS6b.pdf", width=9, height=3)
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


dev.off()
