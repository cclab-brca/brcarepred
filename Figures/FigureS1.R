rm(list=ls())
cols.er <- c("olivedrab", "red")
Clinical <- read.table(file="../../TableS5.txt", header=T, sep="\t", quote="", comment.char="")
library(cmprsk)

## We remove Samples with no follow-up time

ids <- which(Clinical$T==0)
if (length(ids)>0) Clinical <- Clinical[-ids,]
## We remove Samples with no follow-up time or death known

## We remove samples with stage 4
Clinical <- Clinical[-which(Clinical$Stage==4),]

## We remove Samples with no follow-up time or death known
Clinical <- Clinical[which(!is.na(Clinical$T)),]
Clinical <- Clinical[which(!is.na(Clinical$Death)),]
## We remove benign, DCIS or PHYL
bad.hist <- which(Clinical$Histological.Type %in% c("BENIGN", "PHYL"))
if (length(bad.hist)>0) Clinical <- Clinical[-bad.hist,]


pdf("FigureS1.pdf", width=12, height=5)
par(mfrow=c(1,3))
plot(survfit(Surv(T, DeathBreast) ~ ER.Status, data=Clinical),
     fun="event", xscale=365.25, cex.lab=1.3, lwd=2,
             xmax=20 *365.25, xlab="Years", ylab="Cumulative Incidence",
             col=cols.er, ylim=c(0,1), mark.time=F, main="a")
annots <- survfit(Surv(T, DeathBreast) ~ ER.Status, data=Clinical)
N <- annots$n
legend(0, 1, lwd=2, col=cols.er, legend=paste0(c("ER- 1-KM Estimate (n=",
                                          "ER+ 1-KM Estimate (n="),N, c(")", ")")) , bty="n")

library(cmprsk)
event <- rep(NA, nrow(Clinical))
event[which(Clinical$Last.Followup.Status=="a")] <- "ALIVE"
event[which(Clinical$Last.Followup.Status=="d")] <- "DEATH/UNKNOWN"
event[which(Clinical$Last.Followup.Status=="d-d.s.")] <- "DEATH/CANCER"
event[which(Clinical$Last.Followup.Status=="d-o.c.")] <- "DEATH/OTHER"
X <- factor(Clinical$ER.Status, levels=c("neg", "pos"), labels=c("ER-", "ER+"))
m1 <- cuminc(ftime=I(Clinical$T/365.25), fstatus=event,
             group=X, cencode="ALIVE")
cols.er <- c("olivedrab", "red")

plot(m1, ylab="Cumulative Incidence", cex.axis=1.3,
     color=rep(cols.er, 3), cex.lab=1.3, cex.axis=1.3, xlim=c(0, 20),
     lwd=2, xlab="Years", lty=rep(c(1, 2, 3), rep(2, 3)), main="b",
     curvlab=paste0(names(m1)[1:6], rep("(n=", 6), as.vector(t(table(event, X)[-1,])), rep(")", 6)))
box()


boxplot(Age.At.Diagnosis ~ ER.Status, data=Clinical, col=cols.er,
        names=paste0(c("ER-","ER+"), "\n(n=", table(Clinical$ER.Status), ")"), ylab="Age at diagnosis", cex.lab=1.3,
        cex.axis=1.3, main="c", axes=F)
axis(2)
axis(1, cex.axis=1.3, at=c(1, 2), labels=paste0(c("ER-","ER+"), "\n(n=", table(Clinical$ER.Status), ")"), line=2, tick=F)
box()
t.test(Age.At.Diagnosis ~ ER.Status, data=Clinical)
dev.off()


## Extra LRR

Clinical$Mast <- 1 * (Clinical$Breast.Surgery=="MASTECTOMY")
Clinical$Mast[which(Clinical$Breast.Surgery=="null")] <- NA
Clinical$Mast <- factor(Clinical$Mast, levels=c(0,1), labels=c("YES", "NO"))
Clinical$Radio <- 1 * I(1*(Clinical$RT!="null"))
Clinical$Radio <- factor(Clinical$Radio, levels=c(0,1), labels=c("YES", "NO"))
feo <- summary(survfit(Surv(TLR, LR) ~ interaction(Mast, Radio), data=Clinical), xscale=365.25)
