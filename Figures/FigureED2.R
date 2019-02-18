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

## Panel a
par(mfrow=c(1,3))

plot(survfit(Surv(T, DeathBreast) ~ ER.Status, data=Clinical),
     fun="event", xscale=365.25, cex.lab=1.3, lwd=2,
             xmax=20 *365.25, xlab="Years", ylab="Cumulative Incidence",
             col=cols.er, ylim=c(0,1), mark.time=F, main="b")
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
     lwd=2, xlab="Years", lty=rep(c(1, 2, 3), rep(2, 3)), main="c",
     curvlab=paste0(names(m1)[1:6], rep("(n=", 6), as.vector(t(table(event, X)[-1,])), rep(")", 6)))
box()

boxplot(Age.At.Diagnosis ~ ER.Status, data=Clinical, col=cols.er,
        names=paste0(c("ER-","ER+"), "\n(n=", table(Clinical$ER.Status), ")"), ylab="Age at diagnosis", cex.lab=1.3,
        cex.axis=1.3, main="a", axes=F)
axis(2)
axis(1, cex.axis=1.3, at=c(1, 2), labels=paste0(c("ER-","ER+"), "\n(n=", table(Clinical$ER.Status), ")"), line=2, tick=F)
box()
t.test(Age.At.Diagnosis ~ ER.Status, data=Clinical)

## Second panel


library(mstate)
library(lattice)
load(file="../Models/ERM.RData")
fm <- fm2

model.cofs <- data.frame(Name=names(coef(fm)), HR=summary(fm)$conf.int[,1], LI=summary(fm)$conf.int[,3], UI=summary(fm)$conf.int[,4])
model.cofs$Effect <- sapply(strsplit(as.character(model.cofs$Name), ".", fixed=TRUE), function(x) x[1])
model.cofs$Time <- sapply(strsplit(as.character(model.cofs$Name), ".", fixed=TRUE), function(x) x[2])
model.cofs$ER <- sapply(strsplit(as.character(model.cofs$Name), ".", fixed=TRUE), function(x) x[3])
model.cofs$Time <- factor(model.cofs$Time, levels=c("PS", "LR", "DR"))
model.cofs$ER <- factor(model.cofs$ER, levels=c("NEG", "POS"), labels=c("ER-", "ER+"))
model.cofs <- model.cofs[which(!is.na(model.cofs$ER)),]
model.cofs$Effect[which(model.cofs$Effect=="TLastSurgery")] <- "Time from Surgery"
model.cofs$Effect[which(model.cofs$Effect=="TLastLocal")] <- "Time from LR"

## First plot: effect on time of each clinical parameter?
prepanel.ci <- function(x, y, lx, ux, subscripts, ...)
{
    x <- as.numeric(y)
    lx <- as.numeric(lx[subscripts])
    ux <- as.numeric(ux[subscripts])
    print(range(x, ux, lx, finite = TRUE))
    list(ylim = range(x, ux, lx, finite = TRUE))
}


panel.ci <- function(x, y, lx, ux, subscripts, pch = 16, ...)
{
    x <- as.numeric(x)
    y <- as.numeric(y)
    lx <- as.numeric(lx[subscripts])
    ux <- as.numeric(ux[subscripts])
    panel.dotplot(x, y, pch = pch, col.line="white", ...)
    panel.abline(h=0, col="grey")
    panel.arrows(x, lx, x, ux, col = 'black',
                 length = 0.25, unit = "native",
                 angle = 90, code = 3, lwd=0.5)

}
model.cofs$Effect <- factor(model.cofs$Effect, levels=c("GRADE", "LN", "SIZE", "Time from Surgery", "Time from LR"))
print(with(model.cofs,
     dotplot(log(HR) ~ Time | Effect*ER,
            lx = log(LI), ux = log(UI),
            prepanel = prepanel.ci,
             panel = panel.ci, ylab="Log-hazard Ratio")))



