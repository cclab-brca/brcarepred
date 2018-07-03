
library(survival)
library(RColorBrewer)
####################################################################
####################################################################
## ALL EVENTS (When available)
####################################################################
####################################################################

library(survival)
X <- read.table(file="TableS7.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
X <- X[which(X$Stage!=4),]

## We remove benign, DCIS or PHYL
bad.hist <- which(X$Histological.Type %in% c("BENIGN", "PHYL"))
if (length(bad.hist)>0) X <- X[-bad.hist,]


X$TDR[which(X$TDR==0)] <- 0.1
X$TLR[which(X$TLR==0)] <- 0.1
X$TIME.RELAPSE[which(X$TIME.RELAPSE==0)] <- 0.1
X$T <- X$T/365.25
X$TLR <- X$TLR/365.25
X$TDR <- X$TDR/365.25
X$TIME.RELAPSE <- X$TIME.RELAPSE/365.25
X <- X[which(X$TYPE.RELAPSE=="DISTANT"),]
X <- X[order(X$METABRIC.ID, X$TIME.RELAPSE),]
X <- X[which(X$T>0),]
X <- X[which(!is.na(X$TIME.RELAPSE)),]
X$LN <- X$Lymph.Nodes.Positive
X$LN[which(X$LN>=10)] <- 10
X$AGE <- X$Age.At.Diagnosis
X$GRADE <- as.numeric(as.character(X$Grade))
X$SIZE <- as.numeric(as.character(X$Size))
X$SITE <- factor(X$SITE)
levels(X$SITE) <- list("BRAIN.MENINGEAL"=c("BRAIN", "MENINGEAL"),
                       "PULMONARY"=c("LUNG", "PLEURA"),
                       "LIVER"="LIVER",
                       "BONE"="BONE",
                       "OTHER"=c("ABDOMEN", "ADRENALS", "ASCITES","BREAST",
                           "BLADDER", "LNS", "MEDIASTINAL",
                           "EYE", "GI_TRACT", "KIDNEY", "UNSPECIFIED",
                           "OTHER", "OVARY", "PANCREAS","PERITONEUM", "SKIN",
                           "SOFT_TISSUES",
                           "PERICARDIAL"))




####################################################################
####################################################################
## Conditional (PWP) model
####################################################################
####################################################################
ID <- unique(X$METABRIC.ID)
AG <- NULL
for (i in ID) {
    tmp <- X[which(X$METABRIC.ID==i),]
    tmp <- tmp[which(!is.na(tmp$TIME.RELAPSE)),]
    time1 <- 0
    j <- 1
    while(j <= nrow(tmp)) {
        time2 <- tmp$TIME.RELAPSE[j]
        if (time1>=time2) time2 <- time1 + 0.01
        feo <- data.frame(ID=i, time1=time1, time2=time2,
                          status=1*(!is.na(tmp$SITE[j])),
                          ER=tmp$ER.Status[1],
                              Pam50=tmp$Pam50Subtype[1], iC10=tmp$iC10[1],
                          enum=j)
        AG <- rbind(AG, feo)
        time1 <- time2
        j <- j + 1
    }
    if (time1 < tmp$T[j-1]) {
        feo <- data.frame(ID=i, time1=time1, time2=tmp$T[j-1],
                          status=0,
                          ER=tmp$ER.Status[1],
                              Pam50=tmp$Pam50Subtype[1], iC10=tmp$iC10[1],
                          enum=j)
        AG <- rbind(AG, feo)
    }
}
AG$T <- (AG$time2 - AG$time1)/365.25



pdf("Figure6.pdf", width=10, height=10)
layout(matrix(1:6, ncol=2, byrow=T), height=c(1, 5, 5))
par(mar=c(5, 5, 5, 2), oma=c(3, 5, 3, 1))
medt <- matrix(NA, 4, 2)
feo <- survfit(Surv(time2, status) ~ strata(enum), data=AG, subset=ER=="pos")
medt[,1] <- survival:::survmean(feo, scale=1, "common")$matrix[1:4,'median']

feo <- survfit(Surv(time2, status) ~ strata(enum), data=AG, subset=ER=="neg")
medt[,2] <- survival:::survmean(feo, scale=1, "common")$matrix[1:4,'median']
rownames(medt) <- c("1st r", "2nd r",
                            "3rd r", "4th r")
Gr <- rev(brewer.pal(9, "Greens"))
par(mar=c(2, 2, 1, 1))
tt <- barplot(matrix(diff(c(0, medt[,1,drop=F])), ncol=1), beside=F, las=2, ylab="Years until event", horiz=T, col=Gr, xlim=c(0, 16), cex.main=2,
              cex.lab=1.5, cex.axis=1.5)
mtext("ER+", side=3, line=1, cex=1.4)
start <- c(0, cumsum(diff(c(0, medt[,1]))))[1:4]
end <- cumsum(diff(c(0, medt[,1])))[1:4]
text(start + (end-start)/2, 0.65, labels=rownames(medt), srt=90, cex=1.5,
     col="grey")
tt <- barplot(matrix(diff(c(0, medt[,2,drop=F])), ncol=1), beside=F, las=2, ylab="Years until event", horiz=T, main="" ,col=Gr, xlim=c(0, 16), cex.main=2,
              cex.lab=1.5, cex.axis=1.5)
mtext("ER-", side=3, line=1, cex=1.4)
start <- c(0, diff(c(0, cumsum(medt[,2]))))[1:4]
end <- diff(c(0, cumsum(medt[,2])))[1:4]
text(start + (end-start)/2, 0.65, labels=rownames(medt), srt=90, cex=1.5,
     col="grey")

mtext("a", side=3, outer=T, line=1, cex=2)

Gr <- rev(brewer.pal(9, "Greens"))
plot(survfit(Surv(time2, status) ~ strata(enum), data=AG, subset=ER=="pos")[1:4], col=Gr[1:5], mark.time=FALSE, main="", conf.int=FALSE, xmax=20, cex.lab=1.5, cex.axis=1.5,     xlab="Years", ylab="Subsequent Relapse-free survival", lwd=2)
mtext(side=2, "Subsequent Relapse-free survival", cex=1.4, line=3)
plot(survfit(Surv(time2, status) ~ strata(enum), data=AG, subset=ER=="neg")[1:4], col=Gr[1:5], mark.time=FALSE, main="", conf.int=FALSE, xmax=20,cex.lab=1.5, cex.axis=1.5,
     xlab="Years", ylab="", lwd=2)

legend("topright", col=Gr[1:4], lwd=2,
       legend=c("1st recurrence", "2nd recurrence",
           "3rd recurrence", "4th recurrence"),
           bty="n", cex=1.5)

## Time-dependent Cox model
## Prognosis

par(mar=c(4, 10, 2, 2))

CR2 <- NULL
Met <- data.frame("LOCAL.REGIONAL"=0, "BRAIN/MENINGEAL"=0, "PULMONARY"=0,
                  "LIVER"=0,
                  "BONE"=0, "OTHER"=0)

for (i in ID) {
    tmp <- X[which(X$METABRIC.ID==i),]
    time1 <- 0
    enum <- 0
    tts <- unique(tmp$TIME.RELAPSE)
    for (j in tts) {
            if (j==0) {
                events <- as.character(tmp$SITE[which(tmp$TIME.RELAPSE==j)])
                for (ev in events)  Met[,ev] <- 1
            } else {
                feo <- data.frame(ID=i, time1=time1,
                              time2=j,
                              Outcome="CANCERDEATH",
                              status=0, iC10=tmp$iC10[1],
                              ER=tmp$ER.Status[1], LN=tmp$LN[1],
                              AGE=tmp$AGE[1], GRADE=tmp$GRADE[1],
                              SIZE=tmp$SIZE[1])
            feo <- cbind(feo, Met)
            feo$enum <- enum
            CR2 <- rbind(CR2, feo)
            feo <- data.frame(ID=i, time1=time1,
                              time2=j,
                              Outcome="OTHERDEATH",
                              status=0,iC10=tmp$iC10[1],
                              ER=tmp$ER.Status[1], LN=tmp$LN[1],
                              AGE=tmp$AGE[1], GRADE=tmp$GRADE[1],
                              SIZE=tmp$SIZE[1])
            feo <- cbind(feo, Met)
            time1 <- j
            events <- as.character(tmp$SITE[which(tmp$TIME.RELAPSE==j)])
            for (ev in events)  Met[,ev] <- 1
            feo$enum <- enum
            CR2 <- rbind(CR2, feo)
            enum <- enum + length(events)
            }
            feo <- data.frame(ID=i, time1=time1,
                          time2=tmp$T[1],
                          Outcome="CANCERDEATH",
                              status=1 * (tmp$DeathBreast[1]==1),
                              iC10=tmp$iC10[1],
                          ER=tmp$ER.Status[1], LN=tmp$LN[1],
                          AGE=tmp$AGE[1], GRADE=tmp$GRADE[1],
                          SIZE=tmp$SIZE[1])
            feo <- cbind(feo, Met)
            feo$enum <- enum
            CR2 <- rbind(CR2, feo)
            feo <- data.frame(ID=i, time1=time1,
                          time2=tmp$T[1],
                          Outcome="OTHERDEATH",
                              status=1 * (tmp$DeathBreast[1]==0 & tmp$Death[1]==1),
                              iC10=tmp$iC10[1],
                          ER=tmp$ER.Status[1], LN=tmp$LN[1],
                          AGE=tmp$AGE[1], GRADE=tmp$GRADE[1],
                          SIZE=tmp$SIZE[1])
            feo <- cbind(feo, Met)
            feo$enum <- enum
            CR2 <- rbind(CR2, feo)
            Met <- data.frame("LOCAL.REGIONAL"=0, "BRAIN.MENINGEAL"=0,
                              "PULMONARY"=0,
                              "LIVER"=0,
                              "BONE"=0, "OTHER"=0)
        }
}




CR2$AGE.OT <- CR2$AGE * (CR2$Outcome=="OTHERDEATH")
CR2$AGE.CD <- CR2$AGE * (CR2$Outcome=="CANCERDEATH")
CR2$ER.CD <- (CR2$ER=="neg") * 1*(CR2$Outcome=="CANCERDEATH")
CR2$ER.OT <- (CR2$ER=="neg") * 1*(CR2$Outcome=="OTHERDEATH")
CR2$LN.CD <- CR2$LN * (CR2$Outcome=="CANCERDEATH")
CR2$GRADE.CD <- CR2$GRADE * (CR2$Outcome=="CANCERDEATH")
CR2$SIZE.CD <- CR2$SIZE * (CR2$Outcome=="CANCERDEATH")
CR2$LOCAL.REGIONAL.CD <- CR2$LOCAL.REGIONAL * (CR2$Outcome=="CANCERDEATH")
CR2$BRAIN.MENINGEAL.CD <- CR2$BRAIN.MENINGEAL * (CR2$Outcome=="CANCERDEATH")
CR2$PULMONARY.CD <- CR2$PULMONARY * (CR2$Outcome=="CANCERDEATH")
CR2$LIVER.CD <- CR2$LIVER * (CR2$Outcome=="CANCERDEATH")
CR2$BONE.CD <- CR2$BONE * (CR2$Outcome=="CANCERDEATH")
CR2$OTHER.CD <- CR2$OTHER * (CR2$Outcome=="CANCERDEATH")
CR2$enum.CD <- CR2$enum * (CR2$Outcome=="CANCERDEATH")
CR2$iC10.CD <- CR2$iC10
CR2$iC10.CD[which(CR2$Outcome!="CANCERDEATH")] <- NA
CR2$LN.OT <- CR2$LN * (CR2$Outcome=="OTHERDEATH")
CR2$GRADE.OT <- CR2$GRADE * (CR2$Outcome=="OTHERDEATH")
CR2$SIZE.OT <- CR2$SIZE * (CR2$Outcome=="OTHERDEATH")
CR2$LOCAL.REGIONAL.OT <- CR2$LOCAL.REGIONAL * (CR2$Outcome=="OTHERDEATH")
CR2$BRAIN.MENINGEAL.OT <- CR2$BRAIN.MENINGEAL * (CR2$Outcome=="OTHERDEATH")
CR2$PULMONARY.OT <- CR2$PULMONARY * (CR2$Outcome=="OTHERDEATH")
CR2$LIVER.OT <- CR2$LIVER * (CR2$Outcome=="OTHERDEATH")
CR2$BONE.OT <- CR2$BONE * (CR2$Outcome=="OTHERDEATH")
CR2$OTHER.OT <- CR2$OTHER * (CR2$Outcome=="OTHERDEATH")
CR2$enum.OT <- CR2$enum * (CR2$Outcome=="OTHERDEATH")

rownames(CR2) <- paste("ID", 1:nrow(CR2), sep="")

m3 <- coxph(Surv(time1, time2, status) ~ LN.CD + GRADE.CD + SIZE.CD +
                        strata(Outcome) + cluster(ID), data=CR2)



feo <- subset(CR2, CR2$ER == "pos")


m3 <- coxph(Surv(time1, time2, status) ~ LN.CD + GRADE.CD + SIZE.CD + enum.CD +
                BRAIN.MENINGEAL.CD + PULMONARY.CD +
                    LIVER.CD + BONE.CD + OTHER.CD +
                        strata(Outcome) + cluster(ID), data=feo)
cox.zph(m3)

ids <- order(coef(m3))
coef.labels <- names(coef(m3))
coef.labels <- sub(".CD","", coef.labels, fixed=TRUE)
coef.labels <- sub(".OT","", coef.labels, fixed=TRUE)
coef.labels[which(coef.labels=="enum")] <- "#RELAPSES"
coef.labels[which(coef.labels=="LOCAL.REGIONAL")] <- "L/R"
coef.labels[which(coef.labels=="BRAIN.MENINGEAL")] <- "BRAIN/MENINGEAL"
K <- length(coef(m3))
YLIM <- c(1, K)
## XLIM <- c(0, 10.5)
XLIM <- c(-0.5, 3)
plot(0,0, type="n", axes=F, xlab="", ylab="", xlim=XLIM,
     ylim=YLIM, xaxs="i", main="", cex.axis=1.7, cex.lab=1.7)
mtext("Log-Hazard Ratio\nER+", side=1, line=4)
## points(exp(coef(m3))[ids], 1:K, pch=19)
points(coef(m3)[ids], 1:K, pch=19, cex=1.5)
## confs <- exp(confint(m3))[ids,]
confs <- confint(m3)[ids,]
points(confs[,1], 1:K, pch="[", cex=1.5)
points(confs[,2], 1:K, pch="]", cex=1.5)
for (i in 1:K) lines(confs[i,], c(i,i))
abline(h=1:K, lty=2, col="grey")
abline(v=0, lty=2, col="grey")
over <- which(confs[,2]>XLIM[2])
if (length(over)>0) {
    arrows(XLIM[2], over, x1=XLIM[2]-0.1, code=1, length=0.1, cex=1.5)
    }
axis(1, cex=1.5)
axis(2, at=1:K, coef.labels[ids], las=2, cex.axis=1.5)

feo <- subset(CR2, CR2$ER == "neg")


m3 <- coxph(Surv(time1, time2, status) ~ LN.CD + GRADE.CD + SIZE.CD + enum.CD +
                BRAIN.MENINGEAL.CD + PULMONARY.CD +
                    LIVER.CD + BONE.CD + OTHER.CD +
                        strata(Outcome) + cluster(ID), data=feo)

## tmp <- cox.zph(m3)
## plot(tmp[2])
## plot(tmp[4])

coef.labels <- names(coef(m3))
coef.labels <- sub(".CD","", coef.labels, fixed=TRUE)
coef.labels <- sub(".OT","", coef.labels, fixed=TRUE)
coef.labels[which(coef.labels=="enum")] <- "#RELAPSES"
coef.labels[which(coef.labels=="LOCAL.REGIONAL")] <- "L/R"
K <- length(coef(m3))
YLIM <- c(1, K)

par(mar=c(4, 2, 2, 10))

plot(0,0, type="n", axes=F, xlab="", ylab="", xlim=XLIM,
     ylim=YLIM, xaxs="i", main="", cex.lab=1.7, cex.axis=1.7)
mtext("Log-Hazard Ratio\nER-", side=1, line=4)
## points(exp(coef(m3))[ids], 1:K, pch=19)
points(coef(m3)[ids], 1:K, pch=19, cex=1.5)
## confs <- exp(confint(m3))[ids,]
confs <- confint(m3)[ids,]
points(confs[,1], 1:K, pch="[", cex=1.5)
points(confs[,2], 1:K, pch="]", cex=1.5)
for (i in 1:K) lines(confs[i,], c(i,i))
abline(h=1:K, lty=2, col="grey")
abline(v=0, lty=2, col="grey")
over <- which(confs[,2]>XLIM[2])
if (length(over)>0) {
    arrows(XLIM[2], over, x1=XLIM[2]-0.1, code=1, length=0.1, cex=1.5)
    }
axis(1, cex=1.5)
axis(2, at=1:K, rep("", K), las=2, cex.axis=1.5)

mtext("b", side=1, outer=T, line=1, cex=2)

dev.off()



m3 <- coxph(Surv(time1, time2, status) ~ ER.CD + ER.OT + (LN.CD + GRADE.CD + SIZE.CD + enum.CD +
                BRAIN.MENINGEAL.CD + PULMONARY.CD +
                    LIVER.CD + BONE.CD + OTHER.CD) +
                        strata(Outcome) + cluster(ID), data=CR2)
cox.zph(m3)

m3 <- coxph(Surv(time1, time2, status) ~ ER.CD*(LN.CD + GRADE.CD + SIZE.CD + enum.CD +
                BRAIN.MENINGEAL.CD + PULMONARY.CD +
                    LIVER.CD + BONE.CD + OTHER.CD) +
                        strata(Outcome) + cluster(ID), data=CR2)
cox.zph(m3)
