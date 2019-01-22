rm(list=ls())
library(RColorBrewer)

cols <- c("black", "lightblue", "darkblue", "red", "olivedrab")
cols <- rep(cols, 2)
all.X <- list()
coliClust <- c('#FF5500', '#00EE76', '#CD3278','#00C5CD', '#B5D0D2', '#8B0000',
               '#FFFF40', '#0000CD', '#FFAA00', '#EE82EE', '#7D26CD')
cols.er <- c("olivedrab", "red")
library(mstate)
library(survival)
library(lattice)
##################################
##################################
## Common Model for all subtypes
##################################
##################################


##################################
##################################
## FD
##################################
##################################

Clinical <- read.table(file="../../TableS5.txt", header=T, sep="\t", quote="", comment.char="", stringsAsFactors=FALSE)
Clinical$NaturalDeath <- 1 * (Clinical$Last.Followup.Status %in% c("d", "d-o.c."))

## We remove Samples with no follow-up time

ids <- which(Clinical$T==0)
if (length(ids)>0) Clinical <- Clinical[-ids,]

## We remove samples with stage 4
Clinical <- Clinical[-which(Clinical$Stage==4),]


## We remove Samples with no follow-up time or death known

Clinical <- Clinical[which(!is.na(Clinical$T)),]
Clinical <- Clinical[which(!is.na(Clinical$Death)),]

## We remove benign, DCIS or PHYL
bad.hist <- which(Clinical$Histological.Type %in% c("BENIGN", "PHYL"))
if (length(bad.hist)>0) Clinical <- Clinical[-bad.hist,]

dim(na.omit(Clinical[,c('ER.Status', 'Size', 'Grade', 'Lymph.Nodes.Positive', 'T', 'LR', 'DR', 'TLR', 'TDR', 'Death', 'DeathBreast')]))

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
Clinical$HT <- 1 * (Clinical$HT!="NO/NA")
Clinical$HT <- factor(Clinical$HT, levels=c(0, 1), labels=c("NO", "YES"))
Clinical$CT <- 1 * (Clinical$CT!="NO/NA")
Clinical$CT <- factor(Clinical$CT, levels=c(0, 1), labels=c("NO", "YES"))
Clinical$RT[which(Clinical$RT %in% c("NO/NA", "NONE RECORDED IN LANTIS", "Nnne"))] <- 0
Clinical$RT[which(Clinical$RT!=0)] <- 1
Clinical$RT <- factor(Clinical$RT, levels=c(0, 1), labels=c("NO", "YES"))
Clinical$BS <- factor(Clinical$Breast.Surgery, levels=c("BREAST CONSERVING", "MASTECTOMY"), labels=c("BC", "M"))
Clinical$ER.Status <- factor(Clinical$ER.Status,
                             levels=c("neg", "pos"), labels=c("ER-", "ER+"))

sC <- survfit(Surv(TLR, LR) ~ ER.Status, data=Clinical)
sC2 <- subset(Clinical, (Clinical$LR==1 | Clinical$T >= 15))
round(prop.table(table(sC2$LR), ), 2)
round(prop.table(table(sC2$LR, sC2$ER.Status), 2), 2)
tapply(sC2$TLR[which(sC2$LR==1)], sC2$ER.Status[which(sC2$LR==1)], summary)
sC3 <- subset(Clinical, (Clinical$DR==1 | Clinical$DeathBreast==1 | Clinical$T >= 15))
sC3$Event <- 1 * (sC3$DR==1 | sC3$DeathBreast==1)
round(prop.table(table(sC3$Event), ), 2)
round(prop.table(table(sC3$Event, sC3$ER.Status), 2), 2)
tapply(pmin(sC3$TDR[which(sC3$Event==1)], sC3$T[which(sC3$Event==1)]), sC3$ER.Status[which(sC3$Event==1)], summary)

sC2 <- subset(Clinical, Clinical$LR==1)
tapply(sC2$TLR, sC2$ER.Status, summary)
round(prop.table(table(sC2$DR, sC2$ER.Status), 2), 2)
sC2 <- subset(Clinical, (Clinical$DR==1))
tapply(sC2$TDR, sC2$ER.Status, summary)

sC2 <- subset(Clinical, (Clinical$LR==1 & (Clinical$DeathBreast==1 | Clinical$DR==1)))

tapply(pmin(I(sC2$T-sC2$TLR), I(sC2$TDR-sC2$TLR)), sC2$ER.Status, summary)

## Transition matrix
tra <- matrix(NA, 5*2, 5*2)

state.names <- paste(c("PostSurgery",
                            "LocalRelapse", "DistantRelapse",
                       "CancerDeath", "NaturalDeath"),
                     rep(levels(Clinical$ER), rep(5, 2)), sep="_")
rownames(tra) <- state.names
colnames(tra) <- rownames(tra)
counter <- 1
for (i in seq(from=1, by=5, length=2)) {
    tra[i,(i+1):(i+4)] <- counter:(counter+3)
    counter <- counter + 4
    tra[i+1,(i+2):(i+4)] <- counter:(counter+2)
    counter <- counter + 3
    tra[i+2,(i+3):(i+4)] <- counter:(counter+1)
    counter <- counter + 2
}


covs <- c("LN", "AGE", "GRADE", "SIZE", "BS")

ALL.CLINICAL <- Clinical


set.seed(12324)
all.X <- list()
for (ic in levels(ALL.CLINICAL$ER.Status)) {
    Clinical <- ALL.CLINICAL[which(ALL.CLINICAL$ER.Status==ic),]
    X <- NULL
    for (i in 1:nrow(Clinical)) {
    ## Surgery to Local
    tmp <- data.frame(id=Clinical$METABRIC.ID[i],
                      from=1,
                      to=2,
                      trans=1,
                      Tstart=0,
                      Tstop=min(Clinical$TLR[i], Clinical$TDR[i]),
                      status=Clinical$LR[i],
                      LN=Clinical$LN[i],
                      AGE=Clinical$AGE[i],
                      GRADE=Clinical$GRADE[i],
                      SIZE=Clinical$SIZE[i],
                      BS=Clinical$BS[i],
                      ER=Clinical$ER.Status[i],
                      TLastSurgery=0,
                      TLastLocal=0)

    X <- rbind(X, tmp)
    ## Surgery to Distant
    tmp <- data.frame(id=Clinical$METABRIC.ID[i],
                      from=1,
                      to=3,
                      trans=2,
                      Tstart=0,
                      Tstop=min(Clinical$TLR[i], Clinical$TDR[i]),
                      status=1 * (Clinical$LR[i]==0) * Clinical$DR[i],
                      LN=Clinical$LN[i],
                      AGE=Clinical$AGE[i],
                      GRADE=Clinical$GRADE[i],
                      SIZE=Clinical$SIZE[i],
                      BS=Clinical$BS[i],
                      ER=Clinical$ER.Status[i],
                      TLastSurgery=0,
                      TLastLocal=0)
    X <- rbind(X, tmp)
    ## Surgery to Cancer Death
    tmp <- data.frame(id=Clinical$METABRIC.ID[i],
                      from=1,
                      to=4,
                      trans=3,
                      Tstart=0,
                      Tstop=min(Clinical$TLR[i], Clinical$TDR[i]),
                      status=1 * (Clinical$LR[i]==0 & Clinical$DR[i]==0) *
                          Clinical$DeathBreast[i],
                      LN=Clinical$LN[i],
                      AGE=Clinical$AGE[i],
                      GRADE=Clinical$GRADE[i],
                      SIZE=Clinical$SIZE[i],
                      BS=Clinical$BS[i],
                      ER=Clinical$ER.Status[i],
                      TLastSurgery=0,
                      TLastLocal=0)
    X <- rbind(X, tmp)
    ## Surgery to Natural Death
    tmp <- data.frame(id=Clinical$METABRIC.ID[i],
                      from=1,
                      to=5,
                      trans=4,
                      Tstart=0,
                      Tstop=min(Clinical$TLR[i], Clinical$TDR[i]),
                      status=1 * (Clinical$LR[i]==0 & Clinical$DR[i]==0) *
                          Clinical$NaturalDeath[i],
                      LN=Clinical$LN[i],
                      AGE=Clinical$AGE[i],
                      GRADE=Clinical$GRADE[i],
                      SIZE=Clinical$SIZE[i],
                      BS=Clinical$BS[i],
                      ER=Clinical$ER.Status[i],
                      TLastSurgery=0,
                      TLastLocal=0)
    X <- rbind(X, tmp)
    ## Local to distant
    if (Clinical$LR[i]==1) {
        tmp <- data.frame(id=Clinical$METABRIC.ID[i],
                          from=2,
                          to=3,
                          trans=5,
                          Tstart=Clinical$TLR[i],
                          Tstop=min(Clinical$TDR[i], Clinical$T[i]),
                          status=Clinical$DR[i],
                          LN=Clinical$LN[i],
                          AGE=Clinical$AGE[i] + Clinical$TLR[i],
                          GRADE=Clinical$GRADE[i],
                          SIZE=Clinical$SIZE[i],
                          BS=Clinical$BS[i],
                          ER=Clinical$ER.Status[i],
                          TLastSurgery=Clinical$TLR[i],
                          TLastLocal=0)
        X <- rbind(X, tmp)
        ## Local to cancer Death
        tmp <- data.frame(id=Clinical$METABRIC.ID[i],
                          from=2,
                          to=4,
                          trans=6,
                          Tstart=Clinical$TLR[i],
                          Tstop=min(Clinical$TDR[i], Clinical$T[i]),
                          status=1 * (Clinical$DR[i]==0) * Clinical$DeathBreast[i],
                          LN=Clinical$LN[i],
                          AGE=Clinical$AGE[i] + Clinical$TLR[i],
                          GRADE=Clinical$GRADE[i],
                          SIZE=Clinical$SIZE[i],
                          BS=Clinical$BS[i],
                          ER=Clinical$ER.Status[i],
                          TLastSurgery=Clinical$TLR[i],
                          TLastLocal=0)
        X <- rbind(X, tmp)
        ## Local to Natural Death
        tmp <- data.frame(id=Clinical$METABRIC.ID[i],
                          from=2,
                          to=5,
                          trans=7,
                          Tstart=Clinical$TLR[i],
                          Tstop=min(Clinical$TDR[i], Clinical$T[i]),
                          status=1 * (Clinical$DR[i]==0) * Clinical$NaturalDeath[i],
                          LN=Clinical$LN[i],
                          AGE=Clinical$AGE[i] + Clinical$TLR[i],
                          GRADE=Clinical$GRADE[i],
                          SIZE=Clinical$SIZE[i],
                          BS=Clinical$BS[i],
                          ER=Clinical$ER.Status[i],
                          TLastSurgery=0,
                          TLastLocal=0)
        X <- rbind(X, tmp)
    }
    if (Clinical$DR[i]==1) {
        if (Clinical$LR[i]==1) {
            TLastSurgery <- Clinical$TDR[i]
            TLastLocal <- Clinical$TDR[i] - Clinical$TLR[i]
        } else {
            TLastSurgery <- Clinical$TDR[i]
            TLastLocal <- 0
        }
        ## Distant to Cancer Death
        tmp <- data.frame(id=Clinical$METABRIC.ID[i],
                          from=3,
                          to=4,
                          trans=8,
                          Tstart=Clinical$TDR[i],
                          Tstop=Clinical$T[i],
                      status=Clinical$DeathBreast[i],
                          LN=Clinical$LN[i],
                          AGE=Clinical$AGE[i] + Clinical$TDR[i],
                          GRADE=Clinical$GRADE[i],
                          SIZE=Clinical$SIZE[i],
                          BS=Clinical$BS[i],
                      ER=Clinical$ER.Status[i],
                          TLastSurgery=TLastSurgery,
                          TLastLocal=TLastLocal)
        X <- rbind(X, tmp)
        ## Distant to Natural Death
        tmp <- data.frame(id=Clinical$METABRIC.ID[i],
                          from=3,
                          to=5,
                          trans=9,
                          Tstart=Clinical$TDR[i],
                          Tstop=Clinical$T[i],
                          status=Clinical$NaturalDeath[i],
                          LN=Clinical$LN[i],
                          AGE=Clinical$AGE[i] + Clinical$TDR[i],
                          GRADE=Clinical$GRADE[i],
                          SIZE=Clinical$SIZE[i],
                          BS=Clinical$BS[i],
                      ER=Clinical$ER.Status[i],
                          TLastSurgery=0,
                          TLastLocal=0)
        X <- rbind(X, tmp)
    }
    cat("Sample ", Clinical$METABRIC.ID[i], " done\n")
}
    all.X[[ic]] <- X
}

res <- all.X

for (i in 1:length(res)) {
    res[[i]]$from <- res[[i]]$from + 5 * (i-1)
    res[[i]]$to <- res[[i]]$to + 5 * (i-1)
    res[[i]]$trans <- res[[i]]$trans + 9 * (i-1)
}

res <- do.call("rbind", res)

res$time <- res$Tstop-res$Tstart
 res <- res[,c(1:6, 16,7:15)]

class(res) <- c("msdata", "data.frame")
attr(res, "trans") <- tra
library(mstate)
tmp <- events(res)

nat.death <- as.numeric((tra[,grep("NaturalDeath", colnames(tra))]))
nat.death <- unique(na.omit(nat.death))
res$AGE[which(!res$trans %in% nat.death)] <- 0
res$GRADE[which(res$trans %in% nat.death)] <- 0
res$SIZE[which(res$trans %in% nat.death)] <- 0
res$LN[which(res$trans %in% nat.death)] <- 0

res$AGE.NEG <- res$AGE * (res$ER=="ER-")
res$AGE.POS <- res$AGE * (res$ER=="ER+")
res$GRADE.NEG <- res$GRADE * (res$ER=="ER-")
res$GRADE.POS <- res$GRADE * (res$ER=="ER+")
res$SIZE.NEG <- res$SIZE * (res$ER=="ER-")
res$SIZE.POS <- res$SIZE * (res$ER=="ER+")
res$LN.NEG <- res$LN * (res$ER=="ER-")
res$LN.POS <- res$LN * (res$ER=="ER+")
res$TLastSurgery.NEG <- res$TLastSurgery * (res$ER=="ER-")
res$TLastSurgery.POS <- res$TLastSurgery * (res$ER=="ER+")
res$TLastLocal.DR.NEG <- res$TLastLocal * (res$ER=="ER-")
res$TLastLocal.DR.POS <- res$TLastLocal * (res$ER=="ER+")

res$AGE.PS <- res$AGE * (res$from %in% grep("PostSurgery", colnames(tra)))
res$AGE.LR <- res$AGE * (res$from %in% grep("LocalRelapse", colnames(tra)))
res$AGE.DR <- res$AGE * (res$from %in% grep("DistantRelapse", colnames(tra)))
res$GRADE.PS <- res$GRADE * (res$from %in% grep("PostSurgery", colnames(tra)))
res$GRADE.LR <- res$GRADE * (res$from %in% grep("LocalRelapse", colnames(tra)))
res$GRADE.DR <- res$GRADE * (res$from %in% grep("DistantRelapse", colnames(tra)))
res$SIZE.PS <- res$SIZE * (res$from %in% grep("PostSurgery", colnames(tra)))
res$SIZE.LR <- res$SIZE * (res$from %in% grep("LocalRelapse", colnames(tra)))
res$SIZE.DR <- res$SIZE * (res$from %in% grep("DistantRelapse", colnames(tra)))
res$LN.PS <- res$LN * (res$from %in% grep("PostSurgery", colnames(tra)))
res$LN.LR <- res$LN * (res$from %in% grep("LocalRelapse", colnames(tra)))
res$LN.DR <- res$LN * (res$from %in% grep("DistantRelapse", colnames(tra)))
res$TLastSurgery.LR <- res$TLastSurgery * (res$from %in% grep("LocalRelapse", colnames(tra)))
res$TLastSurgery.DR <- res$TLastSurgery * (res$from %in% grep("DistantRelapse", colnames(tra)))

res$AGE.PS.NEG <- res$AGE.NEG * (res$from %in% grep("PostSurgery", colnames(tra)))
res$AGE.PS.POS <- res$AGE.POS * (res$from %in% grep("PostSurgery", colnames(tra)))
res$AGE.LR.NEG <- res$AGE.NEG * (res$from %in% grep("LocalRelapse", colnames(tra)))
res$AGE.LR.POS <- res$AGE.POS * (res$from %in% grep("LocalRelapse", colnames(tra)))
res$AGE.DR.NEG <- res$AGE.NEG * (res$from %in% grep("DistantRelapse", colnames(tra)))
res$AGE.DR.POS <- res$AGE.POS * (res$from %in% grep("DistantRelapse", colnames(tra)))

res$AGE.PS.LR.NEG <- res$AGE.PS.NEG + res$AGE.LR.NEG
res$AGE.PS.LR.POS <- res$AGE.PS.POS + res$AGE.LR.POS


res$GRADE.PS.NEG <- res$GRADE.NEG * (res$from %in% grep("PostSurgery", colnames(tra)))
res$GRADE.PS.POS <- res$GRADE.POS * (res$from %in% grep("PostSurgery", colnames(tra)))
res$GRADE.LR.NEG <- res$GRADE.NEG * (res$from %in% grep("LocalRelapse", colnames(tra)))
res$GRADE.LR.POS <- res$GRADE.POS * (res$from %in% grep("LocalRelapse", colnames(tra)))
res$GRADE.DR.NEG <- res$GRADE.NEG * (res$from %in% grep("DistantRelapse", colnames(tra)))
res$GRADE.DR.POS <- res$GRADE.POS * (res$from %in% grep("DistantRelapse", colnames(tra)))
res$SIZE.PS.NEG <- res$SIZE.NEG * (res$from %in% grep("PostSurgery", colnames(tra)))
res$SIZE.PS.POS <- res$SIZE.POS * (res$from %in% grep("PostSurgery", colnames(tra)))
res$SIZE.LR.NEG <- res$SIZE.NEG * (res$from %in% grep("LocalRelapse", colnames(tra)))
res$SIZE.LR.POS <- res$SIZE.POS * (res$from %in% grep("LocalRelapse", colnames(tra)))
res$SIZE.DR.NEG <- res$SIZE.NEG * (res$from %in% grep("DistantRelapse", colnames(tra)))
res$SIZE.DR.POS <- res$SIZE.POS * (res$from %in% grep("DistantRelapse", colnames(tra)))
res$LN.PS.NEG <- res$LN.NEG * (res$from %in% grep("PostSurgery", colnames(tra)))
res$LN.PS.POS <- res$LN.POS * (res$from %in% grep("PostSurgery", colnames(tra)))
res$LN.LR.NEG <- res$LN.NEG * (res$from %in% grep("LocalRelapse", colnames(tra)))
res$LN.LR.POS <- res$LN.POS * (res$from %in% grep("LocalRelapse", colnames(tra)))
res$LN.DR.NEG <- res$LN.NEG * (res$from %in% grep("DistantRelapse", colnames(tra)))
res$LN.DR.POS <- res$LN.POS * (res$from %in% grep("DistantRelapse", colnames(tra)))
res$TLastSurgery.LR.NEG <- res$TLastSurgery.NEG * (res$from %in% grep("LocalRelapse", colnames(tra)))
res$TLastSurgery.LR.POS <- res$TLastSurgery.POS * (res$from %in% grep("LocalRelapse", colnames(tra)))
res$TLastSurgery.DR.NEG <- res$TLastSurgery.NEG * (res$from %in% grep("DistantRelapse", colnames(tra)))
res$TLastSurgery.DR.POS <- res$TLastSurgery.POS * (res$from %in% grep("DistantRelapse", colnames(tra)))

## SIMPLEST MODEL

m0 <- coxph(Surv(time, status) ~ strata(trans) + AGE + SIZE + GRADE + LN + TLastSurgery + TLastLocal, data=res)

## AGE stays common

## GRADE. DOES IT CHANGE?

m1 <- coxph(Surv(time, status) ~ strata(trans) + AGE + SIZE + GRADE.NEG + GRADE.POS + LN + TLastSurgery + TLastLocal, data=res)
m2 <- coxph(Surv(time, status) ~ strata(trans)  + AGE + SIZE + GRADE.PS + GRADE.LR + GRADE.DR + LN + TLastSurgery + TLastLocal, data=res)
m3 <- coxph(Surv(time, status) ~ strata(trans) + AGE + SIZE + GRADE.PS.NEG + GRADE.LR.NEG + GRADE.DR.NEG + GRADE.PS.POS + GRADE.LR.POS + GRADE.DR.POS + LN +
                TLastSurgery + TLastLocal , data=res)

## SIZE

m1 <- coxph(Surv(time, status) ~ strata(trans) + AGE + GRADE + SIZE.NEG + SIZE.POS + LN + TLastSurgery + TLastLocal, data=res)
m2 <- coxph(Surv(time, status) ~ strata(trans) + AGE + GRADE + SIZE.PS + SIZE.LR + SIZE.DR + LN + TLastSurgery + TLastLocal, data=res)
m3 <- coxph(Surv(time, status) ~ strata(trans) + AGE + GRADE + SIZE.PS.NEG + SIZE.LR.NEG + SIZE.DR.NEG + SIZE.PS.POS + SIZE.LR.POS + SIZE.DR.POS + LN +
                TLastSurgery + TLastLocal, data=res)

## LN

m1 <- coxph(Surv(time, status) ~ strata(trans) + AGE + GRADE + SIZE + LN.NEG + LN.POS + TLastSurgery + TLastLocal, data=res)
m2 <- coxph(Surv(time, status) ~ strata(trans) + AGE + GRADE + SIZE + LN.PS + LN.LR + LN.DR + TLastSurgery + TLastLocal, data=res)
m3 <- coxph(Surv(time, status) ~ strata(trans) + AGE + GRADE + SIZE + LN.PS.NEG + LN.LR.NEG + LN.DR.NEG + LN.PS.POS + LN.LR.POS + LN.DR.POS +
                TLastSurgery + TLastLocal, data=res)

## TLast

m2 <- coxph(Surv(time, status) ~ strata(trans) + AGE + GRADE + SIZE + LN + TLastSurgery.LR + TLastSurgery.DR + TLastLocal, data=res)
m3 <- coxph(Surv(time, status) ~ strata(trans) + AGE + GRADE + SIZE + LN + TLastSurgery.LR.NEG + TLastSurgery.DR.NEG +
                TLastSurgery.LR.POS + TLastSurgery.DR.POS + TLastLocal.DR.NEG + TLastLocal.DR.POS + cluster(id), data=res)

## Full model

fm <- coxph(Surv(time, status) ~ strata(trans) + AGE.PS.NEG + AGE.LR.NEG + AGE.DR.NEG +
                AGE.PS.POS + AGE.LR.POS + AGE.DR.POS + GRADE.PS.NEG + GRADE.LR.NEG + GRADE.DR.NEG + GRADE.PS.POS + GRADE.LR.POS + GRADE.DR.POS +
                SIZE.PS.NEG + SIZE.LR.NEG + SIZE.DR.NEG + SIZE.PS.POS + SIZE.LR.POS + SIZE.DR.POS +
                    LN.PS.NEG + LN.LR.NEG + LN.DR.NEG + LN.PS.POS + LN.LR.POS + LN.DR.POS +
                        TLastSurgery.LR.NEG + TLastSurgery.DR.NEG +
                            TLastSurgery.LR.POS + TLastSurgery.DR.POS + TLastLocal.DR.NEG + TLastLocal.DR.POS+
                          cluster(id), data=res)
## Final model
## WE ARE HERE

model.cofs <- data.frame(Name=names(coef(fm)), HR=summary(fm)$conf.int[,1], LI=summary(fm)$conf.int[,3], UI=summary(fm)$conf.int[,4])
model.cofs$Effect <- sapply(strsplit(as.character(model.cofs$Name), ".", fixed=TRUE), function(x) x[1])
model.cofs$Time <- sapply(strsplit(as.character(model.cofs$Name), ".", fixed=TRUE), function(x) x[2])
model.cofs$ER <- sapply(strsplit(as.character(model.cofs$Name), ".", fixed=TRUE), function(x) x[3])
model.cofs$Time <- factor(model.cofs$Time, levels=c("PS", "LR", "DR"))
model.cofs$ER <- factor(model.cofs$ER, levels=c("NEG", "POS"), labels=c("ER-", "ER+"))
model.cofs <- model.cofs[which(!is.na(model.cofs$ER)),]
model.cofs$Effect[which(model.cofs$Effect=="TLastSurgery")] <- "Time from Surgery"
model.cofs$Effect[which(model.cofs$Effect=="TLastLocal")] <- "Time from L Relapse"



fm <- coxph(Surv(time, status) ~ strata(trans) + AGE.PS.NEG + AGE.LR.NEG + AGE.DR.NEG +
                AGE.PS.POS + AGE.LR.POS + AGE.DR.POS  + GRADE.PS.NEG + GRADE.LR.NEG + GRADE.DR.NEG + GRADE.PS.POS + GRADE.LR.POS + GRADE.DR.POS +
                SIZE.PS.NEG + SIZE.LR.NEG + SIZE.DR.NEG + SIZE.PS.POS + SIZE.LR.POS + SIZE.DR.POS +
                    LN.PS.NEG + LN.LR.NEG + LN.DR.NEG + LN.PS.POS + LN.LR.POS + LN.DR.POS +
                        TLastSurgery.LR.NEG + TLastSurgery.DR.NEG +
                            TLastSurgery.LR.POS + TLastSurgery.DR.POS, data=res)

fm2 <- coxph(Surv(time, status) ~ strata(trans) + AGE.PS.NEG + AGE.LR.NEG + AGE.DR.NEG +
                AGE.PS.POS + AGE.LR.POS + AGE.DR.POS + GRADE.PS.NEG + GRADE.LR.NEG + GRADE.DR.NEG + GRADE.PS.POS + GRADE.LR.POS + GRADE.DR.POS +
                SIZE.PS.NEG + SIZE.LR.NEG + SIZE.DR.NEG + SIZE.PS.POS + SIZE.LR.POS + SIZE.DR.POS +
                    LN.PS.NEG + LN.LR.NEG + LN.DR.NEG + LN.PS.POS + LN.LR.POS + LN.DR.POS +
                        TLastSurgery.LR.NEG + TLastSurgery.DR.NEG +
                            TLastSurgery.LR.POS + TLastSurgery.DR.POS + TLastLocal.DR.NEG + TLastLocal.DR.POS+
                          cluster(id), data=res)


Clinical <- ALL.CLINICAL
ids <- which(Clinical$ER.Status=="ER-")
newdata <- data.frame(strata=1:9,
                          AGE=mean(Clinical$AGE[ids], na.rm=T),
                          LN=mean(Clinical$LN[ids], na.rm=T),
                          GRADE=mean(Clinical$GRADE[ids], na.rm=T),
                          SIZE=mean(Clinical$SIZE[ids], na.rm=T),
                      TLastSurgery.LR=mean(Clinical$TLR[ids][which(Clinical$LR[ids]==1)]),
                      TLastSurgery.DR=mean(Clinical$TDR[ids][which(Clinical$DR[ids]==1)]))
ids <- which(Clinical$ER.Status=="ER+")
newdata <- rbind(newdata, data.frame(strata=10:18,
                          AGE=mean(Clinical$AGE[ids], na.rm=T),
                          LN=mean(Clinical$LN[ids], na.rm=T),
                          GRADE=mean(Clinical$GRADE[ids], na.rm=T),
                          SIZE=mean(Clinical$SIZE[ids], na.rm=T),
                      TLastSurgery.LR=mean(Clinical$TLR[ids][which(Clinical$LR[ids]==1)]),
                      TLastSurgery.DR=mean(Clinical$TDR[ids][which(Clinical$DR[ids]==1)])))

newdata$AGE[which(!newdata$strata %in% nat.death)] <- 0
newdata$TLastSurgery.LR[which(newdata$strata %in% c(1:4, 7:13, 16:18))] <- 0
newdata$TLastSurgery.DR[which(newdata$strata %in% c(1:7, 9:16, 18))] <- 0
newdata$TLastSurgery.LR.NEG <- newdata$TLastSurgery.LR * (newdata$strata %in% 1:8)
newdata$TLastSurgery.LR.POS <- newdata$TLastSurgery.LR * (newdata$strata %in% 10:17)
newdata$TLastSurgery.DR.NEG <- newdata$TLastSurgery.DR * (newdata$strata %in% 1:8)
newdata$TLastSurgery.DR.POS <- newdata$TLastSurgery.DR * (newdata$strata %in% 10:17)

newdata$AGE.PS.NEG <- newdata$AGE * (newdata$strata %in% 4)
newdata$AGE.PS.POS <- newdata$AGE * (newdata$strata %in% 13)
newdata$AGE.LR.NEG <- newdata$AGE * (newdata$strata %in% 7)
newdata$AGE.LR.POS <- newdata$AGE * (newdata$strata %in% 16)
newdata$AGE.DR.NEG <- newdata$AGE * (newdata$strata %in% 9)
newdata$AGE.DR.POS <- newdata$AGE * (newdata$strata %in% 18)

newdata$AGE.PS.LR.NEG <- newdata$AGE.PS.NEG + newdata$AGE.LR.NEG
newdata$AGE.PS.LR.POS <- newdata$AGE.PS.POS + newdata$AGE.LR.POS

newdata$GRADE.PS.NEG <- newdata$GRADE * (newdata$strata %in% 1:3)
newdata$GRADE.PS.POS <- newdata$GRADE * (newdata$strata %in% 10:12)
newdata$GRADE.LR.NEG <- newdata$GRADE * (newdata$strata %in% 5:6)
newdata$GRADE.LR.POS <- newdata$GRADE * (newdata$strata %in% 14:15)
newdata$GRADE.DR.NEG <- newdata$GRADE * (newdata$strata %in% 8)
newdata$GRADE.DR.POS <- newdata$GRADE * (newdata$strata %in% 17)

newdata$SIZE.PS.NEG <- newdata$SIZE * (newdata$strata %in% 1:3)
newdata$SIZE.PS.POS <- newdata$SIZE * (newdata$strata %in% 10:12)
newdata$SIZE.LR.NEG <- newdata$SIZE * (newdata$strata %in% 5:6)
newdata$SIZE.LR.POS <- newdata$SIZE * (newdata$strata %in% 14:15)
newdata$SIZE.DR.NEG <- newdata$SIZE * (newdata$strata %in% 8)
newdata$SIZE.DR.POS <- newdata$SIZE * (newdata$strata %in% 17)

newdata$LN.PS.NEG <- newdata$LN * (newdata$strata %in% 1:3)
newdata$LN.PS.POS <- newdata$LN * (newdata$strata %in% 10:12)
newdata$LN.LR.NEG <- newdata$LN * (newdata$strata %in% 5:6)
newdata$LN.LR.POS <- newdata$LN * (newdata$strata %in% 14:15)
newdata$LN.DR.NEG <- newdata$LN * (newdata$strata %in% 8)
newdata$LN.DR.POS <- newdata$LN * (newdata$strata %in% 17)

## ## Adjust AGE with relapse
## newdata$AGE[which(newdata$strata %in% c(7,16))] <-
##     newdata$AGE[which(newdata$strata %in% c(7,16))] +
##         newdata$TLastSurgery.LR[which(newdata$strata %in% c(6,15))]
## newdata$AGE[which(newdata$strata %in% c(9,18))] <-
##     newdata$AGE[which(newdata$strata %in% c(9,18))] +
##         newdata$TLastSurgery.DR[which(newdata$strata %in% c(8,17))]

class(newdata) <- c("msdata", "data.frame")
attr(newdata, "trans") <- tra
Clinical <- ALL.CLINICAL
Clinical$C <- 1 * (Clinical$Death==0)
cens <- survfit(Surv(T, C) ~ 1, data=Clinical)
cens <- data.frame(time=cens$time, surv=cens$surv, Haz=-log(cens$surv))
cens <- cens[which(is.finite(cens$Haz)),]

save(fm, fm2, res, Clinical, newdata, tra, cens, file="../Models/ERM.RData")
