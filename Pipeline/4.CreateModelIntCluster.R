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
## MD
##################################
##################################

Clinical <- read.table(file="TableS6.txt", header=T, sep="\t", quote="", comment.char="", stringsAsFactors=FALSE)
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

## Transition matrix
tra <- matrix(NA, 5*11, 5*11)

state.names <- paste(c("PostSurgery",
                            "LocoregionalRelapse", "DistantRelapse",
                       "CancerDeath", "NaturalDeath"),
                     rep(levels(Clinical$iC10), rep(5, 11)), sep="_")
rownames(tra) <- state.names
colnames(tra) <- rownames(tra)
counter <- 1
for (i in seq(from=1, by=5, length=11)) {
    tra[i,(i+1):(i+4)] <- counter:(counter+3)
    counter <- counter + 4
    tra[i+1,(i+2):(i+4)] <- counter:(counter+2)
    counter <- counter + 3
    tra[i+2,(i+3):(i+4)] <- counter:(counter+1)
    counter <- counter + 2
}


covs <- c("LN", "AGE", "GRADE", "SIZE", "BS")

Clinical$Group <- Clinical$iC10
ALL.CLINICAL <- Clinical


set.seed(12324)
all.X <- list()
for (ic in levels(ALL.CLINICAL$Group)) {
    Clinical <- ALL.CLINICAL[which(ALL.CLINICAL$Group==ic),]
    X <- NULL
    for (i in 1:nrow(Clinical)) {
    ## Surgery to Local
    tmp <- data.frame(id=Clinical$METABRIC.ID[i],
                      from=1,
                      to=2,
                      trans=1,
                      Tstart=0,
                      Tstop=Clinical$TLR[i],
                      status=Clinical$LR[i],
                      LN=Clinical$LN[i],
                      AGE=Clinical$AGE[i],
                      GRADE=Clinical$GRADE[i],
                      SIZE=Clinical$SIZE[i],
                      BS=Clinical$BS[i],
                      ER=Clinical$ER.Status[i],
                      Group=ic,
                      TLastSurgery=0,
                      TLastLocal=0)
    X <- rbind(X, tmp)
    ## Surgery to Distant
    tmp <- data.frame(id=Clinical$METABRIC.ID[i],
                      from=1,
                      to=3,
                      trans=2,
                      Tstart=0,
                      Tstop=Clinical$TDR[i],
                      status=1 * (Clinical$LR[i]==0) * Clinical$DR[i],
                      LN=Clinical$LN[i],
                      AGE=Clinical$AGE[i],
                      GRADE=Clinical$GRADE[i],
                      SIZE=Clinical$SIZE[i],
                      BS=Clinical$BS[i],
                      ER=Clinical$ER.Status[i],
                      Group=ic,
                      TLastSurgery=0,
                      TLastLocal=0)
    X <- rbind(X, tmp)
    ## Surgery to Cancer Death
    tmp <- data.frame(id=Clinical$METABRIC.ID[i],
                      from=1,
                      to=4,
                      trans=3,
                      Tstart=0,
                      Tstop=Clinical$T[i],
                      status=1 * (Clinical$LR[i]==0 & Clinical$DR[i]==0) *
                          Clinical$DeathBreast[i],
                      LN=Clinical$LN[i],
                      AGE=Clinical$AGE[i],
                      GRADE=Clinical$GRADE[i],
                      SIZE=Clinical$SIZE[i],
                      BS=Clinical$BS[i],
                      ER=Clinical$ER.Status[i],
                      Group=ic,
                      TLastSurgery=0,
                      TLastLocal=0)
    X <- rbind(X, tmp)
    ## Surgery to Natural Death
    tmp <- data.frame(id=Clinical$METABRIC.ID[i],
                      from=1,
                      to=5,
                      trans=4,
                      Tstart=0,
                      Tstop=Clinical$T[i],
                      status=1 * (Clinical$LR[i]==0 & Clinical$DR[i]==0) *
                          Clinical$NaturalDeath[i],
                      LN=Clinical$LN[i],
                      AGE=Clinical$AGE[i],
                      GRADE=Clinical$GRADE[i],
                      SIZE=Clinical$SIZE[i],
                      BS=Clinical$BS[i],
                      ER=Clinical$ER.Status[i],
                      Group=ic,
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
                          Tstop=Clinical$TDR[i],
                          status=Clinical$DR[i],
                          LN=Clinical$LN[i],
                          AGE=Clinical$AGE[i],
                          GRADE=Clinical$GRADE[i],
                          SIZE=Clinical$SIZE[i],
                          BS=Clinical$BS[i],
                          ER=Clinical$ER.Status[i],
                      Group=ic,
                      TLastSurgery=Clinical$TLR[i],
                      TLastLocal=0)
        X <- rbind(X, tmp)
        ## Local to cancer Death
        tmp <- data.frame(id=Clinical$METABRIC.ID[i],
                          from=2,
                          to=4,
                          trans=6,
                          Tstart=Clinical$TLR[i],
                          Tstop=Clinical$T[i],
                          status=1 * (Clinical$DR[i]==0) * Clinical$DeathBreast[i],
                          LN=Clinical$LN[i],
                          AGE=Clinical$AGE[i],
                          GRADE=Clinical$GRADE[i],
                          SIZE=Clinical$SIZE[i],
                          BS=Clinical$BS[i],
                      ER=Clinical$ER.Status[i],
                      Group=ic,
                      TLastSurgery=Clinical$TLR[i],
                      TLastLocal=0)
        X <- rbind(X, tmp)
        ## Local to Natural Death
        tmp <- data.frame(id=Clinical$METABRIC.ID[i],
                          from=2,
                          to=5,
                          trans=7,
                          Tstart=Clinical$TLR[i],
                          Tstop=Clinical$T[i],
                          status=1 * (Clinical$DR[i]==0) * Clinical$NaturalDeath[i],
                          LN=Clinical$LN[i],
                          AGE=Clinical$AGE[i],
                          GRADE=Clinical$GRADE[i],
                          SIZE=Clinical$SIZE[i],
                          BS=Clinical$BS[i],
                      ER=Clinical$ER.Status[i],
                      Group=ic,
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
                          AGE=Clinical$AGE[i],
                          GRADE=Clinical$GRADE[i],
                          SIZE=Clinical$SIZE[i],
                          BS=Clinical$BS[i],
                      ER=Clinical$ER.Status[i],
                      Group=ic,
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
                          AGE=Clinical$AGE[i],
                          GRADE=Clinical$GRADE[i],
                          SIZE=Clinical$SIZE[i],
                          BS=Clinical$BS[i],
                      ER=Clinical$ER.Status[i],
                      Group=ic,
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
res <- res[,c(1:6, 17,7:15)]

class(res) <- c("msdata", "data.frame")
attr(res, "trans") <- tra
library(mstate)
tmp <- events(res)


tmp <- events(res)


nat.death <- as.numeric((tra[,grep("NaturalDeath", colnames(tra))]))
nat.death <- unique(na.omit(nat.death))
res$AGE[which(!res$trans %in% nat.death)] <- 0
res$GRADE[which(res$trans %in% nat.death)] <- 0
res$SIZE[which(res$trans %in% nat.death)] <- 0
res$LN[which(res$trans %in% nat.death)] <- 0

res$GRADE.PS <- res$GRADE * (res$from %in% grep("PostSurgery", colnames(tra)))
res$GRADE.LR <- res$GRADE * (res$from %in% grep("LocoregionalRelapse", colnames(tra)))
res$GRADE.DR <- res$GRADE * (res$from %in% grep("DistantRelapse", colnames(tra)))
res$GRADE.R <- res$GRADE * (res$from %in% c(grep("LocoregionalRelapse", colnames(tra)), grep("DistantRelapse", colnames(tra))))
res$SIZE.PS <- res$SIZE * (res$from %in% grep("PostSurgery", colnames(tra)))
res$SIZE.LR <- res$SIZE * (res$from %in% grep("LocoregionalRelapse", colnames(tra)))
res$SIZE.DR <- res$SIZE * (res$from %in% grep("DistantRelapse", colnames(tra)))
res$SIZE.R <- res$SIZE * (res$from %in% c(grep("LocoregionalRelapse", colnames(tra)), grep("DistantRelapse", colnames(tra))))
res$LN.PS <- res$LN * (res$from %in% grep("PostSurgery", colnames(tra)))
res$LN.LR <- res$LN * (res$from %in% grep("LocoregionalRelapse", colnames(tra)))
res$LN.DR <- res$LN * (res$from %in% grep("DistantRelapse", colnames(tra)))
res$LN.R <- res$LN * (res$from %in% c(grep("LocoregionalRelapse", colnames(tra)), grep("DistantRelapse", colnames(tra))))

res$TLastSurgery.LR <- res$TLastSurgery * (res$from %in% grep("LocoregionalRelapse", colnames(tra)))
res$TLastSurgery.DR <- res$TLastSurgery * (res$from %in% grep("DistantRelapse", colnames(tra)))


m0 <- coxph(Surv(time, status) ~ strata(trans)  + AGE + SIZE + GRADE + LN + TLastSurgery, data=res)

fm <- coxph(Surv(time, status) ~ strata(trans)  + AGE + GRADE.PS + GRADE.LR + GRADE.DR + SIZE.PS + SIZE.LR + SIZE.DR +
                LN.PS + LN.LR + LN.DR +
                    TLastSurgery.LR + TLastSurgery.DR, data=res)

im <- coxph(Surv(time, status) ~ strata(trans)  + AGE + GRADE.PS + SIZE.PS +
                LN.PS +
                    TLastSurgery.LR + TLastSurgery.DR, data=res)

m <- coxph(Surv(time, status) ~ strata(trans)  + AGE + GRADE.PS + GRADE.R +
               SIZE.PS + SIZE.LR + LN.PS + LN.R +
                    TLastSurgery + cluster(id), data=res)

m <- coxph(Surv(time, status) ~ strata(trans)  + AGE + GRADE.PS + GRADE.R +
               SIZE.PS + SIZE.LR + LN.PS + LN.R +
                    TLastSurgery, data=res)

## Test effects that might need different IntClusts
res$Risk <- factor(res$Group)
levels(res$Risk) <- list("ER+Bad"=c("1", "2", "6", "9"),
                         "ER+Good"=c("3", "4ER+", "7", "8"),
                         "5"=c("5"), "ER-"=c("4ER-", "10"))
mAlt <- coxph(Surv(time, status) ~ strata(trans)  + AGE + GRADE.PS:Risk + GRADE.R +
               SIZE.PS + SIZE.LR + LN.PS + LN.R +
                    TLastSurgery, data=res)
anova(m, mAlt)
mAlt <- coxph(Surv(time, status) ~ strata(trans)  + AGE + GRADE.PS + GRADE.R:Risk +
               SIZE.PS + SIZE.LR + LN.PS + LN.R +
                    TLastSurgery, data=res)
anova(m, mAlt)

mAlt <- coxph(Surv(time, status) ~ strata(trans)  + AGE + GRADE.PS + GRADE.R +
               SIZE.PS:Risk + SIZE.LR + LN.PS + LN.R +
                    TLastSurgery, data=res)
anova(mAlt, m)
mAlt <- coxph(Surv(time, status) ~ strata(trans)  + AGE + GRADE.PS + GRADE.R +
               SIZE.PS + SIZE.LR:Risk + LN.PS + LN.R +
                    TLastSurgery, data=res)
anova(mAlt, m)
mAlt <- coxph(Surv(time, status) ~ strata(trans)  + AGE + GRADE.PS + GRADE.R +
               SIZE.PS + SIZE.LR + LN.PS:Risk + LN.R +
                    TLastSurgery, data=res)
anova(mAlt, m)
mAlt <- coxph(Surv(time, status) ~ strata(trans)  + AGE + GRADE.PS + GRADE.R +
               SIZE.PS + SIZE.LR + LN.PS + LN.R:Risk +
                    TLastSurgery, data=res)
anova(mAlt, m)

mAlt <- coxph(Surv(time, status) ~ strata(trans)  + AGE + GRADE.PS + GRADE.R +
               SIZE.PS + SIZE.LR + LN.PS + LN.R  +
                    TLastSurgery:Risk, data=res)
anova(mAlt, m)

## m <- coxph(Surv(time, status) ~ strata(trans)  + AGE + GRADE.PS + GRADE.R +
##                SIZE.PS + SIZE.LR + LN.PS + LN.R +
##                     TLast.LR + TLast.DR, data=res)

Clinical <- ALL.CLINICAL
Clinical$C <- 1 * (Clinical$Death==0)
cens <- survfit(Surv(T, C) ~ 1, data=Clinical)
cens <- data.frame(time=cens$time, surv=cens$surv, Haz=-log(cens$surv))
cens <- cens[which(is.finite(cens$Haz)),]

IC <- levels(Clinical$Group)
newdata <- NULL
for (ic in 1:length(IC)) {
    ids <- which(Clinical$Group==IC[ic])
    newdata <- rbind(newdata,
                     data.frame(strata=(9*(ic-1)+1):(9*(ic-1)+9),
                          AGE=mean(Clinical$AGE[ids], na.rm=T),
                          LN=mean(Clinical$LN[ids], na.rm=T),
                          GRADE=mean(Clinical$GRADE[ids], na.rm=T),
                          SIZE=mean(Clinical$SIZE[ids], na.rm=T),
                          ## BS=names(which.max(table(Clinical$BS))),
                      TLastSurgery=c(rep(0, 4),
                          rep(mean(Clinical$TLR[ids][which(Clinical$LR[ids]==1)]), 3), rep(mean(Clinical$TDR[ids][which(Clinical$DR[ids]==1)]), 2))))
}

rel.can <- c(1:4, 7, 9, 10:13, 16, 18, 19:22, 25, 27, 28:31, 34, 36, 37:40,
             43, 45, 46:49, 52, 54, 55:58, 61, 63, 64:67, 70, 72, 73:76,
             79, 81, 82:85, 88, 90, 91:94, 97, 99)
newdata$AGE[which(!newdata$strata %in% nat.death)] <- 0
newdata$TLastSurgery[which(newdata$strata %in% rel.can)] <- 0

newdata$GRADE.PS <- newdata$GRADE * (newdata$strata %in%
                                         as.vector(tra[grep("Post", rownames(tra)),
                                                       -grep("Nat", colnames(tra))]))
newdata$GRADE.R <- newdata$GRADE * (newdata$strata %in%
                                        as.vector(tra[c(grep("Loco", rownames(tra)),
                                                        grep("Distant", rownames(tra))),
                                                       -grep("Nat", colnames(tra))]))

newdata$SIZE.PS <- newdata$SIZE * (newdata$strata %in%
                                         as.vector(tra[grep("Post", rownames(tra)),
                                                       -grep("Nat", colnames(tra))]))
newdata$SIZE.LR <- newdata$SIZE * (newdata$strata %in%
                                         as.vector(tra[grep("Loco", rownames(tra)),
                                                       -grep("Nat", colnames(tra))]))

newdata$LN.PS <- newdata$LN * (newdata$strata %in%
                                         as.vector(tra[grep("Post", rownames(tra)),
                                                       -grep("Nat", colnames(tra))]))
newdata$LN.R <- newdata$LN * (newdata$strata %in%
                                  as.vector(tra[c(grep("Loco", rownames(tra)),
                                                  grep("Distant", rownames(tra))),
                                                       -grep("Nat", colnames(tra))]))




class(newdata) <- c("msdata", "data.frame")
attr(newdata, "trans") <- tra
save(m, res, Clinical, newdata, tra, cens, file="ICM.RData")
