rm(list=ls())


getCindex <- function(probs, model, Clinical) {
    require(rms)
    ## Another way of computing
    load(probs)
    pslc <- do.call("cbind", pslc)
    pslc <- pslc[,-which(colnames(pslc)=="Times")]
    pslo <- do.call("cbind", pslo)
    pslo <- pslo[,-which(colnames(pslo)=="Times")]
    psld <- do.call("cbind", psld)
    psld <- psld[,-which(colnames(psld)=="Times")]
    psldc <- do.call("cbind", psldc)
    psldc <- psldc[,-which(colnames(psldc)=="Times")]
    psldo <- do.call("cbind", psldo)
    psldo <- psldo[,-which(colnames(psldo)=="Times")]
    psl <- do.call("cbind", psl)
    psl <- psl[,-which(colnames(psl)=="Times")]
    psdc <- do.call("cbind", psdc)
    psdc <- psdc[,-which(colnames(psdc)=="Times")]
    psdo <- do.call("cbind", psdo)
    psdo <- psdo[,-which(colnames(psdo)=="Times")]
    psd <- do.call("cbind", psd)
    psd <- psd[,-which(colnames(psd)=="Times")]
    psc <- do.call("cbind", psc)
    psc <- psc[,-which(colnames(psc)=="Times")]
    pso <- do.call("cbind", pso)
    pso <- pso[,-which(colnames(pso)=="Times")]
    pld <- do.call("cbind", pld)
    pld <- pld[,-which(colnames(pld)=="Times")]
    plo <- do.call("cbind", plo)
    plo <- plo[,-which(colnames(plo)=="Times")]
    pldc <- do.call("cbind", pldc)
    pldc <- pldc[,-which(colnames(pldc)=="Times")]
    pldo <- do.call("cbind", pldo)
    pldo <- pldo[,-which(colnames(pldo)=="Times")]
    plc <- do.call("cbind", plc)
    plc <- plc[,-which(colnames(plc)=="Times")]
    pdc <- do.call("cbind", pdc)
    pdc <- pdc[,-which(colnames(pdc)=="Times")]
    pdo <- do.call("cbind", pdo)
    pdo <- pdo[,-which(colnames(pdo)=="Times")]

    ALLCI <- NULL

    timepoints <- c(2, 5, 10, 15, 20)

    ## Cancer death
    pred1 <- pslc + psldc + psdc + psc
    pred1 <- t(pred1)
    pred1 <- data.frame(METABRIC.ID=rownames(pred1), pred1)
    colnames(pred1)[-1] <- paste("pred", c(2, 5, 10, 15 ,20), sep=".")

    cindex <- numeric(5)
    cindex.se <- numeric(5)
    tmp <- merge(Clinical, pred1)
    id <- which(colnames(tmp)=="pred.2")-1
    for (i in 1:5) {
        TRUNC <- 365.25 * timepoints[i]
        tmp$TIME <- tmp$T
        tmp$EVENT <- tmp$DeathBreast
        tmp$EVENT[which(tmp$TIME>TRUNC)] <- 0
        tmp$TIME[which(tmp$TIME>TRUNC)] <- TRUNC
        res <- survConcordance(Surv(tmp$TIME, tmp$EVENT) ~ tmp[,id+i])
        cindex[i] <- res$concordance
        cindex.se[i] <- res$std.err
    }

    ALLCI <- rbind(ALLCI, data.frame(Model=model, Type="C/D", Year=c(2, 5, 10, 15, 20), cindex=cindex, cindex.se=cindex.se))

    ## Other death
    pred1 <- pslo + psldo + psdo + pso
    pred1 <- t(pred1)
    pred1 <- data.frame(METABRIC.ID=rownames(pred1), pred1)
    colnames(pred1)[-1] <- paste("pred", c(2, 5, 10, 15 ,20), sep=".")

    cindex <- numeric(5)
    cindex.se <- numeric(5)
    tmp <- merge(Clinical, pred1)
    for (i in 1:5) {
        TRUNC <- 365.25 * timepoints[i]
        tmp$TIME <- tmp$T
        tmp$EVENT <- 1 * (tmp$Death==1 & tmp$DeathBreast==0)
        tmp$EVENT[which(tmp$TIME>TRUNC)] <- 0
        tmp$TIME[which(tmp$TIME>TRUNC)] <- TRUNC
        res <- survConcordance(Surv(tmp$TIME, tmp$EVENT) ~ tmp[,id+i])
        cindex[i] <- res$concordance
        cindex.se[i] <- res$std.err
    }

    ALLCI <- rbind(ALLCI, data.frame(Model=model, Type="O/D", Year=c(2, 5, 10, 15, 20), cindex=cindex, cindex.se=cindex.se))

    ## Overall death
    pred1 <- pslc + psldc + psdc + psc + pslo + psldo + psdo + pso
    pred1 <- t(pred1)
    pred1 <- data.frame(METABRIC.ID=rownames(pred1), pred1)
    colnames(pred1)[-1] <- paste("pred", c(2, 5, 10, 15 ,20), sep=".")

    cindex <- numeric(5)
    cindex.se <- numeric(5)
    tmp <- merge(Clinical, pred1)
    for (i in 1:5) {
        TRUNC <- 365.25 * timepoints[i]
        tmp$TIME <- tmp$T
        tmp$EVENT <- tmp$Death
        tmp$EVENT[which(tmp$TIME>TRUNC)] <- 0
        tmp$TIME[which(tmp$TIME>TRUNC)] <- TRUNC
        res <- survConcordance(Surv(tmp$TIME, tmp$EVENT) ~ tmp[,id+i])
        cindex[i] <- res$concordance
        cindex.se[i] <- res$std.err
    }

    ALLCI <- rbind(ALLCI, data.frame(Model=model, Type="Death", Year=c(2, 5, 10, 15, 20), cindex=cindex, cindex.se=cindex.se))

    ## Local relapse
    pred1 <- psl + psldo + psldc + psld + pslc + pslo
    pred1 <- t(pred1)
    pred1 <- data.frame(METABRIC.ID=rownames(pred1), pred1)
    colnames(pred1)[-1] <- paste("pred", c(2, 5, 10, 15 ,20), sep=".")

    cindex <- numeric(5)
    cindex.se <- numeric(5)
    tmp <- merge(Clinical, pred1)
    for (i in 1:5) {
        TRUNC <- 365.25 * timepoints[i]
        tmp$TIME <- tmp$TLR
        tmp$EVENT <- tmp$LR
        tmp$EVENT[which(tmp$TIME>TRUNC)] <- 0
        tmp$TIME[which(tmp$TIME>TRUNC)] <- TRUNC
        res <- survConcordance(Surv(tmp$TIME, tmp$EVENT) ~ tmp[,id+i])
        cindex[i] <- res$concordance
        cindex.se[i] <- res$std.err
    }

    ALLCI <- rbind(ALLCI, data.frame(Model=model, Type="LR-R", Year=c(2, 5, 10, 15, 20), cindex=cindex, cindex.se=cindex.se))

## Distant relapse
    pred1 <- psd + psdc + psdo + psld + psldc + psldo
    pred1 <- t(pred1)
    pred1 <- data.frame(METABRIC.ID=rownames(pred1), pred1)
    colnames(pred1)[-1] <- paste("pred", c(2, 5, 10, 15 ,20), sep=".")

    cindex <- numeric(5)
    cindex.se <- numeric(5)
    tmp <- merge(Clinical, pred1)
    for (i in 1:5) {
        TRUNC <- 365.25 * timepoints[i]
        tmp$TIME <- tmp$TDR
        tmp$EVENT <- tmp$DR
        tmp$EVENT[which(tmp$TIME>TRUNC)] <- 0
        tmp$TIME[which(tmp$TIME>TRUNC)] <- TRUNC
        res <- survConcordance(Surv(tmp$TIME, tmp$EVENT) ~ tmp[,id+i])
        cindex[i] <- res$concordance
        cindex.se[i] <- res$std.err
}

    ALLCI <- rbind(ALLCI, data.frame(Model=model, Type="D-R", Year=c(2, 5, 10, 15, 20), cindex=cindex, cindex.se=cindex.se))

## Relapse
    pred1 <- pslc + psldc + psdc + psc + psl + psldo + pslo + psld + psd + psdo
    pred1 <- t(pred1)
    pred1 <- data.frame(METABRIC.ID=rownames(pred1), pred1)
    colnames(pred1)[-1] <- paste("pred", c(2, 5, 10, 15 ,20), sep=".")

    cindex <- numeric(5)
    cindex.se <- numeric(5)
    tmp <- merge(Clinical, pred1)
    for (i in 1:5) {
        TRUNC <- 365.25 * timepoints[i]
        tmp$TIME <- pmin(tmp$TLR, tmp$TDR)
        tmp$EVENT <- 1 * (tmp$LR==1 | tmp$DR==1 | tmp$DeathBreast==1)
        tmp$EVENT[which(tmp$TIME>TRUNC)] <- 0
        tmp$TIME[which(tmp$TIME>TRUNC)] <- TRUNC
        res <- survConcordance(Surv(tmp$TIME, tmp$EVENT) ~ tmp[,id+i])
        cindex[i] <- res$concordance
        cindex.se[i] <- res$std.err
}


ALLCI <- rbind(ALLCI, data.frame(Model=model, Type="Relapse", Year=c(2, 5, 10, 15, 20), cindex=cindex, cindex.se=cindex.se))



## Distant relapse after Loco-regional
pred1 <- pld + pldc + pldo
pred1 <- t(pred1)

pred1 <- data.frame(METABRIC.ID=rownames(pred1), pred1)
colnames(pred1)[-1] <- paste("pred", seq(from=0, to=20, by=0.25), sep=".")
pred1 <- pred1[,c('METABRIC.ID', 'pred.2', 'pred.5', 'pred.10', 'pred.15', 'pred.20')]
    cindex <- numeric(5)
        cindex.se <- numeric(5)
tmp <- merge(Clinical, pred1)
tmp <- tmp[which(tmp$LR==1),]
for (i in 1:5) {
    TRUNC <- 365.25 * timepoints[i]
    tmp$TIME <- I(tmp$TDR-tmp$TLR)
    tmp$EVENT <- tmp$DR
    tmp$EVENT[which(tmp$TIME>TRUNC)] <- 0
    tmp$TIME[which(tmp$TIME>TRUNC)] <- TRUNC
    res <- survConcordance(Surv(tmp$TIME, tmp$EVENT) ~ tmp[,id+i])
        cindex[i] <- res$concordance
        cindex.se[i] <- res$std.err
}

ALLCI <- rbind(ALLCI, data.frame(Model=model, Type="D-R After LR-R", Year=c(2, 5, 10, 15, 20), cindex=cindex, cindex.se=cindex.se))

## Relapse (+cancer death) after Loco-regional
pred1 <- pld + pldc + pldo + plc
pred1 <- t(pred1)

pred1 <- data.frame(METABRIC.ID=rownames(pred1), pred1)
colnames(pred1)[-1] <- paste("pred", seq(from=0, to=20, by=0.25), sep=".")
pred1 <- pred1[,c('METABRIC.ID', 'pred.2', 'pred.5', 'pred.10', 'pred.15', 'pred.20')]
    cindex <- numeric(5)
        cindex.se <- numeric(5)
tmp <- merge(Clinical, pred1)
tmp <- tmp[which(tmp$LR==1),]
for (i in 1:5) {
    TRUNC <- 365.25 * timepoints[i]
    tmp$TIME <- I(tmp$TDR-tmp$TLR)
    tmp$EVENT <- 1 * (tmp$DR==1 | tmp$DeathBreast==1)
    tmp$EVENT[which(tmp$TIME>TRUNC)] <- 0
    tmp$TIME[which(tmp$TIME>TRUNC)] <- TRUNC
    res <- survConcordance(Surv(tmp$TIME, tmp$EVENT) ~ tmp[,id+i])
        cindex[i] <- res$concordance
        cindex.se[i] <- res$std.err
}

ALLCI <- rbind(ALLCI, data.frame(Model=model, Type="R.C/D After LR-R", Year=c(2, 5, 10, 15, 20), cindex=cindex, cindex.se=cindex.se))

## Other death after Loco-regional
pred1 <- plo
pred1 <- t(pred1)

pred1 <- data.frame(METABRIC.ID=rownames(pred1), pred1)
colnames(pred1)[-1] <- paste("pred", seq(from=0, to=20, by=0.25), sep=".")
pred1 <- pred1[,c('METABRIC.ID', 'pred.2', 'pred.5', 'pred.10', 'pred.15', 'pred.20')]
    cindex <- numeric(5)
    cindex.se <- numeric(5)
tmp <- merge(Clinical, pred1)
tmp <- tmp[which(tmp$LR==1),]
for (i in 1:5) {
    TRUNC <- 365.25 * timepoints[i]
    tmp$TIME <- I(tmp$T-tmp$TLR)
    tmp$EVENT <- 1 * (tmp$Death==1 & tmp$DeathBreast==0)
    tmp$EVENT[which(tmp$TIME>TRUNC)] <- 0
    tmp$TIME[which(tmp$TIME>TRUNC)] <- TRUNC
    res <- survConcordance(Surv(tmp$TIME, tmp$EVENT) ~ tmp[,id+i])
        cindex[i] <- res$concordance
        cindex.se[i] <- res$std.err
}

ALLCI <- rbind(ALLCI, data.frame(Model=model, Type="O/D After LR-R", Year=c(2, 5, 10, 15, 20), cindex=cindex, cindex.se=cindex.se))


## Cancer Death After Distant Relapse
pred1 <- pdc
pred1 <- t(pred1)

pred1 <- data.frame(METABRIC.ID=rownames(pred1), pred1)
colnames(pred1)[-1] <- paste("pred", seq(from=0, to=20, by=0.25), sep=".")
pred1 <- pred1[,c('METABRIC.ID', 'pred.2', 'pred.5', 'pred.10', 'pred.15', 'pred.20')]
    cindex <- numeric(5)
    cindex.se <- numeric(5)
tmp <- merge(Clinical, pred1)
tmp <- tmp[which(tmp$DR==1),]
for (i in 1:5) {
    TRUNC <- 365.25 * timepoints[i]
    tmp$TIME <- I(tmp$T-tmp$TDR)
    tmp$EVENT <- tmp$DeathBreast
    tmp$EVENT[which(tmp$TIME>TRUNC)] <- 0
    tmp$TIME[which(tmp$TIME>TRUNC)] <- TRUNC
    res <- survConcordance(Surv(tmp$TIME, tmp$EVENT) ~ tmp[,id+i])
        cindex[i] <- res$concordance
        cindex.se[i] <- res$std.err
}

ALLCI <- rbind(ALLCI, data.frame(Model=model, Type="C/D After D-R", Year=c(2, 5, 10, 15, 20), cindex=cindex, cindex.se=cindex.se))

## Other Death After Distant Relapse
pred1 <- pdo
pred1 <- t(pred1)

pred1 <- data.frame(METABRIC.ID=rownames(pred1), pred1)
colnames(pred1)[-1] <- paste("pred", seq(from=0, to=20, by=0.25), sep=".")
pred1 <- pred1[,c('METABRIC.ID', 'pred.2', 'pred.5', 'pred.10', 'pred.15', 'pred.20')]
    cindex <- numeric(5)
    cindex.se <- numeric(5)
tmp <- merge(Clinical, pred1)
tmp <- tmp[which(tmp$DR==1),]
for (i in 1:5) {
    TRUNC <- 365.25 * timepoints[i]
    tmp$TIME <- I(tmp$T-tmp$TDR)
    tmp$EVENT <- 1 * (tmp$Death==1 & tmp$DeathBreast==0)
    tmp$EVENT[which(tmp$TIME>TRUNC)] <- 0
    tmp$TIME[which(tmp$TIME>TRUNC)] <- TRUNC
    res <- survConcordance(Surv(tmp$TIME, tmp$EVENT) ~ tmp[,id+i])
        cindex[i] <- res$concordance
        cindex.se[i] <- res$std.err
}

ALLCI <- rbind(ALLCI, data.frame(Model=model, Type="O/D After D-R", Year=c(2, 5, 10, 15, 20), cindex=cindex, cindex.se=cindex.se))
    ALLCI
}



X <- NULL
load("../Models/ERM.RData")
X <- rbind(X, getCindex("../Models/ER_AllProbs.RData", "ERM", Clinical=Clinical))
load("../Models/ICM.RData")
X <- rbind(X, getCindex("../Models/IntClust_AllProbs.RData", "IntClustM", Clinical=Clinical))
load("../Models/Pam50M.RData")
X <- rbind(X, getCindex("../Models/Pam50_AllProbs.RData", "Pam50M", Clinical=Clinical))
load("../Models/FOURGROUPSM.RData")
X <- rbind(X, getCindex("../Models/FourGroups_AllProbs.RData", "IHCM", Clinical=Clinical))


## Now we add predict
load(file="predictCindex.RData")
X <- rbind(X, predict.cindex)



pdf("cindex_Models.pdf", width=12, height=9)
print(xyplot(cindex ~ Year|Type, group=Model, data=X, auto.key=T, pch=21))
X <- X[order(X$Type, X$Year, X$Model),]
X$li <- X$cindex - 1.96 * X$cindex.se
X$ui <- X$cindex + 1.96 * X$cindex.se
print(xYplot(Cbind(cindex, li, ui) ~ Year|Type, group=Model, data=X, auto.key=T, pch=21))
dev.off()
pdf("cindex.pdf", width=12, height=9)
for (i in c(2, 5, 10, 15, 20)) {
    print(Dotplot(Model ~ Cbind(cindex, li, ui)| Type, data=subset(X, Year==i), pch=19, main=paste(i, "years predictions"), xlab="c-index"))
}
for (i in unique(X$Model)) {
    print(xYplot(Cbind(cindex, li, ui) ~ Year| Type, data=subset(X, Model==i), pch=19, main=i, xlab="Years", ylim=c(0.47, 0.95)))
}
dev.off()
save(X, file="Cindex.RData")

####################################################
####################################################
####################################################
## ER+ Samples
####################################################
####################################################
####################################################
X <- NULL
load("~/Documents/Projects/Metastasis/PaperRevision/VersionAGE/Scripts/Models/ERM.RData")
X <- rbind(X, getCindex("~/Documents/Projects/Metastasis/PaperRevision/Scripts/Models/ER_AllProbs.RData", "ERM", Clinical=subset(Clinical, ER.Status=="ER+")))
load("~/Documents/Projects/Metastasis/PaperRevision/VersionAGE/Scripts/Models/ICM.RData")
X <- rbind(X, getCindex("~/Documents/Projects/Metastasis/PaperRevision/Scripts/Models/IntClust_AllProbs.RData", "IntClustM", Clinical=subset(Clinical, ER.Status=="pos")))
load("~/Documents/Projects/Metastasis/PaperRevision/VersionAGE/Scripts/Models/Pam50M.RData")
X <- rbind(X, getCindex("~/Documents/Projects/Metastasis/PaperRevision/Scripts/Models/Pam50_AllProbs.RData", "Pam50M", Clinical=subset(Clinical, ER.Status=="ER+")))
load("~/Documents/Projects/Metastasis/PaperRevision/VersionAGE/Scripts/Models/FOURGROUPSM.RData")
X <- rbind(X, getCindex("~/Documents/Projects/Metastasis/PaperRevision/Scripts/Models/FourGroups_AllProbs.RData", "ClinicalM", Clinical=subset(Clinical, ER.Status=="ER+")))


## Now we add predict
load(file="/Users/rueda01/Documents/Projects/Metastasis/PaperRevision/VersionAGE/Validation/predictCindexERPOS.RData")
X <- rbind(X, predict.cindex)



pdf("cindex_ModelsER+.pdf", width=12, height=9)
xyplot(cindex ~ Year|Type, group=Model, data=X, auto.key=T, pch=21)
X <- X[order(X$Type, X$Year, X$Model),]
X$li <- X$cindex - 1.96 * X$cindex.se
X$ui <- X$cindex + 1.96 * X$cindex.se
xYplot(Cbind(cindex, li, ui) ~ Year|Type, group=Model, data=X, auto.key=T, pch=21)
dev.off()
pdf("cindexER+.pdf", width=12, height=9)
for (i in c(2, 5, 10, 15, 20)) {
    print(Dotplot(Model ~ Cbind(cindex, li, ui)| Type, data=subset(X, Year==i), pch=19, main=paste(i, "years predictions"), xlab="c-index"))
}
for (i in unique(X$Model)) {
    print(xYplot(Cbind(cindex, li, ui) ~ Year| Type, data=subset(X, Model==i), pch=19, main=i, xlab="Years", ylim=c(0.47, 0.95)))
}
dev.off()
save(X, file="~/Documents/Projects/Metastasis/PaperRevision/VersionAGE/Validation/CindexERPOS.RData")

####################################################
####################################################
####################################################
## ER- Samples
####################################################
####################################################
####################################################
X <- NULL
load("~/Documents/Projects/Metastasis/PaperRevision/VersionAGE/Scripts/Models/ERM.RData")
X <- rbind(X, getCindex("~/Documents/Projects/Metastasis/PaperRevision/Scripts/Models/ER_AllProbs.RData", "ERM", Clinical=subset(Clinical, ER.Status=="ER-")))
load("~/Documents/Projects/Metastasis/PaperRevision/VersionAGE/Scripts/Models/ICM.RData")
X <- rbind(X, getCindex("~/Documents/Projects/Metastasis/PaperRevision/Scripts/Models/IntClust_AllProbs.RData", "IntClustM", Clinical=subset(Clinical, ER.Status=="neg")))
load("~/Documents/Projects/Metastasis/PaperRevision/VersionAGE/Scripts/Models/Pam50M.RData")
X <- rbind(X, getCindex("~/Documents/Projects/Metastasis/PaperRevision/Scripts/Models/Pam50_AllProbs.RData", "Pam50M", Clinical=subset(Clinical, ER.Status=="ER-")))
load("~/Documents/Projects/Metastasis/PaperRevision/VersionAGE/Scripts/Models/FOURGROUPSM.RData")
X <- rbind(X, getCindex("~/Documents/Projects/Metastasis/PaperRevision/Scripts/Models/FourGroups_AllProbs.RData", "ClinicalM", Clinical=subset(Clinical, ER.Status=="ER-")))


## Now we add predict
load(file="/Users/rueda01/Documents/Projects/Metastasis/PaperRevision/VersionAGE/Validation/predictCindexERNEG.RData")
X <- rbind(X, predict.cindex)



pdf("cindex_ModelsER-.pdf", width=12, height=9)
xyplot(cindex ~ Year|Type, group=Model, data=X, auto.key=T, pch=21)
X <- X[order(X$Type, X$Year, X$Model),]
X$li <- X$cindex - 1.96 * X$cindex.se
X$ui <- X$cindex + 1.96 * X$cindex.se
xYplot(Cbind(cindex, li, ui) ~ Year|Type, group=Model, data=X, auto.key=T, pch=21)
dev.off()
pdf("cindexER-.pdf", width=12, height=9)
for (i in c(2, 5, 10, 15, 20)) {
    print(Dotplot(Model ~ Cbind(cindex, li, ui)| Type, data=subset(X, Year==i), pch=19, main=paste(i, "years predictions"), xlab="c-index"))
}
for (i in unique(X$Model)) {
    print(xYplot(Cbind(cindex, li, ui) ~ Year| Type, data=subset(X, Model==i), pch=19, main=i, xlab="Years", ylim=c(0.47, 0.95)))
}
dev.off()
save(X, file="~/Documents/Projects/Metastasis/PaperRevision/VersionAGE/Validation/CindexERNEG.RData")


########################################################################
########################################################################
## cindex: New Metabric samples
########################################################################
########################################################################
########################################################################




X <- NULL
tmp1 <- read.table("~/Documents/Projects/Metastasis/PaperRevision/VersionAGE/Validation/ValidationMetabricExp.txt", header=T, sep="\t")
tmp2 <- read.table("~/Documents/Projects/Metastasis/PaperRevision/VersionAGE/Validation/ValidationMetabricCN.txt", header=T, sep="\t")
tmp2 <- tmp2[-which(tmp2$METABRIC.ID %in% tmp1$METABRIC.ID),]

Clinical <- rbind(tmp1, tmp2)
X <- rbind(X, getCindex("../Models/IntClust_AllProbs.RData", "IntClustM", Clinical=Clinical))

