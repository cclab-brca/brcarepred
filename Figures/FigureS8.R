rm(list=ls())
library(RColorBrewer)

cols.er <- c("olivedrab", "red")
cols <- c("lightskyblue", "lightsalmon", "plum2", "wheat4")
coliPam <- c("#E41A1C", "#FB9A99", "#1F78B4", "#A6CEE3", "#66A61E")
coliClust <- c('#FF5500', '#00EE76', '#CD3278','#00C5CD', '#B5D0D2', '#8B0000',
                '#FFFF40', '#0000CD', '#FFAA00', '#EE82EE', '#7D26CD')

Clinical <- read.table(file="../../TableS6.txt", header=T, sep="\t", quote="", comment.char="", stringsAsFactors=FALSE)
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
Clinical$HT <- 1 * (Clinical$HT!="null")
Clinical$HT <- factor(Clinical$HT, levels=c(0, 1), labels=c("NO", "YES"))
Clinical$CT <- 1 * (Clinical$CT!="null")
Clinical$CT <- factor(Clinical$CT, levels=c(0, 1), labels=c("NO", "YES"))
Clinical$RT[which(Clinical$RT %in% c("null", "NONE RECORDED IN LANTIS", "Nnne"))] <- 0
Clinical$RT[which(Clinical$RT!=0)] <- 1
Clinical$RT <- factor(Clinical$RT, levels=c(0, 1), labels=c("NO", "YES"))
Clinical$BS <- factor(Clinical$Breast.Surgery, levels=c("BREAST CONSERVING", "MASTECTOMY"), labels=c("BC", "M"))
Clinical$ER.Status <- factor(Clinical$ER.Status,
                             levels=c("neg", "pos"), labels=c("ER-", "ER+"))

Clinical$NaturalDeath <- 1 * (Clinical$Last.Followup.Status %in% c("d", "d-o.c."))


load(file="../Models/IntClust_AllProbs.RData")
pld <- do.call("cbind", pld)
pld <- pld[,-which(colnames(pld)=="Times")]
plc <- do.call("cbind", plc)
plc <- plc[,-which(colnames(plc)=="Times")]
pldc <- do.call("cbind", pldc)
pldc <- pldc[,-which(colnames(pldc)=="Times")]
pldo <- do.call("cbind", pldo)
pldo <- pldo[,-which(colnames(pldo)=="Times")]
RL <- pld + plc + pldc + pldo

RL <- RL[3,]
RL <- data.frame(METABRIC.ID=names(RL), RL=RL)
Clinical <- merge(Clinical, RL)

pdf("FigureS8.pdf", width=12, height=6)
par(mfrow=c(2, 4), oma=c(3, 3, 1, 2), cex.lab=1.4, cex.axis=1.3)
boxplot(RL ~ GRADE, data=Clinical, xlab="Tumor Grade", col=brewer.pal(3, "Blues"), ylim=c(0,1),
        axes=F)
axis(1, at=c(1, 2, 3))
axis(2, las=2)
box()
boxplot(RL ~ LN, data=Clinical, xlab="Number of Lymph Nodes", col=rev(brewer.pal(11, "RdBu")),
        ylim=c(0,1), axes=F)
axis(1)
axis(2, las=2)
box()

plot(RL ~ SIZE, data=Clinical, xlab="Tumor Size", pch=".", ylim=c(0,1), axes=F, ylab="")
axis(1)
axis(2, las=2)
box()

ids <- order(Clinical$SIZE)
lines(Clinical$SIZE[ids], predict(loess(RL~SIZE, data=Clinical))[ids], col=2, ylim=c(0,1))
plot(RL ~ TLR, data=Clinical, xlab="Time of relapse", ylim=c(0,1), pch=".", axes=F, ylab="")
axis(1)
axis(2, las=2)
box()

ids <- order(Clinical$TLR)
lines(Clinical$TLR[ids], predict(loess(RL~TLR, data=Clinical))[ids], col=2, ylim=c(0,1))
ids <- order(tapply(Clinical$RL, Clinical$ER.Status, median))
Clinical$ER.Status <- factor(Clinical$ER.Status, levels=levels(Clinical$ER.Status)[ids])
boxplot(RL ~ ER.Status, data=Clinical, col=cols.er[ids], ylim=c(0,1))
Clinical$Group <- NA
Clinical$Group[which(Clinical$ER.Status=="ER-" & Clinical$Her2.Expr=="-")] <- "ER-/HER2-"
Clinical$Group[which(Clinical$ER.Status=="ER-" & Clinical$Her2.Expr=="+")] <- "ER-/HER2+"
Clinical$Group[which(Clinical$ER.Status=="ER+" & Clinical$Her2.Expr=="+")] <- "ER+/HER2+"
Clinical$Group[which(Clinical$ER.Status=="ER+" & Clinical$Her2.Expr=="-")] <- "ER+/HER2-"
Clinical$Group <- factor(Clinical$Group)
Clinical$iC10 <- factor(Clinical$iC10, levels=c(1:3, "4ER+", "4ER-", 5:10))
ids <- order(tapply(Clinical$RL, Clinical$Group, median))
Clinical$Group <- factor(Clinical$Group, levels=levels(Clinical$Group)[ids])
boxplot(RL ~ Group, data=Clinical, las=2, col=cols[ids], ylim=c(0,1))
Clinical$Group <- Clinical$Pam50Subtype
Clinical$Group <- factor(Clinical$Group, levels=c("Basal", "Her2", "LumA", "LumB", "Normal"))
ids <- order(tapply(Clinical$RL, Clinical$Group, median))
Clinical$Group <- factor(Clinical$Group, levels=levels(Clinical$Group)[ids])
boxplot(RL ~ Group, data=Clinical, las=2, col=coliPam[ids], ylim=c(0,1))
ids <- order(tapply(Clinical$RL, Clinical$iC10, median))
Clinical$iC10 <- factor(Clinical$iC10, levels=levels(Clinical$iC10)[ids])
boxplot(RL ~ iC10, data=Clinical, las=2, col=coliClust[ids], ylim=c(0,1))
mtext("Probability of DR/Death 10 years after LRR", side=2, outer=T, line=1, cex=1.5)
dev.off()

