rm(list=ls())
source("predict.R")
source("getbrcapredProbs.R")

plotCalibrationPredict <- function(df, model) {
all.res <- NULL
    pos.x <- 0.3
        pos.y <- 0.9
    par(mfrow=c(2,4), oma=c(0, 0, 2, 0))
res <- cor(df$BDy02,df$CD.2)
all.res <- c(all.res, res)
plot(BDy02 ~ CD.2, data=df, ylab="predict", xlab="brcarepred",
     main="D/C at 2 years", pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,1))
text(pos.x, pos.y, paste("r=", round(res, 2)))
abline(a=0, b=1)
res <- cor(df$BDy05,df$CD.5)
all.res <- c(all.res, res)
plot(BDy05 ~ CD.5, data=df, ylab="predict", xlab="brcarepred",
     main="C/D at 5 years", pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,1))
text(pos.x, pos.y, paste("r=", round(res, 2)))
abline(a=0, b=1)
res <- cor(df$BDy10,df$CD.10)
all.res <- c(all.res, res)
plot(BDy10 ~ CD.10, data=df, ylab="predict", xlab="brcarepred",
     main="C/D at 10 years", pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,1))
text(pos.x, pos.y, paste("r=", round(res, 2)))
abline(a=0, b=1)
res <- cor(df$BDy15,df$CD.15)
all.res <- c(all.res, res)
plot(BDy15 ~ CD.15, data=df, ylab="predict", xlab="brcarepred",
     main="C/D at 15 years", pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,1))
text(pos.x, pos.y, paste("r=", round(res, 2)))
abline(a=0, b=1)
res <- cor(df$ODy02,df$OD.2)
all.res <- c(all.res, res)
plot(ODy02 ~ OD.2, data=df, ylab="predict", xlab="brcarepred",
     main="O/D at 2 years", pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,1))
text(pos.x, pos.y, paste("r=", round(res, 2)))
abline(a=0, b=1)
res <- cor(df$ODy05,df$OD.5)
all.res <- c(all.res, res)
plot(ODy05 ~ OD.5, data=df, ylab="predict", xlab="brcarepred",
     main="O/D at 5 years", pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,1))
text(pos.x, pos.y, paste("r=", round(res, 2)))
abline(a=0, b=1)
res <- cor(df$ODy10,df$OD.10)
all.res <- c(all.res, res)
plot(ODy10 ~ OD.10, data=df, ylab="predict", xlab="brcarepred",
     main="O/D at 10 years", pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,1))
text(pos.x, pos.y, paste("r=", round(res, 2)))
abline(a=0, b=1)
res <- cor(df$ODy15,df$OD.15)
all.res <- c(all.res, res)
plot(ODy15 ~ OD.15, data=df, ylab="predict", xlab="brcarepred",
     main="O/D at 15 years", pch=19, cex=0.7, xlim=c(0,1), ylim=c(0,1))
text(pos.x, pos.y, paste("r=", round(res, 2)))
abline(a=0, b=1)
mtext(model, outer=T)
all.res <- data.frame(Year=rep(c(2, 5, 10, 15), 2), Event=rep(c("C/D", "O/D"), c(4, 4)), cor=all.res)
all.res
}

########################################################################
########################################################################
########################################################################
## Calibration predict vs brcarepred
########################################################################
########################################################################
########################################################################

Clinical <- read.table("../../TableS6.txt", header=T, sep="\t")
df <- data.frame(METABRIC.ID=Clinical$METABRIC.ID, age.start=Clinical$Age.At.Diagnosis, size=Clinical$Size, grade=Clinical$Grade,
                 nodes=Clinical$Lymph.Nodes.Positive, er=1 * (Clinical$ER.Status=="pos"), her2=1 * (Clinical$Her2.Expr=="+"),
                 ki67=9, horm= 1 * (Clinical$HT!="NO/NA"), traz=0, bis=0, generation=1, screen=0, T=Clinical$T, DeathBreast=Clinical$DeathBreast, Death=Clinical$Death)
df <- na.omit(df)

df <- PPredict(df)

library(rms)
all.bb <-  NULL
pdf("predictCalibration.pdf", width=12, height=9)
pred <- getbrcarepredProbs("../Models/ER_AllProbs.RData")
X <- merge(df, pred$pred.S)
bb <-  plotCalibrationPredict(X, model="ER Model")
all.bb <- rbind(all.bb, data.frame(Model="ER", bb))
pred <- getbrcarepredProbs("../Models/FourGroups_AllProbs.RData")
X <- merge(df, pred$pred.S)
bb <- plotCalibrationPredict(X, model="IHC Model")
all.bb <- rbind(all.bb, data.frame(Model="IHC", bb))
pred <- getbrcarepredProbs("../Models/Pam50_AllProbs.RData")
X <- merge(df, pred$pred.S)
bb <- plotCalibrationPredict(X, model="Pam50 Model")
all.bb <- rbind(all.bb, data.frame(Model="Pam50", bb))
pred <- getbrcarepredProbs("../Models/IntClust_AllProbs.RData")
X <- merge(df, pred$pred.S)
bb <- plotCalibrationPredict(X, model="IntClust Model")
all.bb <- rbind(all.bb, data.frame(Model="IntClust", bb))
## plotCalibrationPredict(subset(X, er==1), model="IntClust Model (ER+ Patients)")
## plotCalibrationPredict(subset(X, er==0), model="IntClust Model (ER- Patients)")
dev.off()
write.table(all.bb, file="ExtCalibration.txt", row.names=F, sep="\t", quote=F)
## Now cindex of predict

predict.cindex <- NULL

tmp <- survConcordance(Surv(T, DeathBreast) ~ BDy02, data=df)
predict.cindex <- rbind(predict.cindex,
                        data.frame(Model="Predict",
                                   Type="C/D",
                                   Year=2,
                                   cindex=tmp$concordance,
                                   cindex.se=tmp$std.err))
tmp <- survConcordance(Surv(T, Death) ~ I(ODy02 + BDy02), data=df)
predict.cindex <- rbind(predict.cindex,
                        data.frame(Model="Predict",
                                   Type="Death",
                                   Year=2,
                                   cindex=tmp$concordance,
                                   cindex.se=tmp$std.err))
tmp <- survConcordance(Surv(T, 1 * (Death==1 & DeathBreast==0)) ~ ODy02, data=df)
predict.cindex <- rbind(predict.cindex,
                        data.frame(Model="Predict",
                                   Type="O/D",
                                   Year=2,
                                   cindex=tmp$concordance,
                                   cindex.se=tmp$std.err))
tmp <- survConcordance(Surv(T, DeathBreast) ~ BDy05, data=df)
predict.cindex <- rbind(predict.cindex,
                        data.frame(Model="Predict",
                                   Type="C/D",
                                   Year=5,
                                   cindex=tmp$concordance,
                                   cindex.se=tmp$std.err))
tmp <- survConcordance(Surv(T, Death) ~ I(ODy05 + BDy05), data=df)
predict.cindex <- rbind(predict.cindex,
                        data.frame(Model="Predict",
                                   Type="Death",
                                   Year=5,
                                   cindex=tmp$concordance,
                                   cindex.se=tmp$std.err))
tmp <- survConcordance(Surv(T, 1 * (Death==1 & DeathBreast==0)) ~ ODy05, data=df)
predict.cindex <- rbind(predict.cindex,
                        data.frame(Model="Predict",
                                   Type="O/D",
                                   Year=5,
                                   cindex=tmp$concordance,
                                   cindex.se=tmp$std.err))
tmp <- survConcordance(Surv(T, DeathBreast) ~ BDy10, data=df)
predict.cindex <- rbind(predict.cindex,
                        data.frame(Model="Predict",
                                   Type="C/D",
                                   Year=10,
                                   cindex=tmp$concordance,
                                   cindex.se=tmp$std.err))
tmp <- survConcordance(Surv(T, Death) ~ I(ODy10 + BDy10), data=df)
predict.cindex <- rbind(predict.cindex,
                        data.frame(Model="Predict",
                                   Type="Death",
                                   Year=10,
                                   cindex=tmp$concordance,
                                   cindex.se=tmp$std.err))
tmp <- survConcordance(Surv(T, 1 * (Death==1 & DeathBreast==0)) ~ ODy10, data=df)
predict.cindex <- rbind(predict.cindex,
                        data.frame(Model="Predict",
                                   Type="O/D",
                                   Year=10,
                                   cindex=tmp$concordance,
                                   cindex.se=tmp$std.err))
tmp <- survConcordance(Surv(T, DeathBreast) ~ BDy15, data=df)
predict.cindex <- rbind(predict.cindex,
                        data.frame(Model="Predict",
                                   Type="C/D",
                                   Year=15,
                                   cindex=tmp$concordance,
                                   cindex.se=tmp$std.err))
tmp <- survConcordance(Surv(T, Death) ~ I(ODy15 + BDy15), data=df)
predict.cindex <- rbind(predict.cindex,
                        data.frame(Model="Predict",
                                   Type="Death",
                                   Year=15,
                                   cindex=tmp$concordance,
                                   cindex.se=tmp$std.err))
tmp <- survConcordance(Surv(T, 1 * (Death==1 & DeathBreast==0)) ~ ODy15, data=df)
predict.cindex <- rbind(predict.cindex,
                        data.frame(Model="Predict",
                                   Type="O/D",
                                   Year=15,
                                   cindex=tmp$concordance,
                                   cindex.se=tmp$std.err))
save(predict.cindex, file="predictCindex.RData")


########################################################################
########################################################################
## Calibration predict vs brcarepred: New Metabric samples
########################################################################
########################################################################
########################################################################

pdf("predictCalibrationNewMetabricSamples.pdf", width=12, height=9)
tmp1 <- read.table("~/Documents/Projects/Metastasis/PaperRevision/VersionAGE/Validation/ValidationMetabricExp.txt", header=T, sep="\t")
tmp2 <- read.table("~/Documents/Projects/Metastasis/PaperRevision/VersionAGE/Validation/ValidationMetabricCN.txt", header=T, sep="\t")
tmp2 <- tmp2[-which(tmp2$METABRIC.ID %in% tmp1$METABRIC.ID),]

Clinical <- tmp1
df <- data.frame(METABRIC.ID=Clinical$METABRIC.ID, age.start=Clinical$Age.At.Diagnosis, size=Clinical$Size, grade=Clinical$Grade,
                 nodes=Clinical$Lymph.Nodes.Positive, er=1 * (Clinical$ER.Status=="pos"), her2=9,
                 ki67=9, horm= 1 * (Clinical$HT!="NO/NA"), traz=0, bis=0, generation=1, screen=0, T=Clinical$T, DeathBreast=Clinical$DeathBreast, Death=Clinical$Death)
df <- na.omit(df)

df <- PPredict(df)

pred <- getbrcarepredProbs("../Models/Validation/ValidationEXPMETABRIC.RData")
X <- merge(df, pred$pred.S)
plotCalibrationPredict(X, model="IntClust Model New Metabric Expression Samples")
pred <- getbrcarepredProbs("../Models/Validation/ValidationCNMETABRIC.RData")
Clinical <- tmp2
df2 <- data.frame(METABRIC.ID=Clinical$METABRIC.ID, age.start=Clinical$Age.At.Diagnosis, size=Clinical$Size, grade=Clinical$Grade,
                 nodes=Clinical$Lymph.Nodes.Positive, er=1 * (Clinical$ER.Status=="pos"), her2=9,
                 ki67=9, horm= 1 * (Clinical$HT!="NO/NA"), traz=0, bis=0, generation=1, screen=0, T=Clinical$T, DeathBreast=Clinical$DeathBreast, Death=Clinical$Death)
df2 <- na.omit(df2)

df2 <- PPredict(df2)

X <- merge(df2, pred$pred.S)
plotCalibrationPredict(X, model="IntClust Model New Metabric Copy Number Samples")
dev.off()


df <- rbind(df, df2)


predict.cindex <- NULL

tmp <- survConcordance(Surv(T, DeathBreast) ~ BDy02, data=df)
predict.cindex <- rbind(predict.cindex,
                        data.frame(Model="Predict",
                                   Type="C/D",
                                   Year=2,
                                   cindex=tmp$concordance,
                                   cindex.se=tmp$std.err))
tmp <- survConcordance(Surv(T, Death) ~ I(ODy02 + BDy02), data=df)
predict.cindex <- rbind(predict.cindex,
                        data.frame(Model="Predict",
                                   Type="Death",
                                   Year=2,
                                   cindex=tmp$concordance,
                                   cindex.se=tmp$std.err))
tmp <- survConcordance(Surv(T, 1 * (Death==1 & DeathBreast==0)) ~ ODy02, data=df)
predict.cindex <- rbind(predict.cindex,
                        data.frame(Model="Predict",
                                   Type="O/D",
                                   Year=2,
                                   cindex=tmp$concordance,
                                   cindex.se=tmp$std.err))
tmp <- survConcordance(Surv(T, DeathBreast) ~ BDy05, data=df)
predict.cindex <- rbind(predict.cindex,
                        data.frame(Model="Predict",
                                   Type="C/D",
                                   Year=5,
                                   cindex=tmp$concordance,
                                   cindex.se=tmp$std.err))
tmp <- survConcordance(Surv(T, Death) ~ I(ODy05 + BDy05), data=df)
predict.cindex <- rbind(predict.cindex,
                        data.frame(Model="Predict",
                                   Type="Death",
                                   Year=5,
                                   cindex=tmp$concordance,
                                   cindex.se=tmp$std.err))
tmp <- survConcordance(Surv(T, 1 * (Death==1 & DeathBreast==0)) ~ ODy05, data=df)
predict.cindex <- rbind(predict.cindex,
                        data.frame(Model="Predict",
                                   Type="O/D",
                                   Year=5,
                                   cindex=tmp$concordance,
                                   cindex.se=tmp$std.err))
tmp <- survConcordance(Surv(T, DeathBreast) ~ BDy10, data=df)
predict.cindex <- rbind(predict.cindex,
                        data.frame(Model="Predict",
                                   Type="C/D",
                                   Year=10,
                                   cindex=tmp$concordance,
                                   cindex.se=tmp$std.err))
tmp <- survConcordance(Surv(T, Death) ~ I(ODy10 + BDy10), data=df)
predict.cindex <- rbind(predict.cindex,
                        data.frame(Model="Predict",
                                   Type="Death",
                                   Year=10,
                                   cindex=tmp$concordance,
                                   cindex.se=tmp$std.err))
tmp <- survConcordance(Surv(T, 1 * (Death==1 & DeathBreast==0)) ~ ODy10, data=df)
predict.cindex <- rbind(predict.cindex,
                        data.frame(Model="Predict",
                                   Type="O/D",
                                   Year=10,
                                   cindex=tmp$concordance,
                                   cindex.se=tmp$std.err))
tmp <- survConcordance(Surv(T, DeathBreast) ~ BDy15, data=df)
predict.cindex <- rbind(predict.cindex,
                        data.frame(Model="Predict",
                                   Type="C/D",
                                   Year=15,
                                   cindex=tmp$concordance,
                                   cindex.se=tmp$std.err))
tmp <- survConcordance(Surv(T, Death) ~ I(ODy15 + BDy15), data=df)
predict.cindex <- rbind(predict.cindex,
                        data.frame(Model="Predict",
                                   Type="Death",
                                   Year=15,
                                   cindex=tmp$concordance,
                                   cindex.se=tmp$std.err))
tmp <- survConcordance(Surv(T, 1 * (Death==1 & DeathBreast==0)) ~ ODy15, data=df)
predict.cindex <- rbind(predict.cindex,
                        data.frame(Model="Predict",
                                   Type="O/D",
                                   Year=15,
                                   cindex=tmp$concordance,
                                   cindex.se=tmp$std.err))
save(predict.cindex, file="predictCindexNewMetabric.RData")

