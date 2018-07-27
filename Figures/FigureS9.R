rm(list=ls())
cols.er <- c("red", "olivedrab")

coliClust <- c('#FF5500', '#00EE76', '#CD3278','#00C5CD', '#B5D0D2', '#8B0000',
               '#FFFF40', '#0000CD', '#FFAA00', '#EE82EE', '#7D26CD')
coliPam <- c("#E41A1C", "#FB9A99", "#1F78B4", "#A6CEE3", "#66A61E")
coliGroup <- c("lightskyblue", "lightsalmon", "plum2", "wheat4")
Rec <- read.table("../../TableS7.txt", header=TRUE, sep="\t", stringsAsFactors=F)
Rec <- Rec[-which(Rec$Stage==4),]
## We remove benign, DCIS or PHYL
bad.hist <- which(Rec$Histological.Type %in% c("BENIGN", "PHYL"))
if (length(bad.hist)>0) Rec <- Rec[-bad.hist,]
Rec$TDR[which(Rec$TDR==0)] <- 0.1
Rec$TLR[which(Rec$TLR==0)] <- 0.1

Clinical <- read.table(file="../../TableS6.txt", header=T, sep="\t", quote="", comment.char="", stringsAsFactors=FALSE)
Clinical <- Clinical[,c('METABRIC.ID', 'ER.Expr', 'Her2.Expr', 'Pam50Subtype', 'iC10')]
Rec <- merge(Rec, Clinical, all.x=T, sort=F)

Rec$ER.Status <- factor(Rec$ER.Status,
                             levels=c("neg", "pos"), labels=c("ER-", "ER+"))
Rec$ER.Status[which(is.na(Rec$ER.Status) & Rec$ER.Expr=="+")] <- "ER+"
Rec$ER.Status[which(is.na(Rec$ER.Status) & Rec$ER.Expr=="-")] <- "ER-"
Rec$Group <- NA
Rec$Group[which(Rec$ER.Status=="ER-" & Rec$Her2.Expr=="-")] <- "ER-/HER2-"
Rec$Group[which(Rec$ER.Status=="ER-" & Rec$Her2.Expr=="+")] <- "ER-/HER2+"
Rec$Group[which(Rec$ER.Status=="ER+" & Rec$Her2.Expr=="+")] <- "ER+/HER2+"
Rec$Group[which(Rec$ER.Status=="ER+" & Rec$Her2.Expr=="-")] <- "ER+/HER2-"
Rec$Group <- factor(Rec$Group)



Rec <- Rec[which(Rec$TYPE.RELAPSE=="DISTANT"),]
Rec <- Rec[,c('METABRIC.ID', 'SITE', 'TIME.RELAPSE', 'ER.Status', 'Group',
              'Pam50Subtype', 'iC10')]
Times <- sapply(split(Rec, Rec$METABRIC.ID),
       function(x) nrow(x))
Times <- data.frame(METABRIC.ID=names(Times), N=Times)
Times <- merge(Times, unique(Rec[,c(1, 4, 5, 6, 7)]))
tapply(Times$N, Times$ER.Status, mean)
tapply(Times$N, Times$Group, mean)
tapply(Times$N, Times$Pam50Subtype, mean)
tapply(Times$N, Times$iC10, mean)
Times$N <- factor(Times$N, levels=1:max(Times$N))
Times$ER.Status <- factor(Times$ER.Status,
                          levels=c("ER+", "ER-"))
names(cols.er) <- levels(Times$ER.Status)
Times$Group <- factor(Times$Group,
                      levels=c("ER+/HER2+", "ER+/HER2-" ,
                          "ER-/HER2+", "ER-/HER2-"))
names(coliGroup) <- levels(Times$Group)
Times$Pam50Subtype <- factor(Times$Pam50Subtype,
                             levels=c("Basal", "Her2",
                                 "LumA", "LumB", "Normal"))
names(coliPam) <- levels(Times$Pam50Subtype)
Times$iC10 <- factor(Times$iC10,
                     levels=c(1:3, "4ER+", "4ER-", 5:10))
names(coliClust) <- levels(Times$iC10)
pdf("FigureS9.pdf", width=12, height=9)
par(oma=c(2, 3, 2, 2), mar=c(3, 1, 5, 2), mfrow=c(4, 6), xaxs="i", yaxs="i")
for (i in levels(Times$ER.Status)) {

    barplot(prop.table(table(Times$N[which(Times$ER.Status==i)])), col=cols.er[i], main=i, ylim=c(0, 0.5))
    title(main="Clinical/histopathologic subgroups", line=-1, outer=T, adj=0, cex.main=1.5)
}
for (i in levels(Times$Group)[c(2, 4, 1, 3)]) {
    barplot(prop.table(table(Times$N[which(Times$Group==i)])), col=coliGroup[i], main=i, ylim=c(0, 0.5))
}
library(rms)
x1 <- c(tapply(Times$N, Times$ER.Status, function(x) mean(as.numeric(x))),
        tapply(Times$N, Times$Group, function(x) mean(as.numeric(x))),
        tapply(Times$N, Times$Pam50Subtype, function(x) mean(as.numeric(x))))
x2 <- c(tapply(Times$N, Times$ER.Status, function(x) sd(as.numeric(x))/sqrt(length(x))),
        tapply(Times$N, Times$Group, function(x) sd(as.numeric(x))/sqrt(length(x))),
        tapply(Times$N, Times$Pam50Subtype, function(x) sd(as.numeric(x))/sqrt(length(x))))
all.cols <- c(cols.er, coliGroup, coliPam)
errbar(1:11, x1, x1+x2, x1-x2, add=F, col=all.cols,
       errbar.col=all.cols, axes=F, ylab="Number of metastases",
       cex=1.5, xlab="", ylim=c(0, 4), lwd=1.5, xlim=c(0.5, 11.5))
axis(2)
axis(1, at=1:11, c(levels(Times$ER.Status), levels(Times$Group),
                   levels(Times$Pam50Subtype)), cex.axis=0.5, las=2)
title("Number of metastases", cex.main=0.8, line=0)
title(main="Pam 50 subtypes", line=3, outer=F, adj=0, cex.main=1.5)

for (i in levels(Times$Pam50Subtype)[c(1, 3, 4, 5, 2)]) {
    barplot(prop.table(table(Times$N[which(Times$Pam50Subtype==i)])), col=coliPam[i], main=i, ylim=c(0, 0.5))
}


for (i in levels(Times$iC10)) {
    barplot(prop.table(table(Times$N[which(Times$iC10==i)])), col=coliClust[i], main="", ylim=c(0, 0.5))
     title(i, line=0.5)
    if (i=="1") title(main="Integrative clusters", line=3, outer=F, adj=0, cex.main=1.5)
}

library(rms)
errbar(1:11, tapply(Times$N, Times$iC10, function(x) mean(as.numeric(x))),
       tapply(Times$N, Times$iC10, function(x) mean(as.numeric(x))) +
           tapply(Times$N, Times$iC10, function(x) sd(as.numeric(x))/sqrt(length(x))), tapply(Times$N, Times$iC10, function(x) mean(as.numeric(x))) -
           tapply(Times$N, Times$iC10, function(x) sd(as.numeric(x))/sqrt(length(x))), add=F, col=coliClust,
       errbar.col=coliClust, axes=F, ylab="Number of metastases",
       cex=1.5, xlab="", ylim=c(0, 4), lwd=1.5, xlim=c(0.5, 11.5))
title("Number of metastases", cex.main=0.8)
axis(2)
axis(1, at=1:11, levels(Times$iC10), cex.axis=0.9, las=2)

dev.off()


