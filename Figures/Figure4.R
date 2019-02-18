## First panel

pdf("Figure4a.pdf", width=6, height=12)

library(multcomp)
library(png)
library(cmprsk)
library(gridBase)
library(grid)
cols.er <- c("red", "olivedrab")

coliClust <- c('#FF5500', '#00EE76', '#CD3278','#00C5CD', '#B5D0D2', '#8B0000',
               '#FFFF40', '#0000CD', '#FFAA00', '#EE82EE', '#7D26CD')
coliPam <- c("#E41A1C", "#FB9A99", "#1F78B4", "#A6CEE3", "#66A61E")
coliClaudin <- c("#E41A1C", "yellow", "#FB9A99", "#1F78B4", "#A6CEE3", "#66A61E")





library(png)
library(survival)
maxpiesize <- unit(1, "inches")
####################################
####################################
####################################
## ALL SAMPLES
####################################
####################################
####################################

Rec <- read.table("../../TableS7.txt", header=TRUE, sep="\t", stringsAsFactors=F)
Rec <- Rec[-which(Rec$Stage==4),]
## We remove benign, DCIS or PHYL
bad.hist <- which(Rec$Histological.Type %in% c("BENIGN", "PHYL"))
if (length(bad.hist)>0) Rec <- Rec[-bad.hist,]
Rec$TDR[which(Rec$TDR==0)] <- 0.1
Rec$TLR[which(Rec$TLR==0)] <- 0.1
Rec$TIME.RELAPSE[which(Rec$TIME.RELAPSE==0)] <- 0.1
Rec$T <- Rec$T/365.25
Rec$TLR <- Rec$TLR/365.25
Rec$TDR <- Rec$TDR/365.25
Rec$TIME.RELAPSE <- Rec$TIME.RELAPSE/365.25
Rec <- Rec[which(Rec$TYPE.RELAPSE=="DISTANT"),]

Clinical <- read.table("../../TableS6.txt", header=TRUE, sep="\t", stringsAsFactors=F)
Clinical <- Clinical[,c('METABRIC.ID', 'iC10')]
Rec <- merge(Rec, Clinical, all.x=T)
Rec <- Rec[,c('METABRIC.ID', 'SITE', 'TIME.RELAPSE', 'iC10', 'ER.Status')]
Rec$SITE <- factor(Rec$SITE)
levels(Rec$SITE) <- list(
    LIVER="LIVER",
    "PULMONARY"=c("LUNG", "PLEURA"),
    "BRAIN/MENINGEAL"=c("BRAIN", "MENINGEAL"),
    "BONE"="BONE",
    "SKIN"="SKIN",
    "MEDIASTINAL"="MEDIASTINAL",
    "PERITONEUM"="PERITONEUM",
    "LNS"="LNS",
    "OTHER"=c("ADRENALS", "PANCREAS", "GI_TRACT", "OVARY",
        "OTHER", "ABDOMEN", "PERICARDIAL", "ASCITES", "EYE",
        "KIDNEY", "BLADDER", "BREAST", "UNSPECIFIED"),
    "SOFT_TISSUES"="SOFT_TISSUES")
Clinical <- Rec[,c('METABRIC.ID', 'iC10', 'ER.Status')]
Clinical <- unique(Clinical)
all.sites <- levels(Rec$SITE)
all.sites <- sort(all.sites)

## color.Sites <-  c("#8DD3C7",  "#FFFFB3",  "#BEBADA",
##                   "#FB8072", "#80B1D3", "#FDB462",
##                   "#B3DE69", "#FCCDE5", "#D9D9D9",
##                   "#BC80BD", "#CCEBC5")




par(mar=c(4, 5, 4, 4), oma=c(1, 2, 1, 1))
oldpar <- par(no.readonly = TRUE)

start.x <- c(-0.9, -0.2,
             -0.9,
             -0.9,
             -0.2,
             -0.2,
             -0.9,
             -0.2,
             -0.2,
             -0.9)

start.y <- c(-0.45, 0.95,
             0.3, 0.8,
             0.45,
             -0.9,
             0.1,
             0.7, -0.1, -0.2)

names(start.x) <- all.sites
names(start.y) <- all.sites
plot(start.x, start.y, xlim = c(-1, 0),
     ylim = c(-1, 1), type = "n", axes=F, xlab="", ylab="")
library(gridBase)
library(lattice)
vps <- baseViewports()
par(new=T)
## get from open source image
silh <- readPNG("Silhouette_Female.png", FALSE, )
silh2 <- as.raster(silh)
silh2[which(silh2=="#E8E9E6FF")] <-
adjustcolor("#FFB6C1", alpha.f=0.6)
rasterImage(silh2, -0.8,-1, -0.2, 1)
size.text <- 1.5
mid.text <- c(-0.6, -0.5,
              -0.55, -0.5, -0.5, -0.5, -0.3,
              -0.5,
              -0.3, -0.3, 0.6)
text(-0.5, 0.9, "Brain/\nMeningeal", cex=size.text, adj=0.5)
text(-0.5, 0.7, "Lymph Nodes", cex=size.text, adj=0.5)
text(-0.5, 0.6, "Pulmonary", cex=size.text, adj=0.5)
text(-0.5, 0.4, "Mediastinal", cex=size.text, adj=0.5)
text(-0.55, 0.3, "Liver", cex=size.text, adj=0.5)
text(-0.6, -0.5, "Bone", cex=size.text, adj=0.5)
text(-0.3, 0, "Skin", cex=size.text, adj=0.5)
text(-0.5, 0.1, "Peritoneum", cex=size.text, adj=0.5)
text(-0.5, -0.2, "Soft Tissues", cex=size.text, adj=0.5)
text(-0.4, -1, "Other", cex=size.text, adj=0.5)

pushViewport(vps$inner, vps$figure,
     vps$plot)


totals <- NULL
res <- list()
all.pvals <- list()
all.coefs <- list()
for (site in all.sites) {
    tmp <- sapply(split(Rec, Rec$METABRIC.ID), function(x) 1 * any(x$SITE==site))
    tmp <-data.frame(METABRIC.ID=names(tmp), tmp)
    colnames(tmp)[2] <- site
    tmp <- merge(Clinical, tmp, all.x=T)
    tmp[,site] <- factor(tmp[,site], levels=c(1, 0))
    tmp[,'ER.Status'] <- factor(tmp[,'ER.Status'],
                                levels=c("pos", "neg"),
                                labels=c("+", "-"))
    K <- contrMat(table(tmp$ER.Status), type="Tukey")
    m <- data.frame(x=tmp$ER.Status, y=tmp[,site])
    m <- glm(y ~ x - 1, data=m, family="binomial")
    m <- glht(m, linfct=K)
    all.pvals[[site]] <- summary(m)$test$pvalues
    all.coefs[[site]] <- summary(m)$test$coef

    tmp <- table(tmp[,site], tmp$ER.Status)
    totals <- c(totals, tmp['1',])
    res[[site]] <- tmp
}

lines.start.x <- c(-0.82,-0.5,
                   -0.82, -0.82,
                   -0.5,
                   -0.5,
                   -0.3,
                   -0.82,
                 -0.3, -0.82)
lines.end.x <- c(-0.6, -0.25,
                 -0.55, -0.5,
                 -0.25,
                 -0.25,
                 -0.25,
                 -0.5, -0.25,
               -0.6)
lines.start.y <- c(-0.45,0.9,
                   0.3, 0.8,
                   0.6,
                   0.4,
                   -1,
                   0.1, 0, -0.2)
lines.end.y <- c(-0.5, 0.95,
                 0.3, 0.7,
                 0.7,
                 0.4,
                 -0.9,
                 0.1, -0.1, -0.2)
for (i in 1:length(all.sites)) {
    pushViewport(viewport(x = unit(start.x[i],
    "native"), y = unit(start.y[i],
    "native"), width = 0.5 *
    maxpiesize, height = 0.7 *
    maxpiesize))
    par(plt = gridPLT(), new = TRUE, xpd=NA)
    feo <- barplot(cbind(prop.table(res[[i]], 2)[,1,drop=F],
                     matrix(c(NA, NA), nrow=2)), axes=F, names.arg=c("", ""),
               col=c(cols.er[1], "white"), cex.axis=1.5, cex.lab=1.5, cex.names=1.5)
    feo <- barplot(cbind(matrix(c(NA, NA), nrow=2),
                         prop.table(res[[i]], 2)[,2,drop=F]),
                   axes=F, names.arg=c("", ""),cex.axis=1.5, cex.lab=1.5, cex.names=1.5,
               col=c(cols.er[2], "white"), add=T)

    text(feo, c(0.9, 0.9), c("+", "-"), cex=1.5)
    if (i==4) {
         text(mean(feo), 1.2, "ER", cex=1.5)
         text(2.8, 0.4, "# Cases", srt=90)
     }
    if (all.pvals[[i]] < 0.05) {
        pch <- ifelse(all.coefs[[i]]>0, 24, 25)
        points(feo[1], 1.1, pch=pch, col="black", bg="black", cex=1)
    }
    axis(1, at=feo, res[[i]]['1',], line=-1, cex.axis=0.7, tick=FALSE)
    popViewport()
    grid.segments(x0=unit(lines.start.x[i], "native"),
                  x1=unit(lines.end.x[i], "native"),
                  y0=unit(lines.start.y[i], "native"),
                  y1=unit(lines.end.y[i], "native"),
                   gp=gpar(lty="dashed", col="grey"))
}

start.x <- c(-1, -0.1,
             -1,
             -1,
             -0.1, -0.1, -1, -0.1,
             -0.1,
             -1)

start.y <- c(-0.45, 0.95,
             0.3, 0.8,
             0.45, -0.9, 0.1,
             0.7, -0.1, -0.2)
names(start.x) <- all.sites
names(start.y) <- all.sites
library(gridBase)
library(lattice)


library(survival)
X <- read.table(file="../../TableS7.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
X <- X[-which(X$Stage==4),]
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
X$SITE <- factor(X$SITE)
levels(X$SITE) <- list("BRAIN/MENINGEAL"=c("BRAIN", "MENINGEAL"), "LNS"="LNS",
                       "PULMONARY"=c("LUNG", "PLEURA"),
                       "MEDIASTINAL"="MEDIASTINAL",
                       "LIVER"="LIVER",
                       "BONE"="BONE",
                       ## "BREAST"="BREAST",
                       "PERITONEUM"="PERITONEUM",
                       "SKIN"="SKIN",
                       "SOFT_TISSUES"="SOFT TISSUES",
                       "OTHER"=c("ABDOMEN", "ADRENALS", "ASCITES","BREAST",
                           "BLADDER",
                           "EYE", "GI_TRACT", "KIDNEY", "UNSPECIFIED",
                           "OTHER", "OVARY", "PANCREAS",
                           "PERICARDIAL"))


## Now we only keep the first relapse of each site

ID <- unique(X$METABRIC.ID)
TOT <- levels(X$SITE)

CR <- NULL
for (i in ID) {
    tmp <- X[which(X$METABRIC.ID==i),]
    for (j in TOT) {
        id <- which(tmp$SITE==j)
        if (length(id)>0) {
            min.time <- tmp$TIME.RELAPSE[id]
            if(length(!is.na(min.time))>1) {
                min.time <- min(min.time, na.rm=T)
            } else {
                min.time <- min.time[1]
            }
            feo <- data.frame(ID=i, time=min.time,
                          status=1,
                              ER=tmp$ER.Status[1],
                              SITE=j)
        } else {
            feo <- data.frame(ID=i, time=tmp$T[1],
                          status=0,
                          ER=tmp$ER.Status[1],
                              SITE=j)
        }
        CR <- rbind(CR, feo)
    }

    feo <- data.frame(ID=i, time=tmp$T[1],
                      status=1 * (tmp$Death[1]==1),
                      ER=tmp$ER.Status[1],
                      SITE="DEATH")
    CR <- rbind(CR, feo)
}
CR$ER <- factor(CR$ER, levels=c("pos", "neg"),
                                labels=c("+", "-"))

m.neg <- survfit(Surv(time, status) ~ strata(SITE), data=CR, subset=ER=="-")
m.pos <- survfit(Surv(time, status) ~ strata(SITE), data=CR, subset=ER=="+")

extra <- ifelse(start.x> -0.4, 0.1, -0.1)
for (i in 1:length(all.sites)) {
    pushViewport(viewport(x = unit(start.x[i]+extra[i],
    "native"), y = unit(start.y[i],
    "native"), width = 1 *
    maxpiesize, height = 0.7 *
    maxpiesize))
    par(plt = gridPLT(), new = TRUE, xpd=NA, mgp=c(3,0.25,0))
    plot(m.neg[which(names(m.neg$strata)==paste0("strata(SITE)=", all.sites[i]))],
         col=cols.er[2], axes=F, lwd=2, fun="event",
         ylim=c(0,1), conf.int=F, xmax=15, cex.axes=1.5, cex.lab=1.5)
    lines(m.pos[which(names(m.neg$strata)==paste0("strata(SITE)=", all.sites[i]))],          col=cols.er[1], axes=F, lwd=2, fun="event",
          ylim=c(0,1), conf.int=F, xmax=15)
    axis(1, at=c(0, 5, 10, 15), labels=c(0, 5, 10, 15),
         cex.axis=0.8, tick=TRUE)
    if (start.x[i] < -0.4) {
            axis(2, cex.axis=0.8, at=c(0,0.5, 1))
        } else {
            axis(4, cex.axis=0.8, at=c(0,0.5, 1))
        }
    if (i==4) {
        title(ylab="Probability", line=1, cex.lab=1.2)
        title(xlab="Years", line=1, cex.lab=1.2)
    }
    box()
    popViewport()
    grid.segments(x0=unit(lines.start.x[i], "native"),
                  x1=unit(lines.end.x[i], "native"),
                  y0=unit(lines.start.y[i], "native"),
                  y1=unit(lines.end.y[i], "native"),
                   gp=gpar(lty="dashed", col="grey"))
}
popViewport(3)
par(oldpar)
dev.off()



rm(list=ls())

## Panel b

library(survival)
library(RColorBrewer)
####################################################################
####################################################################
## ALL EVENTS (When available)
####################################################################
####################################################################

library(survival)
X <- read.table(file="../../TableS7.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
X <- X[-which(X$Stage==4),]

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
                          enum=j)
        AG <- rbind(AG, feo)
        time1 <- time2
        j <- j + 1
    }
    if (time1 < tmp$T[j-1]) {
        feo <- data.frame(ID=i, time1=time1, time2=tmp$T[j-1],
                          status=0,
                          ER=tmp$ER.Status[1],
                          enum=j)
        AG <- rbind(AG, feo)
    }
}
AG$T <- (AG$time2 - AG$time1)/365.25



pdf("Figure4bc.pdf", width=10, height=10)
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
## Panel c
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
                events <- na.omit(events)
                if (length(events)>0) for (ev in events)  Met[,ev] <- 1
            } else {
                feo <- data.frame(ID=i, time1=time1,
                              time2=j,
                              Outcome="CANCERDEATH",
                              status=0,
                              ER=tmp$ER.Status[1], LN=tmp$LN[1],
                              AGE=tmp$AGE[1], GRADE=tmp$GRADE[1],
                              SIZE=tmp$SIZE[1])
            feo <- cbind(feo, Met)
            feo$enum <- enum
            CR2 <- rbind(CR2, feo)
            feo <- data.frame(ID=i, time1=time1,
                              time2=j,
                              Outcome="OTHERDEATH",
                              status=0,
                              ER=tmp$ER.Status[1], LN=tmp$LN[1],
                              AGE=tmp$AGE[1], GRADE=tmp$GRADE[1],
                              SIZE=tmp$SIZE[1])
            feo <- cbind(feo, Met)
            time1 <- j
            events <- as.character(tmp$SITE[which(tmp$TIME.RELAPSE==j)])
                events <- na.omit(events)
                if (length(events)>0) for (ev in events)  Met[,ev] <- 1
            feo$enum <- enum
            CR2 <- rbind(CR2, feo)
            enum <- enum + length(events)
            }
            feo <- data.frame(ID=i, time1=time1,
                          time2=tmp$T[1],
                          Outcome="CANCERDEATH",
                              status=1 * (tmp$DeathBreast[1]==1),
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
XLIM <- c(-0.5, 3)
plot(0,0, type="n", axes=F, xlab="", ylab="", xlim=XLIM,
     ylim=YLIM, xaxs="i", main="", cex.axis=1.7, cex.lab=1.7)
mtext("Log-Hazard Ratio\nER+", side=1, line=4)
points(coef(m3)[ids], 1:K, pch=19, cex=1.5)
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
points(coef(m3)[ids], 1:K, pch=19, cex=1.5)
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


dev.off()



## m3 <- coxph(Surv(time1, time2, status) ~ ER.CD + ER.OT + (LN.CD + GRADE.CD + SIZE.CD + enum.CD +
##                 BRAIN.MENINGEAL.CD + PULMONARY.CD +
##                     LIVER.CD + BONE.CD + OTHER.CD) +
##                         strata(Outcome) + cluster(ID), data=CR2)
## cox.zph(m3)

## m3 <- coxph(Surv(time1, time2, status) ~ ER.CD*(LN.CD + GRADE.CD + SIZE.CD + enum.CD +
##                 BRAIN.MENINGEAL.CD + PULMONARY.CD +
##                     LIVER.CD + BONE.CD + OTHER.CD) +
##                         strata(Outcome) + cluster(ID), data=CR2)
## cox.zph(m3)
