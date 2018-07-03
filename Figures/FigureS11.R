library(grid)
library(multcomp)
library(png)
library(cmprsk)
cols.er <- c("olivedrab", "red")

coliClust <- c('#FF5500', '#00EE76', '#CD3278','#00C5CD', '#B5D0D2', '#8B0000',
               '#FFFF40', '#0000CD', '#FFAA00', '#EE82EE', '#7D26CD')
coliPam <- c("#E41A1C", "#FB9A99", "#1F78B4", "#A6CEE3", "#66A61E")





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


Rec <- read.table("TableS7", header=TRUE, sep="\t", stringsAsFactors=F)
Rec <- Rec[which(Rec$Stage!=4),]
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
Rec <- Rec[,c('METABRIC.ID', 'SITE', 'TIME.RELAPSE', 'Pam50Subtype', 'ER.Status')]
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
Clinical <- Rec[,c('METABRIC.ID', 'Pam50Subtype', 'ER.Status')]
Clinical <- unique(Clinical)
all.sites <- levels(Rec$SITE)
all.sites <- sort(all.sites)

color.Sites <-  c("#8DD3C7",  "#FFFFB3",  "#BEBADA",
                  "#FB8072", "#80B1D3", "#FDB462",
                  "#B3DE69", "#FCCDE5", "#D9D9D9",
                  "#BC80BD", "#CCEBC5")


pdf("FigureS11.pdf", width=12, height=12)
par(mfrow=c(1,2))
oldpar <- par(no.readonly = TRUE)

start.x <- c(-0.9, -0.2,
             ## -0.9,
             -0.9,
             -0.9,
             -0.2,
             -0.2,
             -0.9,
             -0.2,
             -0.2,
             -0.9)

start.y <- c(-0.45, 0.95,
             ## 0.55,
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
silh <- readPNG("Silhouette_Female.png", FALSE, )
silh2 <- as.raster(silh)
silh2[which(silh2=="#E8E9E6FF")] <-
adjustcolor("#FFB6C1", alpha.f=0.6)
rasterImage(silh2, -0.8,-1, -0.2, 1)
size.text <- 1.5
mid.text <- c(-0.6, -0.5,
              ## -0.5,
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


maxpiesize <- unit(1, "inches")
totals <- NULL
all.pvals <- list()
all.coefs <- list()
res <- list()
for (site in all.sites) {
    tmp <- sapply(split(Rec, Rec$METABRIC.ID), function(x) 1 * any(x$SITE==site))
    tmp <-data.frame(METABRIC.ID=names(tmp), tmp)
    colnames(tmp)[2] <- site
    tmp <- merge(Clinical, tmp, all.x=T)
    tmp[,site] <- factor(tmp[,site], levels=c(1, 0))
    tmp[,'iC10'] <- factor(tmp[,'Pam50Subtype'],
                                levels=c("Basal", "Her2", "LumA", "LumB", "Normal"))
    K <- contrMat(table(tmp$iC10), type="GrandMean")
    m <- data.frame(x=tmp$iC10, y=tmp[,site])
    m$y <- as.numeric(as.character(m$y))
    m <- glm(y ~ x - 1, data=m, family="binomial")
    m <- glht(m, linfct=K)
    all.pvals[[site]] <- summary(m)$test$pvalues
    all.coefs[[site]] <- summary(m)$test$coef
    tmp <- table(tmp[,site], tmp$iC10)
    totals <- c(totals, tmp['1',])
    res[[site]] <- tmp
}

lines.start.x <- c(-0.82,-0.5,
                   ## -0.82,
                   -0.82, -0.82,
                   -0.5,
                   -0.5,
                   -0.3,
                   -0.82,
                 -0.3, -0.82)
lines.end.x <- c(-0.6, -0.25,
                 ## -0.5,
                 -0.55, -0.5,
                 -0.25,
                 -0.25,
                 -0.25,
                 -0.5, -0.25,
               -0.6)
lines.start.y <- c(-0.45,0.9,
                   ## 0.55,
                   0.3, 0.8,
                   0.6,
                   0.4,
                   -1,
                   0.1, 0, -0.2)
lines.end.y <- c(-0.5, 0.95,
                 ## 0.5,
                 0.3, 0.7,
                 0.7,
                 0.4,
                 -0.9,
                 0.1, -0.1, -0.2)
bar.cols <- rep("white", 10)
bar.cols[seq(from=1, by=2, length=5)] <- coliPam
for (i in 1:length(all.sites)) {
    pushViewport(viewport(x = unit(start.x[i],
    "native"), y = unit(start.y[i],
    "native"), width = 1.8 *
    maxpiesize, height = 0.5 *
    maxpiesize))
    par(plt = gridPLT(), new = TRUE, xpd=NA)
    feo <- barplot(cbind(prop.table(res[[i]], 2)[,1,drop=F], matrix(rep(NA, 2 * (5-1)), nrow=2)), axes=F, names.arg=rep("", 5),
             col=c(coliPam[1], "white"), cex.axis=1.5, cex.lab=1.5, cex.names=1.5)
    for (j in 2:5) {
        barplot(cbind(matrix(rep(NA, 2 * (j-1)), nrow=2), prop.table(res[[i]], 2)[,j,drop=F], matrix(rep(NA, 2 * (5-j)), nrow=2)),
                cex.axis=1.5, cex.lab=1.5, cex.names=1.5,
                axes=F, names.arg=rep("", 5),
             col=c(coliPam[j], "white"), add=T)
    }
    if (i==4) mtext(at=feo, side=3, text=c("Basal", "Her2", "LumA", "LumB", "Normal"), cex=1, las=2, line=1)
    if (any(all.pvals[[i]] < 0.05)) {
        pch <- ifelse(all.coefs[[i]][all.pvals[[i]]<0.05]>0, 24, 25)
        points(feo[which(all.pvals[[i]]<0.05)],
               rep(1.1, sum(all.pvals[[i]]<0.05)), pch=pch, col="black", bg="black", cex=1)
    }
    axis(1, at=feo, res[[i]]['1',], line=-1, cex.axis=0.7, tick=FALSE)
    popViewport()
    grid.segments(x0=unit(lines.start.x[i], "native"),
                  x1=unit(lines.end.x[i], "native"),
                  y0=unit(lines.start.y[i], "native"),
                  y1=unit(lines.end.y[i], "native"),
                   gp=gpar(lty="dashed", col="grey"))
}
popViewport(3)
par(oldpar)
par(mfg=c(1,2))



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
plot(start.x, start.y, xlim = c(-1, 0),
     ylim = c(-1, 1), type = "n", axes=F, xlab="", ylab="")
library(gridBase)
library(lattice)
vps <- baseViewports()
par(new = FALSE)

## get from open source images
silh <- readPNG("Silhouette_Female.png", FALSE, )
silh2 <- as.raster(silh)
silh2[which(silh2=="#E8E9E6FF")] <-
adjustcolor("#FFB6C1", alpha.f=0.6)
rasterImage(silh2, -0.8,-1, -0.2, 1)
size.text <- 1.5
mid.text <- c(-0.6, -0.5, -0.5, -0.55, -0.5, -0.5, -0.5, -0.3, -0.5, -0.3, -0.3, 0.6)
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


library(survival)
X <- read.table(file="TableS7.txt", header=TRUE, sep="\t", stringsAsFactors=F)
X <- X[which(X$Stage!=4),]
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
                       "PERITONEUM"="PERITONEUM",
                       "SKIN"="SKIN",
                       "SOFT_TISSUES"="SOFT TISSUES",
                       "OTHER"=c("ABDOMEN", "ADRENALS", "ASCITES", "BREAST",
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
                              Pam50=tmp$Pam50Subtype[1],
                              SITE=j)
        } else {
            feo <- data.frame(ID=i, time=tmp$T[1],
                          status=0,
                              Pam50=tmp$Pam50Subtype[1],
                              SITE=j)
        }
        CR <- rbind(CR, feo)
    }

    feo <- data.frame(ID=i, time=tmp$T[1],
                      status=1 * (tmp$Death[1]==1),
                      Pam50=tmp$Pam50Subtype[1],
                      SITE="DEATH")
    CR <- rbind(CR, feo)
}

m <- list()
for (i in c("Basal", "Her2", "LumA", "LumB", "Normal")) {
m[[i]] <- survfit(Surv(time, status) ~ strata(SITE), data=CR, subset=Pam50==i)
}

for (i in 1:length(all.sites)) {
    pushViewport(viewport(x = unit(start.x[i],
    "native"), y = unit(start.y[i],
    "native"), width = 1 *
    maxpiesize, height = 0.7 *
    maxpiesize))
    par(plt = gridPLT(), new = TRUE, xpd=NA, mgp=c(3,0.25,0))
    plot(m[[1]][which(names(m[[1]]$strata)==paste0("strata(SITE)=", all.sites[i]))],
         col=coliPam[1], axes=F, lwd=2, fun="event",
         ylim=c(0,1), conf.int=F, xmax=15)
    for (j in 2:5) {
    lines(m[[j]][which(names(m[[j]]$strata)==paste0("strata(SITE)=", all.sites[i]))],
          col=coliPam[j], axes=F, lwd=2, fun="event",
          ylim=c(0,1), conf.int=F, xmax=15)
}
    axis(2, cex.axis=0.8, at=c(0,0.5, 1))
    axis(1, at=c(0, 5, 10, 15), labels=c(0, 5, 10, 15),
         cex.axis=0.8, tick=TRUE)
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

