rm(list=ls())
par(oma=c(2, 3, 2, 2), mar=c(3, 1, 4, 2), mfrow=c(4,6), xaxs="i", yaxs="i")

library(survival)
TP <- list()
for (i in c("NEG", "POS") ) {
    load(paste0("../Models/TP/ERTP_", i, ".RData"))
    times <- pt$DR$Times
    TP[[i]] <- list(pt$S$p.s, pt$S$p.s.l, I(pt$S$p.s.l.d +
                        pt$S$p.s.d), I(pt$S$p.s.l.c +
                            pt$S$p.s.l.d.c + pt$S$p.s.d.c +
                                pt$S$p.s.c),
                    I(pt$S$p.s.l.o + pt$S$p.s.l.d.o +
                          pt$S$p.s.d.o + pt$S$p.s.o),
                    pt$LR$p.l, pt$LR$p.l.d,
                    I(pt$LR$p.l.c + pt$LR$p.l.d.c),
                 I(pt$LR$p.l.o + pt$LR$p.l.d.o),
                    pt$DR$p.d, pt$DR$p.d.c,
                 pt$DR$p.d.o)
}
pt <-  TP
rm(TP)

load(file="../Models/Bootstraps/BootstrapModel_ER_Results.RData")


boot.TP <- list()
for (i in 1:2) {
    boot.TP[[i]] <- list(ER.probs[[i]]$S$p.s, ER.probs[[i]]$S$p.s.l, I(ER.probs[[i]]$S$p.s.l.d +
                        ER.probs[[i]]$S$p.s.d), I(ER.probs[[i]]$S$p.s.l.c +
                            ER.probs[[i]]$S$p.s.l.d.c + ER.probs[[i]]$S$p.s.d.c +
                                ER.probs[[i]]$S$p.s.c),
                    I(ER.probs[[i]]$S$p.s.l.o + ER.probs[[i]]$S$p.s.l.d.o +
                          ER.probs[[i]]$S$p.s.d.o + ER.probs[[i]]$S$p.s.o),
                    ER.probs[[i]]$LR$p.l, ER.probs[[i]]$LR$p.l.d,
                    I(ER.probs[[i]]$LR$p.l.c + ER.probs[[i]]$LR$p.l.d.c),
                 I(ER.probs[[i]]$LR$p.l.o + ER.probs[[i]]$LR$p.l.d.o),
                    ER.probs[[i]]$DR$p.d, ER.probs[[i]]$DR$p.d.c,
                 ER.probs[[i]]$DR$p.d.o)
}

ER.mean <- lapply(boot.TP, function(x) lapply(x, function(y) apply(y, 1, mean)))
ER.sd <- lapply(boot.TP, function(x) lapply(x, function(y) apply(y, 1, sd)))

names(ER.sd) <- names(pt)


cols <- c("darkblue", "red", "olivedrab")
cols <- rep(cols, 2)
plot(times, rep(0, length(times)), type="n",
     ylim=c(0,1), xlim=c(0, 20), xlab="", ylab="")
title(main="ER+", line=1)
title(main="a. IHC subtypes", line=-1, outer=T, adj=0, cex.main=1.5)
for (k in 10:12) {
    lines(times, pt$POS[[k]], lwd=2, col=cols[k-9], cex.axis=3)
    polygon(c(times, rev(times), times[1]),
            c(pt$POS[[k]]- 1.96*ER.sd$POS[[k]],
              rev(pt$POS[[k]]+ 1.96*ER.sd$POS[[k]]),
              pt$POS[[k]][1]- 1.96*ER.sd$POS[[k]][1]),
            col=adjustcolor(cols[k-9], alpha=0.2), border=NA)
    }
plot(times, rep(0, length(times)), type="n",
     ylim=c(0,1), xlim=c(0, 20), xlab="", ylab="")
title(main="ER-", line=1)
for (k in 10:12) {
    lines(times, pt$NEG[[k]], lwd=2, col=cols[k-9])
    polygon(c(times, rev(times), times[1]),
            c(pt$NEG[[k]]- 1.96*ER.sd$NEG[[k]],
              rev(pt$NEG[[k]]+ 1.96*ER.sd$NEG[[k]]),
              pt$NEG[[k]][1]- 1.96*ER.sd$NEG[[k]][1]),
            col=adjustcolor(cols[k-9], alpha=0.2), border=NA)
    }

S <- c(2, 5, 10, 15, 20)
TP <- list()
for (i in 1:4) {
    load(paste0("../Models/TP/ERTP_4GROUPS_", i, ".RData"))
    times <- pt$DR$Times
    TP[[i]] <- list(pt$S$p.s, pt$S$p.s.l, I(pt$S$p.s.l.d +
                        pt$S$p.s.d), I(pt$S$p.s.l.c +
                            pt$S$p.s.l.d.c + pt$S$p.s.d.c +
                                pt$S$p.s.c),
                    I(pt$S$p.s.l.o + pt$S$p.s.l.d.o +
                          pt$S$p.s.d.o + pt$S$p.s.o),
                    pt$LR$p.l, pt$LR$p.l.d,
                    I(pt$LR$p.l.c + pt$LR$p.l.d.c),
                 I(pt$LR$p.l.o + pt$LR$p.l.d.o),
                    pt$DR$p.d, pt$DR$p.d.c,
                 pt$DR$p.d.o)
}
pt <-  TP
rm(TP)
library(survival)
load(file="../Models/Bootstraps/BootstrapModel_4Groups_Results.RData")
boot.TP <- list()
for (i in 1:4) {
    boot.TP[[i]] <- list(FOURGROUPS.probs[[i]]$S$p.s, FOURGROUPS.probs[[i]]$S$p.s.l, I(FOURGROUPS.probs[[i]]$S$p.s.l.d +
                        FOURGROUPS.probs[[i]]$S$p.s.d), I(FOURGROUPS.probs[[i]]$S$p.s.l.c +
                            FOURGROUPS.probs[[i]]$S$p.s.l.d.c + FOURGROUPS.probs[[i]]$S$p.s.d.c +
                                FOURGROUPS.probs[[i]]$S$p.s.c),
                    I(FOURGROUPS.probs[[i]]$S$p.s.l.o + FOURGROUPS.probs[[i]]$S$p.s.l.d.o +
                          FOURGROUPS.probs[[i]]$S$p.s.d.o + FOURGROUPS.probs[[i]]$S$p.s.o),
                    FOURGROUPS.probs[[i]]$LR$p.l, FOURGROUPS.probs[[i]]$LR$p.l.d,
                    I(FOURGROUPS.probs[[i]]$LR$p.l.c + FOURGROUPS.probs[[i]]$LR$p.l.d.c),
                 I(FOURGROUPS.probs[[i]]$LR$p.l.o + FOURGROUPS.probs[[i]]$LR$p.l.d.o),
                    FOURGROUPS.probs[[i]]$DR$p.d, FOURGROUPS.probs[[i]]$DR$p.d.c,
                 FOURGROUPS.probs[[i]]$DR$p.d.o)
}

BOOT.mean <- lapply(boot.TP, function(x) lapply(x, function(y) apply(y, 1, mean)))
BOOT.sd <- lapply(boot.TP, function(x) lapply(x, function(y) apply(y, 1, sd)))

names(pt) <- c("ER+/HER2+", "ER+/HER2-", "ER-/HER2+", "ER-/HER2-")
labs <- names(pt)

cols <- c("darkblue", "red", "olivedrab")
cols <- rep(cols, 2)
for (id in c(2, 4, 1, 3)) {
    plot(times, rep(0, length(times)), type="n",
         ylim=c(0,1), xlim=c(0, 20), xlab="", ylab="")
    title(main=paste(labs[id]), line=1)
for (k in 10:12) {
    lines(times, pt[[id]][[k]], lwd=2, col=cols[k-6], cex.axis=3)
    polygon(c(times, rev(times), times[1]),
            c(pt[[id]][[k]]- 1.96*BOOT.sd[[id]][[k]],
              rev(pt[[id]][[k]]+ 1.96*BOOT.sd[[id]][[k]]),
              pt[[id]][[k]][1]- 1.96*BOOT.sd[[id]][[k]][1]),
            col=adjustcolor(cols[k-9], alpha=0.2), border=NA)
    }
}

TP <- list()
for (i in 1:5) {
    load(paste0("../Models/TP/TP_Pam50_", i, ".RData"))
    times <- pt$DR$Times
    TP[[i]] <- list(pt$S$p.s, pt$S$p.s.l, I(pt$S$p.s.l.d +
                        pt$S$p.s.d), I(pt$S$p.s.l.c +
                            pt$S$p.s.l.d.c + pt$S$p.s.d.c +
                                pt$S$p.s.c),
                    I(pt$S$p.s.l.o + pt$S$p.s.l.d.o +
                          pt$S$p.s.d.o + pt$S$p.s.o),
                    pt$LR$p.l, pt$LR$p.l.d,
                    I(pt$LR$p.l.c + pt$LR$p.l.d.c),
                 I(pt$LR$p.l.o + pt$LR$p.l.d.o),
                    pt$DR$p.d, pt$DR$p.d.c,
                 pt$DR$p.d.o)
}
pt <-  TP
rm(TP)

library(survival)
load(file="../Models/Bootstraps/BootstrapModel_Pam50_Results.RData")

boot.TP <- list()
for (i in 1:5) {
    boot.TP[[i]] <- list(PAM50.probs[[i]]$S$p.s, PAM50.probs[[i]]$S$p.s.l, I(PAM50.probs[[i]]$S$p.s.l.d +
                        PAM50.probs[[i]]$S$p.s.d), I(PAM50.probs[[i]]$S$p.s.l.c +
                            PAM50.probs[[i]]$S$p.s.l.d.c + PAM50.probs[[i]]$S$p.s.d.c +
                                PAM50.probs[[i]]$S$p.s.c),
                    I(PAM50.probs[[i]]$S$p.s.l.o + PAM50.probs[[i]]$S$p.s.l.d.o +
                          PAM50.probs[[i]]$S$p.s.d.o + PAM50.probs[[i]]$S$p.s.o),
                    PAM50.probs[[i]]$LR$p.l, PAM50.probs[[i]]$LR$p.l.d,
                    I(PAM50.probs[[i]]$LR$p.l.c + PAM50.probs[[i]]$LR$p.l.d.c),
                 I(PAM50.probs[[i]]$LR$p.l.o + PAM50.probs[[i]]$LR$p.l.d.o),
                    PAM50.probs[[i]]$DR$p.d, PAM50.probs[[i]]$DR$p.d.c,
                 PAM50.probs[[i]]$DR$p.d.o)
}


BOOT.mean <- lapply(boot.TP, function(x) lapply(x, function(y) apply(y, 1, mean)))
BOOT.sd <- lapply(boot.TP, function(x) lapply(x, function(y) apply(y, 1, sd)))
names(pt) <- c("Basal", "Her2", "LumA", "LumB", "Normal")
labs <- names(pt)



cols <- c("darkblue", "red", "olivedrab")
cols <- rep(cols, 2)
plot(0,0, type="n", xlab="", ylab="", axes=F)
title(main="b. PAM50 subtypes", line=3, outer=F, adj=0, cex.main=1.5)
for (id in c(5, 3,4,1, 2)) {
    plot(times, rep(0, length(times)), type="n",
     ylim=c(0,1), xlim=c(0, 20), xlab="", ylab="")
    title(main=paste(labs[id]), line=1)
for (k in 10:12) {
    lines(times, pt[[id]][[k]], lwd=2, col=cols[k-9], cex.axis=3)
    polygon(c(times, rev(times), times[1]),
            c(pt[[id]][[k]]- 1.96*BOOT.sd[[id]][[k]],
              rev(pt[[id]][[k]]+ 1.96*BOOT.sd[[id]][[k]]),
              pt[[id]][[k]][1]- 1.96*BOOT.sd[[id]][[k]][1]),
            col=adjustcolor(cols[k-9], alpha=0.2), border=NA)
    }
}



TP <- list()
for (i in 1:11) {
    load(paste0("../Models/TP/TP_IntClust_", i, ".RData"))
    times <- pt$DR$Times
    TP[[i]] <- list(pt$S$p.s, pt$S$p.s.l, I(pt$S$p.s.l.d +
                        pt$S$p.s.d), I(pt$S$p.s.l.c +
                            pt$S$p.s.l.d.c + pt$S$p.s.d.c +
                                pt$S$p.s.c),
                    I(pt$S$p.s.l.o + pt$S$p.s.l.d.o +
                          pt$S$p.s.d.o + pt$S$p.s.o),
                    pt$LR$p.l, pt$LR$p.l.d,
                    I(pt$LR$p.l.c + pt$LR$p.l.d.c),
                 I(pt$LR$p.l.o + pt$LR$p.l.d.o),
                    pt$DR$p.d, pt$DR$p.d.c,
                 pt$DR$p.d.o)
}
pt <-  TP
rm(TP)

load(file="../Models/Bootstraps/BootstrapModel_INTCLUST_Results.RData")

boot.TP <- list()
for (i in 1:11) {
    boot.TP[[i]] <- list(INTCLUST.probs[[i]]$S$p.s, INTCLUST.probs[[i]]$S$p.s.l, I(INTCLUST.probs[[i]]$S$p.s.l.d +
                        INTCLUST.probs[[i]]$S$p.s.d), I(INTCLUST.probs[[i]]$S$p.s.l.c +
                            INTCLUST.probs[[i]]$S$p.s.l.d.c + INTCLUST.probs[[i]]$S$p.s.d.c +
                                INTCLUST.probs[[i]]$S$p.s.c),
                    I(INTCLUST.probs[[i]]$S$p.s.l.o + INTCLUST.probs[[i]]$S$p.s.l.d.o +
                          INTCLUST.probs[[i]]$S$p.s.d.o + INTCLUST.probs[[i]]$S$p.s.o),
                    INTCLUST.probs[[i]]$LR$p.l, INTCLUST.probs[[i]]$LR$p.l.d,
                    I(INTCLUST.probs[[i]]$LR$p.l.c + INTCLUST.probs[[i]]$LR$p.l.d.c),
                 I(INTCLUST.probs[[i]]$LR$p.l.o + INTCLUST.probs[[i]]$LR$p.l.d.o),
                    INTCLUST.probs[[i]]$DR$p.d, INTCLUST.probs[[i]]$DR$p.d.c,
                 INTCLUST.probs[[i]]$DR$p.d.o)
}



BOOT.mean <- lapply(boot.TP, function(x) lapply(x, function(y) apply(y, 1, mean)))
BOOT.sd <- lapply(boot.TP, function(x) lapply(x, function(y) apply(y, 1, sd)))
names(pt) <- c(1:3, "4ER+", "4ER-", 5:10)
labs <- names(pt)


cols <- c("darkblue", "red", "olivedrab")
cols <- rep(cols, 2)
for (id in 1:11) {
    plot(times, rep(0, length(times)), type="n",
     ylim=c(0,1), xlim=c(0, 20), xlab="", ylab="")
    title(main=paste(labs[id]), line=1)
    if (id==1) mtext("c. Integrative subtypes", line=3, outer=F, adj=0, cex.main=1.5, font=2)
for (k in 10:12) {
    lines(times, pt[[id]][[k]], lwd=2, col=cols[k-9], cex.axis=3)
    polygon(c(times, rev(times), times[1]),
            c(pt[[id]][[k]]- 1.96*BOOT.sd[[id]][[k]],
              rev(pt[[id]][[k]]+ 1.96*BOOT.sd[[id]][[k]]),
              pt[[id]][[k]][1]- 1.96*BOOT.sd[[id]][[k]][1]),
            col=adjustcolor(cols[k-9], alpha=0.2), border=NA)
    }
}


par(mar=rep(0, 4))
plot(0,0, type="n", axes=F, xlab="", ylab="")
legend("left", fill=cols[1:3], legend=c("DR", "D/C", "D/O"), bty="n", ncol=2,
       cex=1.2)
mtext("Years after distant relapse",1,outer=TRUE, at=0.5, cex=1.3)
mtext("Transition Probability",2,outer=TRUE, line=1, at=0.5, cex=1.3)




