library(survival)
pdf("FigureS4.pdf", width=8, height=8)

coliClust <- c('#FF5500', '#00EE76', '#CD3278','#00C5CD', '#B5D0D2', '#8B0000',
               '#FFFF40', '#0000CD', '#FFAA00', '#EE82EE', '#7D26CD')

TP <- list()
for (i in 1:11) {
    load(paste0("./Models/TP/TP_IntClust_", i, ".RData"))
    times <- pt$S$Times
    TP[[i]] <- I(pt$S$p.s.l.d +
                        pt$S$p.s.d + pt$S$p.s.l.c +
                            pt$S$p.s.l.d.c + pt$S$p.s.d.c +
                                pt$S$p.s.c)
}
pt <-  TP
rm(TP)


plot(times, rep(0, length(times)), type="n",
     ylim=c(0,1), xlim=c(0, 20), xlab="Years",
     ylab="Probability of relapse")
for (id in 1:11) {
    lines(times, pt[[id]],
          lwd=2, col=coliClust[id])
    }
legend("topleft", col=coliClust, lwd=2, legend=paste0("IntClust", c(1:3, "4ER+",
                                            "4ER-", 5:10)), bty="n")
dev.off()
