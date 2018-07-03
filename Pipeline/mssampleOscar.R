oscar.msfit <- function (object, newdata, variance = TRUE, vartype = c("aalen",
    "greenwood"), trans)
{
    if (!is.null((object$call)$weights) || !is.null(object$weights))
        stop("msfit cannot (yet) compute the result for a weighted model")
    Terms <- terms(object)
    strat <- attr(Terms, "specials")$strata
    cluster <- attr(Terms, "specials")$cluster
    if (length(cluster))
        stop("cluster terms are not supported")
    if (!is.null(attr(object$terms, "specials")$tt))
        stop("msfit cannot yet process coxph models with a tt term")
    resp <- attr(Terms, "variables")[attr(Terms, "response")]
    nvar <- length(object$coefficients)
    score <- exp(object$linear.predictors)
    vartype <- match.arg(vartype)
    if (is.na(vartype))
        stop("Invalid variance type specified")
    has.strata <- !is.null(attr(object$terms, "specials")$strata)
    if (is.null(object$y) || is.null(object[["x"]]) || !is.null(object$call$weights) ||
        (has.strata && is.null(object$strata)) || !is.null(attr(object$terms,
        "offset"))) {
        mf <- model.frame(object)
    }
    else mf <- NULL
    if (is.null(mf))
        y <- object[["y"]]
    else {
        y <- model.response(mf)
        y2 <- object[["y"]]
        print(length(y))
        print(length(y2))
        if (!is.null(y2) && any(as.matrix(y2) != as.matrix(y)))
            stop("Could not reconstruct the y vector")
    }
    if (is.null(object[["x"]]))
        x <- model.matrix(object, data = mf)
    else x <- object[["x"]]
    n <- nrow(y)
    if (n != object$n[1] || nrow(x) != n)
        stop("Failed to reconstruct the original data set")
    type <- attr(y, "type")
    if (type == "counting")
        lasty <- max(y[, 2])
    else if (type == "right")
        lasty <- max(y[, 1])
    else stop("Cannot handle \"", type, "\" type survival data")
    if (is.null(mf))
        offset <- 0
    else {
        offset <- model.offset(mf)
        if (is.null(offset))
            offset <- 0
    }
    Terms <- object$terms
    if (!has.strata)
        strata <- rep(0L, n)
    else {
        stangle <- untangle.specials(Terms, "strata")
        strata <- object$strata
        if (is.null(strata)) {
            if (length(stangle$vars) == 1)
                strata <- mf[[stangle$vars]]
            else strata <- strata(mf[, stangle$vars], shortlabel = TRUE)
        }
    }
    if (has.strata) {
        temp <- attr(Terms, "specials")$strata
        factors <- attr(Terms, "factors")[temp, ]
        strata.interaction <- any(t(factors) * attr(Terms, "order") >
            1)
        if (strata.interaction)
            stop("Interaction terms with strata not supported")
    }
    if (vartype == "greenwood") {
        if (missing(trans))
            stop("argument trans missing; needed for vartype=\"greenwood\"")
        labels <- attr(Terms, "term.labels")
        if (length(labels) != 1)
            stop("Invalid formula for greenwood, ~strata(trans) needed, no covariates allowed")
        if (attr(Terms, "term.labels") != "strata(trans)")
            stop("Invalid formula for greenwood, ~strata(trans) needed, no covariates allowed")
        sf0 <- summary(survfit(object))
        norisk <- sf0$n.risk
        noevent <- sf0$n.event
        sf0 <- data.frame(time = sf0$time, Haz = -log(sf0$surv),
            norisk = norisk, noevent = noevent, trans = as.numeric(sf0$strata))
        allt <- sort(unique(c(sf0$time, lasty)))
        nt <- length(allt)
        K <- nrow(mstate:::to.trans2(trans))
        Haz <- data.frame(time = rep(allt, K), Haz = NA, trans = rep(1:K,
            rep(nt, K)))
        if (variance) {
            tr12 <- data.frame(trans1 = rep(1:K, rep(K, K)),
                trans2 = rep(1:K, K))
            tr12 <- tr12[tr12$trans1 <= tr12$trans2, ]
            varHaz <- data.frame(time = rep(allt, K * (K + 1)/2),
                varHaz = 0, trans1 = rep(tr12$trans1, rep(nt,
                  K * (K + 1)/2)), trans2 = rep(tr12$trans2,
                  rep(nt, K * (K + 1)/2)))
        }
        S <- nrow(trans)
        for (s in 1:S) {
            trs <- trans[s, ]
            trs <- trs[!is.na(trs)]
            ntrs <- length(trs)
            if (ntrs > 0) {
                for (i in 1:ntrs) {
                  trans1 <- trs[i]
                  sf1 <- sf0[sf0$trans == trans1, ]
                  Haz$Haz[(trans1 - 1) * nt + match(sf1$time,
                    allt)] <- sf1$Haz
                  Haz$Haz[(trans1 - 1) * nt + 1:nt] <- mstate:::NAfix(Haz$Haz[(trans1 -
                    1) * nt + 1:nt], subst = 0)
                  if (variance) {
                    varHaz1 <- cumsum((sf1$norisk - sf1$noevent) *
                      sf1$noevent/sf1$norisk^3)
                    varHaz11 <- varHaz[varHaz$trans1 == trans1 &
                      varHaz$trans2 == trans1, ]
                    varHaz11$varHaz <- NA
                    varHaz11$varHaz[match(sf1$time, allt)] <- varHaz1
                    varHaz11$varHaz <- mstate:::NAfix(varHaz11$varHaz,
                      subst = 0)
                    varHaz[varHaz$trans1 == trans1 & varHaz$trans2 ==
                      trans1, ] <- varHaz11
                    if (i < ntrs) {
                      for (j in ((i + 1):ntrs)) {
                        trans2 <- trs[j]
                        sf2 <- sf0[sf0$trans == trans2, ]
                        jointt <- intersect(sf1$time, sf2$time)
                        if (length(jointt) > 0) {
                          varHazij <- rep(NA, length(jointt))
                          ik <- match(jointt, sf1$time)
                          jk <- match(jointt, sf2$time)
                          varHazij <- cumsum(-sf1$noevent[ik] *
                            sf2$noevent[jk]/sf1$norisk[ik]^3)
                          varHaz12 <- varHaz[varHaz$trans1 ==
                            trans1 & varHaz$trans2 == trans2,
                            ]
                          varHaz12$varHaz <- NA
                          varHaz12$varHaz[match(jointt, allt)] <- varHazij
                          varHaz12$varHaz <- mstate:::NAfix(varHaz12$varHaz,
                            subst = 0)
                          varHaz[varHaz$trans1 == trans1 & varHaz$trans2 ==
                            trans2, ] <- varHaz12
                        }
                      }
                    }
                  }
                }
            }
        }
    }
    else {
        labels <- attr(Terms, "term.labels")
        if (length(labels) == 1) {
            if (labels == "strata(trans)") {
                sf0 <- summary(survfit(object))
                norisk <- sf0$n.risk
                noevent <- sf0$n.event
                sf0 <- data.frame(time = sf0$time, Haz = -log(sf0$surv),
                  norisk = norisk, noevent = noevent, var = sf0$std.err^2/(sf0$surv)^2,
                  trans = as.numeric(sf0$strata))
                allt <- sort(unique(c(sf0$time, lasty)))
                nt <- length(allt)
                K <- max(sf0$trans)
                Haz <- data.frame(time = rep(allt, K), Haz = NA,
                  trans = rep(1:K, rep(nt, K)))
                if (variance) {
                  tr12 <- data.frame(trans1 = rep(1:K, rep(K,
                    K)), trans2 = rep(1:K, K))
                  tr12 <- tr12[tr12$trans1 <= tr12$trans2, ]
                  varHaz <- data.frame(time = rep(allt, K * (K +
                    1)/2), varHaz = 0, trans1 = rep(tr12$trans1,
                    rep(nt, K * (K + 1)/2)), trans2 = rep(tr12$trans2,
                    rep(nt, K * (K + 1)/2)))
                }
                for (k in 1:K) {
                  sfk <- sf0[sf0$trans == k, ]
                  wht <- match(sfk$time, allt)
                  Hazk <- Haz[Haz$trans == k, ]
                  Hazk$Haz[wht] <- sfk$Haz
                 if (k==6) browser()
                  Hazk$Haz <- mstate:::NAfix(Hazk$Haz, subst = 0)
                  Haz[Haz$trans == k, ] <- Hazk
                  if (variance) {
                    varHazkk <- varHaz[varHaz$trans1 == k & varHaz$trans2 ==
                      k, ]
                    varHazkk$varHaz <- NA
                    varHazkk$varHaz[wht] <- sfk$var

                    varHazkk$varHaz <- mstate:::NAfix(varHazkk$varHaz,
                      subst = 0)
                    varHaz[varHaz$trans1 == k & varHaz$trans2 ==
                      k, ] <- varHazkk
                  }
                }
            }
        }
        else {
            method <- object$method
            if (method == "breslow")
                method <- 1
            else if (method == "efron")
                method <- 2
            else stop("Only \"efron\" and \"breslow\" methods for ties supported")
            type <- attr(y, "type")
            if (type == "counting") {
                if (has.strata)
                  ord <- order(mf[, strat], y[, 2], -y[, 3])
                else ord <- order(y[, 2], -y[, 3])
            }
            else if (type == "right") {
                if (has.strata)
                  ord <- order(mf[, strat], y[, 1], -y[, 2])
                else ord <- order(y[, 1], -y[, 2])
                miny <- min(y[, 1])
                if (miny < 0)
                  y <- cbind(2 * miny - 1, y)
                else y <- cbind(-1, y)
            }
            else stop("Cannot handle \"", type, "\" type survival data")
            if (variance)
                x <- x[ord, ]
            else x <- 0
            if (has.strata)
                newstrat <- (as.numeric(mf[, strat]))[ord]
            else newstrat <- rep(1, n)
            newstrat <- cumsum(table(newstrat))
            H <- length(newstrat)
            subterms <- function(tt, i) {
                dataClasses <- attr(tt, "dataClasses")
                predvars <- attr(tt, "predvars")
                oldnames <- dimnames(attr(tt, "factors"))[[1]]
                tt <- tt[i]
                index <- match(dimnames(attr(tt, "factors"))[[1]],
                  oldnames)
                if (length(index) > 0) {
                  if (!is.null(predvars))
                    attr(tt, "predvars") <- predvars[c(1, index +
                      1)]
                  if (!is.null(dataClasses))
                    attr(tt, "dataClasses") <- dataClasses[index]
                }
                tt
            }
            if (has.strata) {
                temp <- untangle.specials(Terms, "strata")
                if (length(temp$vars))
                  Terms <- subterms(Terms, -temp$terms)
            }
            Terms2 <- delete.response(Terms)
            if (has.strata) {
                if (length(attr(Terms2, "specials")$strata))
                  Terms2 <- subterms(Terms2, -attr(Terms2, "specials")$strata)
                if (!is.null(object$xlevels)) {
                  myxlev <- object$xlevels[match(attr(Terms2,
                    "term.labels"), names(object$xlevels), nomatch = 0)]
                  if (length(myxlev) == 0)
                    myxlev <- NULL
                }
                else myxlev <- NULL
                mf2 <- model.frame(Terms2, data = newdata, xlev = myxlev)
            }
            else mf2 <- model.frame(Terms2, data = newdata, xlev = object$xlevels)
            offset2 <- 0
            if (!missing(newdata)) {
                offset2 <- model.offset(mf2)
                if (length(offset2) > 0)
                  offset2 <- offset2 - mean(offset)
                else offset2 <- 0
                x2 <- model.matrix(Terms2, mf2)[, -1, drop = FALSE]
            }
            else stop("newdata missing")
            if (has.strata & is.null(newdata$strata))
                stop("no \"strata\" column present in newdata")
            n2 <- nrow(x2)
            coef <- ifelse(is.na(object$coefficients), 0, object$coefficients)
            newrisk <- exp(c(x2 %*% coef) + offset2 - sum(coef *
                object$means))
            dimnames(y) <- NULL
            storage.mode(y) <- "double"
            ndead <- sum(y[, 3])
            untimes <- sort(unique(y[, 2][y[, 3] == 1]))
            nt <- length(untimes)
            surv <- .C("agmssurv", sn = as.integer(n), sp = as.integer(nvar),
                svar = as.integer(variance), smethod = as.integer(method),
                sH = as.integer(H), sK = as.integer(n2), snt = as.integer(nt),
                y = y[ord, ], score = as.double(score[ord]),
                xmat = as.double(x), varcov = as.double(object$var),
                strata = as.integer(c(0, newstrat)), kstrata = as.integer(newdata$strata),
                unt = as.double(untimes), newx = as.double(x2),
                newrisk = as.double(newrisk), Haz = double(nt *
                  n2), varHaz = double(nt * n2 * (n2 + 1)/2),
                d = double(3 * nvar), work = double(nt * n2 *
                  (nvar + 1)))
            Haz <- data.frame(time = rep(untimes, n2), Haz = surv$Haz,
                trans = rep(1:n2, rep(nt, n2)))
            varHaz <- as.vector(t(matrix(surv$varHaz, ncol = nt)))
            hlp <- matrix(c(rep(1:n2, rep(n2, n2)), rep(1:n2,
                n2)), n2^2, 2)
            hlp <- hlp[hlp[, 1] <= hlp[, 2], ]
            varHaz <- data.frame(time = rep(untimes, n2 * (n2 +
                1)/2), varHaz = varHaz, trans1 = rep(hlp[, 1],
                rep(nt, n2 * (n2 + 1)/2)), trans2 = rep(hlp[,
                2], rep(nt, n2 * (n2 + 1)/2)))
            if (lasty > max(untimes)) {
                Hmat <- matrix(Haz$Haz, nrow = nt)
                Hmat <- rbind(Hmat, Hmat[nt, ])
                vHmat <- matrix(varHaz$varHaz, nrow = nt)
                vHmat <- rbind(vHmat, vHmat[nt, ])
                untimes <- c(untimes, lasty)
                nt <- nt + 1
                Haz <- data.frame(time = rep(untimes, n2), Haz = as.vector(Hmat),
                  trans = rep(1:n2, rep(nt, n2)))
                varHaz <- data.frame(time = rep(untimes, n2 *
                  (n2 + 1)/2), varHaz = as.vector(vHmat), trans1 = rep(hlp[,
                  1], rep(nt, n2 * (n2 + 1)/2)), trans2 = rep(hlp[,
                  2], rep(nt, n2 * (n2 + 1)/2)))
            }
        }
    }
    if (variance)
        res <- list(Haz = Haz, varHaz = varHaz, trans = trans)
    else res <- list(Haz = Haz, trans = trans)
    class(res) <- "msfit"
    return(res)
}



oscar.mssample <-
function (Haz, trans, history = list(state = 1, time = 0, tstate = NULL),
    beta.state = NULL, clock = c("forward", "reset"), output = c("state",
        "path", "data"), tvec, cens = NULL, M = 10, do.trace = NULL)
{
    output <- match.arg(output)
    clock <- match.arg(clock)
    K <- dim(trans)[1]
    trans2 <- mstate:::to.trans2(trans)
    ntrans <- nrow(trans2)
    if (length(history$state) == 1)
        history$state <- rep(history$state, M)
    if (length(history$time) == 1)
        history$time <- rep(history$time, M)
    if (length(history$state) != length(history$time))
        stop("lengths of history$state and history$time differ")
    if (!is.null(history$tstate)) {
        if (is.vector(history$tstate))
            if (length(history$tstate) != K)
                stop("length of history$tstate should equal no of states")
            else history$tstate <- matrix(history$tstate, K,
                M)
        if (is.null(beta.state))
            stop("beta.state should be specified when history$tstate not null")
    }
    if (!is.null(beta.state))
        if (any(dim(beta.state) != c(K, ntrans)))
            stop("incorrect dimension of beta.state")
    if (output == "state")
        res <- matrix(0, length(tvec), K)
    else if (output == "path") {
        thepaths <- paths(trans, start=history$state[1])
        L <- nrow(thepaths)
        res <- matrix(0, length(tvec), L)
    }
    else res <- NULL
    for (m in 1:M) {
        if (!is.null(history$tstate))
            res1 <- oscar.mssample1(Haz, trans, history = list(state = history$state[m],                time = history$time[m], tstate = history$tstate[,
                  m]), beta.state = beta.state, clock = clock,
                output = output, tvec = tvec, cens = cens)
        else res1 <- oscar.mssample1(Haz, trans, history = list(state = history$state[m],
            time = history$time[m], tstate = rep(0, K)), beta.state = beta.state,
            clock = clock, output = output, tvec = tvec, cens = cens)
        if (output == "data") {
            res1[, 1] <- m
            res <- rbind(res, res1)
        }
        else {
            res <- res + res1
        }
        if (!is.null(do.trace))
            if (m%%do.trace == 0) {
                cat("Replication", m, "finished at", date(),
                  "\n")
                flush.console()
            }
    }
    if (output == "state") {
        res <- data.frame(cbind(tvec, res/M))
        names(res) <- c("time", paste("pstate", 1:K, sep = ""))
    }
    else if (output == "path") {
        res <- data.frame(cbind(tvec, res/M))
        names(res) <- c("time", paste("ppath", 1:L, sep = ""))
    }
    else if (output == "data") {
        res <- data.frame(res)
        names(res) <- c("id", "Tstart", "Tstop", "duration",
            "from", "to", "status", "trans")
        attr(res, "trans") <- trans
        class(res) <- "msdata"
    }
    return(res)
}


 oscar.mssample1 <-
function (Haz, trans, history, beta.state, clock, output, tvec,
    cens)
{
    if (!is.null(cens)) {
        pcens <- diff(c(0, 1 - cens$surv))
        idx <- sample(1:length(cens$time), size = 1, prob = pcens)
        fut <- cens$time[idx]
        censtime <- list(time = fut, jump = ifelse(idx > 1, cens$Haz[idx] -
            cens$Haz[idx - 1], cens$Haz[idx]))
    }
    else censtime <- NULL
    K <- dim(trans)[1]
    trans2 <- mstate:::to.trans2(trans)
    from <- to <- history$state
    tcond <- t0 <- Tstart <- history$time
    if (output == "state")
        res <- matrix(0, length(tvec), K)
    else if (output == "path") {
        thepaths <- paths(trans, history$state[1])
        path <- c(to, rep(NA, ncol(thepaths) - 1))
        res <- matrix(0, length(tvec), nrow(thepaths))
    }
    else res <- NULL
    tstates <- history$tstate
    while (!is.na(to)) {
        from <- to
        nstates <- trans[from, ]
        transs <- nstates[!is.na(nstates)]
        allto <- which(!is.na(nstates))
        ntr <- length(transs)
        if (ntr != 0) {
            transnos <- transs
            for (tr in 1:ntr) Haz$Haz[Haz$trans == transnos[tr]] <- exp(sum(beta.state[,
                transnos[tr]] * tstates)) * Haz$Haz[Haz$trans ==
                transnos[tr]]
            whh <- which(!is.na(match(Haz$trans, transnos)))
            if (clock == "forward") {
                crs <- crsample(Haz[whh, ], tcond, censtime)
                tcond <- Tstop <- crs$t
            }
            else {
                crs <- oscar.crsample(Haz[whh, ], t0, censtime)
                t0 <- 0
                tcond <- Tstop <- crs$t + tcond
            }
            transno <- crs$trans
            if (is.na(transno))
                to <- NA
            else {
                to <- trans2$to[transno]
                tstates[to] <- Tstop
            }
            if (output == "state") {
                res[((tvec >= Tstart) & (tvec < Tstop)), from] <- 1
                Tstart <- Tstop
            }
            else if (output == "path") {
                idx <- which(apply(thepaths, 1, function(x) identical(x,
                  path)))
                res[((tvec >= Tstart) & (tvec < Tstop)), idx] <- 1
                path[which(is.na(path))[1]] <- to
                Tstart <- Tstop
            }
            else {
                res1 <- matrix(c(rep(NA, ntr), rep(Tstart, ntr),
                  rep(Tstop, ntr), rep(Tstop - Tstart, ntr),
                  rep(from, ntr), allto, rep(0, 2 * ntr)), ntr,
                  8)
                res1[res1[, 6] == to, 7] <- 1
                res1[, 8] <- trans[from, allto]
                Tstart <- Tstop
                res <- rbind(res, res1)
            }
        }
        else {
            to <- NA
            if (output == "state") {
                res[tvec >= Tstart, from] <- 1
            }
            else if (output == "path") {
                idx <- which(apply(thepaths, 1, function(x) identical(x,
                  path)))
                res[tvec >= Tstart, idx] <- 1
                path[which(is.na(path))[1]] <- to
            }
            else {
                res1 <- matrix(c(rep(NA, ntr), rep(Tstart, ntr),
                  rep(Tstop, ntr), rep(Tstop - Tstart, ntr),
                  rep(from, ntr), allto, rep(0, 2 * ntr)), ntr,
                  8)
                res1[res1[, 6] == to, 7] <- 1
                res1[, 8] <- trans[from, allto]
                res <- rbind(res, res1)
            }
        }
    }
    return(res)
}

oscar.crsample <- function (Haz, tcond = 0, censtime = NULL)
{
    if (is.null(censtime))
        fut <- Inf
    else fut <- censtime$time
    transs <- Haz$trans
    transun <- unique(transs)
    K <- length(transun)
    tt <- sort(unique(Haz$time))
    n <- length(tt)
    cim <- matrix(NA, n, 3 * K + 4)
    ci <- as.data.frame(cim)
    names(ci)[1] <- "time"
    names(ci)[2:(K + 1)] <- paste("Haz", as.character(1:K), sep = "")
    names(ci)[(K + 2):(2 * K + 1)] <- paste("haz", as.character(1:K),
        sep = "")
    names(ci)[(2 * K + 2):(3 * K + 1)] <- paste("CI", as.character(1:K),
        sep = "")
    names(ci)[3 * K + 2] <- "hazsum"
    names(ci)[3 * K + 3] <- "Hazsum"
    names(ci)[3 * K + 4] <- "S0"
    ci$time <- tt
    for (k in 1:K) {
        wh <- which(Haz$trans == transun[k])
        idx <- match(Haz$time[wh], tt)
        ci[, k + 1][idx] <- Haz$Haz[wh]
        ci[, k + 1] <- mstate:::NAfix(ci[, k + 1], subst = 0)
        ci[, K + 1 + k] <- diff(c(0, ci[, k + 1]))
    }
    ci <- ci[ci$time > tcond, ]
    n <- nrow(ci)
    for (k in 1:K) ci[, k + 1] <- cumsum(ci[, K + 1 + k])
    if (K == 1)
        ci$hazsum <- ci[, 3]
    else ci$hazsum <- apply(ci[, ((K + 2):(2 * K + 1))], 1, sum)
    ci$S0 <- cumprod(1 - ci$hazsum)
    ci$Hazsum <- -log(ci$S0)
    nci <- nrow(ci)
    k <- NA
    tmp <- data.frame(time = ci$time, Haz = ci$Hazsum)
    tmp <- tmp[which(!is.nan(tmp$Haz)),]
    tsample <- oscar.Hazsample(tmp)
    if (fut < tsample)
        crt <- fut
    else {
        crt <- tsample
        if (fut > tsample) {
            k <- sample(1:K, size = 1, prob = ci[which(ci$time ==
                tsample), (K + 2):(2 * K + 1)])
        }
        else if (crt != Inf) {
            k <- sample(c(1:K, NA), size = 1, prob = c(ci[which(ci$time ==
                tsample), (K + 2):(2 * K + 1)], censtime$jump))
        }
    }
    if (!is.na(k))
        trans <- unique(Haz$trans)[k]
    else trans <- NA
    return(list(t = crt, trans = trans))
}

oscar.Hazsample <-
function (Haz, size = 1, replace = TRUE)
{
    p <- diff(c(0, 1 - exp(-Haz$Haz)))
    p <- c(p, exp(-Haz$Haz[nrow(Haz)]))
    ## I've added this
    p <- pmax(0, p)
    ## I've added this
    return(sample(c(Haz$time, Inf), size = size, prob = p, replace = replace))
}

oscar.probtrans <-
function (object, predt, direction = c("forward", "fixedhorizon"),
    method = c("aalen", "greenwood"), variance = TRUE, covariance = FALSE)
{
    if (!inherits(object, "msfit"))
        stop("'object' must be a 'msfit' object")
    method <- match.arg(method)
    direction <- match.arg(direction)
    trans <- object$trans
    transit <- mstate:::to.trans2(trans)
    numtrans <- nrow(transit)
    stackhaz <- object$Haz
    stackvarhaz <- object$varHaz
    for (i in 1:numtrans) stackhaz$dhaz[stackhaz$trans == i] <- diff(c(0,
                                            stackhaz$Haz[stackhaz$trans == i]))
    return(stackhaz$dhaz)
    if (direction == "forward")
        stackhaz <- stackhaz[stackhaz$time > predt, ]
    else stackhaz <- stackhaz[stackhaz$time <= predt, ]
    untimes <- sort(unique(stackhaz$time))
    TT <- length(untimes)
    S <- nrow(trans)
    if (covariance)
        variance <- TRUE
    if (direction == "forward") {
        if (variance == TRUE)
            res <- array(0, c(TT + 1, 2 * S + 1, S))
        else res <- array(0, c(TT + 1, S + 1, S))
        res[1, 1, ] <- predt
        for (j in 1:S) res[1, 1 + j, ] <- rep(c(0, 1, 0), c(j -
            1, 1, S - j))
        if (variance)
            res[1, (S + 2):(2 * S + 1), ] <- 0
    }
    else {
        if (predt %in% untimes) {
            if (variance)
                res <- array(0, c(TT + 1, 2 * S + 1, S))
            else res <- array(0, c(TT + 1, S + 1, S))
            res[TT + 1, 1, ] <- predt
            for (j in 1:S) res[TT + 1, 1 + j, ] <- rep(c(0, 1,
                0), c(j - 1, 1, S - j))
            if (variance)
                res[TT + 1, (S + 2):(2 * S + 1), ] <- 0
        }
        else {
            if (variance)
                res <- array(0, c(TT + 2, 2 * S + 1, S))
            else res <- array(0, c(TT + 2, S + 1, S))
            res[TT + 1, 1, ] <- max(untimes)
            for (j in 1:S) res[TT + 1, 1 + j, ] <- rep(c(0, 1,
                0), c(j - 1, 1, S - j))
            if (variance)
                res[TT + 1, (S + 2):(2 * S + 1), ] <- 0
            res[TT + 2, 1, ] <- predt
            for (j in 1:S) res[TT + 2, 1 + j, ] <- rep(c(0, 1,
                0), c(j - 1, 1, S - j))
            if (variance)
                res[TT + 2, (S + 2):(2 * S + 1), ] <- 0
        }
    }
    P <- diag(S)
    if (covariance) {
        varParr <- array(0, c(S^2, S^2, TT + 1))
        if ((direction == "fixedhorizon") & !(predt %in% untimes))
            varParr <- array(0, c(S^2, S^2, TT + 2))
        ffrom <- rep(1:S, S)
        tto <- rep(1:S, rep(S, S))
        fromto <- paste("from", ffrom, "to", tto, sep = "")
        if (direction == "forward")
            dimnames(varParr) <- list(fromto, fromto, c(predt,
                untimes))
        else {
            if (predt %in% untimes)
                dimnames(varParr) <- list(fromto, fromto, c(0,
                  untimes))
            else dimnames(varParr) <- list(fromto, fromto, c(0,
                untimes, predt))
        }
    }
    if (variance) {
        varP <- matrix(0, S^2, S^2)
        if (direction == "forward") {
            varAnew <- array(0, c(S, S, S, S))
            if (predt != 0) {
                tmin <- max(stackvarhaz$time[stackvarhaz$time <=
                  predt])
                varHaz <- stackvarhaz[stackvarhaz$time == tmin,
                  ]
                lHaz <- nrow(varHaz)
                for (j in 1:lHaz) {
                  from1 <- transit$from[transit$transno == varHaz$trans1[j]]
                  to1 <- transit$to[transit$transno == varHaz$trans1[j]]
                  from2 <- transit$from[transit$transno == varHaz$trans2[j]]
                  to2 <- transit$to[transit$transno == varHaz$trans2[j]]
                  varAnew[from1, to1, from2, to2] <- varAnew[from2,
                    to2, from1, to1] <- varHaz$varHaz[j]
                }
            }
        }
        else {
            varA <- array(0, c(S, S, S, S))
            varHaz <- stackvarhaz[stackvarhaz$time == untimes[TT],
                ]
            lHaz <- nrow(varHaz)
            for (j in 1:lHaz) {
                from1 <- transit$from[transit$transno == varHaz$trans1[j]]
                to1 <- transit$to[transit$transno == varHaz$trans1[j]]
                from2 <- transit$from[transit$transno == varHaz$trans2[j]]
                to2 <- transit$to[transit$transno == varHaz$trans2[j]]
                varA[from1, to1, from2, to2] <- varA[from2, to2,
                  from1, to1] <- varHaz$varHaz[j]
            }
        }
    }
    for (i in 1:TT) {
        idx <- ifelse(direction == "forward", i, TT + 1 - i)
        tt <- untimes[idx]
        Haztt <- stackhaz[stackhaz$time == tt, ]
        lHaztt <- nrow(Haztt)
        IplusdA <- diag(S)
        for (j in 1:lHaztt) {
            from <- transit$from[transit$transno == Haztt$trans[j]]
            to <- transit$to[transit$transno == Haztt$trans[j]]
            IplusdA[from, to] <- Haztt$dhaz[j]
            IplusdA[from, from] <- IplusdA[from, from] - Haztt$dhaz[j]
        }
        if (any(diag(IplusdA) < 0))
            warning("Warning! Negative diagonal elements of (I+dA); the estimate may not be meaningful. \n")
        if (variance) {
            if (direction == "forward") {
                varA <- varAnew
                varAnew <- array(0, c(S, S, S, S))
                varHaztt <- stackvarhaz[stackvarhaz$time == tt,
                  ]
                lHaztt <- nrow(varHaztt)
                for (j in 1:lHaztt) {
                  from1 <- transit$from[transit$transno == varHaztt$trans1[j]]
                  to1 <- transit$to[transit$transno == varHaztt$trans1[j]]
                  from2 <- transit$from[transit$transno == varHaztt$trans2[j]]
                  to2 <- transit$to[transit$transno == varHaztt$trans2[j]]
                  varAnew[from1, to1, from2, to2] <- varAnew[from2,
                    to2, from1, to1] <- varHaztt$varHaz[j]
                }
                vardA <- varAnew - varA
            }
            else {
                varAttmin <- array(0, c(S, S, S, S))
                varHazttmin <- stackvarhaz[stackvarhaz$time ==
                  untimes[idx - 1], ]
                lHazttmin <- nrow(varHazttmin)
                for (j in 1:lHazttmin) {
                  from1 <- transit$from[transit$transno == varHazttmin$trans1[j]]
                  to1 <- transit$to[transit$transno == varHazttmin$trans1[j]]
                  from2 <- transit$from[transit$transno == varHazttmin$trans2[j]]
                  to2 <- transit$to[transit$transno == varHazttmin$trans2[j]]
                  varAttmin[from1, to1, from2, to2] <- varAttmin[from2,
                    to2, from1, to1] <- varHazttmin$varHaz[j]
                }
                vardA <- varA - varAttmin
                varA <- varAttmin
            }
            for (from in 1:S) {
                for (from2 in 1:S) {
                  for (to2 in 1:S) {
                    if (to2 != from2)
                      vardA[from, from, from2, to2] <- vardA[from2,
                        to2, from, from] <- -sum(vardA[from,
                        -from, from2, to2])
                  }
                }
            }
            for (from in 1:S) {
                for (from2 in 1:S) vardA[from, from, from2, from2] <- vardA[from2,
                  from2, from, from] <- -sum(vardA[from, from,
                  from2, -from2])
            }
            vardA <- matrix(vardA, S^2, S^2)
        }
        if (method == "aalen") {
            if (direction == "forward") {
                P <- P %*% IplusdA
                if (variance) {
                  tmp1 <- kronecker(t(IplusdA), diag(S)) %*%
                    varP %*% kronecker(IplusdA, diag(S))
                  tmp2 <- kronecker(diag(S), P) %*% vardA %*%
                    kronecker(diag(S), t(P))
                  varP <- tmp1 + tmp2
                }
            }
            else {
                if (variance) {
                  tmp1 <- kronecker(diag(S), IplusdA) %*% varP %*%
                    kronecker(diag(S), t(IplusdA))
                  tmp2 <- kronecker(t(P), IplusdA) %*% vardA %*%
                    kronecker(P, t(IplusdA))
                  varP <- tmp1 + tmp2
                }
                P <- IplusdA %*% P
            }
        }
        if (method == "greenwood") {
            if (direction == "forward") {
                if (variance) {
                  tmp1 <- kronecker(t(IplusdA), diag(S)) %*%
                    varP %*% kronecker(IplusdA, diag(S))
                  tmp2 <- kronecker(diag(S), P) %*% vardA %*%
                    kronecker(diag(S), t(P))
                  varP <- tmp1 + tmp2
                }
                P <- P %*% IplusdA
            }
            else {
                if (variance) {
                  tmp1 <- kronecker(diag(S), IplusdA) %*% varP %*%
                    kronecker(diag(S), t(IplusdA))
                  tmp2 <- kronecker(t(P), diag(S)) %*% vardA %*%
                    kronecker(P, diag(S))
                  varP <- tmp1 + tmp2
                }
                P <- IplusdA %*% P
            }
        }
        if (variance) {
            seP <- sqrt(diag(varP))
            seP <- matrix(seP, S, S)
        }
        if (covariance) {
            if (direction == "forward")
                varParr[, , i + 1] <- varP
            else {
                varParr[, , idx + 1] <- varP
            }
        }
        if (direction == "forward") {
            res[idx + 1, 1, ] <- tt
            res[idx + 1, 2:(S + 1), ] <- t(P)
            if (variance)
                res[idx + 1, (S + 2):(2 * S + 1), ] <- t(seP)
        }
        else {
            res[idx, 1, ] <- ifelse(i == TT, 0, untimes[TT -
                i])
            res[idx, 2:(S + 1), ] <- t(P)
            if (variance)
                res[idx, (S + 2):(2 * S + 1), ] <- t(seP)
        }
    }
    if (covariance & (direction == "fixedhorizon"))
        varParr[, , 1] <- varParr[, , 2]
    res2 <- vector("list", S)
    for (s in 1:S) {
        tmp <- as.data.frame(res[, , s])
        if (min(dim(tmp)) == 1)
            tmp <- res[, , s]
        if (variance)
            names(tmp) <- c("time", paste("pstate", 1:S, sep = ""),
                paste("se", 1:S, sep = ""))
        else names(tmp) <- c("time", paste("pstate", 1:S, sep = ""))
        res2[[s]] <- tmp
    }
    if (covariance)
        res2$varMatrix <- varParr
    res2$trans <- trans
    res2$method <- method
    res2$predt <- predt
    res2$direction <- direction
    class(res2) <- "probtrans"
    return(res2)
}
