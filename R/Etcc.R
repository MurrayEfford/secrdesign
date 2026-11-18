## Etcc.R  DRAFT 2026-11-02 (needs .cuemoments from Echat.R; split of the 2026-11-01 version; earlier version in Etcc_before_split_2026-11-02.R)
## Expected counts and predicted precision for telemetry-calibrated cue counts (secr telemetrytype 'marking').
## Analogue of secrdesign::Enrm (expected counts) and secrdesign::minnrRSE (RSE from the counts):
##   Etcc()    expected counts from lambda0, sigma, D, the camera array, the number of occasions and the number of collars:
##             E(n) animals detected, J1 all cues, U unmarked cues, R identified collar cues, Tm unidentified collar cues,
##             collars seen, c-hat of the unmarked total
##   RSEtcc()  same inputs as Etcc() (as minnrRSE has the inputs of Enrm); returns the counts of Etcc() and the predicted RSE(D-hat)
##
## Model (hazard half-normal 'HHN' or hazard exponential 'HEX', count detectors, AC Poisson with density D(x), one cue class):
##   L(x)  = S sum_k lambda0 exp(-d_k(x)^2 / (2 sigma^2))       expected cues from an animal centred at x (S occasions)
##           (HEX: exp(-d_k(x) / sigma); Vs below then uses 2 sigma^2 in place of sigma^2, the per-axis variance of the AC given the fixes)
##   En = a sum_x D (1 - exp(-L)),  J1 = a sum_x D L,  J2 = a sum_x D L^2,  Nexp = a sum_x D
##   neff = J1^2 / J2   effective number of animals (Kish effective sample size, animals weighted by their expected cues);
##        the relative variance of the pooled cue total is about 1 / J1 + 1 / neff
##   collars: AC distribution w(x) proportional to D(x)              ('uniform', i.e. widespread)
##                                         or to D(x)(1 - exp(-L1(x)))  ('pdot', centred; L1 = L / S, one occasion)
##   Cm = ncollar * E_w[L]            expected cues of the collared animals;  R = q Cm identified, Tm = (1 - q) Cm
##   U = mu = J1 - Cm                 expected unmarked cues
##   Vs = ncollar * E_w[|grad L|^2] sigma^2 / nfix          location uncertainty of the collared animals' cues
##   c-hat = 1 + (J2 + Vs) / mu       (fixed N: J2 - J1^2 / Nexp in place of J2), as in secr_chatmarking
## Approximate precision (RSEtcc):
##   RSE^2 ~ c-hat * mu / J1^2  +  CV^2(Lambda-hat)            RSEunmarked^2 + RSEcalib^2
##   RSEcoarse = sqrt(1 / En + 1 / Cm)                          coarse form, in the spirit of Efford and Boulanger (2019)
##   calib = 'mc' (default): CV^2(Lambda-hat) = E[1 / realised collar cues] over random collar locations (nsimcollar draws);
##        the information depends on where the collars happen to be, and E(1/R) exceeds 1/E(R) when few collars are near the
##        cameras (pnocues = proportion of draws with no cues at all, excluded; the fit would fail there)
##   calib = 'fisher': g' I^-1 g, I the Fisher information for (log lambda0, log sigma) from the collared animals (identified
##        plus unidentified counts, activity centres taken as known, distributed as w(x)) plus the telemetry fixes
##        (information 2 (nfix - 1) per collar for log sigma), g the gradient of log Lambda-bar. No fitted constants.
##        The cross-information between q and (lambda0, sigma) is zero, so q does not change this term.
##   calib = 'poisson': kcal / R (Poisson variance of the identified collar cues; used for the first draft)
## Checked against simsN.RDS and rse_check.RDS only for q = 1, uniform AC, count detectors.
##
## Arguments (D, lambda0, sigma, noccasions, ncollar, q, nfix are recycled to a common length, one row per setting):
##   traps  count-detector traps object;  mask  habitat mask (D is per hectare, a = attr(mask, 'area'))
##   D      scalar, or a vector with one value per mask cell (single setting only)
## warn is TRUE when fewer than `warnseen` collars are expected to be seen (prediction optimistic).
## Value: data.frame with one row per setting.

## common set-up: recycle the settings and compute the distances
.etcc_setup <- function (D, lambda0, sigma, traps, mask, noccasions, ncollar, deployment, q, nfix, detectfn = 'HHN') {
    if (length(D) > 1 && length(D) != nrow(mask)) stop ("vector D must have one value per mask cell")
    perCell <- length(D) > 1
    if (ms(traps)) stop ("Etcc does not accept multi-session traps")
    us <- usage(traps)
    ## noccasions defaults to ncol(usage) when traps has usage (as pdot), otherwise 5
    if (is.null(noccasions)) noccasions <- if (is.null(us)) 5 else ncol(us)
    n <- max(if (perCell) 1 else length(D), length(lambda0), length(sigma), length(noccasions), length(ncollar), length(q), length(nfix))
    rc <- function(x) rep_len(x, n)
    ## mean effort per occasion of each detector, fk = rowMeans(usage) (1 if traps has no usage attribute); usage is any
    ## non-negative effort that multiplies the hazard (secr::usage), not necessarily 0/1 or a proportion;
    ## expected cues from a detector are noccasions * fk * hazard (= sum of usage over occasions when noccasions = ncol(usage))
    fk <- if (is.null(us)) rep(1, nrow(traps)) else rowMeans(us)
    dx <- outer(mask$x, traps$x, '-'); dy <- outer(mask$y, traps$y, '-')
    list(n = n, perCell = perCell, detectfn = detectfn, D = if (perCell) D else rc(D), lambda0 = rc(lambda0), sigma = rc(sigma),
         noccasions = rc(noccasions), ncollar = rc(ncollar), q = rc(q), nfix = rc(nfix), deployment = deployment,
         a = attr(mask, 'area'), nmask = nrow(mask), dx = dx, dy = dy, d2 = dx^2 + dy^2,
         fk = fk)
}

## expected counts and the pieces needed for the RSE, for setting i
.etcc_one <- function (z, i, Ndist) {
    s2 <- z$sigma[i]^2; S <- z$noccasions[i]; m <- z$ncollar[i]; lam <- z$lambda0[i]
    if (z$detectfn == 'HEX') {                              # h = lambda0 exp(-d / sigma)
        d  <- pmax(sqrt(z$d2), 1e-8)
        g  <- exp(-d / z$sigma[i])
        gk <- g / (z$sigma[i] * d)                          # -(dh/dx) / (lambda0 dx)
        vloc <- 2 * s2                                      # per-axis variance of the AC given the fixes: 1 / Fisher information
    }
    else {                                                  # h = lambda0 exp(-d^2 / (2 sigma^2))
        g  <- exp(-z$d2 / (2 * s2))
        gk <- g / s2
        vloc <- s2
    }
    fk <- z$fk
    L1 <- lam * as.vector(g %*% fk)                         # per occasion (mean over occasions if usage varies)
    L  <- S * L1
    Lx <- -S * lam * as.vector((gk * z$dx) %*% fk)          # gradient of L
    Ly <- -S * lam * as.vector((gk * z$dy) %*% fk)
    Dx <- if (z$perCell) z$D else rep(z$D[i], z$nmask)
    mom <- .cuemoments(L, Dx, z$a, Ndist)               # J1, J2, N, E(n), as in Echat()
    J1 <- mom$J1; Nexp <- mom$N; En <- mom$En               # En: expected number of animals detected
    w  <- if (z$deployment == 'pdot') Dx * (1 - exp(-L1)) else Dx
    w  <- w / sum(w)
    Cm <- m * sum(w * L)
    Vs <- m * sum(w * (Lx^2 + Ly^2)) * vloc / z$nfix[i]
    mu <- J1 - Cm
    vex <- mom$vex                                          # J2, or J2 - J1^2 / N if N is fixed
    list(N = Nexp, En = En, neff = J1^2 / mom$J2, J1 = J1, Lambda = z$a * sum(L), U = mu, C = Cm, R = z$q[i] * Cm, Tm = (1 - z$q[i]) * Cm,
         seen = m * sum(w * (1 - exp(-L))), chat = if (mu > 0) 1 + (vex + Vs) / mu else NA,
         w = w, L = L, g = g, Dx = Dx, dlh = if (z$detectfn == 'HEX') d / z$sigma[i] else z$d2 / s2)   # dlh = d log h / d log sigma
}

.etcc_row <- function (z, i, o, warnseen)
    data.frame(D = if (z$perCell) NA else z$D[i], lambda0 = z$lambda0[i], sigma = z$sigma[i], noccasions = z$noccasions[i],
               ncollar = z$ncollar[i], deployment = z$deployment, q = z$q[i],
               N = o$N, En = o$En, neff = o$neff, J1 = o$J1, Lambda = o$Lambda, U = o$U, R = o$R, Tm = o$Tm, seen = o$seen,
               collarshare = o$C / o$J1, chat = o$chat, warn = o$seen < warnseen)

Etcc <- function (D, lambda0, sigma, traps, mask, noccasions = NULL, ncollar = 20,
                  deployment = c('pdot', 'uniform'), q = 1, nfix = 100,
                  Ndist = c('poisson', 'fixed'), warnseen = 3, detectfn = c('HHN', 'HEX')) {
    deployment <- match.arg(deployment)
    Ndist  <- match.arg(Ndist)
    detectfn <- match.arg(detectfn)
    z <- .etcc_setup(D, lambda0, sigma, traps, mask, noccasions, ncollar, deployment, q, nfix, detectfn)
    do.call(rbind, lapply(seq_len(z$n), function(i) .etcc_row(z, i, .etcc_one(z, i, Ndist), warnseen)))
}

RSEtcc <- function (D, lambda0, sigma, traps, mask, noccasions = NULL, ncollar = 20,
                    deployment = c('pdot', 'uniform'), q = 1, nfix = 100,
                    Ndist = c('poisson', 'fixed'), calib = c('mc', 'fisher', 'poisson'), kcal = 1,
                    warnseen = 3, nsimcollar = 2000, detectfn = c('HHN', 'HEX')) {
    calib <- match.arg(calib)
    deployment <- match.arg(deployment)
    Ndist  <- match.arg(Ndist)
    detectfn <- match.arg(detectfn)
    z <- .etcc_setup(D, lambda0, sigma, traps, mask, noccasions, ncollar, deployment, q, nfix, detectfn)
    out <- vector('list', z$n)
    for (i in seq_len(z$n)) {
        o <- .etcc_one(z, i, Ndist)
        m <- z$ncollar[i]; S <- z$noccasions[i]; lam <- z$lambda0[i]; s2 <- z$sigma[i]^2
        unm2 <- if (o$U > 0) o$chat * o$U / o$J1^2 else NA     # = (mu + J2 + Vs) / J1^2
        pzero <- NA
        if (calib == 'poisson') cal2 <- kcal / o$R
        else if (calib == 'mc') {
            ## the information from the collars depends on where the collars happened to be: draw collar activity centres from
            ## w(x), the realised total collar cues is Poisson(sum L(s_i)), and the calibration variance is 1 / (realised cues)
            cells <- matrix(sample.int(length(o$w), m * nsimcollar, replace = TRUE, prob = o$w), m)
            Cr <- rpois(nsimcollar, colSums(matrix(o$L[cells], m)))
            pzero <- mean(Cr == 0)
            cal2 <- mean(1 / Cr[Cr > 0])
        }
        else {
            g <- o$g; w <- o$w; dlh <- o$dlh
            fk <- z$fk
            bx <- S * lam * as.vector((g * dlh) %*% fk)          # dL/dlog(sigma)
            M  <- m * matrix(c(sum(w * o$L), sum(w * bx), sum(w * bx), sum(w * as.vector(S * lam * (g * dlh^2) %*% fk))), 2)
            M[2, 2] <- M[2, 2] + 2 * m * (z$nfix[i] - 1)
            gr <- c(1, sum(o$Dx * bx) / sum(o$Dx * o$L))
            cal2 <- as.numeric(t(gr) %*% solve(M, gr))
        }
        out[[i]] <- cbind(.etcc_row(z, i, o, warnseen),
                          RSEcoarse = sqrt(1 / o$En + 1 / o$C), RSEunmarked = sqrt(unm2), RSEcalib = sqrt(cal2),
                          RSE = sqrt(unm2 + cal2), pnocues = pzero)
    }
    do.call(rbind, out)
}
