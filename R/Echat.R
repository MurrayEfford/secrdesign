## Echat.R  DRAFT 2026-11-02
## Expected variance inflation (c-hat) of the pooled total of sightings of unmarked animals, from the parameters of the
## model, the detector array and the number of occasions alone: the design-stage counterpart of the analytic c-hat that
## secr computes at the fitted parameters (secr_chatsighting in secr/R/chatmarking.R; details = list(chatmethod = 'analytic')).
## Family: Enrm (expected n, r, m), Echat (expected c-hat), Etcc / RSEtcc (telemetry-calibrated cue counts).
##
## Sighting-only models: all n_u animals unmarked, AC independent with distribution u(x) proportional to D(x), n_u Poisson
## (variance = mean) or fixed at N. Cells c = (k, s) are detector-occasions; all detectors are used on all occasions.
## count detectors:   T = total cues, E(T) = J1,  Var(T) = J1 + J2 (Poisson N) or J1 + J2 - J1^2 / N (fixed N);
##                    J1 = a sum D L, J2 = a sum D L^2, L(x) = noccasions sum_k h_k(x);  c-hat = Var(T) / E(T)
## binary (proximity) detectors:  T = number of positive cells, a cell is positive with probability 1 - exp(-X_c), X_c the
##                    summed hazard of the animals present;  with G the probability generating function of n_u,
##                    M1_k = G(sum_x u(x) exp(-h_k(x))),  M2_kk' = G(sum_x u(x) exp(-h_k(x) - h_k'(x)));
##                    E(T) = S sum_k (1 - M1_k),
##                    Var(T) = S sum_k (M1_k - M2_kk) + S^2 sum_kk' M2_kk' - S^2 (sum_k M1_k)^2   (cells in different
##                    occasions share the same hazards, so the sums over cells reduce to K x K sums);
##                    c-hat = Var(T) / (C pbar (1 - pbar)), C = K S cells, pbar = E(T) / C  (naive binomial variance)
## Shared with Etcc(): .cuemoments() (J1, J2, N, E(n) and the count-detector variance term).

## moments of the cue total for a cue-rate surface L(x) (cues of an animal centred at x over the survey)
.cuemoments <- function (L, Dx, a, Ndist = c('poisson', 'fixed')) {
    Ndist <- match.arg(Ndist)
    J1 <- a * sum(Dx * L); J2 <- a * sum(Dx * L^2); N <- a * sum(Dx)
    list(N = N, En = a * sum(Dx * (1 - exp(-L))), J1 = J1, J2 = J2,
         vex = if (Ndist == 'fixed') J2 - J1^2 / N else J2)     # Var(T) = J1 + vex
}

Echat <- function (D, traps, mask, detectpar, noccasions,
                   detectfn = c('HHN', 'HHR', 'HEX', 'HN', 'HR', 'EX'),
                   Ndist = c('poisson', 'fixed')) {
    if (ms(traps) | ms(mask))
        stop ("Echat does not accept multi-session traps or mask")
    Ndist <- match.arg(Ndist)
    detect <- detector(traps)[1]
    if (!detect %in% c('count', 'proximity'))
        stop ("Echat is only for 'count' and 'proximity' detectors")
    if (is.character(detectfn))
        detectfn <- match.arg(detectfn)
    detectfn <- secr:::secr_valid.detectfn(detectfn, valid = c(0,1,2,14,15,16))
    dfc <- dfcast (detectfn, detectpar)            # transforms detectfn 0 to 14, 1 to 15, 2 to 16
    detectfn <- dfc$detectfn
    detectpar <- dfc$detectpar
    if (length(D) == 1) D <- rep(D, nrow(mask))
    if (length(D) != nrow(mask)) stop ("D should be a scalar or have one value per mask cell")
    d <- edist(mask, traps)                        # mask cells x detectors
    h <- switch(as.character(detectfn),            # hazard per occasion
        '14' = detectpar$lambda0 * exp(-d^2 / (2 * detectpar$sigma^2)),
        '15' = detectpar$lambda0 * (1 - exp(-(d / detectpar$sigma)^-detectpar$z)),
        '16' = detectpar$lambda0 * exp(-d / detectpar$sigma))
    a <- attr(mask, 'area')
    S <- noccasions
    if (detect == 'count') {
        mom <- .cuemoments(S * rowSums(h), D, a, Ndist)
        ET <- mom$J1
        VT <- ET + mom$vex
        chat <- if (ET > 0) VT / ET else 1
    }
    else {
        nbar <- a * sum(D)
        u <- D / sum(D)
        G <- if (Ndist == 'poisson') function(s) exp(nbar * (s - 1)) else function(s) s^nbar
        W  <- exp(-h)
        M1 <- G(colSums(W * u))                    # K
        M2 <- G(crossprod(W * sqrt(u)))            # K x K
        ET <- S * sum(1 - M1)
        VT <- S * sum(M1 - diag(M2)) + S^2 * sum(M2) - S^2 * sum(M1)^2
        npc <- ncol(h) * S                         # number of cells
        pbar <- ET / npc
        chat <- if (pbar > 0 && pbar < 1) VT / (npc * pbar * (1 - pbar)) else 1
    }
    c(ET = ET, VT = VT, chat = chat)
}
