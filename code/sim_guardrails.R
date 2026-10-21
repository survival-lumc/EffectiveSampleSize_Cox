####################### ADEMP SIMULATION: WHEN IS N* WELL-BEHAVED #####################

# At which event fraction and from which predicted survival S is the effective sample size
# N* well-behaved (N* <= N)? Designs: "factorial" (C x event rate x administrative censoring),
# "censoring". Start with a pilot: R_REP <- 50 and the first 2 scenarios.
#
# A  Aims: (1) smallest event fraction at which N*-hat > N is rare, per model strength;
#          (2) lowest S above which N*-hat is reliable; (3) is N*-hat > N real, i.e. does the
#          empirical-variance N* also exceed N?
# D  Data: Weibull PH (H0(t) = t^shape); k continuous (N(0,1)) + k binary (p = 0.5) true
#          predictors with equal contribution to the linear predictor, scaled to a target C;
#          optionally k extra pure-noise predictors (alternating continuous/binary) that are
#          fitted but unrelated to the outcome; administrative censoring so that a proportion
#          e of patients has an event.
# E  Estimands: S(t|x) at 100 fixed profiles x 7 horizons (fractions of all events);
#          N*-true = S(1-S)/Var_empirical(S-hat).
# M  Methods: coxph (Efron) with all predictors (true + noise); survfit delta-method SE.
# P  Performance: P(N*-hat > N), P(> 1.25 N), q95 of N*-hat/N, N*-true/N, median N*-hat/N*-true,
#          variance ratio, 95% log-log CI coverage; per scenario x event fraction x S band.

library(survival)
library(future.apply)
library(ggplot2)

####################### ARTICLE PALETTE (SAME HEX CODES AS ARTICLE_FIGURES.R) #####################
PAL <- c( navy = "#37293F", teal = "#2E7691", purple = "#9f84af", red = "#C2666B", gold = "#c6aa2c" )
# sequential ramp for ordered levels (event rate, C, ...): gold -> red -> teal -> navy
ramp_pal <- function( n ) if( n <= 1 ) PAL[["navy"]] else
  colorRampPalette( unname(PAL[c("gold", "red", "teal", "navy")]), space = "Lab" )(n)
theme_art <- function( base_size = 11 )
  theme_minimal(base_size = base_size) +
  theme(strip.text = element_text(hjust = 0), legend.position = "bottom", panel.grid.minor = element_blank())

####################### SETTINGS #####################
FRACS   <- c( 0.05, 0.10, 0.25, 0.40, 0.60, 0.75, 0.90 )   # fraction of all events observed
S_BANDS <- c( 0, 0.01, 0.05, 0.20, 0.50, 0.80, 0.90, 0.95, 0.98, 1 )   # extra resolution at high S
R_REP   <- 1000          # replications per scenario: SET FROM pilot_mc() (see bottom)
N_PROF  <- 200
D_EVENTS <- 500          # event-rate design: N = D_EVENTS / e, i.e. ~500 events in every scenario

# Registry-style censoring (design "censoring"; columns `admin` and `rc` of the scenario table)
#   The end of follow-up is tau = 1 in every scenario. The baseline hazard H0(t) = lambda * t^shape is scaled
#   so that e = 1 - S(tau), with S the marginal (population) survival function: the event rate e is a
#   property of the hazard and identical with and without administrative censoring.
#   admin = FALSE : one common end of follow-up tau
#   admin = TRUE  : staggered entry. Administrative follow-up ~ tau * Unif(ADMIN_MIN, 1), i.e. nobody is
#                   censored before ADMIN_MIN * tau; then censored uniformly until tau at the end of their
#                   potential follow-up. The proportion of patients with an OBSERVED event is then below e
#                   (reported in sc$cens["events"]).
#   rc            : intensity of RANDOM censoring (loss to follow-up, independent exponential).
#                   rc = 0 none; rc = 1 'regular': LTFU_BASE of the patients lost by tau (if no event first);
#                   rc = 5 severe: five times that hazard (about 15-40% lost in the scenarios of the design).
TAU       <- 1
ADMIN_MIN <- 0.2
LTFU_BASE <- 0.10

# Reference scenario; all other scenarios vary ONE factor at a time (see make_scen)
REF   <- list( C = 0.7, N = 1000, e = 0.9, shape = 1, k = 3, noise = 0 )
LEVELS <- list( C = c(0.6, 0.7, 0.8, 0.9), N = 1000, e = c(0.3, 0.9),
                shape = 1, k = c(1, 3, 10), noise = c(0, 1) )   # shape: only varied in design 'shape_only'
# noise = 1: add k pure-noise predictors (k = number of true predictors per type; the noise
#            variables total k, alternating continuous/binary)

# design = "oat":  C x e grid at the reference (shape 1, k, noise), plus k and noise varied
#                  one at a time at every C (e = 0.9)   -> 20 scenarios. (Weibull shape had no effect, see "shape_only".)
# design = "full": full factorial (4*1*2*3*3*2 = 144 scenarios; large)
# design = "event_rate": non-factorial study of the OVERALL event rate e (N = D / e so that the
#                  number of events is the same in every scenario), at C = 0.7 and 0.8; noise predictors
#                  are added at C = 0.7 for e = 0.1 and 0.9 to show (ir)relevance. 12 scenarios.
# design = "shape_only": the reference scenario (C = 0.7, N = 1000, e = 0.9, k = 3, no noise) with
#                  ONLY the Weibull shape changed (0.7, 1, 1.5). 3 scenarios.
# design = "factorial" (MAIN simulation of the article): C = 0.8 (default; other values via `C`) x overall event rate
#                  e in {0.05, 0.1, 0.3, 0.6, 0.9} x administrative censoring off / on (staggered entry),
#                  N = D / e (about D = 500 events in every scenario), k = 3, no noise predictors, no random
#                  censoring (rc = 0). 12 scenarios. e = 1 - S(tau = 1) in all scenarios; with staggered entry
#                  the number of observed events is below D.
# design = "csens": sensitivity of the main simulation to C: C in {0.6, 0.9} x e in {0.05, 0.3, 0.9}, common end of
#                  follow-up only, N = D / e. 6 scenarios; combine with the C = 0.8 rows of "factorial".
# design = "censoring": registry-style censoring at C = 0.7 (N = D / e): overall event rate e in
#                  {0.1, 0.3, 0.6} x administrative censoring (off = common end of follow-up, on = staggered
#                  entry) x random censoring rc in {0, 1, 5}. 18 scenarios; (off, rc = 0) is the old design.
make_scen <- function( design = c("oat", "full", "event_rate", "shape_only", "censoring", "factorial", "csens"),
                       D = D_EVENTS, C = NULL, e = NULL ){
  design <- match.arg(design)
  if( design == "csens" ){                                      # sensitivity of the main simulation to C
    if( is.null(C) ) C <- c( 0.6, 0.9 )
    if( is.null(e) ) e <- c( 0.05, 0.3, 0.9 )
    s <- expand.grid( C = C, e = e, admin = FALSE, shape = 1, k = 3, noise = 0, rc = 0, KEEP.OUT.ATTRS = FALSE )
    s$N <- round(D / s$e)
    s <- s[order(s$e, s$C), c("C", "N", "e", "shape", "k", "noise", "admin", "rc")]
    rownames(s) <- NULL
    return(s)
  }
  if( design == "factorial" ){
    if( is.null(C) ) C <- 0.8                                   # C hardly moves the curves; sensitivity: C = c(0.6, 0.9)
    if( is.null(e) ) e <- c( 0.01, 0.05, 0.1, 0.3, 0.6, 0.9 )   # 0.01: N = 50,000 at D = 500 (heavy; drop it or lower D to save time)
    s <- expand.grid( C = C, e = e, admin = c(FALSE, TRUE), shape = 1, k = 3, noise = 0, rc = 0,
                      KEEP.OUT.ATTRS = FALSE )
    s$N <- round(D / s$e)
    s <- s[order(s$admin, s$e, s$C), c("C", "N", "e", "shape", "k", "noise", "admin", "rc")]
    rownames(s) <- NULL
    return(s)
  }
  if( design == "censoring" ){
    if( is.null(C) ) C <- 0.7
    if( is.null(e) ) e <- c( 0.1, 0.3, 0.6 )
    s <- expand.grid( C = C, e = e, admin = c(FALSE, TRUE), rc = c(0, 1, 5), shape = 1, k = 3, noise = 0,
                      KEEP.OUT.ATTRS = FALSE )
    s$N <- round(D / s$e)
    s <- s[order(s$C, s$e, s$admin, s$rc), c("C", "N", "e", "shape", "k", "noise", "admin", "rc")]
    rownames(s) <- NULL
    return(s)
  }
  if( design == "shape_only" )
    return(data.frame(C = REF$C, N = REF$N, e = REF$e, shape = c(0.7, 1, 1.5), k = REF$k, noise = REF$noise))
  if( design == "event_rate" ){
    base <- expand.grid( C = c(0.7, 0.8), e = c(0.05, 0.1, 0.3, 0.6, 0.9), shape = 1, k = 3, noise = 0,
                         KEEP.OUT.ATTRS = FALSE )
    ex <- expand.grid( C = 0.7, e = c(0.1, 0.9), shape = 1, k = 3, noise = 1, KEEP.OUT.ATTRS = FALSE )
    s <- rbind( base, ex ); s$N <- round(D / s$e)
    s <- s[, c("C", "N", "e", "shape", "k", "noise")]; rownames(s) <- NULL
    return(s)
  }
  if( design == "full" ) return(expand.grid(LEVELS, KEEP.OUT.ATTRS = FALSE))
  base <- expand.grid( C = LEVELS$C, N = LEVELS$N, e = LEVELS$e, shape = REF$shape,
                       k = REF$k, noise = REF$noise, KEEP.OUT.ATTRS = FALSE )
  vary <- function( name ){
    g <- expand.grid( C = LEVELS$C, v = LEVELS[[name]], KEEP.OUT.ATTRS = FALSE )
    out <- data.frame( C = g$C, N = REF$N, e = REF$e, shape = REF$shape, k = REF$k, noise = REF$noise )
    out[[name]] <- g$v
    out
  }
  s <- rbind( base, vary("k"), vary("noise") )
  s <- s[!duplicated(s), ]
  rownames(s) <- NULL
  s
}
SCEN <- make_scen("oat")

####################### DATA GENERATION #####################
# columns: xc1..xck (continuous), xb1..xbk (binary), z1..zm (noise)
rcov <- function( n, k, noise ){
  X <- c( lapply(seq_len(k), function(i) rnorm(n)),
          lapply(seq_len(k), function(i) rbinom(n, 1, 0.5)) )
  names(X) <- c( paste0("xc", seq_len(k)), paste0("xb", seq_len(k)) )
  if( noise > 0 ){
    Z <- lapply( seq_len(k), function(i) if( i %% 2 == 1 ) rnorm(n) else rbinom(n, 1, 0.5) )
    names(Z) <- paste0( "z", seq_len(k) )
    X <- c( X, Z )
  }
  as.data.frame(X)
}
# equal contribution: binary(0.5) has SD 0.5, so coefficient 2 for binary, 1 for continuous
beta_true <- function( k ) c( rep(1, k), rep(2, k) )
lp_raw <- function( X, k ) as.numeric(as.matrix(X[, seq_len(2 * k), drop = FALSE]) %*% beta_true(k))
rtime  <- function( lp, shape ) (-log(runif(length(lp))) / exp(lp))^(1 / shape)

calibrate_kappa <- function( C, e, shape, k, n = 20000 ){
  set.seed(1)
  X <- rcov( n, k, 0 ); raw <- lp_raw( X, k ); u <- runif(n)
  f <- function( kappa ){
    T <- (-log(u) / exp(kappa * raw))^(1 / shape); tau <- quantile( T, e )
    concordance( Surv(pmin(T, tau), T <= tau) ~ I(kappa * raw), reverse = TRUE )$concordance - C
  }
  uniroot( f, c(0.001, 6) )$root
}

# Censoring times: common end of follow-up tau (admin = FALSE) or staggered administrative follow-up
# tau * Unif(ADMIN_MIN, 1) (admin = TRUE), plus optional random loss to follow-up (rate scaled by rc).
# attr "ltfu": TRUE where the censoring time is the random (LTFU) one.
cens_times <- function( n, tau, admin, rc ){
  cadm <- if( admin ) tau * runif( n, ADMIN_MIN, 1 ) else rep( tau, n )
  crnd <- if( rc > 0 ) rexp( n, rate = rc * -log(1 - LTFU_BASE) / tau ) else rep( Inf, n )
  structure( pmin(cadm, crnd), ltfu = crnd < cadm )
}
setup_scenario <- function( C, N, e, shape, k, noise, admin = FALSE, rc = 0, n_pop = 1e5 ){
  kappa <- calibrate_kappa( C, e, shape, k )    # C refers to the model truncated at the e-quantile of T (as before)
  set.seed(2)
  X <- rcov( n_pop, k, noise ); raw <- kappa * lp_raw( X, k )
  loglam <- shape * log( unname(quantile(rtime(raw, shape), e)) )   # H0(t) = lambda t^shape: P(T <= TAU) = e
  lp <- raw + loglam
  T <- rtime( lp, shape )
  tau <- TAU
  set.seed(4); ct <- cens_times( n_pop, tau, admin, rc ); obs <- T <= ct
  lt <- attr( ct, "ltfu" )
  cens <- c( events = mean(obs), admin_cens = mean(!obs & !lt), ltfu = mean(!obs & lt),
             median_fu = median(pmin(T, ct)) )
  t_f <- unname(quantile(T[obs], FRACS))                    # fixed horizons: fractions of the OBSERVED events
  q   <- c( 0.0005, 0.002, seq(0.01, 0.99, length.out = N_PROF - 4), 0.998, 0.9995 )   # incl. extremes
  pos <- order(lp)[pmax(1, round(q * n_pop))]
  prof <- X[pos, , drop = FALSE]; lp_prof <- lp[pos]
  S_true <- outer( exp(lp_prof), t_f^shape, function(h, H) exp(-h * H) )   # profiles x horizons
  list( C = C, N = N, e = e, shape = shape, k = k, noise = noise, admin = admin, rc = rc, kappa = kappa,
        loglam = loglam, tau = tau, cens = cens, t_f = t_f, prof = prof, S_true = t(S_true),            # S_true: horizons x profiles
        frm = reformulate(names(X), response = "Surv(time, status)") )
}

####################### ONE REPLICATION #####################
one_rep <- function( sc ){
  X <- rcov( sc$N, sc$k, sc$noise ); lp <- sc$kappa * lp_raw( X, sc$k ) + sc$loglam
  T <- rtime( lp, sc$shape )
  ct <- cens_times( sc$N, sc$tau, sc$admin, sc$rc )
  d <- cbind( X, time = pmin(T, ct), status = as.integer(T <= ct) )
  # Speed: survfit() returns a row for every distinct time, censoring times included. With staggered entry
  # every patient has a distinct censoring time (up to N = 50,000 times x 200 profiles), with a common end of
  # follow-up there is a single one. A censoring time only matters through who is at risk at the event times, so
  # moving each censoring time down to the last event time before it gives identical estimates and variances.
  ev <- sort( unique(d$time[d$status == 1]) )
  if( length(ev) > 0 ){
    cc <- d$status == 0
    jj <- findInterval( d$time[cc], ev )
    d$time[cc] <- ifelse( jj > 0, ev[pmax(jj, 1)], ev[1] / 2 )   # before the first event: any time before it
  }
  frm <- sc$frm; environment(frm) <- environment()   # so that survfit() finds `d` (formula env = this function)
  fit <- coxph( frm, data = d, ties = "efron", model = TRUE )
  sf  <- survfit( fit, newdata = sc$prof )
  idx <- findInterval( sc$t_f, sf$time )
  S   <- matrix( NA_real_, length(idx), nrow(sc$prof) )
  seH <- S
  ok  <- idx >= 1
  S[ok, ]   <- as.matrix(sf$surv)[idx[ok], , drop = FALSE]
  seH[ok, ] <- as.matrix(sf$std.err)[idx[ok], , drop = FALSE]   # SE of log S (cumulative hazard)
  list( S = S, seH = seH,
        frac = vapply(sc$t_f, function(t) sum(d$status == 1 & d$time <= t) / sum(d$status), numeric(1)) )
}

####################### RUN ONE SCENARIO #####################
sim_scenario <- function( i, scen = SCEN, R = R_REP ){
  s <- scen[i, ]
  sc <- setup_scenario( s$C, s$N, s$e, s$shape, s$k, s$noise,
                        admin = isTRUE(s$admin), rc = if( is.null(s$rc) ) 0 else s$rc )   # old designs: no admin/rc columns
  reps <- future_lapply( seq_len(R), function(r) one_rep(sc), future.seed = TRUE )
  nH <- length(FRACS); nP <- N_PROF
  to_arr <- function( nm ) array( unlist(lapply(reps, `[[`, nm)), dim = c(nH, nP, length(reps)) )
  S   <- to_arr("S")                                    # horizons x profiles x reps
  seH <- to_arr("seH")
  stopifnot( length(dim(S)) == 3, all(dim(S)[1:2] == c(nH, nP)) )
  fr  <- rowMeans(matrix(unlist(lapply(reps, `[[`, "frac")), nrow = nH))
  list( s = s, sc = sc, S = S, seH = seH, fr = fr, R = length(reps) )
}

summarise_reps <- function( x ){
  s <- x$s; S <- x$S; seH <- x$seH; sc <- x$sc
  nH <- length(FRACS); nP <- N_PROF
  seS <- S * seH
  Nhat <- S * (1 - S) / seS^2
  Nhat[!is.finite(Nhat)] <- NA_real_
  Hh  <- -log(S); z <- qnorm(0.975)
  lo  <- array( exp(-Hh * exp( z * seH / Hh)), dim = dim(S) )
  up  <- array( exp(-Hh * exp(-z * seH / Hh)), dim = dim(S) )
  St  <- matrix( sc$S_true, nrow = nH, ncol = nP )
  out <- list()
  for( h in seq_len(nH) ) for( p in seq_len(nP) ){
    sh <- S[h, p, ]; nh <- Nhat[h, p, ]
    v_emp <- var( sh, na.rm = TRUE )
    n_true <- St[h, p] * (1 - St[h, p]) / v_emp
    out[[length(out) + 1]] <- data.frame(
      C = s$C, N = s$N, e = s$e, shape = s$shape, k = s$k, noise = s$noise,
      admin = sc$admin, rc = sc$rc, R = x$R,
      frac = FRACS[h], frac_realised = x$fr[h],
      profile = p, S_true = St[h, p], S_hat_mean = mean(sh, na.rm = TRUE),
      Ntrue_ratio = n_true / s$N,
      Nhat_med_ratio = median(nh, na.rm = TRUE) / s$N,
      Nhat_q95_ratio = unname(quantile(nh, 0.95, na.rm = TRUE)) / s$N,
      P_gt_N = mean(nh > s$N, na.rm = TRUE), P_gt_1.25N = mean(nh > 1.25 * s$N, na.rm = TRUE),
      Nhat_over_Ntrue = median(nh, na.rm = TRUE) / n_true,
      var_ratio = mean(seS[h, p, ]^2, na.rm = TRUE) / v_emp,
      coverage = mean(lo[h, p, ] <= St[h, p] & St[h, p] <= up[h, p, ], na.rm = TRUE) )
  }
  do.call( rbind, out )
}

run_scenario <- function( i, scen = SCEN, R = R_REP ) summarise_reps(sim_scenario(i, scen, R))

####################### PILOT: MONTE CARLO STANDARD ERRORS AND THE NUMBER OF REPLICATIONS NEEDED #####################
# For each (event fraction x S band) cell the MC SE of the three headline quantities is
# estimated from the replications themselves (cells = averages over the profiles in the band):
#   P(N*-hat > N)   and   CI coverage : sd of the per-replication cell mean / sqrt(R)
#   variance ratio                    : bootstrap over replications (B resamples)
# The R needed for a target MC SE follows from SE ~ 1/sqrt(R).
find_row <- function( scen = SCEN, C, e = REF$e, shape = REF$shape, k = REF$k, noise = REF$noise )
  which(scen$C == C & scen$N == REF$N & scen$e == e & scen$shape == shape & scen$k == k & scen$noise == noise)

mc_cells <- function( x, B = 100 ){
  nH <- length(FRACS); nP <- N_PROF; R <- x$R; N <- x$s$N
  S <- x$S; seS <- x$S * x$seH
  Nhat <- S * (1 - S) / seS^2; Nhat[!is.finite(Nhat)] <- NA_real_
  Hh <- -log(S); z <- qnorm(0.975)
  lo <- exp(-Hh * exp(z * x$seH / Hh)); up <- exp(-Hh * exp(-z * x$seH / Hh))
  St <- matrix( x$sc$S_true, nH, nP )
  ind <- (Nhat > N) + 0
  cov <- (lo <= as.vector(St) & as.vector(St) <= up) + 0   # St recycled along horizons x profiles
  band <- matrix( as.integer(cut(St, S_BANDS, include.lowest = TRUE, right = FALSE)), nH, nP )
  boot_idx <- replicate( B, sample.int(R, R, replace = TRUE), simplify = FALSE )
  vr_boot <- vapply( boot_idx, function(ix) {
    v <- apply(S[, , ix, drop = FALSE], c(1, 2), var, na.rm = TRUE)
    mean_se2 <- apply(seS[, , ix, drop = FALSE]^2, c(1, 2), mean, na.rm = TRUE)
    as.vector(mean_se2 / v)
  }, numeric(nH * nP) )                                      # (nH*nP) x B
  rows <- list()
  for( h in seq_len(nH) ) for( b in sort(unique(band[h, ])) ){
    ps <- which(band[h, ] == b)
    if( !length(ps) ) next
    mP  <- colMeans( matrix(ind[h, ps, ], nrow = length(ps)), na.rm = TRUE )
    mC  <- colMeans( matrix(cov[h, ps, ], nrow = length(ps)), na.rm = TRUE )
    cell_idx <- (ps - 1) * nH + h                            # position in as.vector(nH x nP)
    vrb <- colMeans( vr_boot[cell_idx, , drop = FALSE], na.rm = TRUE )
    rows[[length(rows) + 1]] <- data.frame(
      C = x$s$C, e = x$s$e, shape = x$s$shape, k = x$s$k, noise = x$s$noise, frac = FRACS[h],
      band_id = b,
      n_profiles = length(ps), P = mean(mP), se_P = sd(mP) / sqrt(R),
      cover = mean(mC), se_cov = sd(mC) / sqrt(R),
      vr = mean(vrb), se_vr = sd(vrb) )
  }
  do.call( rbind, rows )
}

pilot_mc <- function( rows = NULL, scen = SCEN, R_pilot = 300, target_P = 0.01, target_cov = 0.005,
                      target_vr = 0.05, pct = 0.9, workers = max(1, parallel::detectCores() - 1) ){
  if( is.null(rows) && !identical(scen, SCEN) ) stop("give `rows` (rows of `scen`) for a design other than SCEN")
  if( is.null(rows) ) rows <- c( find_row(C = 0.7), find_row(C = 0.9), find_row(C = 0.7, e = 0.3),
                                 find_row(C = 0.7, k = 10, noise = 1) )
  if( workers > 1 ) plan( multisession, workers = workers ) else plan(sequential)   # workers = 1: no cluster
  cells <- do.call( rbind, lapply(rows, function(i) mc_cells(sim_scenario(i, scen, R_pilot))) )
  plan(sequential)
  lab <- levels(cut(0.5, S_BANDS, include.lowest = TRUE, right = FALSE))
  cells$band <- lab[cells$band_id]
  need <- function( se, target ) ceiling(R_pilot * (se / target)^2)
  cells$R_P   <- need( cells$se_P, target_P )
  cells$R_cov <- need( cells$se_cov, target_cov )
  cells$R_vr  <- need( cells$se_vr, target_vr )
  # only cells with at least 5 profiles and a P that is not trivially 0/1 drive the P requirement
  informative <- cells$n_profiles >= 5
  q <- function( x, ok ) ceiling(as.numeric(quantile(x[ok], pct, na.rm = TRUE)))
  rec <- c( P = q(cells$R_P, informative & cells$P > 0.01 & cells$P < 0.99),
            coverage = q(cells$R_cov, informative), var_ratio = q(cells$R_vr, informative) )
  message(sprintf("Replications needed (%.0fth percentile over cells): P(N*>N) %.0f, coverage %.0f, variance ratio %.0f",
                  100 * pct, rec[["P"]], rec[["coverage"]], rec[["var_ratio"]]))
  message( "Recommended R_REP: ", ceiling(max(rec) / 100) * 100 )
  list( cells = cells, required = rec, R_recommended = ceiling(max(rec) / 100) * 100 )
}

run_all <- function( scen = SCEN, R = R_REP, workers = max(1, parallel::detectCores() - 1) ){
  if( workers > 1 ) plan( multisession, workers = workers ) else plan(sequential)   # workers = 1: no cluster
  res <- do.call( rbind, lapply(seq_len(nrow(scen)), function(i) {
    message(sprintf("Scenario %.0f/%.0f: C=%.1f N=%.0f e=%.2f shape=%.1f k=%.0f noise=%.0f admin=%s rc=%.0f", i,
                    nrow(scen), scen$C[i], scen$N[i], scen$e[i], scen$shape[i], scen$k[i], scen$noise[i],
                    isTRUE(scen$admin[i]), if( is.null(scen$rc) ) 0 else scen$rc[i]))
    run_scenario(i, scen, R)
  }) )
  plan(sequential)
  res
}

####################### SUMMARIES #####################
KEYS <- c( "C", "N", "e", "shape", "k", "noise" )
band_of <- function( S ) cut( S, S_BANDS, include.lowest = TRUE, right = FALSE )

# Profiles sit at quantiles of the true linear predictor. The most extreme ones (below the 1st
# and above the 99th percentile) are extrapolations in a sample of ~1000 patients, so by default
# the summaries use quantiles 0.01-0.99 only (trim = FALSE keeps everything).
PROF_Q <- c( 0.0005, 0.002, seq(0.01, 0.99, length.out = N_PROF - 4), 0.998, 0.9995 )
trim_core <- function( res, trim = TRUE, q = c(0.01, 0.99) ){
  res$lp_q <- PROF_Q[res$profile]
  if( trim ) res[res$lp_q >= q[1] & res$lp_q <= q[2], ] else res
}

summarise_sim <- function( res, trim = TRUE ){
  res <- trim_core( res, trim ); res$S_band <- band_of(res$S_true)
  aggregate( cbind(P_gt_N, P_gt_1.25N, Ntrue_ratio, Nhat_q95_ratio, Nhat_over_Ntrue,
                   var_ratio, coverage) ~ C + N + e + shape + k + noise + frac + S_band,
             data = res, FUN = mean )
}

# A1: per scenario and S band: the smallest event fraction f_min such that P(N*-hat > N) is at
#     most `tol` at f_min AND at all larger event fractions observed for that band.
#     f_first = smallest event fraction at which the band is observed at all. If f_min == f_first
#     the band is never a problem once it exists.
fmin_by_band <- function( res, tol = 0.05, trim = TRUE ){
  res <- trim_core( res, trim ); res$S_band <- band_of(res$S_true)
  keys <- intersect( c(KEYS, "admin", "rc"), names(res) )      # admin, rc: censoring design
  g <- aggregate( res["P_gt_N"], res[c(keys, "S_band", "frac")], FUN = mean )
  g$key <- do.call( paste, c(g[c(keys, "S_band")], sep = "|") )
  do.call( rbind, lapply(split(g, g$key), function(x) {
    x <- x[order(x$frac), ]
    ok_from <- rev(cumsum(rev(x$P_gt_N > tol))) == 0
    data.frame(x[1, c(keys, "S_band")], tol = tol, f_first = min(x$frac),
               f_min = if( any(ok_from) ) min(x$frac[ok_from]) else NA_real_)
  }) )
}

# A2: reliability by S band (restricted to event fractions >= f_min_use)
slow_table <- function( res, f_min_use = 0.25, trim = TRUE ){
  res <- trim_core( res, trim )
  r <- res[res$frac >= f_min_use & res$S_true > 1e-6, ]       # S below 1e-6: variance is 0 in floating point
  r$S_band <- band_of(r$S_true)
  aggregate( cbind(P_gt_N, var_ratio, coverage, Nhat_over_Ntrue, Nhat_q95_ratio) ~ C + S_band,
             data = r, FUN = mean )
}

# Bias of S-hat by band (to understand under-coverage at high C): mean(S_hat) - S_true
bias_table <- function( res, trim = TRUE ){
  res <- trim_core( res, trim ); res$S_band <- band_of(res$S_true)
  res$bias <- res$S_hat_mean - res$S_true
  aggregate( cbind(bias, coverage) ~ C + S_band, data = res, FUN = mean )
}

# Does a design factor (k, noise; shape in design 'shape_only') change P(N*-hat > N)? Compares each level with the
# reference level, cell by cell (C x event fraction x S band; cells with >= min_n profiles).
factor_effect <- function( res, factor = c("shape", "k", "noise"), e = REF$e, min_n = 5, trim = TRUE ){
  factor <- match.arg(factor)
  r <- trim_core( res, trim )
  r <- r[r$e %in% e, ]
  for( o in setdiff(c("shape", "k", "noise"), factor) ) r <- r[r[[o]] == REF[[o]], ]
  r$S_band <- band_of(r$S_true); r$lev <- r[[factor]]
  g <- aggregate( P_gt_N ~ C + e + frac + S_band + lev, data = r, FUN = mean )
  n <- aggregate( profile ~ C + e + frac + S_band + lev, data = r, FUN = length )
  g$n <- n$profile; g <- g[g$n >= min_n, ]
  ref <- g[g$lev == REF[[factor]], c("C", "e", "frac", "S_band", "P_gt_N")]; names(ref)[5] <- "P_ref"
  m <- merge( g[g$lev != REF[[factor]], ], ref, by = c("C", "e", "frac", "S_band") )
  m$absdiff <- abs(m$P_gt_N - m$P_ref)
  out <- aggregate( absdiff ~ lev + e, data = m, FUN = function( x ) c(mean = mean(x), max = max(x), n_cells = length(x)) )
  out <- cbind( out[c("lev", "e")], as.data.frame(out$absdiff) )
  names(out)[1] <- factor; out
}

# Event-rate study: does e still matter after conditioning on S (band) or on the event fraction?
# For each conditioning set, P(N*-hat > N) is computed per cell and e level; the spread (max - min)
# across the e levels in a cell measures how much e still matters. A smaller spread = a better
# organising variable.
e_spread <- function( res, min_n = 5, trim = TRUE ){
  r <- trim_core( res, trim ); r <- r[r$shape == 1 & r$k == 3 & r$noise == 0, ]
  r$S_band <- band_of(r$S_true)
  spread <- function( vars ){
    keys <- r[c("C", "e", vars)]
    g <- aggregate( r["P_gt_N"], keys, FUN = mean )
    n <- aggregate( r["profile"], keys, FUN = length )
    g <- g[n$profile >= min_n, ]
    sp <- aggregate( g["P_gt_N"], g[c("C", vars)], FUN = function( x ) if( length(x) > 1 ) max(x) - min(x) else NA_real_ )
    sp <- sp[!is.na(sp$P_gt_N), ]
    data.frame( conditioning = paste(vars, collapse = " x "), mean_spread = mean(sp$P_gt_N),
                max_spread = max(sp$P_gt_N), n_cells = nrow(sp) )
  }
  rbind( spread("S_band"), spread("frac"), spread(c("S_band", "frac")) )
}

####################### PLOTS #####################
# Scales for the true survival S on the x-axis:
#   logit   : symmetric around S = 0.5 (default)
#   cloglog : -log(-log S)
#   log     : log S
#   linear  : S
S_SCALE <- list(
  logit   = list(f = qlogis, br = c(0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999), width = 0.4,
                 xlim = c(-6.2, 6.2), name = "True survival S(t|x), logit scale"),
  cloglog = list(f = function( S ) -log(-log(S)), br = c(0.01, 0.05, 0.2, 0.5, 0.8, 0.95, 0.99), width = 0.3,
                 xlim = c(-1.7, 5.3), name = "True survival S(t|x), complementary log-log scale"),
  log     = list(f = log, br = c(0.01, 0.05, 0.2, 0.5, 0.8, 1), width = 0.12,
                 xlim = c(log(0.005), 0), name = "True survival S(t|x), log scale"),
  linear  = list(f = identity, br = c(0, 0.25, 0.5, 0.75, 1), width = 0.025,
                 xlim = c(0, 1), name = "True survival S(t|x)") )
scale_x_S <- function( sc )
  scale_x_continuous(breaks = sc$f(sc$br), labels = sc$br, name = sc$name,
                     guide = guide_axis(check.overlap = TRUE))
bin_s <- function( d, vars, groups, sc, width = sc$width, min_n = 5 ){
  d$xs <- sc$f(d$S_true); d <- d[is.finite(d$xs), ]
  d$xbin <- floor(d$xs / width) * width + width / 2
  a <- aggregate( d[vars], d[c(groups, "xbin")], FUN = mean )
  a$n <- aggregate( d[vars[1]], d[c(groups, "xbin")], FUN = length )[[vars[1]]]
  a[a$n >= min_n, ]
}
frac_lab <- function( f ) factor( f, levels = FRACS, labels = paste0(100 * FRACS, "% of events") )

# returns list(p1 = list(logit=, cloglog=, log=, linear=), p2 = ..., p3 = ...); PNGs go to outdir
plot_sim <- function( res, outdir = "results/", scales = names(S_SCALE), trim = TRUE,
                      fracs_show = c(0.10, 0.25, 0.40, 0.60) ){
  dir.create( outdir, showWarnings = FALSE, recursive = TRUE )
  res <- trim_core( res, trim )
  ref <- res[res$N == REF$N & res$e == REF$e & res$shape == REF$shape &
               res$k == REF$k & res$noise == REF$noise, ]
  a3 <- res[res$N == REF$N & res$e == REF$e & res$C %in% c(0.7, 0.8) & res$frac %in% fracs_show, ]
  a3$setting <- sprintf( "k = %.0f%s", a3$k, ifelse(a3$noise == 1, ", +noise", "") )
  Rm <- min(res$R); mc <- 1.96 * sqrt(0.95 * 0.05 / Rm)
  out <- list( p1 = list(), p2 = list(), p3 = list() )
  for( nm in scales ){
    sc <- S_SCALE[[nm]]
    cl <- coord_cartesian(xlim = sc$xlim)
    # (1) P(N*-hat > N) against true S: one line per C, one panel per event fraction
    b1 <- bin_s( ref, "P_gt_N", c("C", "frac"), sc ); b1$f <- frac_lab(b1$frac)
    p1 <- ggplot(b1, aes(xbin, P_gt_N, colour = factor(C))) +
      scale_colour_manual(values = ramp_pal(length(unique(b1$C)))) +
      geom_hline(yintercept = 0.05, linetype = 2, colour = "grey40") +
      geom_line(linewidth = 0.8) + geom_point(size = 1) +
      facet_wrap(~f, ncol = 3) + cl + scale_x_S(sc) +
      scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      labs(y = "Probability that N*-hat > N", colour = "Harrell's C",
           title = "N*-hat > N: reference scenario", subtitle = "dashed line: 5%") +
      theme_art(11)
    ggsave(file.path(outdir, sprintf("sim_guardrails_PgtN_%s.png", nm)), p1, width = 12, height = 8, dpi = 200, bg = "white")
    out$p1[[nm]] <- p1
    
    # (2) CI coverage against true S (all scenarios), one line per C, with Monte Carlo band
    b2 <- bin_s( res, "coverage", "C", sc, width = sc$width * 0.8 )
    rs <- res; rs$xs <- sc$f(rs$S_true)    # computed now: aes() is evaluated lazily, `sc` changes with the loop
    p2 <- ggplot(rs, aes(xs, coverage)) +
      scale_colour_manual(values = ramp_pal(length(unique(b2$C)))) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.95 - mc, ymax = 0.95 + mc, fill = "grey85") +
      geom_point(alpha = 0.04, size = 0.6, colour = "grey30") +
      geom_hline(yintercept = 0.95, linetype = 2) +
      geom_line(data = b2, aes(xbin, coverage, colour = factor(C)), linewidth = 1) +
      coord_cartesian(xlim = sc$xlim, ylim = c(0.88, 1)) + scale_x_S(sc) +
      labs(y = "Coverage of the 95% log-log CI", colour = "Harrell's C",
           subtitle = sprintf("grey band: +/- 1.96 Monte Carlo SE of one coverage estimate (R = %.0f); dots: profiles x horizons", Rm)) +
      theme_art(11)
    ggsave(file.path(outdir, sprintf("sim_guardrails_coverage_%s.png", nm)), p2, width = 10, height = 6, dpi = 200, bg = "white")
    out$p2[[nm]] <- p2
    
    # (3) k and noise (e = 0.9), C = 0.7 and 0.8, selected event fractions
    b3 <- bin_s( a3, "P_gt_N", c("C", "frac", "setting"), sc, width = sc$width * 1.3 )
    b3$f <- frac_lab(b3$frac); b3$Cl <- paste( "C =", b3$C )
    p3 <- ggplot(b3, aes(xbin, P_gt_N, colour = setting, linetype = setting)) +
      scale_colour_manual(values = ramp_pal(length(unique(b3$setting)))) +
      geom_hline(yintercept = 0.05, linetype = 2, colour = "grey40") +
      geom_line(linewidth = 0.7) +
      facet_grid(Cl ~ f) + cl + scale_x_S(sc) +
      scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      labs(y = "Probability that N*-hat > N", colour = NULL, linetype = NULL) +
      theme_art(11)
    ggsave(file.path(outdir, sprintf("sim_guardrails_settings_%s.png", nm)), p3, width = 14, height = 7.5, dpi = 200, bg = "white")
    out$p3[[nm]] <- p3
  }
  invisible(out)
}

# Event-rate study (design "event_rate"). Returns list(p1, p2, p3, p4); p1/p3/p4 per x-axis scale.
plot_erate <- function( res, outdir = "results/", scales = c("logit", "linear"), trim = TRUE,
                        fracs_show = c(0.10, 0.25, 0.40, 0.60, 0.90), C_show = 0.7 ){
  dir.create( outdir, showWarnings = FALSE, recursive = TRUE )
  res <- trim_core( res, trim )
  es <- sort(unique(res$e)); elab <- function( e ) factor( e, levels = es, labels = format(es, drop0trailing = TRUE) )
  r <- res[res$shape == 1 & res$k == 3 & res$noise == 0, ]
  r$elab <- elab(r$e); r$S_band <- band_of(r$S_true)
  out <- list( p1 = list(), p3 = list(), p4 = list() )
  
  # p2: P(N*-hat > N) against the EVENT FRACTION, within bands of true S, one line per overall event rate
  g <- aggregate( P_gt_N ~ C + elab + frac + S_band, data = r, FUN = mean )
  n <- aggregate( profile ~ C + elab + frac + S_band, data = r, FUN = length )
  g <- g[n$profile >= 5 & g$C == C_show, ]
  p2 <- ggplot(g, aes(frac, P_gt_N, colour = elab)) +
    scale_colour_manual(values = ramp_pal(length(es))) +
    geom_hline(yintercept = 0.05, linetype = 2, colour = "grey40") +
    geom_line(linewidth = 0.8) + geom_point(size = 1.2) +
    facet_wrap(~S_band, ncol = 3) +
    scale_x_continuous(breaks = c(0.05, 0.25, 0.5, 0.75, 0.9)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(x = "Fraction of events observed", y = "Probability that N*-hat > N", colour = NULL,
         title = sprintf("Within bands of true S (C = %.1f): does the event fraction still matter?", C_show)) +
    theme_art(11)
  ggsave(file.path(outdir, "sim_erate_vs_fraction.png"), p2, width = 12, height = 8, dpi = 200, bg = "white")
  
  # settings (noise predictors) at C = 0.7 for e = 0.1 and 0.9
  a4 <- res[res$C == 0.7 & res$e %in% c(0.1, 0.9) & res$k == 3 & res$frac %in% fracs_show, ]
  a4$setting <- ifelse( a4$noise == 1, "with noise predictors", "no noise predictors" )
  a4$el <- elab(a4$e)
  
  for( nm in scales ){
    sc <- S_SCALE[[nm]]; cl <- coord_cartesian(xlim = sc$xlim)
    # p1: P(N*-hat > N) against true S, one line per overall event rate, panels = event fraction
    b1 <- bin_s( r[r$C == C_show & r$frac %in% fracs_show, ], "P_gt_N", c("C", "elab", "frac"), sc )
    b1$f <- frac_lab(b1$frac)
    p1 <- ggplot(b1, aes(xbin, P_gt_N, colour = elab)) +
      scale_colour_manual(values = ramp_pal(length(es))) +
      geom_hline(yintercept = 0.05, linetype = 2, colour = "grey40") +
      geom_line(linewidth = 0.8) + geom_point(size = 1) +
      facet_wrap(~f, ncol = 3) + cl + scale_x_S(sc) +
      scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      labs(y = "Probability that N*-hat > N", colour = NULL,
           title = sprintf("N*-hat > N by overall event rate (C = %.1f)", C_show)) +
      theme_art(11)
    ggsave(file.path(outdir, sprintf("sim_erate_vs_S_%s.png", nm)), p1, width = 12, height = 7, dpi = 200, bg = "white")
    out$p1[[nm]] <- p1
    
    # p3: pooled over event fractions, one line per overall event rate and C
    b3 <- bin_s( r, "P_gt_N", c("C", "elab"), sc, width = sc$width * 0.8 )
    p3 <- ggplot(b3, aes(xbin, P_gt_N, colour = elab, linetype = factor(C))) +
      scale_colour_manual(values = ramp_pal(length(es))) +
      geom_hline(yintercept = 0.05, linetype = 2, colour = "grey40") +
      geom_line(linewidth = 0.8) + cl + scale_x_S(sc) +
      scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      labs(y = "Probability that N*-hat > N", colour = NULL, linetype = "C",
           subtitle = "pooled over event fractions") +
      theme_art(11)
    ggsave(file.path(outdir, sprintf("sim_erate_pooled_%s.png", nm)), p3, width = 9, height = 6, dpi = 200, bg = "white")
    out$p3[[nm]] <- p3
    
    # p4: noise predictors (C = 0.7), e = 0.1 and 0.9 (only if the results contain noise scenarios)
    if( any(a4$noise == 1) ){
      b4 <- bin_s( a4, "P_gt_N", c("el", "frac", "setting"), sc, width = sc$width * 1.3 )
      b4$f <- frac_lab(b4$frac)
      p4 <- ggplot(b4, aes(xbin, P_gt_N, colour = setting, linetype = setting)) +
        scale_colour_manual(values = ramp_pal(length(unique(b4$setting)))) +
        geom_hline(yintercept = 0.05, linetype = 2, colour = "grey40") +
        geom_line(linewidth = 0.7) + facet_grid(el ~ f) + cl + scale_x_S(sc) +
        scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
        labs(y = "Probability that N*-hat > N", colour = NULL, linetype = NULL) +
        theme_art(11)
      ggsave(file.path(outdir, sprintf("sim_erate_settings_%s.png", nm)), p4, width = 14, height = 6.5, dpi = 200, bg = "white")
      out$p4[[nm]] <- p4
    }
  }
  out$p2 <- p2
  invisible(out)
}

# Shape-only study (design "shape_only"): everything but the Weibull shape is kept as in REF.
shape_table <- function( res, trim = TRUE ){
  res <- trim_core( res, trim ); res <- res[res$S_true > 1e-6, ]; res$S_band <- band_of(res$S_true)
  aggregate( cbind(P_gt_N, var_ratio, coverage, Nhat_over_Ntrue) ~ shape + S_band, data = res, FUN = mean )
}

plot_shape <- function( res, outdir = "results/", scales = c("logit", "linear"), trim = TRUE,
                        fracs_show = c(0.10, 0.25, 0.40, 0.60, 0.90) ){
  dir.create( outdir, showWarnings = FALSE, recursive = TRUE )
  res <- trim_core( res, trim )
  res$shp <- factor( res$shape, levels = sort(unique(res$shape)), labels = paste("shape", sort(unique(res$shape))) )
  mc <- 1.96 * sqrt(0.95 * 0.05 / min(res$R))
  out <- list( p1 = list(), p2 = list() )
  for( nm in scales ){
    sc <- S_SCALE[[nm]]; cl <- coord_cartesian(xlim = sc$xlim)
    b1 <- bin_s( res[res$frac %in% fracs_show, ], "P_gt_N", c("shp", "frac"), sc ); b1$f <- frac_lab(b1$frac)
    p1 <- ggplot(b1, aes(xbin, P_gt_N, colour = shp, linetype = shp)) +
      scale_colour_manual(values = ramp_pal(nlevels(res$shp))) +
      geom_hline(yintercept = 0.05, linetype = 2, colour = "grey40") +
      geom_line(linewidth = 0.8) + geom_point(size = 1) +
      facet_wrap(~f, ncol = 3) + cl + scale_x_S(sc) +
      scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      labs(y = "Probability that N*-hat > N", colour = NULL, linetype = NULL,
           title = "Only the Weibull shape changes (C = 0.7, 90% events, N = 1000, k = 3)") +
      theme_art(11)
    ggsave(file.path(outdir, sprintf("sim_shape_PgtN_%s.png", nm)), p1, width = 12, height = 7, dpi = 200, bg = "white")
    out$p1[[nm]] <- p1
    b2 <- bin_s( res, "coverage", "shp", sc, width = sc$width * 0.8 )
    p2 <- ggplot(b2, aes(xbin, coverage, colour = shp, linetype = shp)) +
      scale_colour_manual(values = ramp_pal(nlevels(res$shp))) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.95 - mc, ymax = 0.95 + mc, fill = "grey85") +
      geom_hline(yintercept = 0.95, linetype = 2) + geom_line(linewidth = 0.9) +
      coord_cartesian(xlim = sc$xlim, ylim = c(0.9, 1)) + scale_x_S(sc) +
      labs(y = "Coverage of the 95% log-log CI", colour = NULL, linetype = NULL) +
      theme_art(11)
    ggsave(file.path(outdir, sprintf("sim_shape_coverage_%s.png", nm)), p2, width = 9, height = 5.5, dpi = 200, bg = "white")
    out$p2[[nm]] <- p2
  }
  invisible(out)
}

####################### MAIN SIMULATION (DESIGN "FACTORIAL"): C X OVERALL EVENT RATE X ADMINISTRATIVE CENSORING #####################
ADMIN_LAB <- c( "No censoring", "Uniform censoring" )

main_data <- function( res, trim = TRUE ){
  r <- trim_core( res, trim )
  if( is.null(r$admin) ) r$admin <- FALSE                       # results from older designs
  if( !is.null(r$rc) ) r <- r[r$rc == 0, ]
  r <- r[r$shape == 1 & r$k == 3 & r$noise == 0, ]
  r$adm_lab <- factor( ifelse(r$admin, ADMIN_LAB[2], ADMIN_LAB[1]), levels = ADMIN_LAB )
  es <- sort(unique(r$e)); r$el <- factor( r$e, levels = es, labels = format(es, drop0trailing = TRUE) )
  cs <- sort(unique(r$C)); r$Cl <- factor( r$C, levels = cs, labels = sprintf("C = %.1f", cs) )
  r
}

# p1: P(N*-hat > N) against the true survival, pooled over the fractions of events observed (lines);
#     dots = one point per fraction of events observed (dots = FALSE: lines only).
#     Panels: rows = administrative censoring, columns = C; colour = overall proportion with an event.
# p2: coverage of the 95% log-log CI, same layout (pooled over fractions).
plot_main <- function( res, outdir = "results/", scales = c("logit", "linear"), trim = TRUE, dots = TRUE ){
  dir.create( outdir, showWarnings = FALSE, recursive = TRUE )
  r <- main_data( res, trim )
  mc <- 1.96 * sqrt(0.95 * 0.05 / min(r$R))
  cols <- scale_colour_manual(values = ramp_pal(nlevels(r$el)), name = "Event rate, 1 \u2212 S(\u03c4)")
  lab <- labeller(adm_lab = label_wrap_gen(20))
  out <- list( p1 = list(), p2 = list() )
  for( nm in scales ){
    sc <- S_SCALE[[nm]]; cl <- coord_cartesian(xlim = sc$xlim)
    b1 <- bin_s( r, "P_gt_N", c("Cl", "adm_lab", "el"), sc )
    p1 <- ggplot(b1, aes(xbin, P_gt_N, colour = el)) +
      geom_hline(yintercept = 0.05, linetype = 2, colour = "grey40")
    if( dots ){
      b1f <- bin_s( r, "P_gt_N", c("Cl", "adm_lab", "el", "frac"), sc, width = sc$width * 1.5, min_n = 3 )
      p1 <- p1 + geom_point(data = b1f, size = 0.8, alpha = 0.3)
    }
    p1 <- p1 + geom_line(linewidth = 0.8) +
      facet_grid(adm_lab ~ Cl, labeller = lab) + cl + scale_x_S(sc) + cols +
      scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      labs(y = "Predictions with N* > N") +
      guides(colour = guide_legend(nrow = 1, override.aes = list(alpha = 1))) +
      theme_art(11)
    ggsave(file.path(outdir, sprintf("sim_main_PgtN_%s.png", nm)), p1, width = 14, height = 6.5, dpi = 300, bg = "white")
    out$p1[[nm]] <- p1
    
    b2 <- bin_s( r, "coverage", c("Cl", "adm_lab", "el"), sc, width = sc$width * 0.8 )
    p2 <- ggplot(b2, aes(xbin, coverage, colour = el)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.95 - mc, ymax = 0.95 + mc, fill = "grey88") +
      geom_hline(yintercept = 0.95, linetype = 2, colour = "grey30") + geom_line(linewidth = 0.8) +
      facet_grid(adm_lab ~ Cl, labeller = lab) + coord_cartesian(xlim = sc$xlim, ylim = c(0.9, 1)) +
      scale_x_S(sc) + cols +
      labs(y = "Coverage of the 95% log-log CI",
           subtitle = sprintf("grey band: +/- 1.96 Monte Carlo SE of one coverage estimate (R = %.0f)", min(r$R))) +
      guides(colour = guide_legend(nrow = 1)) +
      theme_art(11)
    ggsave(file.path(outdir, sprintf("sim_main_coverage_%s.png", nm)), p2, width = 14, height = 6.5, dpi = 300, bg = "white")
    out$p2[[nm]] <- p2
  }
  invisible(out)
}

####################### MAIN FIGURE, SINGLE PANEL #####################
# What is plotted. Every prediction (one covariate profile, one time point, one simulated dataset) has a TRUE
# survival probability S (known from the data-generating model) and an estimated N*. The predictions are grouped
# by their true S (bins on the logit scale) and, for each event rate, the line gives the percentage of the
# predictions in the bin that have N* > N. Equivalently: the probability that a prediction for a patient with
# true survival S has N* > N. All Harrell's C values, both censoring settings (admin = "pooled") and all
# fractions of events observed are pooled; every prediction has the same weight.
# Dots (dots = TRUE): the same percentage for each fraction of events observed separately (one dot per fraction
# and S bin), i.e. the spread that the lines average over.
# Dashed horizontal line: 5%. Points + labels: the true survival at which the line crosses 5% (upper end, and
# lower end where it exists), interpolated on the logit scale.
S_cross_main <- function( r, tol = 0.05, width = 0.2, min_n = 20, group = "el" ){
  b <- bin_s( r, "P_gt_N", group, S_SCALE$logit, width = width, min_n = min_n )
  first <- function( x, y ){                       # first bin (in the order given) with y > tol, interpolated
    i <- which(y > tol)[1]
    if( is.na(i) || i == 1 ) return(NA_real_)
    plogis(x[i - 1] + (tol - y[i - 1]) / (y[i] - y[i - 1]) * (x[i] - x[i - 1]))
  }
  do.call( rbind, lapply(split(b, b[group], drop = TRUE), function(x) {
    up <- x[x$xbin >= 0, ]; up <- up[order(up$xbin), ]
    lo <- x[x$xbin <  0, ]; lo <- lo[order(-lo$xbin), ]
    data.frame(x[1, group, drop = FALSE], S_up = first(up$xbin, up$P_gt_N), S_low = first(lo$xbin, lo$P_gt_N),
               row.names = NULL)
  }) )
}

plot_main1 <- function( res, outdir = "results/", scales = c("logit", "linear"), trim = TRUE,
                        admin = c("pooled", "common", "staggered"), dots = FALSE, tol = 0.05, label = TRUE ){
  admin <- match.arg(admin)
  dir.create( outdir, showWarnings = FALSE, recursive = TRUE )
  r <- main_data( res, trim )
  if( admin == "common" ) r <- r[!r$admin, ] else if( admin == "staggered" ) r <- r[r$admin, ]
  cols <- scale_colour_manual(values = ramp_pal(nlevels(r$el)), name = "Event rate, 1 \u2212 S(\u03c4)")
  cr <- S_cross_main( r, tol )
  mc <- 1.96 * sqrt(0.95 * 0.05 / min(r$R))
  out <- list( p1 = list(), p2 = list() )
  for( nm in scales ){
    sc <- S_SCALE[[nm]]
    b1 <- bin_s( r, "P_gt_N", "el", sc, min_n = 20 )
    p <- ggplot(b1, aes(xbin, P_gt_N, colour = el)) +
      geom_hline(yintercept = tol, linetype = 2, colour = "grey40")
    if( dots ){
      b1f <- bin_s( r, "P_gt_N", c("el", "frac"), sc, width = sc$width * 1.5, min_n = 5 )
      p <- p + geom_point(data = b1f, size = 0.8, alpha = 0.3)
    }
    p <- p + geom_line(linewidth = 0.9)
    if( label ){
      cu <- cr[is.finite(cr$S_up), ];  cu$x <- sc$f(cu$S_up);  cu$y <- rep( tol, nrow(cu) )    # x computed now (aes is lazy)
      cl_ <- cr[is.finite(cr$S_low), ]; cl_$x <- sc$f(cl_$S_low); cl_$y <- rep( tol, nrow(cl_) )
      dx <- 0.12 * sc$width / 0.4
      if( nrow(cu) ) p <- p + geom_point(data = cu, aes(x = x, y = y), size = 1.8, show.legend = FALSE) +
        geom_text(data = cu, aes(x = x, y = y, label = sprintf("%.3f", S_up)), angle = 90, hjust = 0,
                  nudge_y = 0.025, nudge_x = -dx, size = 3, show.legend = FALSE)
      if( nrow(cl_) ) p <- p + geom_point(data = cl_, aes(x = x, y = y), size = 1.8, show.legend = FALSE) +
        geom_text(data = cl_, aes(x = x, y = y, label = sprintf("%.3f", S_low)), angle = 90, hjust = 0,
                  nudge_y = 0.025, nudge_x = dx, size = 3, show.legend = FALSE)
    }
    p <- p + coord_cartesian(xlim = sc$xlim) + scale_x_S(sc) + cols +
      scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      labs(y = "Predictions with N* > N") +
      guides(colour = guide_legend(nrow = 1, override.aes = list(alpha = 1))) +
      theme_art(11)
    ggsave(file.path(outdir, sprintf("sim_main1_PgtN_%s.png", nm)), p, width = 8.5, height = 5.5, dpi = 300, bg = "white")
    out$p1[[nm]] <- p
    b2 <- bin_s( r, "coverage", "el", sc, width = sc$width * 0.8, min_n = 20 )
    p2 <- ggplot(b2, aes(xbin, coverage, colour = el)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.95 - mc, ymax = 0.95 + mc, fill = "grey88") +
      geom_hline(yintercept = 0.95, linetype = 2, colour = "grey30") + geom_line(linewidth = 0.9) +
      coord_cartesian(xlim = sc$xlim, ylim = c(0.9, 1)) + scale_x_S(sc) + cols +
      labs(y = "Coverage of the 95% log-log CI",
           subtitle = sprintf("grey band: +/- 1.96 Monte Carlo SE of one coverage estimate (R = %.0f)", min(r$R))) +
      guides(colour = guide_legend(nrow = 1)) + theme_art(11)
    ggsave(file.path(outdir, sprintf("sim_main1_coverage_%s.png", nm)), p2, width = 8.5, height = 5.5, dpi = 300, bg = "white")
    out$p2[[nm]] <- p2
  }
  attr( out, "crossings" ) <- cr
  invisible(out)
}

# Supplement: does administrative censoring change anything? One panel, line type = censoring, colour = event
# rate, pooled over C. `admin_check()` gives the largest absolute difference per event rate over the bins.
plot_admin <- function( res, outdir = "results/", scales = "logit", trim = TRUE, min_n = 5 ){
  dir.create( outdir, showWarnings = FALSE, recursive = TRUE )
  r <- main_data( res, trim )
  cols <- scale_colour_manual(values = ramp_pal(nlevels(r$el)), name = "Event rate, 1 \u2212 S(\u03c4)")
  out <- list()
  for( nm in scales ){
    sc <- S_SCALE[[nm]]
    b <- bin_s( r, "P_gt_N", c("el", "adm_lab"), sc, min_n = min_n )   # bins = cells (profile x time) with at least min_n predictions
    p <- ggplot(b, aes(xbin, P_gt_N, colour = el, linetype = adm_lab)) +
      geom_hline(yintercept = 0.05, linetype = 2, colour = "grey40") + geom_line(linewidth = 0.8) +
      coord_cartesian(xlim = sc$xlim) + scale_x_S(modifyList(sc, list(name = "S(t|x), logit scale"))) + cols +
      scale_linetype_manual(values = c("solid", "22"), name = "Censoring",
                            labels = c("No censoring", "Uniform censoring")) +
      scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      labs(y = "Predictions with N* > N") +
      guides(colour = guide_legend(nrow = 1, order = 1, title.position = "left"),
             linetype = guide_legend(nrow = 1, order = 2, title.position = "left",
                                     override.aes = list(colour = PAL[["navy"]]))) +
      theme_art(11) + theme(legend.box = "vertical", legend.box.just = "left", legend.spacing.y = unit(0.1, "cm"))
    ggsave(file.path(outdir, sprintf("sim_main_admin_%s.png", nm)), p, width = 8.5, height = 5.5, dpi = 300, bg = "white")
    out[[nm]] <- p
  }
  invisible(out)
}

# Supplement / alternative main figure: the same percentage against the predicted RISK RELATIVE TO THE EVENT RATE,
# (1 - S) / e on a log axis. If the survival at which N* > N starts is proportional to e (1 - S* = b e), the lines
# of the different event rates fall on top of each other. Colour = event rate, line type = censoring.
plot_relative <- function( res, outdir = "results/", trim = TRUE, width = 0.1, min_n = 20, tol = 0.05 ){
  dir.create( outdir, showWarnings = FALSE, recursive = TRUE )
  r <- main_data( res, trim )
  r$xbin <- round( log10((1 - r$S_true) / r$e) / width ) * width
  r <- r[is.finite(r$xbin) & r$xbin <= 2, ]
  b <- aggregate( r["P_gt_N"], r[c("el", "adm_lab", "xbin")], FUN = mean )
  b$n <- aggregate( r["P_gt_N"], r[c("el", "adm_lab", "xbin")], FUN = length )$P_gt_N
  b <- b[b$n >= min_n, ]
  br <- c( 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100 )
  p <- ggplot(b, aes(xbin, P_gt_N, colour = el, linetype = adm_lab)) +
    geom_hline(yintercept = tol, linetype = 2, colour = "grey40") + geom_line(linewidth = 0.8) +
    scale_x_continuous(breaks = log10(br), labels = br, name = "Predicted risk relative to the event rate, (1 \u2212 S) / e") +
    scale_colour_manual(values = ramp_pal(nlevels(r$el)), name = "Event rate, 1 \u2212 S(\u03c4)") +
    scale_linetype_manual(values = c("solid", "22"), name = "Censoring",
                          labels = c("No censoring", "Uniform censoring")) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(y = "Predictions with N* > N") +
    guides(colour = guide_legend(nrow = 1, order = 1, title.position = "left"),
           linetype = guide_legend(nrow = 1, order = 2, title.position = "left",
                                   override.aes = list(colour = PAL[["navy"]]))) +
    theme_art(11) + theme(legend.box = "vertical", legend.box.just = "left", legend.spacing.y = unit(0.1, "cm"))
  ggsave(file.path(outdir, "sim_main_relative.png"), p, width = 8.5, height = 5.5, dpi = 300, bg = "white")
  p
}

admin_check <- function( res, trim = TRUE, width = 0.4, min_n = 20 ){
  r <- main_data( res, trim )
  b <- bin_s( r, c("P_gt_N", "coverage"), c("el", "admin"), S_SCALE$logit, width = width, min_n = min_n )
  w <- merge( b[!b$admin, ], b[b$admin, ], by = c("el", "xbin"), suffixes = c("_common", "_staggered") )
  w$dP <- w$P_gt_N_staggered - w$P_gt_N_common; w$dCov <- w$coverage_staggered - w$coverage_common
  list( bins = w, max_abs_diff = aggregate(cbind(dP = abs(dP), dCov = abs(dCov)) ~ el, data = w, FUN = max) )
}

# Share of ALL predictions (whole population) with N* > N, against the fraction of events observed; lines = C,
# panels = overall event rate. This is where C matters: it changes how many patients lie in the high-S range.
# (Profiles are equally spaced quantiles of the linear predictor, 1st-99th percentile, so every profile
# stands for the same share of the population.)
plot_share <- function( res, outdir = "results/", trim = TRUE, admin = c("common", "staggered", "pooled") ){
  admin <- match.arg(admin)
  dir.create( outdir, showWarnings = FALSE, recursive = TRUE )
  r <- main_data( res, trim )
  if( admin == "common" ) r <- r[!r$admin, ] else if( admin == "staggered" ) r <- r[r$admin, ]
  g <- aggregate( P_gt_N ~ Cl + el + frac, data = r, FUN = mean )
  p <- ggplot(g, aes(100 * frac, P_gt_N, colour = Cl)) + geom_line(linewidth = 0.8) + geom_point(size = 1.4) +
    facet_wrap(~el, nrow = 1) +
    scale_colour_manual(values = ramp_pal(nlevels(g$Cl)), name = NULL) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(x = "Percentage of events observed", y = "Predictions with N* > N") + theme_art(11)
  ggsave(file.path(outdir, "sim_main_share.png"), p, width = 12, height = 4, dpi = 300, bg = "white")
  p
}

# Functional form of the threshold. S* = true survival at which 5% of the predictions have N* > N (upper
# crossing, S_cross_main) as a function of the overall event proportion e. Candidate: the RISK at the
# threshold, 1 - S*, is a power of e,   1 - S* = a * e^b   (b = 1: proportional to e, i.e. the threshold is a
# fixed multiple of the population's overall risk). Fitted by least squares on the log-log scale to the
# crossings of every (C, admin, e) combination; `b1` is the fit with b fixed at 1 (1 - S* = c * e).
fit_threshold <- function( res, tol = 0.05, trim = TRUE ){
  r <- main_data( res, trim )
  r$grp <- paste( r$C, r$admin, r$e )
  cr <- S_cross_main( r, tol, group = "grp" )
  cr$C <- as.numeric(sub(" .*", "", cr$grp)); cr$admin <- grepl( " TRUE ", cr$grp )
  cr$e <- as.numeric(sub(".* ", "", cr$grp))
  cr <- cr[is.finite(cr$S_up), ]; cr$risk <- 1 - cr$S_up
  fit <- lm( log(risk) ~ log(e), data = cr )
  b1  <- exp(mean(log(cr$risk) - log(cr$e)))
  pooled <- S_cross_main( main_data(res, trim), tol ); pooled$e <- as.numeric(as.character(pooled$el))
  list( crossings = cr, pooled = pooled, fit = fit,
        a = unname(exp(coef(fit)[1])), b = unname(coef(fit)[2]), b_ci = confint(fit)[2, ],
        resid_sd_log = sigma(fit), c_b1 = b1,
        predict_S = function( e ) 1 - exp(coef(fit)[1]) * e^coef(fit)[2] )
}

plot_threshold <- function( ft, outdir = "results/" ){
  dir.create( outdir, showWarnings = FALSE, recursive = TRUE )
  e <- exp(seq(log(0.008), log(1), length.out = 100))
  cr <- ft$crossings
  p <- ggplot(cr, aes(e, 1 - S_up)) +
    geom_line(data = data.frame(e = e, y = ft$c_b1 * e), aes(e, y), colour = "grey60", linetype = 2, inherit.aes = FALSE) +
    geom_line(data = data.frame(e = e, y = ft$a * e^ft$b), aes(e, y), colour = PAL[["navy"]], inherit.aes = FALSE) +
    geom_point(aes(shape = factor(C)), colour = PAL[["teal"]], alpha = 0.8, size = 2) +
    scale_x_log10(breaks = c(0.01, 0.05, 0.1, 0.3, 0.6, 0.9), labels = scales::percent) +
    scale_y_log10(breaks = c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3),
                  labels = function( x ) sprintf("%.3f", 1 - x), name = "True survival S* above which more than 5% have N* > N") +
    labs(x = "Event rate, 1 \u2212 S(\u03c4)", shape = "Harrell's C",
         subtitle = sprintf("solid: 1 - S* = %.2f e^%.2f   dashed: 1 - S* = %.2f e", ft$a, ft$b, ft$c_b1)) +
    theme_art(11)
  ggsave(file.path(outdir, "sim_main_threshold.png"), p, width = 6.5, height = 5, dpi = 300, bg = "white")
  p
}

# Coverage of the 95% CI against the effective sample size relative to N (instead of S), pooled over C,
# censoring and time points; colour = overall event rate. x = N* / N (log scale; right of 1: N* > N).
# nstar = "true": N* of the prediction from the empirical variance over the replications, S(1-S)/var(S-hat)
# (what the variance should be; no conditioning on the estimate); "hat": median of the estimated N* over the
# replications (what a user would see; the median is the same for every replication of a profile and time point).
# x = "N/N*" gives the reciprocal.
plot_coverage_N <- function( res, outdir = "results/", trim = TRUE, nstar = c("true", "hat"),
                             x = c("N*/N", "N/N*"), dots = TRUE, tail = 0.005, h = 0.2, xlim = c(0.025, 40) ){
  nstar <- match.arg(nstar); x <- match.arg(x)
  dir.create( outdir, showWarnings = FALSE, recursive = TRUE )
  r <- main_data( res, trim )
  rat <- if( nstar == "true" ) r$Ntrue_ratio else r$Nhat_med_ratio        # N*/N
  r$lx <- if( x == "N/N*" ) -log10(rat) else log10(rat)
  r <- r[is.finite(r$lx), ]
  # smoothers and dots only within the central 1 - 2 * tail of the predictions of each event rate (no extrapolation)
  lo <- ave( r$lx, r$el, FUN = function(v) quantile(v, tail) )
  hi <- ave( r$lx, r$el, FUN = function(v) quantile(v, 1 - tail) )
  r <- r[r$lx >= lo & r$lx <= hi, ]
  mc <- 1.96 * sqrt(0.95 * 0.05 / min(r$R))
  # kernel smoother: at each grid point the weighted mean of the coverages, with Gaussian weights of
  # standard deviation h on the log10 scale of N*/N (h = 0.2: a factor 10^0.2 = 1.6 in N*/N)
  ker <- do.call( rbind, lapply(split(r, r$el, drop = TRUE), function(d) {
    g <- seq( min(d$lx), max(d$lx), length.out = 100 )
    data.frame( el = d$el[1], lx = g, coverage = sapply(g, function(z) weighted.mean(d$coverage, dnorm(d$lx, z, h))) )
  }) )
  br <- c( 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30, 100 )
  p <- ggplot(r, aes(lx, coverage, colour = el)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.95 - mc, ymax = 0.95 + mc, fill = "grey88") +
    geom_hline(yintercept = 0.95, linetype = 2, colour = "grey30") +
    geom_vline(xintercept = 0, linetype = 3, colour = "grey50") +
    { if( dots ) geom_point(size = 0.5, alpha = 0.08) } +                  # one dot per profile and time point
    geom_line(data = ker, linewidth = 0.9) +
    scale_x_continuous(breaks = log10(br), labels = br, name = if( x == "N/N*" ) "N / N*" else "N* / N",
                       guide = guide_axis(check.overlap = TRUE)) +
    scale_colour_manual(values = ramp_pal(nlevels(r$el)), name = "Event rate, 1 \u2212 S(\u03c4)") +
    coord_cartesian(xlim = log10(xlim), ylim = c(0.9, 1)) +                # fixed x-axis (N*/N)
    labs(y = "Coverage of the 95% log-log CI") +
    guides(colour = guide_legend(nrow = 1, override.aes = list(alpha = 1))) + theme_art(11)
  ggsave(file.path(outdir, sprintf("sim_main_coverage_N_%s.png", nstar)), p, width = 8.5, height = 5.5, dpi = 300, bg = "white")
  p
}

# Range of true survival in which P(N*-hat > N) is at most `tol` (pooled over fractions of events observed),
# per scenario: S_lower / S_upper = lower edge of the first S bin (logit scale, scanning away from S = 0.5)
# in which P exceeds tol; NA = P never exceeds tol within the range of S that occurs. min_cov = lowest
# coverage of the 95% CI in the bins within [S_lower, S_upper].
S_limits <- function( res, tol = 0.05, width = 0.25, trim = TRUE, min_frac = 0, min_n = 5 ){
  r <- main_data( res, trim ); r <- r[r$frac >= min_frac & r$S_true > 1e-6 & r$S_true < 1 - 1e-9, ]
  r$xbin <- floor(qlogis(r$S_true) / width) * width + width / 2
  key <- c( "C", "N", "e", "admin" )
  g <- aggregate( r[c("P_gt_N", "coverage")], r[c(key, "xbin")], FUN = mean )
  n <- aggregate( r["profile"], r[c(key, "xbin")], FUN = length ); g <- g[n$profile >= min_n, ]
  g$k <- do.call( paste, g[key] )
  do.call( rbind, lapply(split(g, g$k), function(x) {
    up <- x[x$xbin >= 0, ]; up <- up[order(up$xbin), ]
    lo <- x[x$xbin < 0, ];  lo <- lo[order(-lo$xbin), ]
    iu <- which(up$P_gt_N > tol)[1]; il <- which(lo$P_gt_N > tol)[1]
    xu <- if( is.na(iu) ) Inf else up$xbin[iu]; xl <- if( is.na(il) ) -Inf else lo$xbin[il]
    ok <- x$xbin > xl & x$xbin < xu
    data.frame(x[1, key], S_lower = if( is.na(il) ) NA_real_ else plogis(lo$xbin[il] + width / 2),
               S_upper = if( is.na(iu) ) NA_real_ else plogis(up$xbin[iu] - width / 2),
               min_cov = if( any(ok) ) min(x$coverage[ok]) else NA_real_, row.names = NULL)
  }) )
}

# Censoring study (design "censoring"): composition of the censoring per scenario (population level)
cens_summary <- function( scen = make_scen("censoring") ){
  do.call( rbind, lapply(seq_len(nrow(scen)), function(i) {
    s <- scen[i, ]
    sc <- setup_scenario(s$C, s$N, s$e, s$shape, s$k, s$noise, admin = isTRUE(s$admin), rc = s$rc)
    data.frame(s[c("C", "N", "e", "admin", "rc")], t(round(sc$cens, 3)), row.names = NULL)
  }) )
}

cens_lab <- function( admin, rc ) factor( paste0(ifelse(admin, "uniform censoring", "no censoring"), ", rc = ", rc),
                                          levels = c(outer(c("no censoring", "uniform censoring"), c(0, 1, 5), function(a, r) paste0(a, ", rc = ", r))) )

# P(N*-hat > N), variance ratio and coverage by censoring setting, e and S band (event fractions >= f_min_use)
cens_table <- function( res, f_min_use = 0.25, trim = TRUE ){
  r <- trim_core( res, trim ); r <- r[r$frac >= f_min_use & r$S_true > 1e-6, ]; r$S_band <- band_of(r$S_true)
  aggregate( cbind(P_gt_N, var_ratio, coverage, Nhat_over_Ntrue) ~ e + admin + rc + S_band, data = r, FUN = mean )
}

# Returns list(p1, p2, p3); p1/p2 per x-axis scale. p1: P(N*-hat > N) vs true S, rows = e, columns = event
# fraction, colour = censoring setting. p2: coverage vs S. p3: P vs event fraction within S bands.
plot_cens <- function( res, outdir = "results/", scales = c("logit", "linear"), trim = TRUE,
                       fracs_show = c(0.10, 0.25, 0.40, 0.60, 0.90), C_show = 0.7 ){
  dir.create( outdir, showWarnings = FALSE, recursive = TRUE )
  res <- trim_core( res, trim ); res <- res[res$C == C_show, ]
  res$cens <- cens_lab( res$admin, res$rc ); res$el <- factor( res$e, labels = format(sort(unique(res$e)), drop0trailing = TRUE) )
  mc <- 1.96 * sqrt(0.95 * 0.05 / min(res$R))
  out <- list( p1 = list(), p2 = list() )
  for( nm in scales ){
    sc <- S_SCALE[[nm]]; cl <- coord_cartesian(xlim = sc$xlim)
    b1 <- bin_s( res[res$frac %in% fracs_show, ], "P_gt_N", c("cens", "el", "frac"), sc ); b1$f <- frac_lab(b1$frac)
    p1 <- ggplot(b1, aes(xbin, P_gt_N, colour = cens, linetype = cens)) +
      geom_hline(yintercept = 0.05, linetype = 2, colour = "grey40") + geom_line(linewidth = 0.7) +
      facet_grid(el ~ f) + cl + scale_x_S(sc) + scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      scale_linetype_manual(values = rep(c("solid", "22"), 3)) +
      scale_colour_manual(values = rep(ramp_pal(3), each = 2)) +
      labs(y = "Probability that N*-hat > N", colour = NULL, linetype = NULL,
           title = sprintf("N*-hat > N under registry-style censoring (C = %.1f)", C_show)) +
      guides(colour = guide_legend(nrow = 2)) +
      theme_art(11)
    ggsave(file.path(outdir, sprintf("sim_cens_PgtN_%s.png", nm)), p1, width = 14, height = 8, dpi = 200, bg = "white")
    out$p1[[nm]] <- p1
    b2 <- bin_s( res, "coverage", c("cens", "el"), sc, width = sc$width * 0.8 )
    p2 <- ggplot(b2, aes(xbin, coverage, colour = cens, linetype = cens)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.95 - mc, ymax = 0.95 + mc, fill = "grey85") +
      geom_hline(yintercept = 0.95, linetype = 2) + geom_line(linewidth = 0.8) + facet_wrap(~el) +
      scale_linetype_manual(values = rep(c("solid", "22"), 3)) +
      scale_colour_manual(values = rep(ramp_pal(3), each = 2)) +
      coord_cartesian(xlim = sc$xlim, ylim = c(0.88, 1)) + scale_x_S(sc) +
      labs(y = "Coverage of the 95% log-log CI", colour = NULL, linetype = NULL) +
      theme_art(11)
    ggsave(file.path(outdir, sprintf("sim_cens_coverage_%s.png", nm)), p2, width = 12, height = 5.5, dpi = 200, bg = "white")
    out$p2[[nm]] <- p2
  }
  res$S_band <- band_of(res$S_true)
  g <- aggregate( P_gt_N ~ el + cens + frac + S_band, data = res, FUN = mean )
  n <- aggregate( profile ~ el + cens + frac + S_band, data = res, FUN = length ); g <- g[n$profile >= 5, ]
  out$p3 <- ggplot(g, aes(frac, P_gt_N, colour = cens, linetype = cens)) +
    geom_hline(yintercept = 0.05, linetype = 2, colour = "grey40") + geom_line(linewidth = 0.7) +
    facet_grid(el ~ S_band) + scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    scale_linetype_manual(values = rep(c("solid", "22"), 3)) +
    scale_colour_manual(values = rep(ramp_pal(3), each = 2)) +
    scale_x_continuous(breaks = c(0.1, 0.5, 0.9)) +
    labs(x = "Fraction of observed events", y = "Probability that N*-hat > N", colour = NULL, linetype = NULL) +
    theme_art(10)
  ggsave(file.path(outdir, "sim_cens_vs_fraction.png"), out$p3, width = 15, height = 8, dpi = 200, bg = "white")
  invisible(out)
}

# Usage:
#   1. Pilot (MC error and number of replications):
#        pm <- pilot_mc(R_pilot = 300); pm$R_recommended
#        R_REP <- pm$R_recommended
#   2. Full run:
#        res <- run_all(); saveRDS(res, "results/sim_guardrails.rds")
#   3. Results:
#        fmin_by_band(res); slow_table(res); bias_table(res)
#        factor_effect(res, "k"); factor_effect(res, "noise")
#        pl <- plot_sim(res)     # all four x-axis scales; pl$p1$logit etc.
# Shape-only study:
#        resS <- run_all(make_scen("shape_only")); saveRDS(resS, "results/sim_guardrails_shape.rds")
#        factor_effect(resS, "shape"); shape_table(resS); plot_shape(resS)
# Event-rate study (non-factorial):
#        SCEN_E <- make_scen("event_rate")          # N = 500 / e; 12 scenarios
#        pe <- pilot_mc(rows = c(1, 3, 9), scen = SCEN_E)   # optional: e.g. e = 0.05, 0.1, 0.9 at C = 0.7
#        resE <- run_all(SCEN_E); saveRDS(resE, "results/sim_guardrails_erate.rds")
#        e_spread(resE); fmin_by_band(resE); factor_effect(resE, "noise", e = c(0.1, 0.9))
#        pe_plots <- plot_erate(resE)
# MAIN simulation (article): C x overall event rate (e = 1 - S(tau = 1)) x administrative censoring, 12 scenarios at C = 0.8 (N = 500 / e)
#        SCEN_M <- make_scen("factorial")                 # fewer: make_scen("factorial", C = c(0.6, 0.8), e = c(0.05, 0.3, 0.9))
#        resM <- run_all(SCEN_M); saveRDS(resM, "results/sim_main.rds")
#        pm <- plot_main(resM)                            # pm$p1$logit (P(N* > N) vs S), pm$p2$logit (coverage)
#        S_limits(resM)                                   # range of S with P(N* > N) <= 5%, per scenario
#        pilot: pilot_mc(rows = c(2, 6, 12), scen = SCEN_M)   # check the rows first: SCEN_M[rows, ]
# Censoring study (registry-style: staggered entry and/or random loss to follow-up), 18 scenarios:
#        SCEN_C <- make_scen("censoring")                 # C = 0.7; e = 0.1, 0.3, 0.6; admin off/on; rc = 0, 1, 5
#        cens_summary(SCEN_C)                             # events / admin-censored / LTFU proportions: check first!
#        resC <- run_all(SCEN_C); saveRDS(resC, "results/sim_guardrails_cens.rds")
#        cens_table(resC); fmin_by_band(resC)             # fmin_by_band keeps admin and rc as keys
#        pc <- plot_cens(resC)                            # pc$p1$logit, pc$p2$logit, pc$p3
# Note: with censoring tied to the time scale the Weibull shape is no longer irrelevant in principle
# (rank invariance only holds for the Cox fit, not for the censoring distribution); shape is fixed at 1 here.
# Fewer scenarios:  SCEN <- SCEN[SCEN$e == 0.9, ]