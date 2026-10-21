####################### EFFECTIVE N FOR COX PREDICTIONS IN OPEN DATASETS #####################

# Effective sample size (N*) of Cox model predictions in several open datasets.
#
# Self-contained: only needs the packages below (no helpers from compile.R).
# N* is computed from the standard error that survival::survfit.coxph returns
# for each patient (delta-method variance, including coefficient uncertainty):
#     N* = S (1 - S) / SE(S)^2
# This is the same quantity as in Figure 3 of the manuscript.
#
# Run validate_colon_fig3() first and compare with the medians in Figure 3.
#
# Usage (from repository root):
#   source("code/ess_realworld.R")
#   res <- run_all()                      # long data.frame, one row per patient x horizon x spec
#   saveRDS(res, "results/ess_realworld.rds")
#   out <- summarise_ess(res); print(out$by_dataset); print(out$rule_table)
#   make_plots(res, "results/")         # ess_vs_S.png (scatter) and ess_by_eventfraction.png (box plots)
#   plot_fig4(res, type = "scatter")      # or "box": candidate manuscript Figure 4, panels ordered by Harrell's C
#   check_predictors()                    # look at the candidate predictors first
#   drv <- surplus_drivers(res); print(driver_correlations(drv))
#   plot_surplus_drivers(drv, x = "cstat_h", y = "pct_above_2N", outfile = "results/ess_drivers.png")

library(survival)
library(ggplot2)

####################### SETTINGS #####################

# Horizons are defined by the fraction of ALL events observed at the horizon
# (dataset-independent, so S and event fraction vary naturally between datasets).
FRACS <- c( 0.01, 0.10, 0.25, 0.40, 0.60, 0.75, 0.90 )

# Model specification: one pre-specified candidate set per dataset, reduced by backward
# elimination on AIC (Efron ties, the coxph default). The delta-method variance ignores the
# selection step, so N* is conditional on the selected model.
# (Efron and Breslow gave practically the same N*; to re-add other variants, extend this table,
#  e.g. spec = "full_efron", ties = "efron", selection = "none".)
MAIN_SPEC <- "backward_efron"
SPECS <- data.frame(
  spec      = "backward_efron",
  ties      = "efron",
  selection = "aic",
  stringsAsFactors = FALSE
)

####################### 1B. COLOURS AND LABELS FOR THE FIGURES #####################
# Same hex codes as the article figures (article_figures.R); FRAC_PAL = sequential ramp for the
# fraction of events observed (gold -> red -> teal -> navy, decreasing lightness).
PAL <- c( navy = "#37293F", teal = "#2E7691", purple = "#9f84af", red = "#C2666B", gold = "#c6aa2c" )
FRAC_PAL <- colorRampPalette( unname(PAL[c("gold", "red", "teal", "navy")]), space = "Lab" )(length(FRACS))
DS_LABELS <- c( colon = "Colon cancer", pbc = "Primary biliary cholangitis",
                veteran = "Veterans' lung cancer", lung = "NCCTG lung cancer",
                rotterdam = "Rotterdam breast cancer", gbsg = "GBSG breast cancer (RFS)",
                flchain = "Serum free light chain", mgus2 = "MGUS" )

####################### 2. DATASETS #####################
# Each builder returns list(data = data.frame with `time`, `status` (0/1) and the
# raw covariates, terms = character vector of model terms (may contain log() etc.)).
# Time units differ per dataset, which does not matter: horizons are event-based.

ds_colon <- function( ){
  d <- survival::colon
  d <- d[d$etype == 2, ]                      # death, as in the manuscript
  d$rx <- factor(d$rx); d$differ <- factor(d$differ); d$extent <- factor(d$extent)
  list( data = d,
        terms = c("rx", "sex", "age", "obstruct", "perfor", "adhere", "nodes",
                  "differ", "extent", "surg") )
}

ds_pbc <- function( ){
  d <- survival::pbc
  d <- d[!is.na(d$trt), ]                     # the 312 randomised patients
  d$status <- as.integer(d$status == 2)       # death; transplant treated as censored
  d$edema <- factor(d$edema)
  list( data = d,
        terms = c("age", "sex", "edema", "ascites", "hepato", "spiders", "log(bili)",
                  "albumin", "log(protime)", "log(copper)", "log(ast)", "stage") )
}

ds_veteran <- function( ){
  d <- survival::veteran
  d$trt <- factor(d$trt); d$prior <- factor(d$prior)
  list( data = d, terms = c("trt", "celltype", "karno", "log(diagtime + 1)", "age", "prior") )
}

ds_lung <- function( ){
  d <- survival::lung
  d$status <- d$status - 1                    # 1/2 -> 0/1
  d$sex <- factor(d$sex)
  list( data = d, terms = c("age", "sex", "ph.ecog", "ph.karno", "pat.karno", "wt.loss") )
}

ds_rotterdam <- function( ){
  d <- survival::rotterdam
  d$time   <- d$dtime; d$status <- d$death            # overall survival (recurrence-free survival: rtime / recur)
  d$meno <- factor(d$meno); d$grade <- factor(d$grade)
  list( data = d,
        terms = c("age", "meno", "size", "grade", "log(nodes + 1)", "log(pgr + 1)",
                  "log(er + 1)", "hormon", "chemo") )
}

ds_gbsg <- function( ){
  d <- survival::gbsg
  d$time <- d$rfstime                         # recurrence-free survival: the only outcome in this dataset (no OS); `status` already 0/1
  d$meno <- factor(d$meno); d$grade <- factor(d$grade)
  list( data = d,
        terms = c("age", "meno", "size", "grade", "log(nodes + 1)", "log(pgr + 1)",
                  "log(er + 1)", "hormon") )
}

ds_flchain <- function( ){
  d <- survival::flchain
  d$time <- d$futime; d$status <- d$death
  d$sex <- factor(d$sex)
  list( data = d,
        terms = c("age", "sex", "log(kappa)", "log(lambda)", "log(creatinine)", "mgus") )
}

ds_mgus2 <- function( ){
  d <- survival::mgus2
  d$time <- d$futime; d$status <- d$death
  d$sex <- factor(d$sex)
  list( data = d, terms = c("age", "sex", "hgb", "log(creat)", "mspike") )
}

DATASETS <- list( colon = ds_colon, pbc = ds_pbc, veteran = ds_veteran, lung = ds_lung,
                  rotterdam = ds_rotterdam, gbsg = ds_gbsg, flchain = ds_flchain,
                  mgus2 = ds_mgus2 )

####################### 3. MODEL FITTING HELPERS #####################

fit_cox <- function( terms, data, ties ){
  f <- as.formula(paste("Surv(time, status) ~", paste(terms, collapse = " + ")))
  # x = TRUE, y = TRUE so survfit() does not need to rebuild the model frame
  coxph( f, data = data, ties = ties, x = TRUE, y = TRUE )
}

backward_aic <- function( terms, data, ties ){
  cur <- terms
  aic <- AIC(fit_cox(cur, data, ties))
  repeat {
    if( length(cur) <= 1 ) break
    cand <- vapply( seq_along(cur), function(i) AIC(fit_cox(cur[-i], data, ties)), numeric(1) )
    best <- which.min(cand)
    if( cand[best] < aic ){ aic <- cand[best]; cur <- cur[-best] } else break
  }
  cur
}

# Horizons: time at which a fraction f of all events has occurred
make_horizons <- function( time, status, fracs = FRACS ){
  et <- sort(time[status == 1])
  d_tot <- length(et)
  t_f <- et[pmax(1, ceiling(fracs * d_tot))]
  data.frame( f_target = fracs, t = t_f,
              d_t = vapply(t_f, function(x) sum(et <= x), numeric(1)),
              event_frac = vapply(t_f, function(x) mean(et <= x), numeric(1)),
              d_total = d_tot )
}

# S(t) and SE(S(t)) for every row of `newdata` at times `tt`.
# Uses the survfit object directly: sf$std.err is the SE of the cumulative hazard
# (= SE of log S), so SE(S) = S * sf$std.err.
predict_S_se <- function( fit, newdata, tt ){
  sf  <- survfit( fit, newdata = newdata )
  idx <- findInterval( tt, sf$time )
  if( any(idx < 1) ) stop("horizon before first event time")
  Sm  <- as.matrix(sf$surv)[idx, , drop = FALSE]
  SEm <- as.matrix(sf$std.err)[idx, , drop = FALSE]
  list( S = Sm, SE_S = Sm * SEm, sf = sf )
}

# One-off check that summary.survfit gives the same S-scale SE as the line above.
check_scale <- function( fit, newdata, tt, label = "" ){
  nd  <- newdata[1:min(3, nrow(newdata)), , drop = FALSE]
  a   <- predict_S_se( fit, nd, tt )
  ss  <- summary( a$sf, times = tt )
  se2 <- matrix( ss$std.err, nrow = length(tt) )
  rel <- max( abs(se2 - a$SE_S) / pmax(a$SE_S, 1e-12), na.rm = TRUE )
  if( is.finite(rel) && rel > 1e-6 )
    warning( sprintf("[%s] SE from summary.survfit differs from S*sf$std.err (max rel. diff %.3g). ",
                     label, rel), "Check which scale std.err is on before trusting N*." )
  invisible(rel)
}


####################### 4B. PREDICTOR SANITY CHECK #####################
# Prints, per dataset, the numeric model terms with the most extreme standardised values and
# skewness, so that clearly problematic predictors (heavy outliers, IDs, ...) are noticed
# BEFORE interpreting N*. Uses the full candidate set.
check_predictors <- function( datasets = DATASETS ){
  out <- lapply( names(datasets), function(nm) {
    b <- tryCatch(datasets[[nm]](), error = function( e ) NULL)
    if( is.null(b) ) return(NULL)
    d <- b$data
    X <- tryCatch(model.matrix(reformulate(b$terms), data = d)[, -1, drop = FALSE],
                  error = function( e ) NULL)
    if( is.null(X) ) return(NULL)
    z <- scale(X)
    data.frame(dataset = nm, term = colnames(X),
               n_unique = apply(X, 2, function(v) length(unique(v))),
               max_abs_z = round(apply(abs(z), 2, max, na.rm = TRUE), 1),
               skew = round(apply(X, 2, function(v) { v <- v[is.finite(v)]
               mean((v - mean(v))^3) / sd(v)^3 }), 2),
               row.names = NULL)
  } )
  do.call( rbind, out )
}

####################### 4. PER-DATASET ANALYSIS #####################

ess_dataset <- function( name, builder, specs = SPECS, fracs = FRACS, chunk = 1000 ){
  b <- builder()
  d <- b$data; terms <- b$terms
  n_raw <- nrow(d)
  vars <- unique(c("time", "status", all.vars(reformulate(terms))))
  stopifnot(all(vars %in% names(d)))
  d <- d[complete.cases(d[, vars]), vars]
  d <- d[d$time > 0, ]
  rownames(d) <- NULL
  d$id <- seq_len(nrow(d))
  n <- nrow(d)
  hz <- make_horizons( d$time, d$status, fracs )
  tt <- sort(unique(hz$t))
  message(sprintf("%s: %d of %d patients complete (%.1f%% dropped), events = %d (%.0f%%), %d candidate predictors",
                  name, n, n_raw, 100 * (1 - n / n_raw), hz$d_total[1], 100 * hz$d_total[1] / n, length(terms)))
  
  out <- list()
  for( k in seq_len(nrow(specs)) ){
    sp <- specs[k, ]
    tr <- if( sp$selection == "aic" ) backward_aic( terms, d, sp$ties ) else terms
    fit <- fit_cox( tr, d, sp$ties )
    
    lp <- as.numeric(predict(fit, type = "lp"))
    rel_hr <- exp(lp) / mean(exp(lp))          # patient's HR relative to the mean HR (Supplement 12)
    maha <- tryCatch( mahalanobis(fit$x, colMeans(fit$x), cov(fit$x)),
                      error = function( e ) rep(NA_real_, n) )
    # model descriptors (constant within dataset x spec)
    cstat  <- unname(concordance(fit)$concordance)   # apparent c-statistic
    sd_lp  <- sd(lp)                                  # spread of the linear predictor (log-HR scale)
    n_par  <- length(coef(fit))
    # Harrell's C restricted to follow-up up to each horizon (concordance() default is Harrell's C)
    cstat_h <- vapply( hz$t, function(th) tryCatch(unname(concordance(fit, ymax = th)$concordance),
                                                   error = function( e ) NA_real_), numeric(1) )
    p_lowhr <- mean(rel_hr < 0.5)                     # share with HR < half the mean HR
    
    ids <- split( seq_len(n), ceiling(seq_len(n) / chunk) )
    for( ii in seq_along(ids) ){
      rows <- ids[[ii]]
      if( ii == 1 ) check_scale( fit, d[rows, , drop = FALSE], tt, paste(name, sp$spec) )
      p <- predict_S_se( fit, d[rows, , drop = FALSE], tt )
      for( h in seq_len(nrow(hz)) ){
        j <- match( hz$t[h], tt )
        S <- p$S[j, ]; se <- p$SE_S[j, ]
        Ns <- ifelse( se > 0 & S > 0 & S < 1, S * (1 - S) / se^2, NA_real_ )
        out[[length(out) + 1]] <- data.frame(
          dataset = name, spec = sp$spec, n_terms = length(tr), n_total = n,
          d_total = hz$d_total[h], f_target = hz$f_target[h], t = hz$t[h],
          d_t = hz$d_t[h], event_frac = hz$event_frac[h],
          id = d$id[rows], S = S, N_star = Ns, N_ratio = Ns / n,
          E_star = Ns * (1 - S),                       # "effective events"
          rel_hr = rel_hr[rows], maha = maha[rows],
          n_raw = n_raw, cstat = cstat, cstat_h = cstat_h[h], sd_lp = sd_lp, n_par = n_par, epv = hz$d_total[h] / n_par,
          p_lowhr = p_lowhr, pct_cens = 100 * (1 - hz$d_total[h] / n),
          stringsAsFactors = FALSE )
      }
    }
    message(sprintf("   %-15s terms kept: %d/%d", sp$spec, length(tr), length(terms)))
  }
  do.call( rbind, out )
}

run_all <- function( datasets = DATASETS, ... ){
  res <- lapply( names(datasets), function(nm) {
    tryCatch(ess_dataset(nm, datasets[[nm]], ...),
             error = function( e ){ message("!! ", nm, " skipped: ", conditionMessage(e)); NULL })
  } )
  do.call( rbind, res )
}

####################### 5. SUMMARIES #####################

summarise_ess <- function( res, main_spec = MAIN_SPEC ){
  q <- function( x, p ) unname(quantile(x, p, na.rm = TRUE))
  
  by_dataset <- do.call( rbind, lapply(
    split(res, list(res$dataset, res$spec, res$f_target), drop = TRUE), function(g) {
      data.frame(dataset = g$dataset[1], spec = g$spec[1], f_target = g$f_target[1],
                 n = g$n_total[1], n_raw = g$n_raw[1], events = g$d_total[1], event_frac = g$event_frac[1],
                 S_median = median(g$S), S_min = min(g$S),
                 ratio_p05 = q(g$N_ratio, .05), ratio_median = q(g$N_ratio, .5),
                 ratio_p95 = q(g$N_ratio, .95),
                 pct_above_N = 100 * mean(g$N_ratio > 1, na.rm = TRUE),
                 Estar_over_events_median = median(g$E_star / g$d_t, na.rm = TRUE))
    }) )
  rownames(by_dataset) <- NULL
  
  # Check of the rule of thumb: S > 0.2 AND > 50% of events observed
  r <- res[res$spec == main_spec, ]
  r$S_band  <- cut( r$S, c(0, .2, .5, .8, .95, 1), include.lowest = TRUE )
  r$ef_band <- factor(ifelse(r$event_frac > 0.5, ">50% events", "<=50% events"))
  rule_table <- do.call( rbind, lapply(split(r, list(r$S_band, r$ef_band), drop = TRUE), function(g)
    data.frame(S_band = g$S_band[1], events = g$ef_band[1], n_obs = nrow(g),
               pct_above_N = 100 * mean(g$N_ratio > 1, na.rm = TRUE),
               ratio_p95 = q(g$N_ratio, .95), ratio_max = max(g$N_ratio, na.rm = TRUE))) )
  rownames(rule_table) <- NULL
  
  # Early-follow-up identity from Supplement 12: N*/N ~ 1 / (HR_i / mean HR) at the first event
  e <- res[res$spec == main_spec & res$f_target == min(res$f_target), ]
  early_check <- do.call( rbind, lapply(split(e, e$dataset), function(g)
    data.frame(dataset = g$dataset[1], cor_log = cor(log(g$N_ratio), -log(g$rel_hr), use = "complete.obs"),
               median_abs_log_diff = median(abs(log(g$N_ratio) + log(g$rel_hr)), na.rm = TRUE))) )
  rownames(early_check) <- NULL
  
  # Efron vs Breslow (only if both specifications were run)
  ties_check <- NULL
  if( all(c("full_efron", "full_breslow") %in% res$spec) ){
    a <- res[res$spec == "full_efron", c("dataset", "f_target", "id", "N_star")]
    b <- res[res$spec == "full_breslow", c("dataset", "f_target", "id", "N_star")]
    m <- merge( a, b, by = c("dataset", "f_target", "id"), suffixes = c("_efron", "_breslow") )
    ties_check <- do.call( rbind, lapply(split(m, m$dataset), function(g)
      data.frame(dataset = g$dataset[1],
                 median_rel_diff_pct = 100 * median(abs(g$N_star_breslow / g$N_star_efron - 1), na.rm = TRUE),
                 max_rel_diff_pct    = 100 * max(abs(g$N_star_breslow / g$N_star_efron - 1), na.rm = TRUE))) )
    rownames(ties_check) <- NULL
  }
  
  list( by_dataset = by_dataset, rule_table = rule_table,
        early_check = early_check, ties_check = ties_check )
}


####################### 5B. WHAT DRIVES THE SURPLUS? #####################
# One row per dataset x model x event fraction: how often N* > N and N* > 2N, together with
# model descriptors (c-statistic, spread of the linear predictor, events per parameter).
# Early identity (Supplement 12): N*/N ~ mean(HR) / HR_i, so N* > 2N when HR_i < mean(HR)/2,
# i.e. the early surplus is governed by the spread of the risk scores (sd_lp).
# At late horizons the surplus comes from patients with S close to 0 instead.
surplus_drivers <- function( res ){
  do.call( rbind, lapply(split(res, list(res$dataset, res$spec, res$f_target), drop = TRUE), function(g) {
    a2 <- g$N_ratio > 2
    data.frame(dataset = g$dataset[1], spec = g$spec[1], f_target = g$f_target[1],
               n = g$n_total[1], events = g$d_total[1], n_par = g$n_par[1], epv = g$epv[1],
               cstat = g$cstat[1], cstat_h = g$cstat_h[1], sd_lp = g$sd_lp[1], p_lowhr = 100 * g$p_lowhr[1],
               pct_above_N  = 100 * mean(g$N_ratio > 1, na.rm = TRUE),
               pct_above_2N = 100 * mean(a2, na.rm = TRUE),
               pct_S_below_0.2 = 100 * mean(g$S < 0.2, na.rm = TRUE),
               share_above2N_with_S_below_0.2 = if( any(a2, na.rm = TRUE) )
                 100 * mean(g$S[a2 & !is.na(a2)] < 0.2) else NA_real_)
  }) )
}

# x = "cstat" or "sd_lp"; y = "pct_above_2N" or "pct_above_N"
plot_surplus_drivers <- function( drv, x = c("cstat_h", "cstat", "sd_lp"), y = c("pct_above_2N", "pct_above_N"),
                                  fracs = c(0.01, 0.10, 0.25, 0.60, 0.90),
                                  outfile = NULL ){
  x <- match.arg(x); y <- match.arg(y)
  d <- drv[drv$f_target %in% fracs, ]
  d$f_lab <- factor( sprintf("%.0f%% of events", 100 * d$f_target),
                     levels = sprintf("%.0f%% of events", 100 * sort(fracs)) )
  xl <- c( cstat_h = "Harrell's C up to the horizon (apparent)", cstat = "Harrell's C, full follow-up (apparent)", sd_lp = "SD of linear predictor (log HR)" )[x]
  yl <- c( pct_above_2N = "% of predictions with N* > 2N",
           pct_above_N = "% of predictions with N* > N" )[y]
  p <- ggplot(d, aes(.data[[x]], .data[[y]], colour = dataset, shape = spec)) +
    geom_point(size = 2.5, alpha = 0.85) +
    facet_wrap(~f_lab, nrow = 1) +
    scale_y_continuous(limits = c(0, 100)) +
    labs(x = xl, y = yl, colour = "Dataset", shape = "Model") +
    theme_bw() + theme(legend.position = "bottom")
  if( !is.null(outfile) ) ggsave(outfile, p, width = 13, height = 4.5, dpi = 200)
  p
}

# Rank correlation across dataset x model combinations, per event fraction
driver_correlations <- function( drv ){
  do.call( rbind, lapply(split(drv, drv$f_target), function(g)
    data.frame(f_target = g$f_target[1], n_models = nrow(g),
               rho_cstat_above2N = suppressWarnings(cor(g$cstat_h, g$pct_above_2N, method = "spearman")),
               rho_sdlp_above2N  = suppressWarnings(cor(g$sd_lp, g$pct_above_2N, method = "spearman")),
               rho_cstat_aboveN  = suppressWarnings(cor(g$cstat_h, g$pct_above_N, method = "spearman")),
               rho_sdlp_aboveN   = suppressWarnings(cor(g$sd_lp, g$pct_above_N, method = "spearman")))) )
}

####################### 6. FIGURES #####################
# Colours: the palette of the article figures (PAL / FRAC_PAL, defined in section 1b).

theme_art <- function( base_size = 11 )
  theme_minimal(base_size = base_size) +
  theme(strip.text = element_text(hjust = 0), legend.position = "bottom",
        panel.grid.minor = element_blank())

# Main-specification rows with plotting helpers: panel label = dataset, N and apparent Harrell's C
# (full follow-up); panels ordered by C (lowest first) so that the link between prognostic strength and
# N* is visible at a glance. Set order_by_C = FALSE for the order of DATASETS.
prep_plot_data <- function( res, order_by_C = TRUE ){
  r <- res[res$spec == MAIN_SPEC & is.finite(res$N_star) & res$N_star > 0, ]
  r$f_lab <- factor( r$f_target, levels = FRACS, labels = format(FRACS, drop0trailing = TRUE) )
  lab <- unname(DS_LABELS[r$dataset]); lab[is.na(lab)] <- r$dataset[is.na(lab)]
  r$ds_lab <- sprintf( "%s\nN = %s, Harrell's C = %.2f", lab, format(r$n_total, big.mark = ","), r$cstat )
  key <- if( order_by_C ) r$cstat else match( r$dataset, names(DATASETS) )
  r$ds_lab <- factor( r$ds_lab, levels = unique(r$ds_lab[order(key)]) )
  r$above_2N <- r$N_star > 2 * r$n_total
  r$above_N  <- r$N_star > r$n_total
  r
}

# Survival probability above which N* > N at the high end of S, per dataset: the predictions are binned on the
# logit scale of S (width 0.25, S > 0.5, bins with >= min_n predictions); going down from the highest bin, the
# crossing is where the bin median of N*/N first falls below 1 (linear interpolation on logit S between bin
# centres). q = 0.5 uses the median; a lower q (e.g. 0.1) gives the S above which N* > N is the rule, not just
# the typical case. NA when N* <= N throughout (no line is drawn). Uses all predictions, not the plotted sample.
S_cross <- function( r, q = 0.5, width = 0.25, min_n = 20 ){
  do.call( rbind, lapply(split(r, r$ds_lab), function(g) {
    g <- g[g$S > 0.5 & g$S < 1, ]
    g$b <- floor(qlogis(g$S) / width)
    b <- aggregate(N_ratio ~ b, g, function(x) c(m = unname(quantile(x, q)), n = length(x)))
    b <- data.frame(b = b$b, m = b$N_ratio[, "m"], n = b$N_ratio[, "n"])
    b <- b[b$n >= min_n, ]; b <- b[order(-b$b), ]
    x <- NA_real_
    if( nrow(b) > 1 && b$m[1] > 1 ){
      i <- which(b$m < 1)[1]
      if( !is.na(i) ){
        l1 <- (b$b[i - 1] + 0.5) * width; l0 <- (b$b[i] + 0.5) * width      # bin centres, logit S
        x <- plogis(l0 + (1 - b$m[i]) / (b$m[i - 1] - b$m[i]) * (l1 - l0))
      }
    }
    data.frame(ds_lab = g$ds_lab[1], S_cross = x)
  }) )
}

# N*/N against predicted survival S, colour = fraction of events observed; points with N*/N > 2 not drawn
# Dashed vertical line + label: S above which N* > N (see S_cross; cross = FALSE removes it)
plot_scatter <- function( r, max_pts = 3000, cross = TRUE, q = 0.5 ){
  set.seed(1)
  shown <- r[r$N_ratio <= 2, ]
  rs <- do.call( rbind, lapply(split(shown, list(shown$dataset, shown$f_target), drop = TRUE),
                               function(g) g[sample(nrow(g), min(nrow(g), max_pts)), ]) )
  ggplot(rs, aes(S, N_ratio, colour = f_lab)) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey30") +
    geom_point(alpha = 0.35, size = 1, stroke = 0) +
    { if( cross ){ sc <- S_cross( r, q ); sc <- sc[is.finite(sc$S_cross), ]
    list( geom_vline(data = sc, aes(xintercept = S_cross), linetype = 2, colour = "grey30"),
          geom_text(data = sc, aes(x = S_cross, y = 1.9, label = sprintf("S = %.2f", S_cross)),
                    hjust = 1.1, size = 3.2, colour = "grey20", inherit.aes = FALSE) ) } } +
    scale_x_continuous(limits = c(0, 1)) +
    coord_cartesian(ylim = c(0, 2)) +
    scale_colour_manual(values = setNames(FRAC_PAL, levels(r$f_lab)), name = "Fraction of events observed, f") +
    guides(colour = guide_legend(nrow = 1, override.aes = list(alpha = 1, size = 3))) +
    facet_wrap(~ds_lab, ncol = if( nlevels(r$ds_lab) > 8 ) 3 else 4) +
    labs(x = "Predicted survival S(t)", y = "N* / N") +
    theme_art()
}

# Candidate manuscript Figure 4. type = "scatter": N*/N against predicted survival, coloured by the fraction
# of events observed (same layout as Figure 3); type = "box": N*/N by fraction of events observed.
# Dashed line: N* = N. Predictions with N*/N > ymax are outside the panels (their share is printed).
plot_fig4 <- function( res, outfile = NULL, type = c("scatter", "box"), ymax = 2 ){
  type <- match.arg(type)
  r <- prep_plot_data(res)
  message(sprintf("%.1f%% of the predictions have N*/N > %g (outside the plotted range)",
                  100 * mean(r$N_star / r$n_total > ymax), ymax))
  p <- if( type == "scatter" ) plot_scatter(r) else
    ggplot(r, aes(f_lab, N_star / n_total, fill = f_lab)) +
    geom_hline(yintercept = 1, linetype = 2, colour = "grey30") +
    geom_boxplot(outlier.shape = NA, linewidth = 0.35, colour = PAL[["navy"]], alpha = 0.9) +
    coord_cartesian(ylim = c(0, ymax)) +
    scale_fill_manual(values = setNames(FRAC_PAL, levels(r$f_lab)), guide = "none") +
    facet_wrap(~ds_lab, ncol = if( nlevels(r$ds_lab) > 8 ) 3 else 4) +
    labs(x = "Fraction of events observed at the horizon, f", y = "N* / N") +
    theme_art()
  if( is.null(outfile) ) outfile <- sprintf( "results/fig4_realworld_%s.png", type )
  dir.create( dirname(outfile), showWarnings = FALSE, recursive = TRUE )
  ggsave(outfile, p, width = 11, height = 6.8, dpi = 300, bg = "white")
  p
}

make_plots <- function( res, outdir = "results/", max_pts = 3000 ){
  dir.create( outdir, showWarnings = FALSE, recursive = TRUE )
  r <- prep_plot_data(res)
  
  ##### Figure A: N*/N against S per dataset (0 to 2), dashed line at N* = N #####
  pA <- plot_scatter( r, max_pts )
  ggsave(file.path(outdir, "ess_vs_S.png"), pA, width = 11, height = 6.8, dpi = 300, bg = "white")
  
  ##### Figure B: how many predictions are above N (and above 2N) #####
  tab <- do.call( rbind, lapply(split(r, list(r$dataset, r$f_target), drop = TRUE), function(g)
    data.frame(dataset = g$dataset[1], ds_lab = g$ds_lab[1], f_lab = g$f_lab[1], f_target = g$f_target[1],
               n = nrow(g), pct_above_N = 100 * mean(g$above_N), pct_above_2N = 100 * mean(g$above_2N))) )
  tab$ds_lab <- factor( tab$ds_lab, levels = rev(levels(r$ds_lab)) )
  heat <- function( var, col, ttl ){
    ggplot(tab, aes(f_lab, ds_lab, fill = .data[[var]])) +
      geom_tile(colour = "white", linewidth = 0.8) +
      geom_text(aes(label = sprintf("%.0f%%", .data[[var]]),
                    colour = .data[[var]] > 60), size = 3.2, show.legend = FALSE) +
      scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = PAL[["navy"]])) +
      scale_fill_gradient(low = "white", high = col, limits = c(0, 100), name = ttl) +
      labs(x = "Fraction of events observed at the horizon, f", y = NULL) +
      theme_art() + theme(panel.grid = element_blank(), legend.position = "right")
  }
  pB1 <- heat( "pct_above_N",  PAL[["teal"]], "% with\nN* > N" )
  pB  <- heat( "pct_above_2N", PAL[["red"]],  "% with\nN* > 2N" )
  ggsave(file.path(outdir, "ess_above_N.png"),  pB1, width = 8, height = 4.5, dpi = 300, bg = "white")
  ggsave(file.path(outdir, "ess_above_2N.png"), pB,  width = 8, height = 4.5, dpi = 300, bg = "white")
  
  ##### Figure C (= candidate manuscript Figure 4): N*/N by event fraction #####
  pC <- plot_fig4( res, file.path(outdir, "ess_by_eventfraction.png"), type = "box" )
  
  print( tab[order(tab$dataset, tab$f_target), c("dataset", "f_target", "n", "pct_above_N", "pct_above_2N")],
         row.names = FALSE )
  invisible(list(scatter = pA, above2N = pB, aboveN = pB1, box = pC, table = tab))
}

####################### 7. SANITY CHECK AGAINST THE MANUSCRIPT (FIGURE 3) #####################
# Fits the four nested colon models and reports the median 5-year N*.
# The manuscript reports 337 / 329 / 325 / 253. Small deviations are expected
# because the manuscript single-imputes `nodes` (mice) and this uses complete cases.
validate_colon_fig3 <- function( ){
  d <- survival::colon; d <- d[d$etype == 2, ]
  d$rx <- factor(d$rx)
  d <- d[complete.cases(d[, c("time", "status", "rx", "adhere", "nodes", "age")]), ]
  tt <- 5 * 365.25
  sets <- list( A = "rx", B = c("rx", "adhere"), C = c("rx", "adhere", "nodes"),
                D = c("rx", "adhere", "nodes", "age") )
  sapply( sets, function(tr) {
    fit <- fit_cox(tr, d, "efron")
    p <- predict_S_se(fit, d, tt)
    median(p$S * (1 - p$S) / p$SE_S^2)
  } )
}