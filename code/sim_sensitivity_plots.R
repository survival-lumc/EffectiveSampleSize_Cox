####################### SUPPLEMENT 7: LIMITS TABLE AND FIGURE-6-LIKE SENSITIVITY PLOTS #####################
# Usage: source("code/sim_guardrails.R") first (main_data, bin_s, S_cross_main, S_SCALE, scale_x_S, ramp_pal, theme_art).
# Not yet run (written without an R session): run supplement.qmd chunk by chunk and check the figures.

# Lower and upper limit of the true survival at which more than `tol` of the predictions have N* > N, per scenario
# (C x censoring x event rate), with the same binning as Figure 6 (plot_admin: bins of at least min_n = 5 cells).
# Why the earlier table only had NA for the lower limit: fit_threshold() used min_n = 20 cells per bin, and the
# bins at low survival contain few cells in a single scenario, so they were dropped before a crossing was found.
threshold_limits <- function( res, tol = 0.05, min_n = 5, trim = TRIM_PROFILES ){
  r <- main_data( res, trim )
  r$grp <- paste( r$C, r$admin, r$e )
  cr <- S_cross_main( r, tol, min_n = min_n, group = "grp" )
  cr$C <- as.numeric(sub(" .*", "", cr$grp)); cr$admin <- grepl( " TRUE ", cr$grp )
  cr$e <- as.numeric(sub(".* ", "", cr$grp))
  cr
}

# Sensitivity to Harrell's C, like Figure 6: percentage of predictions with N* > N against the true survival
# (logit scale); one panel per event rate, colour = Harrell's C; no censoring.
plot_C_sens <- function( res, outdir = tempdir(), min_n = 5, trim = TRIM_PROFILES ){
  dir.create( outdir, showWarnings = FALSE, recursive = TRUE )
  r <- main_data( res, trim )
  r <- r[!r$admin, ]
  sc <- S_SCALE$logit
  b <- bin_s( r, "P_gt_N", c("Cl", "el"), sc, min_n = min_n )
  p <- ggplot(b, aes(xbin, P_gt_N, colour = Cl)) +
    geom_hline(yintercept = 0.05, linetype = 2, colour = "grey40") + geom_line(linewidth = 0.8) +
    facet_wrap(~el, nrow = 1, labeller = labeller(el = function( x ) paste0("Event rate ", 100 * as.numeric(as.character(x)), "%"))) +
    coord_cartesian(xlim = sc$xlim) + scale_x_S(modifyList(sc, list(name = "S(t|x), logit scale"))) +
    scale_colour_manual(values = ramp_pal(nlevels(b$Cl)), name = "Harrell's C", labels = function( x ) sub("C = ", "", x)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(y = "Predictions with N* > N") + theme_art(11)
  ggsave(file.path(outdir, "sim_csens_PgtN.png"), p, width = 10, height = 3.8, dpi = 300, bg = "white")
  p
}

# Sensitivity to the Weibull shape (design "shape_only"), like Figure 6: percentage of predictions with N* > N
# against the true survival, pooled over the time points; colour / line type = shape.
plot_shape_main <- function( res, outdir = tempdir(), min_n = 5, trim = TRIM_PROFILES ){
  dir.create( outdir, showWarnings = FALSE, recursive = TRUE )
  r <- trim_core( res, trim )
  r$shp <- factor( r$shape, levels = sort(unique(r$shape)), labels = paste("shape", sort(unique(r$shape))) )
  sc <- S_SCALE$logit
  b <- bin_s( r, "P_gt_N", "shp", sc, min_n = min_n )
  p <- ggplot(b, aes(xbin, P_gt_N, colour = shp, linetype = shp)) +
    geom_hline(yintercept = 0.05, linetype = 2, colour = "grey40") + geom_line(linewidth = 0.9) +
    coord_cartesian(xlim = sc$xlim) + scale_x_S(modifyList(sc, list(name = "S(t|x), logit scale"))) +
    scale_colour_manual(values = ramp_pal(nlevels(b$shp)), name = "Weibull shape", labels = function( x ) sub("shape ", "", x)) +
    scale_linetype_manual(values = c("solid", "22", "dotted")[seq_len(nlevels(b$shp))], name = "Weibull shape",
                          labels = function( x ) sub("shape ", "", x)) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(y = "Predictions with N* > N") + theme_art(11)
  ggsave(file.path(outdir, "sim_shape_PgtN_main.png"), p, width = 8.5, height = 5, dpi = 300, bg = "white")
  p
}

# Number of predictors and noise predictors (replaces the old 4 x 4 panel figure): reference scenario
# (N, event rate, shape, no censoring), pooled over Harrell's C and the fractions of events observed, like
# Figure S7.4. Colour = number of true predictors k, line type = without / with k added noise predictors.
k_data <- function( res, trim = TRIM_PROFILES ){
  r <- trim_core( res, trim )
  r <- r[r$N == REF$N & r$e == REF$e & r$shape == REF$shape, ]
  if( "rc" %in% names(r) ) r <- r[r$rc == 0, ]
  if( "admin" %in% names(r) ) r <- r[!r$admin, ]
  r$kf <- factor( r$k, levels = sort(unique(r$k)) )
  r$noisef <- factor( r$noise, levels = c(0, 1), labels = c("no noise predictors", "with noise predictors") )
  r$setting <- paste0( "k = ", r$k, ifelse(r$noise == 1, ", +noise", "") )
  r
}

# S*: true survival above which more than `tol` of the predictions have N* > N, per setting (NA: not reached)
threshold_k <- function( res, tol = 0.05, min_n = 5, trim = TRIM_PROFILES ){
  r <- k_data( res, trim )
  cr <- S_cross_main( r, tol, min_n = min_n, group = "setting" )
  cr$k <- as.numeric(sub("k = ([0-9.]+).*", "\\1", cr$setting)); cr$noise <- as.numeric(grepl("noise", cr$setting))
  cr
}

plot_k_sens <- function( res, outdir = tempdir(), min_n = 5, trim = TRIM_PROFILES ){
  dir.create( outdir, showWarnings = FALSE, recursive = TRUE )
  r <- k_data( res, trim )
  sc <- S_SCALE$logit
  b <- bin_s( r, "P_gt_N", c("kf", "noisef"), sc, min_n = min_n )
  p <- ggplot(b, aes(xbin, P_gt_N, colour = kf, linetype = noisef)) +
    geom_hline(yintercept = 0.05, linetype = 2, colour = "grey40") + geom_line(linewidth = 0.9) +
    coord_cartesian(xlim = sc$xlim) + scale_x_S(modifyList(sc, list(name = "S(t|x), logit scale"))) +
    scale_colour_manual(values = ramp_pal(nlevels(b$kf)), name = "True predictors, k") +
    scale_linetype_manual(values = c("solid", "22"), name = NULL) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(y = "Predictions with N* > N") + theme_art(11)
  ggsave(file.path(outdir, "sim_k_PgtN_main.png"), p, width = 8.5, height = 5, dpi = 300, bg = "white")
  p
}

# Censoring (replaces the 3 x 5 panel figure): Harrell's C = 0.7, pooled over the fractions of events observed;
# one panel per event rate, colour = intensity of random censoring, line type = no / uniform censoring.
plot_cens_main <- function( res, outdir = tempdir(), min_n = 5, trim = TRIM_PROFILES, C_show = 0.7 ){
  dir.create( outdir, showWarnings = FALSE, recursive = TRUE )
  r <- trim_core( res, trim ); r <- r[r$C == C_show, ]
  r$el <- factor( r$e, labels = paste0("Event rate ", 100 * sort(unique(r$e)), "%") )
  r$rcf <- factor( r$rc, levels = c(0, 1, 5), labels = c("none", "regular (1)", "severe (5)") )
  r$adminf <- factor( r$admin, levels = c(FALSE, TRUE), labels = c("no uniform censoring", "uniform censoring") )
  sc <- S_SCALE$logit
  b <- bin_s( r, "P_gt_N", c("el", "rcf", "adminf"), sc, min_n = min_n )
  p <- ggplot(b, aes(xbin, P_gt_N, colour = rcf, linetype = adminf)) +
    geom_hline(yintercept = 0.05, linetype = 2, colour = "grey40") + geom_line(linewidth = 0.8) +
    facet_wrap(~el, nrow = 1) +
    coord_cartesian(xlim = sc$xlim) + scale_x_S(modifyList(sc, list(name = "S(t|x), logit scale"))) +
    scale_colour_manual(values = ramp_pal(3), name = "Random censoring") +
    scale_linetype_manual(values = c("solid", "22"), name = NULL) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(y = "Predictions with N* > N") + theme_art(11)
  ggsave(file.path(outdir, "sim_cens_PgtN_main.png"), p, width = 11, height = 4, dpi = 300, bg = "white")
  p
}
