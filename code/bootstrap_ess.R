####################### EFFECTIVE N FROM A BOOTSTRAPPED VARIANCE #####################

# Supplement 17. N*_boot = S(1 - S) / Var_boot(S-hat), where Var_boot is the variance over B
# bootstrap refits of the same Cox model, evaluated for the original patients' covariates;
# S is the estimate from the original fit. Two bootstrap types:
#   "nonparametric": patients resampled with replacement
#   "parametric":    outcomes simulated from the fitted Cox model (Breslow baseline hazard) for the
#                    original covariates; censoring times drawn from the Kaplan-Meier estimate of the
#                    censoring distribution (independent of covariates). Right-censored data only.
# N*_delta uses the delta-method variance of survfit(). Replaces the old bootstrap.R.
# Test on the colon data first with B = 50.

library(survival)
library(future.apply)

boot_ess <- function( dat, formula, years = 5, B = 1000, seed = 1923, model = "model",
                      type = c("nonparametric", "parametric"),
                      workers = max(1, parallel::detectCores() - 1) ){
  type <- match.arg(type)
  times <- years * 365.25                                   # time in the data is in days
  fit <- coxph( formula, data = dat, ties = "efron", model = TRUE, x = TRUE, y = TRUE )
  
  ##### delta-method estimate (same extraction as in the simulations) #####
  sf  <- survfit( fit, newdata = dat )
  idx <- findInterval( times, sf$time )
  stopifnot(all(idx >= 1))
  S   <- as.matrix(sf$surv)[idx, , drop = FALSE]            # horizons x patients
  seH <- as.matrix(sf$std.err)[idx, , drop = FALSE]         # SE of the cumulative hazard (log S)
  se_delta <- S * seH
  
  ##### ingredients of the parametric bootstrap #####
  resp <- all.vars( formula[[2]] )                          # names of the time and status variables
  H_fit <- basehaz( fit, centered = FALSE )                 # Breslow baseline cumulative hazard
  lp_fit <- predict( fit, newdata = dat, type = "lp", reference = "zero" )
  G <- survfit( Surv(dat[[resp[1]]], 1 - dat[[resp[2]]]) ~ 1 )   # KM of the censoring distribution
  
  ##### bootstrap #####
  one <- function( b ){
    out <- tryCatch( {
      if( type == "nonparametric" ){
        db  <- dat[sample.int(nrow(dat), replace = TRUE), ]
      }
      else{
        n   <- nrow(dat)
        target <- -log(runif(n)) / exp(lp_fit)              # H0(T) = target
        j   <- findInterval( target, H_fit$hazard, left.open = TRUE ) + 1
        tev <- c( H_fit$time, Inf )[j]                      # Inf: no event within the observed follow-up
        jc  <- findInterval( -runif(n), -G$surv, left.open = TRUE ) + 1
        cen <- c( G$time, max(dat[[resp[1]]]) )[jc]         # beyond the last censoring time: administrative end
        db  <- dat
        db[[resp[1]]] <- pmin( tev, cen )
        db[[resp[2]]] <- as.numeric( tev <= cen )
      }
      fb_formula <- formula
      environment(fb_formula) <- environment()               # so that coxph finds db and not a global object
      fb  <- coxph( fb_formula, data = db, ties = "efron" )
      H   <- basehaz(fb, centered = FALSE)                  # baseline at x = 0
      H0  <- approx(c(0, H$time), c(0, H$hazard), xout = times, method = "constant", rule = 2)$y
      lpb <- predict(fb, newdata = dat, type = "lp", reference = "zero")
      exp(-outer(H0, exp(lpb)))                             # horizons x patients
    }, error = function( e ){
      m <- matrix(NA_real_, length(times), nrow(dat))
      attr(m, "msg") <- conditionMessage(e)
      m
    } )
    out
  }
  future::plan( future::multisession, workers = workers )
  set.seed(seed)
  reps <- future_lapply( seq_len(B), one, future.seed = seed, future.packages = "survival" )
  future::plan(future::sequential)
  msgs <- unique(unlist(lapply(reps, attr, "msg")))
  if( length(msgs) > 0 ) warning( "bootstrap refits failed: ", paste(msgs, collapse = " | ") )
  arr <- array( unlist(reps), dim = c(length(times), nrow(dat), B) )
  n_fail <- sum(apply(arr, 3, function(a) all(is.na(a))))
  se_boot <- apply( arr, c(1, 2), sd, na.rm = TRUE )
  
  data.frame(
    model = model, type = type, patient = rep(seq_len(nrow(dat)), each = length(years)),
    year = rep(years, times = nrow(dat)),
    S = as.vector(S), se_delta = as.vector(se_delta), se_boot = as.vector(se_boot),
    N_delta = as.vector(S * (1 - S) / se_delta^2),
    N_boot  = as.vector(S * (1 - S) / se_boot^2),
    n_total = nrow(dat), B = B, failed = n_fail )
}

# Colon data, as in Figure 3 of the article (overall survival; nodes single-imputed with mice on the whole
# data frame, as in article_figures.R). Needs colon_df from code/compile.R (etype 2, time in days).
colon_for_bootstrap <- function( seed = 1923 ){
  d <- colon_df
  set.seed(seed)
  d$nodes <- mice::complete(mice::mice(d, m = 1, print = FALSE))$nodes
  d
}

colon_forms <- list(
  "A: treatment"                    = Surv(time, status) ~ rx,
  "B: + adherence"                  = Surv(time, status) ~ rx + adhere,
  "C: + nodes"                      = Surv(time, status) ~ rx + adhere + nodes,
  "D: + age"                        = Surv(time, status) ~ rx + adhere + nodes + age )

run_colon_bootstrap <- function( years = 1:7, B = 1000, seed = 1923,
                                 types = c("nonparametric", "parametric") ){
  d <- colon_for_bootstrap(seed)
  do.call( rbind, lapply(types, function(ty)
    do.call( rbind, lapply(names(colon_forms), function(m)
      boot_ess(d, colon_forms[[m]], years = years, B = B, seed = seed, model = m, type = ty)) )) )
}

# Summary per model and horizon: ratio N*_boot / N*_delta
summarise_boot <- function( res ){
  if( is.null(res$type) ) res$type <- "nonparametric"       # results from before the parametric option
  res <- res[is.finite(res$N_delta) & is.finite(res$N_boot) & res$N_delta > 0 & res$N_boot > 0, ]
  res$ratio <- res$N_boot / res$N_delta
  do.call( rbind, lapply(split(res, list(res$type, res$model, res$year), drop = TRUE), function(g)
    data.frame(type = g$type[1], model = g$model[1], year = g$year[1], n = nrow(g),
               med_delta = median(g$N_delta), med_boot = median(g$N_boot),
               ratio_med = median(g$ratio), ratio_q25 = quantile(g$ratio, 0.25), ratio_q75 = quantile(g$ratio, 0.75),
               within10 = mean(abs(g$ratio - 1) <= 0.10),
               spearman = cor(g$N_delta, g$N_boot, method = "spearman"),
               row.names = NULL)) )
}

PAL_BOOT <- c( "#37293F", "#2E7691", "#9f84af", "#c6aa2c" )

# Figure: N*_boot against N*_delta at one horizon for one bootstrap type, one panel per model
# (2 x 2; log-log; dashed = equality)
plot_boot_scatter <- function( res, year = 5, type = "nonparametric" ){
  if( is.null(res$type) ) res$type <- "nonparametric"
  r <- res[res$type == type & res$year == year & is.finite(res$N_delta) & is.finite(res$N_boot), ]
  r$model <- factor( r$model, levels = names(colon_forms) )
  lim <- range(c(r$N_delta, r$N_boot))
  ggplot2::ggplot(r, ggplot2::aes(N_delta, N_boot, colour = model)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
    ggplot2::geom_point(alpha = 0.25, size = 1.5) +
    ggplot2::scale_x_log10(limits = lim) + ggplot2::scale_y_log10(limits = lim) +
    ggplot2::scale_colour_manual(values = setNames(PAL_BOOT, levels(r$model)), guide = "none") +
    ggplot2::facet_wrap(~model, nrow = 2) +
    ggplot2::labs(x = expression(N[`*`] ~ "from the delta-method variance"),
                  y = expression(N[`*`] ~ "from the bootstrap variance")) +
    ggplot2::theme_minimal()
}