####################### SIMULATION 1: N* EARLY IN FOLLOW-UP #######################

# Binary predictor (5,000 per group, HR 2, no censoring) and continuous predictor (X ~ N(0,1), N = 10,000,
# HR 2 per SD, no censoring), constant baseline hazard. Used for Figure 5 of the article (article_figures.R)
# and for Supplement 8 (cached by code/run_results.R in results/sim_early.rds).
#
# Usage (from repository root, after source("code/compile.R")):
#   source("code/sim_early.R")
#   early <- simulate_early()
#   early$df_bin; early$df_cont           # long data frames: group, time, survival, N*, variance components

simulate_early <- function( seed = 1923 ){
  scale = c(10)
  shape = c(1)                   # shape 1 means exponential
  maxtime = 100
  hr = 2                         # hazard ratio
  
  ##### Part 1: binary predictor #####
  set.seed(seed)
  n = 5000                       # number per group
  events <- simsurv( dist="weibull",
                     gammas=shape,
                     lambdas=scale[1]^(-shape),
                     x=data.frame(tx = c(rep(0, n), rep(1, n))),
                     betas = c("tx" = log(hr)),
                     maxt=maxtime)
  events$tx = c(rep(0, n), rep(1, n))
  events$eventtime[8885] <- 2.997671  # remove ties by hand
  events$eventtime[6298] <- 5.668818  # works because the seed is fixed
  time <- sort(unique(events$eventtime))
  time0 <- sort(unique(events[events$tx == 0, ]$eventtime))
  time1 <- sort(unique(events[events$tx == 1, ]$eventtime))
  
  sf0 <- survfit(Surv(eventtime, status) ~ 1, data = events[events$tx == 0,])
  sf1 <- survfit(Surv(eventtime, status) ~ 1, data = events[events$tx == 1,])
  cox <- coxph(Surv(eventtime, status) ~ tx, data = events)
  
  sf_n0 <- survfit_n(sf0) |> sf_to_df(time = time0)
  sf_n1 <- survfit_n(sf1) |> sf_to_df(time = time1)
  cox_n0 <- survfit_n(survfit(cox, newdata = data.frame(tx = 0)), cox, coef=T, chaz=T)
  cox_n1 <- survfit_n(survfit(cox, newdata = data.frame(tx = 1)), cox, coef=T, chaz=T)
  cox_n0 <- sf_to_df(cox_n0, time = time)
  cox_n1 <- sf_to_df(cox_n1, time = time)
  
  df_bin <- rbind(
    data.frame(grp = "Cox women", time = time, surv = cox_n0$surv*100, n.eff = cox_n0$n.eff,
               v_base = cox_n0$std.chaz^2, v_total = cox_n0$std.err^2),
    data.frame(grp = "Cox men",   time = time, surv = cox_n1$surv*100, n.eff = cox_n1$n.eff,
               v_base = cox_n1$std.chaz^2, v_total = cox_n1$std.err^2),
    data.frame(grp = "KM women",  time = time0, surv = sf_n0$surv*100, n.eff = sf_n0$n.eff,
               v_base = NA_real_, v_total = NA_real_),
    data.frame(grp = "KM men",    time = time1, surv = sf_n1$surv*100, n.eff = sf_n1$n.eff,
               v_base = NA_real_, v_total = NA_real_))
  df_bin$grp <- factor(df_bin$grp, levels = c("Cox women", "Cox men", "KM women", "KM men"))
  df_bin$share <- df_bin$v_base / df_bin$v_total * 100
  
  # First event time: numbers for the heuristic N* = N * mean(HR) / HR_i (Supplement 8)
  hr_hat <- unname(exp(coef(cox)))
  first <- list( hr_hat = hr_hat,
                 n_eff = c( women = cox_n0$n.eff[1], men = cox_n1$n.eff[1] ),
                 n_eff_baseline_only = c( women = (1 - cox_n0$surv[1]) / (cox_n0$surv[1] * cox_n0$std.chaz[1]^2),
                                          men = (1 - cox_n1$surv[1]) / (cox_n1$surv[1] * cox_n1$std.chaz[1]^2) ),
                 heuristic = c( women = n * hr_hat + n, men = (n * hr_hat + n) / hr_hat ) )
  
  ##### Part 2: continuous predictor #####
  set.seed(seed)
  n = 10000
  x <- data.frame(x = rnorm(n, 0, 1))
  events <- simsurv( dist="weibull",
                     gammas=shape,
                     lambdas=scale[1]^(-shape),
                     x=x,
                     betas = c("x" = log(hr)),
                     maxt=maxtime)
  events$x = x$x
  time <- sort(unique(events$eventtime))
  
  cox <- coxph(Surv(eventtime, status) ~ x, data = events)
  
  cox_mean <- survfit_n(survfit(cox, newdata = data.frame(x = 0)), cox, coef=T, chaz=T)
  cox_sd_plus <- survfit_n(survfit(cox, newdata = data.frame(x = 1)), cox, coef=T, chaz=T)
  cox_sd_min <- survfit_n(survfit(cox, newdata = data.frame(x = -1)), cox, coef=T, chaz=T)
  cox_mean <- sf_to_df(cox_mean, time = time)
  cox_sd_plus <- sf_to_df(cox_sd_plus, time = time)
  cox_sd_min <- sf_to_df(cox_sd_min, time = time)
  
  one <- function(d, grp) data.frame(grp = grp, time = time, surv = d$surv*100, n.eff = d$n.eff,
                                     v_base = d$std.chaz^2, v_total = d$std.err^2)
  df_cont <- rbind(one(cox_mean, "Mean"), one(cox_sd_plus, "+1 SD"), one(cox_sd_min, "-1 SD"))
  df_cont$grp <- factor(df_cont$grp, levels = c("Mean", "+1 SD", "-1 SD"))
  df_cont$share <- df_cont$v_base / df_cont$v_total * 100
  
  list( df_bin = df_bin, df_cont = df_cont, first = first )
}
