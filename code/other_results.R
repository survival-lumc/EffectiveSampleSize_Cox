####################### OTHER RESULTS (Supplement 12) #######################

# (1) Linear and logistic regression in the GUSTO data: effective sample size of the mean patient and bounds by subgroup size
# (2) Colon data: the models of Figure 3 stratified by treatment instead of a hazard ratio for treatment
# Cached by code/run_results.R in results/other_results.rds. Needs Hmisc (getHdata downloads the GUSTO data).
#
# Usage (from repository root, after source("code/compile.R")):
#   source("code/other_results.R"); oth <- run_other()

run_other <- function( seed = 1923 ){
  ##### GUSTO subset #####
  Hmisc::getHdata(gusto)
  g <- gusto[which(gusto$grpl %in% c(1, 2, 4) & gusto$grps %in% c(5, 9, 10, 11)), ]
  g$age <- g$age - mean(g$age)
  g$hei <- g$hei - mean(g$hei)
  n_women <- sum(g$sex == "female")
  
  # N* of a new covariate row (intercept first) in a linear model: 1 / x'(X'X)^-1 x
  ols_neff <- function( formula, row ){
    X <- model.matrix(formula, data = g)
    drop(1 / (row %*% solve(t(X) %*% X) %*% t(row)))
  }
  # N* in a logistic model: p(1-p) / var(linear predictor * p(1-p))  =  1 / (x' V x p (1-p))
  glm_neff <- function( fit, newdata ){
    # build the design row with the factor levels of the fit: newdata has a single level of sex
    tt <- delete.response(terms(fit))
    X <- model.matrix(tt, model.frame(tt, newdata, xlev = fit$xlevels))
    p <- predict(fit, newdata = newdata, type = "response")
    drop(1 / ((X %*% vcov(fit) %*% t(X)) * p * (1 - p)))
  }
  
  m_age <- mean(g$age)
  m_age_w <- mean(g$age[g$sex == "female"])
  
  ##### (1a) mean patient in linear regression #####
  mean_patient <- ols_neff(~ age + hei, matrix(c(1, 0, 0), nrow = 1))
  
  ##### (1b) linear regression of height on age and sex #####
  lin <- c( pop_age   = ols_neff(~ age + sex, matrix(c(1, m_age, 1), nrow = 1)),
            women_age = ols_neff(~ age + sex, matrix(c(1, m_age_w, 1), nrow = 1)) )
  
  ##### (1c) logistic regression without and with interaction #####
  fit0 <- glm(day30 ~ age + sex, data = g, family = binomial)
  fit1 <- glm(day30 ~ age + sex + age:sex, data = g, family = binomial)
  nd <- data.frame(age = c(m_age, m_age_w), sex = "female")
  logi <- data.frame(model = c("without interaction", "with interaction"),
                     pop_age = c(glm_neff(fit0, nd[1, ]), glm_neff(fit1, nd[1, ])),
                     women_age = c(glm_neff(fit0, nd[2, ]), glm_neff(fit1, nd[2, ])))
  
  ##### (2) colon data: stratified by treatment #####
  d <- colon_df
  set.seed(seed)
  d$nodes <- mice::complete(mice::mice(d, m = 1, print = FALSE))$nodes
  year <- 5
  forms <- c("rx", "rx + adhere", "rx + adhere + nodes", "rx + adhere + nodes + age")
  pred_at <- function( strat ){
    do.call(rbind, lapply(seq_along(forms), function( i ) {
      f <- if( strat ) gsub("rx", "strata(rx)", forms[i], fixed = TRUE) else forms[i]
      fit <- coxph(as.formula(paste("Surv(time, status) ~", f)), data = d)
      nd <- d[, c("rx", "age", "adhere", "nodes")]
      if( i == 1 && strat ){               # only strata: the prediction is the Kaplan-Meier estimate of the stratum
        sf <- summary(survfit(Surv(time, status) ~ rx, data = d), time = year * 365.25)
        s <- sf$surv[match(nd$rx, levels(d$rx))]; se <- sf$std.err[match(nd$rx, levels(d$rx))]
      } else {
        sf <- summary(survfit(fit, newdata = nd), time = year * 365.25)
        s <- c(sf$surv); se <- c(sf$std.err)
      }
      data.frame(nd, surv = s, n.eff = s * (1 - s) / se^2, model = LETTERS[i], stratified = strat)
    }))
  }
  strat <- rbind(pred_at(FALSE), pred_at(TRUE))
  
  list( n_gusto = nrow(g), n_women = n_women, mean_patient = mean_patient, lin = lin, logi = logi, strat = strat )
}

plot_strat <- function( oth ){
  d <- oth$strat
  d$panel <- factor(paste0("(", d$model, ") ", c("Tx", "+ Adh", "+ Nodes", "+ Age")[match(d$model, LETTERS[1:4])]))
  d$type <- factor(ifelse(d$stratified, "Stratified by treatment", "Treatment as covariate"))
  ggplot(d, aes(surv, n.eff, colour = type)) +
    geom_point(alpha = 0.2, size = 1.5) + facet_wrap(~ panel) +
    scale_colour_manual(values = c("Treatment as covariate" = "#37293F", "Stratified by treatment" = "#c6aa2c"), name = NULL) +
    scale_x_continuous(name = "Predicted 5-year survival probability", limits = c(0, 1)) +
    scale_y_continuous(name = "Effective sample size", limits = c(0, 500)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
    theme_minimal() + theme(legend.position = "bottom")
}