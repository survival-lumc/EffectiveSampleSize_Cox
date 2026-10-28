################################################################################
############################## INTRODUCTION ####################################
################################################################################
# This file contains all code necessary to reproduce the seven figures of the main
# article (numbering as in the manuscript). The figures were edited in Adobe
# Illustrator to make the legends more readable; the contents were not changed.
#
# Run from the repository root. Figures 4, 6 and 7 are based on analyses that take
# long; they are run once by code/run_results.R and saved in results/. Here they
# are only read, so this file runs in about a minute.
#
#   Figure 1  Colon data: Cox vs Kaplan-Meier, treatment and adhesion
#   Figure 2  Colon data: Cox model on age and nodes, common to rare patients
#   Figure 3  Colon data: N* of four nested models at 5 years (single imputation of nodes)
#   Figure 4  Eight open real-world datasets: N*/N by fraction of events observed  [results/ess_realworld.rds]
#   Figure 5  Simulation 1: binary and continuous predictor, no censoring         [computed here, seconds]
#   Figure 6  Simulation 2: % of predictions with N* > N against true survival     [results/sim_main.rds]
#   Figure 7  Simulation 2: coverage of the 95% CI against N*/N                    [results/sim_main.rds]
#
# Output: images/fig1.png ... images/fig7.png, used by manuscript.qmd.

source("code/compile.R")
dir.create("images", showWarnings = FALSE)
save_fig <- function( p, n, width, height ) ggsave(file.path("images", paste0("fig", n, ".png")), p,
                                                   width = width, height = height, dpi = 300, bg = "white")
load_res <- function( file ) {
  path <- file.path("results", file)
  if( !file.exists(path) ) stop("Missing ", path, ". Run code/run_results.R first.")
  readRDS(path)
}

################################################################################
################################ FIGURE 1 ######################################
################################################################################
colon_t <- colon_df
colon_t$time <- colon_t$time / 365.241

# Fit all the models
cox_adhere <- coxph( Surv(time, status) ~ adhere, data = colon_t)
cox_rx <- coxph( Surv(time, status) ~ rx, data = colon_t)

# Predict the effective sample size
km_adhere_0 <- survfit_n( survfit(Surv(time, status) ~ 1, data = colon_t[colon_t$adhere == 0,]) )
km_adhere_1 <- survfit_n( survfit(Surv(time, status) ~ 1, data = colon_t[colon_t$adhere == 1,]) )
km_rx_0 <- survfit_n( survfit(Surv(time, status) ~ 1, data = colon_t[colon_t$rx == "Obs",]) )
km_rx_1 <- survfit_n( survfit(Surv(time, status) ~ 1, data = colon_t[colon_t$rx == "Lev",]) )
km_rx_2 <- survfit_n( survfit(Surv(time, status) ~ 1, data = colon_t[colon_t$rx == "Lev+5FU",]) )
cox_adhere_0 <- survfit_n( survfit(cox_adhere, newdata=data.frame(adhere=0)), cox=cox_adhere)
cox_adhere_1 <- survfit_n( survfit(cox_adhere, newdata=data.frame(adhere=1)), cox=cox_adhere)
cox_rx_0 <- survfit_n( survfit(cox_rx, newdata=data.frame(rx="Obs")), cox=cox_rx)
cox_rx_1 <- survfit_n( survfit(cox_rx, newdata=data.frame(rx="Lev")), cox=cox_rx)
cox_rx_2 <- survfit_n( survfit(cox_rx, newdata=data.frame(rx="Lev+5FU")), cox=cox_rx)

color_values <- c("Obs"= "#2E7691", "Lev" = "#C2666B", "No Adhesion"= "#37293F", "Lev+5FU" = "#c6aa2c", "Adhesion" = "#9f84af")
linetype_values <- c("Cox" = "solid", "Kaplan-Meier" = "dashed")
x_breaks <- seq(0, max(colon_t$time), by = 2)

# Create the plots
plot1 <- ggplot() +
  geom_line(aes(x = cox_rx_0$time, y = cox_rx_0$surv*100, color = "Obs", linetype = "Cox"), lwd = 1.3) +
  geom_line(aes(x = cox_rx_1$time, y = cox_rx_1$surv*100, color = "Lev", linetype = "Cox"), lwd = 1.3) +
  geom_line(aes(x = cox_rx_2$time, y = cox_rx_2$surv*100, color = "Lev+5FU", linetype = "Cox"), lwd = 1.3) +
  geom_line(aes(x = km_rx_0$time, y = km_rx_0$surv*100, color = "Obs", linetype = "Kaplan-Meier"), lwd = 1.3) +
  geom_line(aes(x = km_rx_1$time, y = km_rx_1$surv*100, color = "Lev", linetype = "Kaplan-Meier"), lwd = 1.3) +
  geom_line(aes(x = km_rx_2$time, y = km_rx_2$surv*100, color = "Lev+5FU", linetype = "Kaplan-Meier"), lwd = 1.3) +
  scale_color_manual(values = color_values) +
  scale_linetype_manual(values = linetype_values) +
  theme_minimal() +
  scale_y_continuous(name = "Predicted survival (%)", limits = c(0, 100)) +
  scale_x_continuous(name = NULL, labels = NULL, breaks = x_breaks) +
  ggtitle("Model based on treatment") +
  theme(legend.position = "none")

plot2 <- ggplot() +
  geom_line(aes(x = cox_rx_0$time, y = cox_rx_0$n.eff, color = "Obs", linetype = "Cox"), lwd = 1.3) +
  geom_line(aes(x = cox_rx_1$time, y = cox_rx_1$n.eff, color = "Lev", linetype = "Cox"), lwd = 1.3) +
  geom_line(aes(x = cox_rx_2$time, y = cox_rx_2$n.eff, color = "Lev+5FU", linetype = "Cox"), lwd = 1.3) +
  geom_line(aes(x = km_rx_0$time, y = km_rx_0$n.eff, color = "Obs", linetype = "Kaplan-Meier"), lwd = 1.3) +
  geom_line(aes(x = km_rx_1$time, y = km_rx_1$n.eff, color = "Lev", linetype = "Kaplan-Meier"), lwd = 1.3) +
  geom_line(aes(x = km_rx_2$time, y = km_rx_2$n.eff, color = "Lev+5FU", linetype = "Kaplan-Meier"), lwd = 1.3) +
  geom_hline(yintercept = 929, linetype = "dotted", color = "black", linewidth = 0.6) +
  scale_color_manual(values = color_values) +
  scale_linetype_manual(values = linetype_values) +
  theme_minimal() +
  scale_y_continuous(name = "Effective Sample Size", limits = c(0, 1000)) +
  scale_x_continuous(name = "Time (years)", breaks = x_breaks) +
  theme(legend.position = "none")

plot3 <- ggplot() +
  geom_line(aes(x = cox_adhere_0$time, y = cox_adhere_0$surv*100, color = "No Adhesion", linetype = "Cox"), lwd = 1.3) +
  geom_line(aes(x = cox_adhere_1$time, y = cox_adhere_1$surv*100, color = "Adhesion", linetype = "Cox"), lwd = 1.3) +
  geom_line(aes(x = km_adhere_0$time, y = km_adhere_0$surv*100, color = "No Adhesion", linetype = "Kaplan-Meier"), lwd = 1.3) +
  geom_line(aes(x = km_adhere_1$time, y = km_adhere_1$surv*100, color = "Adhesion", linetype = "Kaplan-Meier"), lwd = 1.3) +
  scale_color_manual(values = color_values) +
  scale_linetype_manual(values = linetype_values) +
  theme_minimal() +
  scale_y_continuous(name = NULL, labels = NULL, limits = c(0, 100)) +
  scale_x_continuous(name = NULL, labels = NULL, breaks = x_breaks) +
  ggtitle("Model based on adhesion") +
  theme(legend.position = "none")

plot4 <- ggplot() +
  geom_line(aes(x = cox_adhere_0$time, y = cox_adhere_0$n.eff, color = "No Adhesion", linetype = "Cox"), lwd = 1.3) +
  geom_line(aes(x = cox_adhere_1$time, y = cox_adhere_1$n.eff, color = "Adhesion", linetype = "Cox"), lwd = 1.3) +
  geom_line(aes(x = km_adhere_0$time, y = km_adhere_0$n.eff, color = "No Adhesion", linetype = "Kaplan-Meier"), lwd = 1.3) +
  geom_line(aes(x = km_adhere_1$time, y = km_adhere_1$n.eff, color = "Adhesion", linetype = "Kaplan-Meier"), lwd = 1.3) +
  geom_hline(yintercept = 929, linetype = "dotted", color = "black", linewidth = 0.6) +
  scale_color_manual(values = color_values) +
  scale_linetype_manual(values = linetype_values) +
  theme_minimal() +
  scale_y_continuous(name = NULL, labels = NULL, limits = c(0, 1000)) +
  scale_x_continuous(name = "Time (years)", breaks = x_breaks) +
  theme(legend.position = "none")

# Combine the plots using patchwork and add a single legend on the right
figure1 <- (plot1 / plot2 | plot3 / plot4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")
print(figure1)
save_fig(figure1, 1, 9, 6)


################################################################################
################################ FIGURE 2 ######################################
################################################################################
# Percentiles of age and nodes
age_percentiles <- quantile(colon_t$age, probs = c(0.5, 0.75, 0.99), na.rm = TRUE)
nodes_percentiles <- round(quantile(colon_t$nodes, probs = c(0.5, 0.75, 0.99), na.rm = TRUE))
nodes_percentiles["99%"] <- 20  # set by hand: more extreme than the 99th percentile

# Fit Cox models for age and nodes
cox_age <- coxph(Surv(time, status) ~ age, data = colon_t)
cox_nodes <- coxph(Surv(time, status) ~ nodes, data = colon_t)

# Effective sample size for age and nodes at the three percentiles
cox_age_median <- survfit_n(survfit(cox_age, newdata = data.frame(age = age_percentiles["50%"])), cox = cox_age)
cox_age_75 <- survfit_n(survfit(cox_age, newdata = data.frame(age = age_percentiles["75%"])), cox = cox_age)
cox_age_99 <- survfit_n(survfit(cox_age, newdata = data.frame(age = age_percentiles["99%"])), cox = cox_age)
cox_nodes_2 <- survfit_n(survfit(cox_nodes, newdata = data.frame(nodes = nodes_percentiles["50%"])), cox = cox_nodes)
cox_nodes_5 <- survfit_n(survfit(cox_nodes, newdata = data.frame(nodes = nodes_percentiles["75%"])), cox = cox_nodes)
cox_nodes_25 <- survfit_n(survfit(cox_nodes, newdata = data.frame(nodes = nodes_percentiles["99%"])), cox = cox_nodes)

# Survival probabilities. The survival curves are those of the same three patients as above
# (the `surv` element of the survfit_n objects would also do)
surv_age_median <- survfit(cox_age, newdata = data.frame(age = age_percentiles["50%"]))
surv_age_75 <- survfit(cox_age, newdata = data.frame(age = age_percentiles["75%"]))
surv_age_99 <- survfit(cox_age, newdata = data.frame(age = age_percentiles["99%"]))
surv_nodes_2 <- survfit(cox_nodes, newdata = data.frame(nodes = nodes_percentiles["50%"]))
surv_nodes_5 <- survfit(cox_nodes, newdata = data.frame(nodes = nodes_percentiles["75%"]))
surv_nodes_25 <- survfit(cox_nodes, newdata = data.frame(nodes = nodes_percentiles["99%"]))

color_values <- c("Common" = "#37293F", "Uncommon" = "#2E7691", "Rare" = "#C2666B")

plot1 <- ggplot() +
  geom_line(aes(x = surv_age_median$time, y = surv_age_median$surv*100, color = "Common"), lwd = 1.3) +
  geom_line(aes(x = surv_age_75$time, y = surv_age_75$surv*100, color = "Uncommon"), lwd = 1.3) +
  geom_line(aes(x = surv_age_99$time, y = surv_age_99$surv*100, color = "Rare"), lwd = 1.3) +
  scale_color_manual(values = color_values) +
  theme_minimal() +
  scale_y_continuous(name = "Predicted survival (%)", limits = c(0, 100)) +
  scale_x_continuous(labels = NULL, name = NULL, breaks = x_breaks) +
  ggtitle("Cox model based on age") +
  theme(legend.position = "none")

plot2 <- ggplot() +
  geom_line(aes(x = surv_nodes_2$time, y = surv_nodes_2$surv*100, color = "Common"), lwd = 1.3) +
  geom_line(aes(x = surv_nodes_5$time, y = surv_nodes_5$surv*100, color = "Uncommon"), lwd = 1.3) +
  geom_line(aes(x = surv_nodes_25$time, y = surv_nodes_25$surv*100, color = "Rare"), lwd = 1.3) +
  scale_color_manual(values = color_values) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 100), labels = NULL, name = NULL) +
  scale_x_continuous(labels = NULL, name = NULL, breaks = x_breaks) +
  ggtitle("Cox model based on number of nodes") +
  theme(legend.position = "none")

plot3 <- ggplot() +
  geom_line(aes(x = cox_age_median$time, y = cox_age_median$n.eff, color = "Common"), lwd = 1.3) +
  geom_line(aes(x = cox_age_75$time, y = cox_age_75$n.eff, color = "Uncommon"), lwd = 1.3) +
  geom_line(aes(x = cox_age_99$time, y = cox_age_99$n.eff, color = "Rare"), lwd = 1.3) +
  geom_hline(yintercept = 929, linetype = "dotted", color = "black", linewidth = 0.6) +
  scale_color_manual(values = color_values) +
  theme_minimal() +
  scale_y_continuous(name = "Effective Sample Size", limits = c(0, 1000)) +
  scale_x_continuous(name = "Time (years)", breaks = x_breaks) +
  theme(legend.position = "none")

plot4 <- ggplot() +
  geom_line(aes(x = cox_nodes_2$time, y = cox_nodes_2$n.eff, color = "Common"), lwd = 1.3) +
  geom_line(aes(x = cox_nodes_5$time, y = cox_nodes_5$n.eff, color = "Uncommon"), lwd = 1.3) +
  geom_line(aes(x = cox_nodes_25$time, y = cox_nodes_25$n.eff, color = "Rare"), lwd = 1.3) +
  geom_hline(yintercept = 929, linetype = "dotted", color = "black", linewidth = 0.6) +
  scale_color_manual(values = color_values) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1000), labels = NULL, name = NULL) +
  scale_x_continuous(name = "Time (years)", breaks = x_breaks) +
  theme(legend.position = "none")

figure2 <- (plot1 / plot3 | plot2 / plot4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")
print(figure2)
save_fig(figure2, 2, 9, 6)


################################################################################
################################ FIGURE 3 ######################################
################################################################################
year <- 5
var1 <- "rx"; var2 <- "rx + adhere"; var3 <- "rx + adhere + nodes"; var4 <- "rx + adhere + nodes + age"

# Single imputation of the missing nodes. The seed makes the figure (and the medians quoted in the
# text and in Supplement 11) reproducible; mice is random.
colon_imp <- colon_df
set.seed(1923)
colon_imp$nodes <- mice::complete(mice::mice(colon_imp, m = 1, print = FALSE))$nodes

# Fit all Cox models
cox_models <- lapply(c(var1, var2, var3, var4), function( v )
  coxph(as.formula(paste0("Surv(time, status) ~ ", v)), data = colon_imp))
model_lab <- paste0(c("(A) ", "(B) ", "(C) ", "(D) "), "Surv( time, status) ~ ", c(var1, var2, var3, var4))

# 5-year predictions for all patients in the data: survival, effective sample size
out_df <- do.call(rbind, lapply(1:4, function( i ) {
  d <- colon_imp[, c("rx", "age", "adhere", "nodes")]
  sf <- summary(survfit(cox_models[[i]], newdata = d), time = year * 365.25)
  d$surv <- c(sf$surv)
  d$n.eff <- d$surv * (1 - d$surv) / c(sf$std.err)^2
  d$model <- model_lab[i]
  d
}))
out_df$model <- factor(out_df$model, levels = model_lab)

colors <- setNames(c("#37293F", "#2E7691", "#9f84af", "#c6aa2c"), model_lab)

figure3 <- ggplot(out_df) +
  geom_point( aes(x = surv, y = n.eff, color = model), alpha = 0.2, size = 2) +
  facet_wrap( ~model ) +
  scale_x_continuous(name = "Predicted 5-Year Survival Probability", limits = c(0, 1)) +
  scale_y_continuous(name = "Effective Sample Size", limits = c(0, 500)) +
  scale_color_manual(values = colors) +
  theme(legend.position = "none")
print(figure3)
save_fig(figure3, 3, 8, 6)

# Increase in Cox-Snell R^2 when adding nodes and age (quoted in the text)
LL_null <- cox_models[[2]]$loglik[1]
n_col <- nrow(colon_imp)
R2_CS <- function( LL_model ) 1 - exp((2 / n_col) * (LL_null - LL_model))
print(paste0("Delta R^2 after addition of nodes is: ", round((R2_CS(cox_models[[3]]$loglik[2]) - R2_CS(cox_models[[2]]$loglik[2])) * 100, 2), "%"))
print(paste0("Delta R^2 after addition of age is: ", round((R2_CS(cox_models[[4]]$loglik[2]) - R2_CS(cox_models[[3]]$loglik[2])) * 100, 2), "%"))


################################################################################
################################ FIGURE 4 ######################################
################################################################################
# N*/N by the fraction of events observed at the horizon, for 8 open datasets (panels ordered by
# Harrell's C). The analysis (code/ess_realworld.R) is run once by code/run_results.R.
source("code/ess_realworld.R")
res_rw <- load_res("ess_realworld.rds")
figure4 <- plot_fig4(res_rw, outfile = file.path("images", "fig4.png"), type = "box")
print(figure4)


################################################################################
################################ FIGURE 5 ######################################
################################################################################
# Simulation 1 (no censoring, N = 10,000, HR 2, constant baseline hazard). Top row: binary predictor,
# bottom row: continuous predictor. Takes seconds; the same function is used for Supplement 8.
source("code/sim_early.R")
early <- simulate_early(seed = 1923)

cols_bin <- c("Cox women" = "#37293F", "Cox men" = "#9f84af", "KM women" = "#2E7691", "KM men" = "#c6aa2c")
cols_cont <- c("Mean" = "#37293F", "+1 SD" = "#C2666B", "-1 SD" = "#2E7691")

fig5_panels <- function( d, cols, tag_row, ylim_n = c(0, 15000), lty_km = TRUE ) {
  d$lt <- ifelse(grepl("^KM", d$grp), "Kaplan-Meier", "Cox")
  base <- list(scale_colour_manual(values = cols, name = NULL), theme_minimal(), theme(legend.position = "right"))
  pA <- ggplot(d, aes(time, surv, colour = grp, linetype = lt)) + geom_line(linewidth = 1.2) + base +
    scale_linetype_manual(values = c(Cox = "solid", "Kaplan-Meier" = "dashed"), name = NULL) +
    scale_x_continuous(name = "Time", limits = c(0, 15)) +
    scale_y_continuous(name = "Predicted survival (%)", limits = c(0, 100)) +
    labs(tag = tag_row[1])
  pB <- ggplot(d, aes(time, n.eff, colour = grp, linetype = lt)) + geom_line(linewidth = 1.2) + base +
    scale_linetype_manual(values = c(Cox = "solid", "Kaplan-Meier" = "dashed"), name = NULL) +
    scale_x_continuous(name = NULL, labels = NULL, limits = c(0, 15)) +
    scale_y_continuous(name = "Effective Sample Size", limits = ylim_n) +
    labs(tag = tag_row[2])
  dv <- d[!is.na(d$share), ]
  pC <- ggplot(dv, aes(time, share, colour = grp)) + geom_line(linewidth = 1.2) + base +
    scale_x_continuous(name = "Time", limits = c(0, 15)) +
    scale_y_continuous(name = "Percentage of variance from h0", limits = c(0, 100)) +
    labs(tag = tag_row[3])
  list(pA, pB, pC)
}
top <- fig5_panels(early$df_bin, cols_bin, c("(A)", "(B)", "(C)"))
bottom <- fig5_panels(early$df_cont, cols_cont, c("(D)", "(E)", "(F)"))
figure5 <- (top[[1]] | top[[2]] | top[[3]]) / (bottom[[1]] | bottom[[2]] | bottom[[3]]) +
  plot_layout(guides = "collect") & theme(legend.position = "right")
print(figure5)
save_fig(figure5, 5, 11, 6)


################################################################################
################################ FIGURE 6 ######################################
################################################################################
# Proportion of predictions with N* > N against the true survival, by event rate and censoring
# (Harrell's C = 0.8, 7 event rates x 2 censoring settings, 1100 replications; code/sim_guardrails.R)
source("code/sim_guardrails.R")
resM <- load_res("sim_main.rds")
figure6 <- plot_admin(resM, outdir = tempdir(), scales = "logit")$logit
print(figure6)
save_fig(figure6, 6, 8.5, 5.5)


################################################################################
################################ FIGURE 7 ######################################
################################################################################
# Coverage of the 95% log-log confidence interval against N*/N (log scale)
figure7 <- plot_coverage_N(resM, outdir = tempdir())
if( is.list(figure7) && !inherits(figure7, "ggplot") ) figure7 <- figure7[[1]]
print(figure7)
save_fig(figure7, 7, 8.5, 5.5)