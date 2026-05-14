################################################################################
############################## INTRODUCTION ####################################
################################################################################
# This file contains all code necessary to reproduce the figures from the main
# article. Note that the figures have been edited in Adobe Illustrator to make
# the legend more readable and visually appealing. None of the contents of the 
# figures was changed. 

################################################################################
################################ FIGURE 1 ######################################
################################################################################
source("code/compile.R")
colon_df$time <- colon_df$time/365.241

# Fit all the models
cox_adhere <- coxph( Surv(time, status) ~ adhere, data = colon_df)
cox_rx <- coxph( Surv(time, status) ~ rx, data = colon_df)

# Predict the effective sample size
km_adhere_0 <- survfit_n( survfit(Surv(time, status) ~ 1, data = colon_df[colon_df$adhere == 0,]) )
km_adhere_1 <- survfit_n( survfit(Surv(time, status) ~ 1, data = colon_df[colon_df$adhere == 1,]) )
km_rx_0 <- survfit_n( survfit(Surv(time, status) ~ 1, data = colon_df[colon_df$rx == "Obs",]) )
km_rx_1 <- survfit_n( survfit(Surv(time, status) ~ 1, data = colon_df[colon_df$rx == "Lev",]) )
km_rx_2 <- survfit_n( survfit(Surv(time, status) ~ 1, data = colon_df[colon_df$rx == "Lev+5FU",]) )
cox_adhere_0 <- survfit_n( survfit(cox_adhere, newdata=data.frame(adhere=0)), cox=cox_adhere)
cox_adhere_1 <- survfit_n( survfit(cox_adhere, newdata=data.frame(adhere=1)), cox=cox_adhere)
cox_rx_0 <- survfit_n( survfit(cox_rx, newdata=data.frame(rx="Obs")), cox=cox_rx)
cox_rx_1 <- survfit_n( survfit(cox_rx, newdata=data.frame(rx="Lev")), cox=cox_rx)
cox_rx_2 <- survfit_n( survfit(cox_rx, newdata=data.frame(rx="Lev+5FU")), cox=cox_rx)

color_values <- c("Obs"= "#2E7691", "Lev" = "#C2666B", "No Adhesion"= "#37293F", "Lev+5FU" = "#c6aa2c", "Adhesion" = "#9f84af")
linetype_values <- c("Cox" = "solid", "Kaplan-Meier" = "dashed")

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
  scale_x_continuous(name = NULL, labels = NULL,
                     breaks = seq(0, max(colon_df$time), by = 2)) +
  ggtitle("Model based on treatment") +
  theme(legend.position = "none") # Remove legend from individual plot

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
  scale_x_continuous(name = "Time (years)",
                     breaks = seq(0, max(colon_df$time), by = 2)) +
  theme(legend.position = "none")  # Remove legend from individual plot

plot3 <- ggplot() +
  geom_line(aes(x = cox_adhere_0$time, y = cox_adhere_0$surv*100, color = "No Adhesion", linetype = "Cox"), lwd = 1.3) +
  geom_line(aes(x = cox_adhere_1$time, y = cox_adhere_1$surv*100, color = "Adhesion", linetype = "Cox"), lwd = 1.3) +
  geom_line(aes(x = km_adhere_0$time, y = km_adhere_0$surv*100, color = "No Adhesion", linetype = "Kaplan-Meier"), lwd = 1.3) +
  geom_line(aes(x = km_adhere_1$time, y = km_adhere_1$surv*100, color = "Adhesion", linetype = "Kaplan-Meier"), lwd = 1.3) +
  scale_color_manual(values = color_values) +
  scale_linetype_manual(values = linetype_values) +
  theme_minimal() +
  scale_y_continuous(name = NULL, labels = NULL, limits = c(0, 100)) +  # Remove y-axis name and tick labels
  scale_x_continuous(name = NULL, labels = NULL,
                     breaks = seq(0, max(colon_df$time), by = 2)) +
  ggtitle("Model based on adhesion") +
  theme(legend.position = "none")  # Remove legend from individual plot

plot4 <- ggplot() +
  geom_line(aes(x = cox_adhere_0$time, y = cox_adhere_0$n.eff, color = "No Adhesion", linetype = "Cox"), lwd = 1.3) +
  geom_line(aes(x = cox_adhere_1$time, y = cox_adhere_1$n.eff, color = "Adhesion", linetype = "Cox"), lwd = 1.3) +
  geom_line(aes(x = km_adhere_0$time, y = km_adhere_0$n.eff, color = "No Adhesion", linetype = "Kaplan-Meier"), lwd = 1.3) +
  geom_line(aes(x = km_adhere_1$time, y = km_adhere_1$n.eff, color = "Adhesion", linetype = "Kaplan-Meier"), lwd = 1.3) +
  geom_hline(yintercept = 929, linetype = "dotted", color = "black", linewidth = 0.6) +
  scale_color_manual(values = color_values) +
  scale_linetype_manual(values = linetype_values) +
  theme_minimal() +
  scale_y_continuous(name = NULL, labels = NULL, limits = c(0, 1000)) +  # Remove y-axis name and tick labels
  scale_x_continuous(name = "Time (years)",
                     breaks = seq(0, max(colon_df$time), by = 2)) +
  theme(legend.position = "none")  # Remove legend from individual plot

# Combine the plots using patchwork and add a single legend on the right
figure1 <- (plot1 / plot2 | plot3 / plot4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# Print the combined plot
print(figure1)


################################################################################
################################ FIGURE 2 ######################################
################################################################################
source("code/compile.R")
colon_df$time <- colon_df$time/365.241

# Calculate percentiles for age
age_percentiles <- quantile(colon_df$age, probs = c(0.5, 0.75, 0.99), na.rm = TRUE)
nodes_percentiles <- round(quantile(colon_df$nodes, probs = c(0.5, 0.75, 0.99), na.rm = TRUE))
nodes_percentiles["99%"] <- 20  # Manually set to 25 as it's more extreme than the 99th percentile

# Fit Cox models for age and nodes
cox_age <- coxph(Surv(time, status) ~ age, data = colon_df)
cox_nodes <- coxph(Surv(time, status) ~ nodes, data = colon_df)

# Predict effective sample size for age percentiles
cox_age_median <- survfit_n(survfit(cox_age, newdata = data.frame(age = age_percentiles["50%"])), cox = cox_age)
cox_age_75 <- survfit_n(survfit(cox_age, newdata = data.frame(age = age_percentiles["75%"])), cox = cox_age)
cox_age_99 <- survfit_n(survfit(cox_age, newdata = data.frame(age = age_percentiles["99%"])), cox = cox_age)

# Predict survival probabilities for age percentiles
surv_age_median <- survfit(cox_age, newdata = data.frame(age = age_percentiles["50%"]))
surv_age_75 <- survfit(cox_age, newdata = data.frame(age = age_percentiles["75%"]))
surv_age_99 <- survfit(cox_age, newdata = data.frame(age = age_percentiles["99%"]))

# Predict effective sample size for nodes with specific values
cox_nodes_2 <- survfit_n(survfit(cox_nodes, newdata = data.frame(nodes = nodes_percentiles["50%"])), cox = cox_nodes)
cox_nodes_5 <- survfit_n(survfit(cox_nodes, newdata = data.frame(nodes = nodes_percentiles["75%"])), cox = cox_nodes)
cox_nodes_25 <- survfit_n(survfit(cox_nodes, newdata = data.frame(nodes = nodes_percentiles["99%"])), cox = cox_nodes)

# Predict survival probabilities for nodes with specific values
surv_nodes_2 <- survfit(cox_nodes, newdata = data.frame(nodes = 2))
surv_nodes_5 <- survfit(cox_nodes, newdata = data.frame(nodes = 5))
surv_nodes_25 <- survfit(cox_nodes, newdata = data.frame(nodes = 25))

# Define color values
color_values <- c("Common" = "#37293F", "Uncommon" = "#2E7691", "Rare" = "#C2666B")
linetype_values <- c("Effective Sample Size" = "solid", "Survival Probability" = "dotted")

# Create the plots for age
plot1 <- ggplot() +
  geom_line(aes(x = surv_age_median$time, y = surv_age_median$surv*100, color = "Common"), lwd = 1.3) +
  geom_line(aes(x = surv_age_75$time, y = surv_age_75$surv*100, color = "Uncommon"), lwd = 1.3) +
  geom_line(aes(x = surv_age_99$time, y = surv_age_99$surv*100, color = "Rare"), lwd = 1.3) +
  scale_color_manual(values = color_values) +
  scale_linetype_manual(values = linetype_values) +
  theme_minimal() +
  scale_y_continuous(
    name = "Predicted survival (%)",
    limits = c(0, 100)
  ) +
  scale_x_continuous(labels = NULL, name = NULL,
                     breaks = seq(0, max(colon_df$time), by = 2)) +
  ggtitle("Cox model based on age") +
  theme(legend.position = "none")  # Remove legend from individual plot

# Create the plots for nodes
plot2 <- ggplot() +
  geom_line(aes(x = surv_nodes_2$time, y = surv_nodes_2$surv*100, color = "Common"), lwd = 1.3) +
  geom_line(aes(x = surv_nodes_5$time, y = surv_nodes_5$surv*100, color = "Uncommon"), lwd = 1.3) +
  geom_line(aes(x = surv_nodes_25$time, y = surv_nodes_25$surv*100, color = "Rare"), lwd = 1.3) +
  scale_color_manual(values = color_values) +
  scale_linetype_manual(values = linetype_values) +
  theme_minimal() +
  scale_y_continuous(
    limits = c(0, 100),
    labels = NULL,
    name = NULL
  ) +  # Remove y-axis name and tick labels
  scale_x_continuous(labels = NULL, name = NULL,
                     breaks = seq(0, max(colon_df$time), by = 2)) +
  ggtitle( "Cox model based on number of nodes") +
  theme(legend.position = "none")  # Remove legend from individual plot

plot3 <- ggplot() +
  geom_line(aes(x = cox_age_median$time, y = cox_age_median$n.eff, color = "Common"), lwd = 1.3) +
  geom_line(aes(x = cox_age_75$time, y = cox_age_75$n.eff, color = "Uncommon"), lwd = 1.3) +
  geom_line(aes(x = cox_age_99$time, y = cox_age_99$n.eff, color = "Rare"), lwd = 1.3) +
  geom_hline(yintercept = 929, linetype = "dotted", color = "black", linewidth = 0.6) +
  scale_color_manual(values = color_values) +
  scale_linetype_manual(values = linetype_values) +
  theme_minimal() +
  scale_y_continuous(
    name = "Effective Sample Size",
    limits = c(0, 1000)
  ) +
  scale_x_continuous(name = "Time (years)",
                     breaks = seq(0, max(colon_df$time), by = 2)) +
  theme(legend.position = "none")  # Remove legend from individual plot

# Create the plots for nodes
plot4 <- ggplot() +
  geom_line(aes(x = cox_nodes_2$time, y = cox_nodes_2$n.eff, color = "Common"), lwd = 1.3) +
  geom_line(aes(x = cox_nodes_5$time, y = cox_nodes_5$n.eff, color = "Uncommon"), lwd = 1.3) +
  geom_line(aes(x = cox_nodes_25$time, y = cox_nodes_25$n.eff, color = "Rare"), lwd = 1.3) +
  geom_hline(yintercept = 929, linetype = "dotted", color = "black", linewidth = 0.6) +
  scale_color_manual(values = color_values) +
  scale_linetype_manual(values = linetype_values) +
  theme_minimal() +
  scale_y_continuous(
    limits = c(0, 1000),
    labels = NULL,
    name = NULL
  ) +  # Remove y-axis name and tick labels
  scale_x_continuous(name = "Time (years)",
                     breaks = seq(0, max(colon_df$time), by = 2)) +
  theme(legend.position = "none")  # Remove legend from individual plot


# Combine the plots using patchwork and add a single legend on the right
figure2 <- (plot1 / plot3 | plot2 / plot4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# Print the combined plot
print(figure2)


################################################################################
################################ FIGURE 3 ######################################
################################################################################
source("code/compile.R")

year <- 5
var1 <- "rx"; var2 <- "rx + adhere"; var3 <- "rx + adhere + nodes"; var4 <- "rx + adhere + nodes + age"
leg1 <- "Tx"; leg2 <- "Tx + Adh"; leg3 <- "Tx + Adh + Nodes"; leg4 <- "Tx + Adh + Nodes + Age"


# Single impute the nodes data
colon_df$nodes <- mice::complete(mice::mice(colon_df, m = 1, print = FALSE))$nodes

# Fit all cox models
cox_model1 <- coxph(as.formula(paste0( "Surv( time, status) ~ ", var1)), data = colon_df)
cox_model2 <- coxph(as.formula(paste0( "Surv( time, status) ~ ", var2)), data = colon_df)
cox_model3 <- coxph(as.formula(paste0( "Surv( time, status) ~ ", var3)), data = colon_df)
cox_model4 <- coxph(as.formula(paste0( "Surv( time, status) ~ ", var4)), data = colon_df)

# Create prediction df
out_df1 = out_df2 = out_df3 = out_df4 <- colon_df |>
  select( rx, age, adhere, nodes )
out_df1$model <- paste0( "(A) Surv( time, status) ~ ", var1)
out_df2$model <- paste0( "(B) Surv( time, status) ~ ", var2)
out_df3$model <- paste0( "(C) Surv( time, status) ~ ", var3)
out_df4$model <- paste0( "(D) Surv( time, status) ~ ", var4)


# Create SF object for correct time
sf1 <- summary( survfit(cox_model1, newdata = out_df1), time = year*365.25 )
sf2 <- summary( survfit(cox_model2, newdata = out_df2), time = year*365.25 )
sf3 <- summary( survfit(cox_model3, newdata = out_df3), time = year*365.25 )
sf4 <- summary( survfit(cox_model4, newdata = out_df4), time = year*365.25 )

# Save survival probabilities
out_df1$surv <- c(sf1$surv)
out_df2$surv <- c(sf2$surv)
out_df3$surv <- c(sf3$surv)
out_df4$surv <- c(sf4$surv)

# Calculate the effective sample size for the 5-year prediction
out_df1$n.eff <- (1-out_df1$surv)*out_df1$surv/(c(sf1$std.err)^2)
out_df2$n.eff <- (1-out_df2$surv)*out_df2$surv/(c(sf2$std.err)^2)
out_df3$n.eff <- (1-out_df3$surv)*out_df3$surv/(c(sf3$std.err)^2)
out_df4$n.eff <- (1-out_df4$surv)*out_df4$surv/(c(sf4$std.err)^2)

out_df <- rbind(out_df1, out_df2, out_df3, out_df4)
out_df$model <- factor(out_df$model, levels = c(paste0( "(A) Surv( time, status) ~ ", var1), 
                                                paste0( "(B) Surv( time, status) ~ ", var2),  
                                                paste0( "(C) Surv( time, status) ~ ", var3), 
                                                paste0( "(D) Surv( time, status) ~ ", var4)))

# Create a scatter plot
colors <- c("#37293F", "#2E7691",  "#9f84af",  "#c6aa2c")
names(colors) <- c(paste0( "(A) Surv( time, status) ~ ", var1), 
                   paste0( "(B) Surv( time, status) ~ ", var2),  
                   paste0( "(C) Surv( time, status) ~ ", var3), 
                   paste0( "(D) Surv( time, status) ~ ", var4))

figure3 <- ggplot(out_df) +
  geom_point( aes(x = surv, y = n.eff, color = model), data = out_df, alpha = 0.2, size = 2) +
  facet_wrap( ~model ) +
  scale_x_continuous(name = "Predicted 5-Year Survival Probability", limits = c(0, 1)) +
  scale_y_continuous(name = "Effective Sample Size", limits = c(0, 500)) +
  #ggtitle("Changes distribution of Eff. N. between models") +
  labs( color = "Legend") +
  scale_color_manual(
    values = colors
  ) +
  #theme_minimal() +
  theme(legend.position = "none") + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

# Print the plot
print(figure3)

LL_null <- cox_model2$loglik[1]   # null model loglik
n <- nrow(colon_df)

R2_CS <- function(LL_model) {
  1 - exp((2/n) * (LL_null - LL_model))
}

R2_base <- R2_CS(cox_model2$loglik[2])
R2_nodes <- R2_CS(cox_model3$loglik[2])
R2_age <- R2_CS(cox_model4$loglik[2])

delta_R2_nodes <- R2_nodes - R2_base
delta_R2_age <- R2_age - R2_nodes
print(paste0("Delta R^2 after addition of nodes is: ", round(delta_R2_nodes*100, 2), "%"))
print(paste0("Delta R^2 after addition of age is: ", round(delta_R2_age*100, 2), "%"))


################################################################################
################################ FIGURE 4 ######################################
################################################################################
source("code/compile.R")
require(surveillance)

# Seed for consistency
set.seed(1923)

scale = c(10)
shape = c(1) # shape 1 means exponential
maxtime = 100 # maxtime
hr = 2 # hazard ratio
n = 5000 # number per group

events <- simsurv( dist="weibull", 
                   gammas=shape, 
                   lambdas=scale[1]^(-shape), 
                   x=data.frame(tx = c(rep(0, n), rep(1, n))),
                   betas = c("tx" = log(hr)),
                   maxt=maxtime)
events$tx = c(rep(0, n), rep(1, n))
events$eventtime[8885] <- 2.997671  # handmatig ties weghalen
events$eventtime[6298] <- 5.668818  # werkt omdat seed vast staats
time <- sort(unique(events$eventtime))
time0 <- sort(unique(events[events$tx == 0, ]$eventtime))
time1 <- sort(unique(events[events$tx == 1, ]$eventtime))

# Fit the models
sf0 <- survfit(Surv(eventtime, status) ~ 1, data = events[events$tx == 0,])
sf1 <- survfit(Surv(eventtime, status) ~ 1, data = events[events$tx == 1,])
cox <- coxph(Surv(eventtime, status) ~ tx, data = events)

# Calculate effective N and stuff
sf_n0 <- survfit_n(sf0) |> sf_to_df(time = time0)
sf_n1 <- survfit_n(sf1) |> sf_to_df(time = time1)
cox_n0 <- survfit_n(survfit(cox, newdata = data.frame(tx = 0)), cox, coef=T, chaz=T)
cox_n1 <- survfit_n(survfit(cox, newdata = data.frame(tx = 1)), cox, coef=T, chaz=T)

cox_n0 <- sf_to_df(cox_n0, time = time)
cox_n1 <- sf_to_df(cox_n1, time = time)

colors <- c("1" = "#37293F", "3" = "#2E7691",  "2" = "#9f84af", "4" = "#c6aa2c")
linetype <- c("1" = "solid", "2" = "dashed")

# plot the effective sample sizes all in one plot
plot1 <- ggplot() +
  geom_line(aes(x = time, y = cox_n0$n.eff, color = "1"), size = 1.2) +
  geom_line(aes(x = time, y = cox_n1$n.eff, color = "2"), size = 1.2) +
  geom_line(aes(x = time1, y = sf_n0$n.eff, color = "3"), size = 1.2) +
  geom_line(aes(x = time1, y = sf_n1$n.eff, color = "4"), size = 1.2) +
  scale_x_continuous(name = NULL, labels = NULL, limits = c(0, 15)) +
  scale_y_continuous(name = "Effective Sample Size (N)", limits = c(0, 15000)) +
  #ggtitle("Effective Sample Size of Coxph and Survfit") +
  labs( color = "Legend", tag = "(B)") +
  scale_color_manual(
    values = colors,
    labels = c("Cox Women", "Cox Men", "KM Women", "KM Men")
  ) +
  theme_minimal() +
  theme(legend.position = "right") + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)))


# plot the two kaplan-meiers overlapping in one plot and the cox survival predictions for the 0 an 1 treatment patients
plot2 <- ggplot() +
  geom_line(aes(x = time0, y = sf_n0$surv*100, color = "3"), size = 1.2) +
  geom_line(aes(x = time1, y = sf_n1$surv*100, color = "4"), size = 1.2) +
  geom_line(aes(x = time, y = cox_n0$surv*100, color = "1"), size = 1.2) +
  geom_line(aes(x = time, y = cox_n1$surv*100, color = "2"), size = 1.2) +
  scale_x_continuous(name = "Time", limits = c(0, 15)) +
  scale_y_continuous(name = "Predicted survival (%)", limits = c(0, 100)) +
  #ggtitle("Survival of Coxph and Survfit") +
  labs( color = "Legend", tag = "(A)") +
  scale_color_manual(
    values = colors,
    labels = c("Cox Women", "Cox Men", "KM Women", "KM Men")
  ) +
  theme_minimal() +
  theme(legend.position = "none") + 
  guides(colour = guide_legend(override.aes = list(alpha = 1))) 

time <- sort(unique(c(cox_n0$time, cox_n1$time)))

cox_incremental_hazvar0 <- (c(cox_n0$std.chaz, 0)-c(0, cox_n0$std.chaz))*cox_n0$surv
cox_incremental_hazvar1 <- (c(cox_n1$std.chaz, 0)-c(0, cox_n1$std.chaz))*cox_n0$surv
cox_incremental_coefvar0 <- (c(cox_n0$std.coef, 0)-c(0, cox_n0$std.coef))*cox_n0$surv
cox_incremental_coefvar1 <- (c(cox_n1$std.coef, 0)-c(0, cox_n1$std.coef))*cox_n0$surv
colors <- c("1" = "#37293F", "3" = "#2E7691",  "2" = "#9f84af", "4" = "#c6aa2c")
linetype <- c("1" = "solid", "2" = "dashed", "3" = "dotted")

scaleFactor = 0.00003/1

plot3 <- ggplot() +
  geom_line(aes(x = time, y = cox_n0$std.chaz^2, color = "1", lty = "1"), size = 1, lwd = 1.2) +
  geom_line(aes(x = time, y = cox_n1$std.chaz^2, color = "2", lty = "1"), size = 1, lwd = 1.2) +
  geom_line(aes(x = time, y = cox_n0$std.err^2, color = "1", lty = "2"), size = 1, lwd = 1.2) +
  geom_line(aes(x = time, y = cox_n1$std.err^2, color = "2", lty = "2"), size = 1, lwd = 1.2) +
  scale_x_continuous(name = "Time", limits = c(0, 15)) +
  scale_y_continuous(name = "Variance of Cox model", limits = c(0, 0.003)) +
  #ggtitle("Increment of standard error") +
  labs( color = "Legend", tag = "(D)") +
  scale_color_manual(
    values = colors,
    labels = c("Women", "Men")
  ) +
  scale_linetype_manual(
    values = linetype,
    labels = c( "Baseline hazard", "Total variance" ) ) +
  theme_minimal() +
  theme(legend.position = "none") + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)))


plot4 <- ggplot() +
  geom_line(aes(x = time, y = cox_n0$std.chaz^2/cox_n0$std.err^2*100, color = "1", lty = "1"), size = 1, lwd = 1.2) +
  geom_line(aes(x = time, y = cox_n1$std.chaz^2/cox_n1$std.err^2*100, color = "2", lty = "1"), size = 1, lwd = 1.2) +
  scale_x_continuous(name = "Time", limits = c(0, 15)) +
  scale_y_continuous(name = "Percentage of variance from h0", limits = c(0, 100)) +
  
  #ggtitle("Proportion of variance from baseline hazard") +
  labs( color = "Legend", tag = "(C)") +
  scale_color_manual(
    values = colors,
    labels = c("Women", "Men")
  ) +
  scale_linetype_manual(
    values = linetype,
    labels = c( "Baseline hazard", "Coefficient" ) ) +
  theme_minimal() +  
  theme(legend.position = "right") + 
  guides(colour = "none")

# Combine the plots using patchwork and add a single legend on the right
figure4 <- ( (plot2 | plot1) / (plot4 | plot3) ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# Print the combined plot
print(figure4)

################################################################################
################################ FIGURE 5 ######################################
################################################################################
source("code/compile.R")
require(surveillance)

# Seed for consistency
set.seed(1923)

scale = c(10)
shape = c(1) # shape 1 means exponential
maxtime = 100 # maxtime
hr = 2 # hazard ratio
n = 10000 # number per group

x <- data.frame(x = rnorm(n, 0, 1))
events <- simsurv( dist="weibull", 
                   gammas=shape, 
                   lambdas=scale[1]^(-shape), 
                   x=x,
                   betas = c("x" = log(hr)),
                   maxt=maxtime)
events$x = x$x
time <- sort(unique(events$eventtime))

# Fit the models
cox <- coxph(Surv(eventtime, status) ~ x, data = events)

# Calculate effective N and stuff
cox_mean <- survfit_n(survfit(cox, newdata = data.frame(x = 0)), cox, coef=T, chaz=T)
cox_sd_plus <- survfit_n(survfit(cox, newdata = data.frame(x = 1)), cox, coef=T, chaz=T)
cox_sd_min <- survfit_n(survfit(cox, newdata = data.frame(x = -1)), cox, coef=T, chaz=T)

cox_mean <- sf_to_df(cox_mean, time = time)
cox_sd_plus <- sf_to_df(cox_sd_plus, time = time)
cox_sd_min <- sf_to_df(cox_sd_min, time = time)

colors <- c("1" = "#37293F", "3" = "#2E7691",  "2" ="#C2666B", "4" = "#c6aa2c")
linetype <- c("1" = "solid", "2" = "dashed")

# plot the effective sample sizes all in one plot
plot1 <- ggplot() +
  geom_line(aes(x = time, y = cox_mean$n.eff, color = "1"), size = 1.2) +
  geom_line(aes(x = time, y = cox_sd_plus$n.eff, color = "2"), size = 1.2) +
  geom_line(aes(x = time, y = cox_sd_min$n.eff, color = "3"), size = 1.2) +
  scale_x_continuous(name = NULL, labels = NULL, limits = c(0, 15)) +
  scale_y_continuous(name = "Effective Sample Size (N)", limits = c(0, 15000)) +
  #ggtitle("Effective Sample Size of Coxph and Survfit") +
  labs( color = "Legend", tag = "(B)") +
  scale_color_manual(
    values = colors,
    labels = c("Mean", "+1 SD", "-1 SD")
  ) +
  theme_minimal() +
  theme(legend.position = "right") + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)))


# plot the two kaplan-meiers overlapping in one plot and the cox survival predictions for the 0 an 1 treatment patients
plot2 <- ggplot() +
  geom_line(aes(x = time, y = cox_mean$surv*100, color = "1"), size = 1.2) +
  geom_line(aes(x = time, y = cox_sd_plus$surv*100, color = "2"), size = 1.2) +
  geom_line(aes(x = time, y = cox_sd_min$surv*100, color = "3"), size = 1.2) +
  scale_x_continuous(name = "Time", limits = c(0, 15)) +
  scale_y_continuous(name = "Predicted survival (%)", limits = c(0, 100)) +
  #ggtitle("Survival of Coxph and Survfit") +
  labs( color = "Legend", tag = "(A)") +
  scale_color_manual(
    values = colors,
    labels = c("Mean", "+1 SD", "-1 SD")
  ) +
  theme_minimal() +
  theme(legend.position = "none") + 
  guides(colour = guide_legend(override.aes = list(alpha = 1))) 

scaleFactor = 0.00003/1

plot3 <- ggplot() +
  geom_line(aes(x = time, y = cox_mean$std.chaz^2, color = "1", lty = "1"), size = 1, lwd = 1.2) +
  geom_line(aes(x = time, y = cox_sd_plus$std.chaz^2, color = "2", lty = "1"), size = 1, lwd = 1.2) +
  geom_line(aes(x = time, y = cox_sd_min$std.chaz^2, color = "3", lty = "1"), size = 1, lwd = 1.2) +
  geom_line(aes(x = time, y = cox_mean$std.err^2, color = "1", lty = "2"), size = 1, lwd = 1.2) +
  geom_line(aes(x = time, y = cox_sd_plus$std.err^2, color = "2", lty = "2"), size = 1, lwd = 1.2) +
  geom_line(aes(x = time, y = cox_sd_min$std.err^2, color = "3", lty = "2"), size = 1, lwd = 1.2) +
  scale_x_continuous(name = "Time", limits = c(0, 15)) +
  scale_y_continuous(name = "Variance of Cox model", limits = c(0, 0.003)) +
  #ggtitle("Increment of standard error") +
  labs( color = "Legend", tag = "(D)") +
  scale_color_manual(
    values = colors,
    labels = c("Mean", "+1 SD", "-1 SD")
  ) +
  scale_linetype_manual(
    values = linetype,
    labels = c( "Baseline hazard", "Total variance" ) ) +
  theme_minimal() +
  theme(legend.position = "none") + 
  guides(colour = guide_legend(override.aes = list(alpha = 1)))


plot4 <- ggplot() +
  geom_line(aes(x = time, y = cox_mean$std.chaz^2/cox_mean$std.err^2*100, color = "1"), size = 1, lwd = 1.2) +
  geom_line(aes(x = time, y = cox_sd_plus$std.chaz^2/cox_sd_plus$std.err^2*100, color = "2"), size = 1, lwd = 1.2) +
  geom_line(aes(x = time, y = cox_sd_min$std.chaz^2/cox_sd_min$std.err^2*100, color = "3"), size = 1, lwd = 1.2) +
  scale_x_continuous(name = "Time", limits = c(0, 15)) +
  scale_y_continuous(name = "Percentage of variance from h0", limits = c(0, 100)) +
  
  #ggtitle("Proportion of variance from baseline hazard") +
  labs( color = "Legend", tag = "(C)") +
  scale_color_manual(
    values = colors,
    labels = c("Mean", "+1 SD", "-1 SD")
  ) +
  theme_minimal() +  
  theme(legend.position = "right") + 
  guides(colour = "none")

# Combine the plots using patchwork and add a single legend on the right
figure5 <- ( (plot2 | plot1) / (plot4 | plot3) ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# Print the combined plot
print(figure5)


################################################################################
################################ FIGURE 6 ######################################
################################################################################
source("code/compile.R")
require(surveillance)
require(future.apply)
require(parallel)

# Seed for consistency
set.seed(1923)

shape = c(0.5) # shape 1 means exponential
scale = c(lambda(0.25, shape), lambda(0.5, shape), lambda(0.75, shape), lambda(0.9, shape), lambda(0.95, shape)) # scale
maxtime = 10 # maxtime
hr = 2 # hazard ratio
n = 1000 # number per group
time_grid <- c(seq(0,1, by = 0.001), seq(1, 5, by = 0.01), seq(5, 10, by = 0.1)) # time grid for predictions
n_sim <- 100  # choose number of repeats

# 1. Define the plan (e.g., use all available cores)
plan(sequential) # Shuts down existing workers
gc()             # Force garbage collection
cores <- detectCores()
plan(multisession, workers = cores-1) # Re-initialize (try fewer workers to save RAM)
df_standard <- future_lapply(1:n_sim, function(r) {
  run_one(r, scale, shape, maxtime, hr, n, time_grid, model_type = "standard")
}, future.seed = TRUE) %>% 
  bind_rows()

# 1. Convert numeric proportions to matching strings
df_standard$prevalence <- paste0(as.numeric(1-df_standard$prevalence) * 100, "%")

# 2. Now apply the factor with matching levels
df_standard$prevalence <- factor(df_standard$prevalence, levels = c("75%", "50%", "25%", "10%", "5%"))

df_standard <- df_standard %>%
  group_by(tx, prevalence) %>%
  mutate(events = max(n.risk) - n.risk)
df_summary_std <- df_standard %>%
  group_by(time, tx, prevalence) %>%
  summarise(
    n.eff = mean(n.eff),
    surv = mean(surv),
    events = mean(events),
    .groups = "drop"
  )  %>%
  ungroup() %>%
  group_by(tx, prevalence) %>%
  mutate(events_rate = events / max(events),
         surv_prop = (1-surv)/max(1-surv)) %>%
  ungroup()

colors <- c("0" = "#37293F", "1" = "#9f84af")

linetypes <- c(
  "75%" = "twodash",
  "50%" = "solid",
  "25%" = "dashed",
  "10%" = "dotted",
  "5%" = "dotdash"
)

p_time <- ggplot(df_summary_std, aes(x = time, y = n.eff,
                                     color = tx,
                                     linetype = prevalence,
                                     group = interaction(tx, prevalence))) +
  geom_hline(yintercept = 1000, linetype = "dotted", color = "grey50") +
  geom_hline(yintercept = 500,  linetype = "dotted", color = "grey50") +
  geom_line(linewidth = 1) +
  scale_color_manual(values = colors, name = "Group") +
  scale_linetype_manual(values = linetypes, name = "Baseline event rate") +
  labs(x = "Time", y = "Effective sample size") +
  scale_x_continuous(limits = c(0, 10)) +
  scale_y_continuous(limits = c(0, 1600)) +
  theme_minimal()

p_event <- ggplot(df_summary_std, aes(x = events, y = n.eff,
                                      color = tx,
                                      linetype = prevalence,
                                      group = interaction(tx, prevalence))) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1000, linetype = "dotted", color = "grey50") +
  geom_hline(yintercept = 500,  linetype = "dotted", color = "grey50") +
  scale_color_manual(values = colors, name = "Group") +
  scale_linetype_manual(values = linetypes, name = "Baseline event rate") +
  labs(x = "Number of events") +
  scale_x_continuous(limits = c(0, 1000)) +
  scale_y_continuous(limits = c(0, 1600)) +
  theme_minimal() +
  theme( axis.text.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.title.y = element_blank())


p_event_prop <- ggplot(df_summary_std, aes(x = events_rate, y = n.eff,
                                           color = tx,
                                           linetype = prevalence,
                                           group = interaction(tx, prevalence))) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1000, linetype = "dotted", color = "grey50") +
  geom_hline(yintercept = 500,  linetype = "dotted", color = "grey50") +
  scale_color_manual(values = colors, name = "Group") +
  scale_linetype_manual(values = linetypes, name = "Baseline event rate") +
  labs(x = "Proportion of total events", y = "Effective sample size") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1600)) +
  theme_minimal()

# Standardize themes and remove redundant y-labels/titles if stacking
p_time <- p_time + labs(title = "A: Absolute Time")
p_event <- p_event + labs(title = "B: Absolute Event Count")
p_event_prop <- p_event_prop + labs(title = "C: Proportion of total events")

# 1. Define the layout: A and B on top, C on bottom left, Legend on bottom right
design <- "
  AABB
  #CCL
"

# 2. Combine using the 'guides = "collect"' argument
# guide_area() acts as a placeholder for the collected legend
figure6 <- p_time + p_event + p_event_prop + guide_area() + 
  plot_layout(design = design, guides = "collect")

print(figure6)


################################################################################
################################ FIGURE 7 ######################################
################################################################################
source("code/compile.R")
require(future.apply)
require(parallel)
set.seed(1923)

shape = c(0.5) # shape 1 means exponential
# make a vector of all lambda's from 0.01 to 0.99 in the scale parameter, so do lambda(x, shape) for all x in that range
scale =  sapply(c(seq(0.001, 0.10, by = 0.001), seq(0.10, 0.5, by = 0.01), seq(0.5, 1, 0.02)), function(x) lambda(x, shape)) # scale
maxtime = 10 # maxtime
hr = 2 # hazard ratio
n = 1000 # number per group
time_grid <- seq(0, 10, by = 0.1)
n_sim <- 100  # choose number of repeats

# 1. Define the plan (e.g., use all available cores)
plan(sequential) # Shuts down existing workers
gc()             # Force garbage collection
cores <- detectCores()
plan(multisession, workers = cores-1) # Re-initialize (try fewer workers to save RAM)
df <- future_lapply(1:n_sim, function(r) {
  run_one(r, scale, shape, maxtime, hr, n, time_grid, model_type = "standard")
}, future.seed = TRUE) %>% 
  bind_rows()

# Levels need to be from 1% to 99% in the prevalence variable as character
#df$prevalence <- factor(df$prevalence, levels = seq(0.01, 0.99, 0.01) )
df <- df %>%
  group_by(tx, prevalence) %>%
  mutate(events = max(n.risk) - n.risk)
df_summary <- df %>%
  group_by(time, tx, prevalence) %>%
  summarise(
    n.eff = mean(n.eff),
    surv = mean(surv),
    events = mean(events),
    .groups = "drop"
  )  %>%
  ungroup() %>%
  group_by(tx, prevalence) %>%
  mutate(events_rate = events / max(events),
         surv_prop = (1-surv)/max(1-surv)) %>%
  ungroup()
df_summary$prevalence = 1-df_summary$prevalence

max_points <- df_summary %>%
  group_by(tx, prevalence) %>%
  slice_max(order_by = 1-surv, n = 1) %>%
  ungroup()

library(ggnewscale) # For multiple color scales

figure7 <- ggplot(df_summary, aes(x = 1-surv, y = n.eff)) +
  # For tx == 0
  geom_line(data = subset(df_summary, tx == 0),
            aes(color = prevalence, group = interaction(tx, prevalence)),
            linewidth = 0.5, alpha = 0.2) +
  scale_color_gradient(high = "#a8d3e3", low = "#153744", name = "Prevalence (tx=0)") +
  # For tx == 1
  new_scale_color() +
  geom_line(data = subset(df_summary, tx == 1),
            aes(color = prevalence, group = interaction(tx, prevalence)),
            linewidth = 1, alpha = 0.2) +
  scale_color_gradient(high = "#dfafb2", low = "#8b373c", name = "Prevalence (tx=1)") +
  # Add solid line for max points (new color scale)
  new_scale_color() +
  geom_line(data = max_points,
            aes(x = 1-surv, y = n.eff, group = tx, color = factor(tx)),
            linewidth = 1, linetype = "solid") +
  scale_color_manual(values = c("0" = "#2E7691", "1" = "#C2666B"),
                     name = "Treatment") +
  labs(x = "Event incidence", y = "Effective sample size") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1600)) +
  theme_minimal()

print(figure7)