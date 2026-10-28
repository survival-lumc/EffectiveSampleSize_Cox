####################### CONVERGENCE OF N* TOWARDS THE KAPLAN-MEIER N* #######################

# Formerly Figure 6 of the article, now Supplement 9.
# Cox model with a single binary predictor (HR 2), N = 1000 (500 per group), Weibull shape 0.5, no censoring,
# for five baseline event rates. N* of the Cox model prediction over time for both groups, against absolute time,
# the number of events and the proportion of all events observed.
#
# Usage (from repository root, after source("code/compile.R")):
#   source("code/sim_convergence.R")
#   conv <- run_convergence()             # takes a while; saved by code/run_results.R
#   plot_convergence( conv )

library(future.apply)

run_convergence <- function( n_sim = 100, seed = 1923, workers = max(1, parallel::detectCores() - 1) ){
  set.seed(seed)
  shape    <- 0.5                                   # shape 1 means exponential
  surv_end <- c( 0.25, 0.5, 0.75, 0.9, 0.95 )       # baseline (X = 0) survival at the end of follow-up, t = 10
  scale    <- sapply( surv_end, lambda, shape = shape )
  maxtime  <- 10
  hr       <- 2
  n        <- 1000                                  # 500 per group
  assign( "time_grid", c(seq(0, 1, by = 0.001), seq(1, 5, by = 0.01), seq(5, 10, by = 0.1)), envir = globalenv() )   # used by extract_sf()
  
  if( workers > 1 ) future::plan( future::multisession, workers = workers ) else future::plan( future::sequential )
  df <- future_lapply( seq_len(n_sim), function( r ) {
    run_one( r, scale, shape, maxtime, hr, n, time_grid, model_type = "standard" )
  }, future.seed = TRUE, future.packages = c("survival", "simsurv", "dplyr", "tidyr") ) |>
    dplyr::bind_rows()
  future::plan( future::sequential )
  
  df$prevalence <- paste0( round((1 - as.numeric(df$prevalence)) * 100), "%" )
  df$prevalence <- factor( df$prevalence, levels = c("75%", "50%", "25%", "10%", "5%") )
  
  df <- df |>
    dplyr::group_by( tx, prevalence ) |>
    dplyr::mutate( events = max(n.risk) - n.risk )
  
  df |>
    dplyr::group_by( time, tx, prevalence ) |>
    dplyr::summarise( n.eff = mean(n.eff), surv = mean(surv), events = mean(events), .groups = "drop" ) |>
    dplyr::group_by( tx, prevalence ) |>
    dplyr::mutate( events_rate = events / max(events) ) |>
    dplyr::ungroup()
}

plot_convergence <- function( conv ){
  colors    <- c( "0" = "#37293F", "1" = "#9f84af" )
  linetypes <- c( "75%" = "twodash", "50%" = "solid", "25%" = "dashed", "10%" = "dotted", "5%" = "dotdash" )
  
  panel <- function( x, xlab, xlim, title ){
    ggplot( conv, aes(x = .data[[x]], y = n.eff, color = tx, linetype = prevalence,
                      group = interaction(tx, prevalence)) ) +
      geom_hline( yintercept = 1000, linetype = "dotted", color = "grey50" ) +
      geom_hline( yintercept = 500,  linetype = "dotted", color = "grey50" ) +
      geom_line( linewidth = 1 ) +
      scale_color_manual( values = colors, name = "Group", labels = c("X = 0", "X = 1") ) +
      scale_linetype_manual( values = linetypes, name = "Baseline event rate" ) +
      labs( x = xlab, y = "Effective sample size", title = title ) +
      scale_x_continuous( limits = xlim ) +
      scale_y_continuous( limits = c(0, 1600) ) +
      theme_minimal()
  }
  pA <- panel( "time", "Time", c(0, 10), "A: Absolute time" )
  pB <- panel( "events", "Number of events", c(0, 1000), "B: Absolute event count" ) +
    theme( axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank() )
  pC <- panel( "events_rate", "Proportion of total events", c(0, 1), "C: Proportion of total events" )
  
  design <- "
    AABB
    #CCL
  "
  pA + pB + pC + guide_area() + plot_layout( design = design, guides = "collect" )
}
