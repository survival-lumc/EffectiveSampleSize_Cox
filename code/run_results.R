################################################################################
############################## RUN ALL RESULTS #################################
################################################################################
# Run this file once (from the repository root). It runs every analysis that takes
# more than a few seconds and saves the output in results/. Both article_figures.R
# and supplement.qmd only read these files, so they compile in seconds.
# A step is skipped if its output file already exists: delete the file to rerun it.
# Expect hours for the simulation studies (Supplement 7); the bootstrap takes minutes.
#
# Output file                      Used in
# ess_realworld.rds                Figure 4, Supplement 6
# sim_main.rds                     Figures 6 and 7, Supplement 7
# sim_csens.rds                    Supplement 7 (sensitivity to Harrell's C)
# sim_guardrails.rds               Supplement 7 (number of predictors, noise predictors)
# sim_guardrails_shape.rds         Supplement 7 (Weibull shape)
# sim_guardrails_cens.rds          Supplement 7 (random censoring)
# sim_main_design.rds              Supplement 7 (scenario table)
# sim_early.rds                    Supplement 8
# sim_convergence.rds              Supplement 9
# other_results.rds                Supplement 12
# ess_bootstrap_colon.rds          Supplement 11
#
# If a cluster fails to start ("workers failed to connect"), set WORKERS <- 1 (sequential, slow)
# or lower the number of workers.

WORKERS <- max( 1, parallel::detectCores() - 1 )

source( "code/compile.R" )
dir.create( here::here("results"), showWarnings = FALSE, recursive = TRUE )

cached <- function( file, expr ){
  path <- here::here( "results", file )
  if( file.exists(path) ){
    message( "Skipped (exists): ", file )
    return( invisible(readRDS(path)) )
  }
  out <- expr
  saveRDS( out, path )
  message( "Saved: ", file )
  invisible( out )
}

################################################################################
# Real-world datasets (Figure 4, Supplement 6)
################################################################################
source( here::here("code", "ess_realworld.R") )
cached( "ess_realworld.rds", run_all() )

################################################################################
# Simulation studies (Figures 6 and 7, Supplement 7)
################################################################################
source( here::here("code", "sim_guardrails.R") )

# Main simulation: 7 event rates (1-99%) x 2 censoring settings, Harrell's C = 0.8 (14 scenarios)
cached( "sim_main.rds", run_all(make_scen("factorial"), workers = WORKERS) )

cached( "sim_main_design.rds", cens_summary(make_scen("factorial")) )


# Sensitivity analyses
cached( "sim_csens.rds",            run_all(make_scen("csens"),      workers = WORKERS) )
cached( "sim_guardrails.rds",       run_all(make_scen("oat"),        workers = WORKERS) )
cached( "sim_guardrails_shape.rds", run_all(make_scen("shape_only"), workers = WORKERS) )
cached( "sim_guardrails_cens.rds",  run_all(make_scen("censoring"),  workers = WORKERS) )

################################################################################
# Simulation 1: N* early in follow-up (Supplement 8; Figure 5 recomputes it in article_figures.R)
################################################################################
source( here::here("code", "sim_early.R") )
cached( "sim_early.rds", simulate_early() )

################################################################################
# Convergence of N* towards the Kaplan-Meier N* (Supplement 9)
################################################################################
source( here::here("code", "sim_convergence.R") )
cached( "sim_convergence.rds", run_convergence(workers = WORKERS) )

################################################################################
# Bootstrap (Supplement 11)
################################################################################
source( here::here("code", "bootstrap_ess.R") )
cached( "ess_bootstrap_colon.rds", run_colon_bootstrap(years = 1:7, B = 1000) )

################################################################################
# Other results (Supplement 12): GUSTO linear/logistic regression, stratified Cox models
################################################################################
source( here::here("code", "other_results.R") )
cached( "other_results.rds", run_other() )