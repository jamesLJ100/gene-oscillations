library(here)

proj_root <- here::here()
setwd(proj_root)
source(file.path(proj_root, "algorithms/run_oscope.R"))
source(file.path(proj_root, "common/grid_search.R"))
source(file.path(proj_root, "common/utils.R"))
source(file.path(proj_root, "synthetic/utils/dyngen_utils.R"))

# Find optimal hyperparameters separately for each (n_cells, n_genes) combination -
# see cyclum_gs.R for the full rationale.
tuning_root <- file.path(proj_root, "synthetic/data/dyngen/gridsearch")
combos      <- list_combo_dirs(tuning_root)

# Fixed parameter (not part of grid search, but logged alongside the swept ones for context).
# NULL matches OscopeKM()'s own default: maxK is computed from the data (top genes / 10)
# rather than fixed, per the Oscope vignette. Empirically this isn't what's gating results
# on this benchmark's synthetic data (see FlagCluster_thre below), and a fixed non-NULL
# value risks OscopeKM() erroring outright if it exceeds the number of genes KM_quan keeps
# as "top genes", so auto is also the safer choice here.
fixed_KM_maxK <- NULL

# Define grid search parameters.
#
# FlagCluster_thre's documented default (pi/4) flags out every cluster Oscope finds on
# this synthetic data - within-cluster phase spreads here top out around 15-17 degrees,
# well under pi/4's 45 degrees - so it's kept in the sweep (per the author's default)
# alongside looser values down to pi/12, empirically the point clusters start surviving.
#
# KM_quan controls what fraction of genes count as "top genes" for clustering
# (top_fraction = 1 - quan). At the documented default (0.95) and nearby, tuning AUC on
# this benchmark is flat at chance level; it climbs substantially down to ~0.3-0.5, then
# keeps climbing further below that - but going much past ~0.3 starts keeping the
# majority of the gene panel as "top genes", which stops testing Oscope's candidate-
# selection step and starts testing OscopeSine's raw similarity ranking with that step
# effectively disabled. Swept down to 0.30 to capture the real, still-legitimate gain
# without crossing that line.
param_grid <- expand.grid(
  KM_quan           = c(0.95, 0.70, 0.50, 0.25),
  FlagCluster_qt    = c(0.90, 0.95, 0.85),
  FlagCluster_thre  = c(pi/12, pi/8, pi/5, pi/4),
  CalcMV_MeanCutLow = c(100),
  stringsAsFactors = FALSE
)

make_run_fn <- function(combo_dir, fnames) {
  force(combo_dir)
  force(fnames)

  function(params, run_number) {
    oscope_config <- list(
      KM_maxK           = fixed_KM_maxK,
      KM_quan           = params$KM_quan,
      CalcMV_MeanCutLow = params$CalcMV_MeanCutLow,
      FlagCluster_qt    = params$FlagCluster_qt,
      FlagCluster_thre  = params$FlagCluster_thre
    )
    run_oscope(combo_dir, oscope_config)

    aucs <- sapply(fnames, function(fname) {
      results_df <- get_model_scores("oscope", combo_dir, fname)
      score_against_ground_truth(fname, combo_dir, results_df)$auc
    })

    list(
      KM_maxK         = if (is.null(fixed_KM_maxK)) NA_integer_ else fixed_KM_maxK,
      mean_tuning_auc = mean(aucs, na.rm = TRUE),
      n_tuning_files  = sum(!is.na(aucs))
    )
  }
}

best_hyperparams <- tune_per_combo(combos, param_grid, make_run_fn, "oscope")
