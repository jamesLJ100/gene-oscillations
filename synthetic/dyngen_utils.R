library(dplyr)
library(ROCR)

#' Discover per-(n_cells, n_genes)-combination subdirectories
#'
#' generate_datasets.R writes each combination's replicates to its own
#' "c<n_cells>g<n_genes>" subdirectory (under synthetic/data/dyngen_new/ for the
#' evaluation set, or .../gridsearch/ for the tuning set) rather than a flat
#' directory - required because run_cyclum.py/run_scPrisma.py each apply one
#' fixed hyperparameter setting to every .h5 file in whatever directory they're
#' pointed at, and different combinations need different (tuned) hyperparameters.
#' This discovers those subdirectories from disk instead of hardcoding the list
#' of combinations in every script that needs it.
#'
#' @param base_dir Directory containing "c<n_cells>g<n_genes>" subdirectories
#' @return Data frame with columns combo (subdirectory name), path, n_cells, n_genes
list_combo_dirs <- function(base_dir) {
  dirs <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
  dirs <- dirs[grepl("^c\\d+g\\d+$", dirs)]

  if (length(dirs) == 0) {
    return(data.frame(combo = character(0), path = character(0),
                       n_cells = integer(0), n_genes = integer(0)))
  }

  data.frame(
    combo     = dirs,
    path      = file.path(base_dir, dirs),
    n_cells   = as.integer(sub("^c(\\d+)g\\d+$", "\\1", dirs)),
    n_genes   = as.integer(sub("^c\\d+g(\\d+)$", "\\1", dirs)),
    stringsAsFactors = FALSE
  )
}

#' Look up the best hyperparameters found for one (n_cells, n_genes) combination
#'
#' Falls back to `defaults` if no grid search has been run yet for this
#' combination (or at all) - lets run_cyclum.R/run_scPrisma.R/run_oscope.R work
#' before any tuning has happened, using the same fixed values they always used to.
#'
#' @param algorithm Algorithm name ("cyclum", "scPrisma", "oscope") - looks up
#'   <gridsearch_root>/best_hyperparams_<algorithm>.csv (written by tune_per_combo())
#' @param gridsearch_root Directory containing that CSV
#'   (synthetic/data/dyngen_new/gridsearch/)
#' @param n_cells,n_genes The combination to look up
#' @param defaults Named list of fallback hyperparameter values, one element per
#'   hyperparameter this algorithm needs
#' @return Named list of hyperparameter values (from the grid search if a match is
#'   found, otherwise `defaults` unchanged)
get_best_hyperparams <- function(algorithm, gridsearch_root, n_cells, n_genes, defaults) {

  best_csv <- file.path(gridsearch_root, paste0("best_hyperparams_", algorithm, ".csv"))

  if (!file.exists(best_csv)) {
    cat("No grid search results found at", best_csv, "- using default hyperparameters.\n")
    return(defaults)
  }

  best_df   <- read.csv(best_csv, stringsAsFactors = FALSE)
  match_row <- best_df[best_df$n_cells == n_cells & best_df$n_genes == n_genes, ]

  if (nrow(match_row) == 0) {
    cat(sprintf("No grid search result for c%dg%d in %s - using default hyperparameters.\n",
                n_cells, n_genes, best_csv))
    return(defaults)
  }

  as.list(match_row[1, names(defaults), drop = FALSE])
}

#' Score an algorithm's cycling scores against dyngen's ground-truth cycling
#' gene labels for one simulated dataset
#'
#' Same label-building + AUC logic evaluate.R uses, factored out so the grid
#' search scripts can score each hyperparameter setting against the tuning-set
#' ground truth to select the best one per (n_cells, n_genes) combination,
#' rather than only logging exit codes.
#'
#' @param fname Dataset basename (e.g. "c1000g200_1"), used to load the
#'   matching "<fname>_sim.rds"
#' @param sim_dir Directory containing that _sim.rds file
#' @param results_df Data frame with `symbol` and `score` columns - the
#'   algorithm's output for this dataset
#' @return List with `auc` (NA_real_ if the sim file is missing or scoring fails)
score_against_ground_truth <- function(fname, sim_dir, results_df) {

  sim_file <- file.path(sim_dir, paste0(fname, "_sim.rds"))
  if (!file.exists(sim_file) || is.null(results_df)) return(list(auc = NA_real_))

  sim <- readRDS(sim_file)
  gene_module_assignments <- propagate_module_assignments(sim, context = fname)

  # Only evaluate on genes at most one hop from a TF (hops <= 1); this also
  # drops cycle genes with NA hops, since NA <= 1 is NA (filtered out).
  cycling_genes <- gene_module_assignments %>%
    filter(module %in% c("B", "C", "D"), hops <= 1) %>%
    pull(gene)

  gene_symbols <- rownames(sim$counts)
  eval_genes   <- gene_module_assignments %>% filter(hops <= 1) %>% pull(gene)
  eval_symbols <- gene_symbols[gene_symbols %in% eval_genes]

  label_df <- data.frame(symbol = eval_symbols, label = 0L)
  label_df$label[eval_symbols %in% cycling_genes] <- 1L

  merge_df <- merge(label_df, results_df, by = "symbol", all.x = TRUE)
  # If the algorithm filtered out a gene, treat it like a non-oscillating prediction
  merge_df$score[is.na(merge_df$score)] <- 0

  auc_val <- tryCatch({
    pred_obj <- prediction(merge_df$score, merge_df$label)
    performance(pred_obj, "auc")@y.values[[1]]
  }, error = function(e) NA_real_)

  list(auc = auc_val)
}

propagate_module_assignments <- function(sim, context = NULL) {

  feature_net <- sim$model$feature_network
  
  tf_modules <- bind_rows(
    feature_net %>%
      filter(!is.na(from_module)) %>%
      dplyr::select(gene = from, module = from_module),
    feature_net %>%
      filter(!is.na(to_module)) %>%
      dplyr::select(gene = to, module = to_module)
  ) %>%
    distinct()
  
  target_edges <- feature_net %>%
    filter(grepl("^Target", to)) %>%
    dplyr::select(from, to)

  target_indegree <- table(target_edges$to)
  overloaded <- names(target_indegree)[target_indegree > 1]
  if (length(overloaded) > 0) {
    prefix <- if (!is.null(context)) sprintf("[%s] ", context) else ""
    stop(sprintf(
      "%s%d target gene(s) have in-degree > 1, which this function assumes cannot happen (dyngen's max_in_degree should be 1): %s",
      prefix,
      length(overloaded),
      paste(overloaded, collapse = ", ")
    ))
  }

  gene_module <- setNames(tf_modules$module, tf_modules$gene)
  hops <- setNames(rep(0L, nrow(tf_modules)), tf_modules$gene)
  
  repeat {
    unresolved <- target_edges$to[!(target_edges$to %in% names(gene_module))]
    if (length(unresolved) == 0) break
    
    newly_resolved <- target_edges %>%
      filter(to %in% unresolved, from %in% names(gene_module)) %>%
      mutate(
        module = gene_module[from],
        hops   = hops[from] + 1L
      )
    
    if (nrow(newly_resolved) == 0) {
      unresolved <- unique(unresolved)
      prefix <- if (!is.null(context)) sprintf("[%s] ", context) else ""
      max_show <- 20L
      shown <- head(unresolved, max_show)
      target_list <- paste(shown, collapse = ", ")
      if (length(unresolved) > max_show) {
        target_list <- sprintf("%s, ... and %d more",
                               target_list, length(unresolved) - max_show)
      }
      warning(sprintf(
        "%s%d target(s) not traceable to a TF (cycle with no TF ancestor); labelling module/hops as NA: %s",
        prefix,
        length(unresolved),
        target_list
      ))
      gene_module <- c(gene_module, setNames(rep(NA_character_, length(unresolved)), unresolved))
      hops        <- c(hops,        setNames(rep(NA_integer_,   length(unresolved)), unresolved))
      break
    }

    gene_module <- c(gene_module, setNames(newly_resolved$module, newly_resolved$to))
    hops        <- c(hops,        setNames(newly_resolved$hops,   newly_resolved$to))
  }
  
  tibble(
    gene   = names(gene_module),
    module = unname(gene_module),
    hops   = unname(hops)
  )
}
# sim_file <- file.path(proj_root, "synthetic/data/dyngen", paste0("c50g207_1", "_sim.rds"))
# sim <- readRDS(sim_file)
# result <- propagate_module_assignments(sim)