library(dplyr)
library(ROCR)

#' Discover per-(n_cells, n_genes)-combination subdirectories
#'
#' generate_datasets.R writes each combination's replicates to its own
#' "c<n_cells>g<n_genes>" subdirectory (under synthetic/data/dyngen/ for the
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
#'   (synthetic/data/dyngen/gridsearch/)
#' @param n_cells,n_genes The combination to look up
#' @param defaults Named list of fallback hyperparameter values, one element per
#'   hyperparameter this algorithm needs
#' @return Named list of hyperparameter values (from the grid search if a match is
#'   found, otherwise `defaults` unchanged)
get_best_hyperparams <- function(algorithm, gridsearch_root, n_cells, n_genes, defaults) {

  best_csv <- file.path(gridsearch_root, paste0("best_hyperparams_", algorithm, ".csv"))
  print_hp <- function(hp) {
    cat(sprintf("Hyperparameters for c%dg%d: %s\n", n_cells, n_genes,
                 paste(sprintf("%s=%s", names(hp), vapply(hp, paste, character(1), collapse = ",")),
                       collapse = ", ")))
  }

  if (!file.exists(best_csv)) {
    warning(sprintf("No grid search results found at %s - using default hyperparameters for c%dg%d.",
                     best_csv, n_cells, n_genes), call. = FALSE, immediate. = TRUE)
    print_hp(defaults)
    return(defaults)
  }

  best_df   <- read.csv(best_csv, stringsAsFactors = FALSE)
  match_row <- best_df[best_df$n_cells == n_cells & best_df$n_genes == n_genes, ]

  if (nrow(match_row) == 0) {
    warning(sprintf("No grid search result for c%dg%d in %s - using default hyperparameters.",
                     n_cells, n_genes, best_csv), call. = FALSE, immediate. = TRUE)
    print_hp(defaults)
    return(defaults)
  }

  # A hyperparameter might not be a logged column (e.g. a fixed, non-swept value the
  # grid search didn't record) - fall back to its own default individually rather than
  # failing the whole lookup.
  missing_cols <- setdiff(names(defaults), names(match_row))
  if (length(missing_cols) > 0) {
    warning(sprintf("%s not logged in %s for c%dg%d - using default value(s) for: %s.",
                     algorithm, best_csv, n_cells, n_genes, paste(missing_cols, collapse = ", ")),
            call. = FALSE, immediate. = TRUE)
  }

  hp <- setNames(lapply(names(defaults), function(param) {
    if (param %in% names(match_row)) match_row[[param]][1] else defaults[[param]]
  }), names(defaults))

  print_hp(hp)
  hp
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
    dplyr::filter(module %in% c("B", "C", "D"), hops <= 1) %>%
    pull(gene)

  gene_symbols <- rownames(sim$counts)
  eval_genes   <- gene_module_assignments %>% dplyr::filter(hops <= 1) %>% pull(gene)
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

#' Compute accuracy/precision/recall/F1 at a single operating point
#'
#' Binary (0/1) scores (e.g. Oscope's gene classification) are evaluated directly at
#' that decision, since there's no cutoff left to search once the algorithm has
#' already made a hard call - searching anyway can land on the degenerate "predict
#' everyone positive" cutoff and look artificially good under class imbalance.
#' Continuous scores use each (algorithm, file)'s own F1-maximising cutoff via ROCR,
#' since a fixed cutoff (e.g. 0.5) isn't comparable across algorithms whose scores
#' live on different scales.
#'
#' @param scores Numeric vector of predicted scores
#' @param labels Integer vector of ground-truth 0/1 labels, same length as scores
#' @return List with accuracy, precision, recall, f1 (NA_real_ where undefined) and
#'   `predicted` - the 0/1 predicted label at the cutoff actually used, same length
#'   and order as `scores`/`labels` (lets a caller classify each gene as a
#'   TP/FP/FN/TN using the exact same cutoff the metrics above are based on)
compute_threshold_metrics <- function(scores, labels) {
  if (all(scores %in% c(0, 1))) {
    pred_pos   <- scores == 1
    actual_pos <- labels == 1
    tp <- sum(pred_pos & actual_pos)
    fp <- sum(pred_pos & !actual_pos)
    fn <- sum(!pred_pos & actual_pos)
    tn <- sum(!pred_pos & !actual_pos)

    accuracy  <- (tp + tn) / length(scores)
    precision <- if (tp + fp > 0) tp / (tp + fp) else NA_real_
    recall    <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
    f1        <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0) {
      2 * precision * recall / (precision + recall)
    } else NA_real_

    return(list(accuracy = accuracy, precision = precision, recall = recall, f1 = f1,
                predicted = as.integer(pred_pos)))
  }

  pred_obj  <- prediction(scores, labels)
  f1_perf   <- performance(pred_obj, "f")
  acc_perf  <- performance(pred_obj, "acc")
  prec_perf <- performance(pred_obj, "prec")
  rec_perf  <- performance(pred_obj, "rec")

  best_idx <- which.max(f1_perf@y.values[[1]])

  predicted <- if (length(best_idx) > 0) {
    cutoff <- f1_perf@x.values[[1]][best_idx]
    as.integer(scores >= cutoff)
  } else {
    rep(NA_integer_, length(scores))
  }

  list(
    f1        = if (length(best_idx) > 0) f1_perf@y.values[[1]][best_idx]   else NA_real_,
    accuracy  = if (length(best_idx) > 0) acc_perf@y.values[[1]][best_idx]  else NA_real_,
    precision = if (length(best_idx) > 0) prec_perf@y.values[[1]][best_idx] else NA_real_,
    recall    = if (length(best_idx) > 0) rec_perf@y.values[[1]][best_idx]  else NA_real_,
    predicted = predicted
  )
}

propagate_module_assignments <- function(sim, context = NULL) {

  feature_net <- sim$model$feature_network
  
  tf_modules <- bind_rows(
    feature_net %>%
      dplyr::filter(!is.na(from_module)) %>%
      dplyr::select(gene = from, module = from_module),
    feature_net %>%
      dplyr::filter(!is.na(to_module)) %>%
      dplyr::select(gene = to, module = to_module)
  ) %>%
    distinct()
  
  target_edges <- feature_net %>%
    dplyr::filter(grepl("^Target", to)) %>%
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
      dplyr::filter(to %in% unresolved, from %in% names(gene_module)) %>%
      dplyr::mutate(
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