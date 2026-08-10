source(file.path(proj_root, "synthetic/utils/dyngen_utils.R"))

# ============================================================================
# list_combo_dirs()
# ============================================================================

test_that("list_combo_dirs discovers c<n_cells>g<n_genes> subdirectories and parses them", {
  base_dir <- tempfile()
  dir.create(file.path(base_dir, "c50g200"), recursive = TRUE)
  dir.create(file.path(base_dir, "c1000g200"), recursive = TRUE)
  dir.create(file.path(base_dir, "notacombo"), recursive = TRUE)
  on.exit(unlink(base_dir, recursive = TRUE))

  combos <- list_combo_dirs(base_dir)

  expect_equal(nrow(combos), 2)
  expect_setequal(combos$combo, c("c50g200", "c1000g200"))
  expect_equal(combos$n_cells[combos$combo == "c50g200"], 50)
  expect_equal(combos$n_genes[combos$combo == "c50g200"], 200)
  expect_equal(combos$n_cells[combos$combo == "c1000g200"], 1000)
})

test_that("list_combo_dirs returns a zero-row data frame when there's nothing to find", {
  base_dir <- tempfile()
  dir.create(base_dir)
  on.exit(unlink(base_dir, recursive = TRUE))

  combos <- list_combo_dirs(base_dir)
  expect_equal(nrow(combos), 0)
  expect_named(combos, c("combo", "path", "n_cells", "n_genes"))
})

# ============================================================================
# get_best_hyperparams()
# ============================================================================

test_that("get_best_hyperparams falls back to defaults and warns when no grid search CSV exists", {
  gs_root <- tempfile()
  dir.create(gs_root)
  on.exit(unlink(gs_root, recursive = TRUE))

  defaults <- list(encoder_width = c(30, 20), epochs = 500, learning_rate = 2e-4)
  hp <- expect_warning(get_best_hyperparams("cyclum", gs_root, 1000, 200, defaults),
                        "No grid search results found")
  expect_identical(hp, defaults)
})

test_that("get_best_hyperparams falls back to defaults and warns when the combination isn't in the CSV", {
  gs_root <- tempfile()
  dir.create(gs_root)
  on.exit(unlink(gs_root, recursive = TRUE))

  defaults <- list(encoder_width = c(30, 20), epochs = 500, learning_rate = 2e-4)
  write.csv(data.frame(n_cells = 50, n_genes = 200, encoder_width = "40, 30",
                        epochs = 300, learning_rate = 1e-4),
            file.path(gs_root, "best_hyperparams_cyclum.csv"), row.names = FALSE)

  hp <- expect_warning(get_best_hyperparams("cyclum", gs_root, 1000, 200, defaults),
                        "No grid search result for c1000g200")
  expect_identical(hp, defaults)
})

test_that("get_best_hyperparams returns the tuned values, with no warning, when the combination is found", {
  gs_root <- tempfile()
  dir.create(gs_root)
  on.exit(unlink(gs_root, recursive = TRUE))

  defaults <- list(encoder_width = c(30, 20), epochs = 500, learning_rate = 2e-4)
  write.csv(data.frame(n_cells = 50, n_genes = 200, encoder_width = "40, 30",
                        epochs = 300, learning_rate = 1e-4),
            file.path(gs_root, "best_hyperparams_cyclum.csv"), row.names = FALSE)

  hp <- expect_no_warning(get_best_hyperparams("cyclum", gs_root, 50, 200, defaults))
  expect_equal(hp$encoder_width, "40, 30")
  expect_equal(hp$epochs, 300)
  expect_equal(hp$learning_rate, 1e-4)
})

test_that("get_best_hyperparams falls back to the default and warns for a hyperparameter missing from the CSV, without failing the whole lookup", {
  gs_root <- tempfile()
  dir.create(gs_root)
  on.exit(unlink(gs_root, recursive = TRUE))

  # regularisation_strength was swept and logged; iternum is a fixed value that
  # wasn't logged as its own column - this is the exact shape that used to crash
  # get_best_hyperparams() with "undefined columns selected".
  defaults <- list(regularisation_strength = 0.1, iternum = 100)
  write.csv(data.frame(n_cells = 50, n_genes = 200, regularisation_strength = 0.3),
            file.path(gs_root, "best_hyperparams_scPrisma.csv"), row.names = FALSE)

  hp <- expect_warning(get_best_hyperparams("scPrisma", gs_root, 50, 200, defaults),
                        "using default value.*iternum")
  expect_equal(hp$regularisation_strength, 0.3)
  expect_equal(hp$iternum, 100)
})

# ============================================================================
# propagate_module_assignments() and score_against_ground_truth()
#
# These both need dplyr/ROCR *attached* (they call filter()/pull()/bind_rows()/
# prediction()/performance() unqualified), which testthat::test_file() doesn't
# carry over from this file's top-level source() call into individual
# test_that() blocks - see test-gsea_functions.R for the fuller explanation.
# Each block below re-attaches what it needs.
#
# Fixture: a small fake dyngen-shaped feature network -
#   A -> B -> C -> D           (TF chain, modules A/B/C/D, hops 0)
#   B -> Target1                (module B, hops 1)
#   Orphan -> Target2           (Orphan is never a TF - unresolvable)
# ============================================================================

make_fake_feature_network <- function() {
  data.frame(
    from        = c("A", "B", "C", "B", "Orphan"),
    to          = c("B", "C", "D", "Target1", "Target2"),
    from_module = c("A", "B", "C", "B", NA),
    to_module   = c("B", "C", "D", NA, NA),
    stringsAsFactors = FALSE
  )
}

test_that("propagate_module_assignments assigns modules/hops through a TF chain", {
  library(dplyr)

  sim <- list(model = list(feature_network = make_fake_feature_network()))
  result <- suppressWarnings(propagate_module_assignments(sim))

  by_gene <- setNames(result$module, result$gene)
  hops_by_gene <- setNames(result$hops, result$gene)

  expect_equal(by_gene[["A"]], "A")
  expect_equal(by_gene[["B"]], "B")
  expect_equal(by_gene[["C"]], "C")
  expect_equal(by_gene[["D"]], "D")
  expect_equal(hops_by_gene[["A"]], 0L)
  expect_equal(by_gene[["Target1"]], "B")     # downstream of B, one hop
  expect_equal(hops_by_gene[["Target1"]], 1L)
})

test_that("propagate_module_assignments warns and assigns NA for genes with no traceable TF ancestor", {
  library(dplyr)

  sim <- list(model = list(feature_network = make_fake_feature_network()))
  expect_warning(
    result <- propagate_module_assignments(sim, context = "test-context"),
    "not traceable to a TF"
  )

  target2 <- result[result$gene == "Target2", ]
  expect_true(is.na(target2$module))
  expect_true(is.na(target2$hops))
})

test_that("propagate_module_assignments accumulates hops through a Target->Target chain", {
  library(dplyr)

  # A -> B -> Target1 -> Target2 -> Target3, each target reached via exactly
  # one parent (matches dyngen's max_in_degree = 1 for target/HK genes - see
  # run_dyngen.R).
  feature_net <- data.frame(
    from        = c("A", "B", "Target1", "Target2"),
    to          = c("B", "Target1", "Target2", "Target3"),
    from_module = c("A", "B", NA, NA),
    to_module   = c("B", NA, NA, NA),
    stringsAsFactors = FALSE
  )
  sim <- list(model = list(feature_network = feature_net))
  result <- propagate_module_assignments(sim)

  by_gene <- setNames(result$module, result$gene)
  hops_by_gene <- setNames(result$hops, result$gene)

  expect_equal(by_gene[["Target1"]], "B")
  expect_equal(hops_by_gene[["Target1"]], 1L)
  expect_equal(by_gene[["Target2"]], "B")
  expect_equal(hops_by_gene[["Target2"]], 2L)
  expect_equal(by_gene[["Target3"]], "B")
  expect_equal(hops_by_gene[["Target3"]], 3L)
})

test_that("propagate_module_assignments errors if a target gene has in-degree > 1", {
  library(dplyr)

  # TargetX has two incoming edges (from A and from B) - shouldn't be
  # reachable given dyngen's max_in_degree = 1 for target genes, but the
  # function should refuse to silently pick/merge one if it ever occurs.
  feature_net <- data.frame(
    from        = c("A", "B", "A", "B"),
    to          = c("B", "C", "TargetX", "TargetX"),
    from_module = c("A", "B", "A", "B"),
    to_module   = c("B", "C", NA, NA),
    stringsAsFactors = FALSE
  )
  sim <- list(model = list(feature_network = feature_net))

  expect_error(
    propagate_module_assignments(sim, context = "test-context"),
    "in-degree > 1"
  )
})

test_that("propagate_module_assignments handles a cyclic TF module network (matches backbone_cycle_simple()'s actual A->B->C->D->B shape)", {
  library(dplyr)

  # run_dyngen.R's real backbone has a feedback edge D->B, forming a cycle among
  # the TF modules themselves - unlike every other test here, which uses a plain
  # linear chain. This function identifies TF genes directly from module-labelled
  # edges rather than walking the TF graph, so the cycle shouldn't matter, but
  # that assumption was never actually exercised by a test until now.
  feature_net <- data.frame(
    from        = c("A", "B", "C", "D", "B"),
    to          = c("B", "C", "D", "B", "Target1"),
    from_module = c("A", "B", "C", "D", NA),
    to_module   = c("B", "C", "D", "B", NA),
    stringsAsFactors = FALSE
  )
  sim <- list(model = list(feature_network = feature_net))
  result <- propagate_module_assignments(sim)

  by_gene      <- setNames(result$module, result$gene)
  hops_by_gene <- setNames(result$hops, result$gene)

  expect_equal(by_gene[["A"]], "A")
  expect_equal(by_gene[["B"]], "B")
  expect_equal(by_gene[["C"]], "C")
  expect_equal(by_gene[["D"]], "D")
  expect_true(all(hops_by_gene[c("A", "B", "C", "D")] == 0L))
  expect_equal(by_gene[["Target1"]], "B")
  expect_equal(hops_by_gene[["Target1"]], 1L)
})

test_that("propagate_module_assignments produces sane output on a real dyngen simulation file", {
  library(dplyr)

  sim_path <- here::here("synthetic", "data", "dyngen", "c1000g200", "c1000g204_1_sim.rds")
  skip_if_not(file.exists(sim_path), "Real simulation fixture not available in this environment")

  sim <- readRDS(sim_path)
  result <- propagate_module_assignments(sim, context = "real-data-check")

  tf_rows <- result[result$gene %in% c("A_TF1", "B_TF1", "C_TF1", "D_TF1"), ]
  expect_equal(nrow(tf_rows), 4)
  expect_true(all(tf_rows$hops == 0))
  expect_setequal(tf_rows$module, c("A", "B", "C", "D"))

  # target_resampling = Inf / max_in_degree = 1 (see run_dyngen.R's make_config())
  # should mean every gene in the simulated count matrix is resolvable in practice.
  gene_symbols <- rownames(sim$counts)
  expect_true(all(gene_symbols %in% result$gene))
  expect_true(all(!is.na(result$module[result$gene %in% gene_symbols])))
})

test_that("score_against_ground_truth returns a high AUC when cycling genes score higher", {
  library(dplyr)
  library(ROCR)

  sim_dir <- tempfile()
  dir.create(sim_dir)
  on.exit(unlink(sim_dir, recursive = TRUE))

  # A clean subset of the shared fixture (A -> B -> C, B -> Target1), deliberately
  # excluding the "Orphan -> Target2" edge used elsewhere to test the
  # unresolvable-gene warning path - not relevant here and would just add noise.
  feature_network <- make_fake_feature_network()
  feature_network <- feature_network[feature_network$from != "Orphan", ]

  sim <- list(
    model = list(feature_network = feature_network),
    counts = matrix(1:20, nrow = 5,
                     dimnames = list(c("A", "B", "C", "Target1", "noise1"), paste0("cell", 1:4)))
  )
  saveRDS(sim, file.path(sim_dir, "c10g5_1_sim.rds"))

  # B, C, Target1 are all module B/C (cycling); A and noise1 aren't
  results_df <- data.frame(symbol = c("A", "B", "C", "Target1", "noise1"),
                            score  = c(0.1, 0.9, 0.9, 0.9, 0.05))
  res <- score_against_ground_truth("c10g5_1", sim_dir, results_df)

  expect_equal(res$auc, 1)
})

test_that("score_against_ground_truth returns NA when the sim file is missing or results are NULL", {
  sim_dir <- tempfile()
  dir.create(sim_dir)
  on.exit(unlink(sim_dir, recursive = TRUE))

  expect_true(is.na(score_against_ground_truth("nonexistent", sim_dir, data.frame(symbol = "a", score = 1))$auc))
  expect_true(is.na(score_against_ground_truth("c10g5_1", sim_dir, NULL)$auc))
})

# ============================================================================
# compute_threshold_metrics()
# ============================================================================

test_that("compute_threshold_metrics evaluates a binary score at its own decision, not a searched-for best cutoff", {
  library(ROCR)

  # Reproduces the exact bug found in production: an all-zero (binary) score against
  # a heavily imbalanced label set. Before the fix, searching for the F1-maximising
  # cutoff would land on "predict everyone positive" (since scores are tied at 0, one
  # of the tested cutoffs treats everyone as positive) and report near-perfect
  # accuracy/precision/recall/F1 despite the algorithm having predicted nothing.
  scores <- rep(0, 123)
  labels <- c(rep(1, 122), 0)

  m <- compute_threshold_metrics(scores, labels)

  # Nothing was predicted positive, so recall/precision/F1 are all zero (or NA for
  # precision, since there are no positive predictions to be precise about) - not the
  # ~0.99 a "predict everyone positive" cutoff would have produced.
  expect_equal(m$recall, 0)
  expect_true(is.na(m$precision))
  expect_true(is.na(m$f1))
  expect_equal(m$accuracy, 1 / 123)
})

test_that("compute_threshold_metrics evaluates a binary score correctly when it does predict some positives", {
  scores <- c(1, 1, 0, 0, 0)
  labels <- c(1, 0, 1, 0, 0)

  m <- compute_threshold_metrics(scores, labels)

  expect_equal(m$accuracy, 3 / 5)
  expect_equal(m$precision, 1 / 2)
  expect_equal(m$recall, 1 / 2)
  expect_equal(m$f1, 1 / 2)
  expect_equal(m$predicted, c(1, 1, 0, 0, 0))
})

test_that("compute_threshold_metrics searches for the best cutoff when scores are continuous", {
  library(ROCR)

  scores <- c(0.9, 0.8, 0.3, 0.2, 0.1)
  labels <- c(1, 1, 0, 0, 0)

  m <- compute_threshold_metrics(scores, labels)

  # A cutoff between 0.8 and 0.3 perfectly separates the classes.
  expect_equal(m$accuracy, 1)
  expect_equal(m$precision, 1)
  expect_equal(m$recall, 1)
  expect_equal(m$f1, 1)
  expect_equal(m$predicted, c(1, 1, 0, 0, 0))
})
