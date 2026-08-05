source(file.path(proj_root, "synthetic/dyngen_utils.R"))

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

test_that("get_best_hyperparams falls back to defaults when no grid search CSV exists", {
  gs_root <- tempfile()
  dir.create(gs_root)
  on.exit(unlink(gs_root, recursive = TRUE))

  defaults <- list(encoder_width = c(30, 20), epochs = 500, learning_rate = 2e-4)
  hp <- get_best_hyperparams("cyclum", gs_root, 1000, 200, defaults)
  expect_identical(hp, defaults)
})

test_that("get_best_hyperparams falls back to defaults when the combination isn't in the CSV", {
  gs_root <- tempfile()
  dir.create(gs_root)
  on.exit(unlink(gs_root, recursive = TRUE))

  defaults <- list(encoder_width = c(30, 20), epochs = 500, learning_rate = 2e-4)
  write.csv(data.frame(n_cells = 50, n_genes = 200, encoder_width = "40, 30",
                        epochs = 300, learning_rate = 1e-4),
            file.path(gs_root, "best_hyperparams_cyclum.csv"), row.names = FALSE)

  hp <- get_best_hyperparams("cyclum", gs_root, 1000, 200, defaults)
  expect_identical(hp, defaults)
})

test_that("get_best_hyperparams returns the tuned values when the combination is found", {
  gs_root <- tempfile()
  dir.create(gs_root)
  on.exit(unlink(gs_root, recursive = TRUE))

  defaults <- list(encoder_width = c(30, 20), epochs = 500, learning_rate = 2e-4)
  write.csv(data.frame(n_cells = 50, n_genes = 200, encoder_width = "40, 30",
                        epochs = 300, learning_rate = 1e-4),
            file.path(gs_root, "best_hyperparams_cyclum.csv"), row.names = FALSE)

  hp <- get_best_hyperparams("cyclum", gs_root, 50, 200, defaults)
  expect_equal(hp$encoder_width, "40, 30")
  expect_equal(hp$epochs, 300)
  expect_equal(hp$learning_rate, 1e-4)
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
