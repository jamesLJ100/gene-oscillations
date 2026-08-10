# common/gsea_functions.R itself library()s clusterProfiler/enrichplot/msigdbr/
# dorothea/TFEA.ChIP/ExperimentHub - not installed in every environment (see
# helper-setup.R). dplyr/ggplot2 are attached for real first since the functions
# under test genuinely need them at runtime; library() is then stubbed to a
# no-op just long enough to source the file without erroring on the missing
# packages, then restored.
#
# The stub and the source() call both have to touch .GlobalEnv explicitly
# (not just assign in this file's own top-level environment): testthat's
# test_file() runs a test file's top-level code in an isolated per-file
# environment, but source()'s default (local = FALSE) always evaluates the
# sourced code in literal .GlobalEnv - so an ordinary `library <- function(...)`
# assignment here wouldn't be visible to gsea_functions.R's own library() calls.
real_library <- base::library
assign("library", function(pkg, ...) invisible(NULL), envir = .GlobalEnv)
source(file.path(proj_root, "common/gsea_functions.R"))
assign("library", real_library, envir = .GlobalEnv)
rm(real_library)

test_that("resolve_org_db rejects an invalid species before touching any annotation package", {
  expect_error(resolve_org_db("cat"), "species must be 'human' or 'mouse'")
})

test_that("resolve_org_db returns the correct OrgDb for a valid species", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("org.Mm.eg.db")
  expect_s4_class(resolve_org_db("human"), "OrgDb")
  expect_s4_class(resolve_org_db("mouse"), "OrgDb")
})

test_that("build_ranked_gene_list converts symbols to Entrez IDs and sorts descending", {
  # bitr() comes from clusterProfiler, which isn't installed here. build_ranked_
  # gene_list() was source()'d into .GlobalEnv, so (as with run_pathway_analyses
  # below) the mock has to be assigned there directly, not just defined inside
  # this test_that() block, to be visible to it.
  bitr_mock <- function(x, fromType, toType, OrgDb) {
    data.frame(SYMBOL = x, ENTREZID = paste0("EID", seq_along(x)), stringsAsFactors = FALSE)
  }
  had_bitr <- exists("bitr", envir = .GlobalEnv, inherits = FALSE)
  old_bitr <- if (had_bitr) get("bitr", envir = .GlobalEnv) else NULL
  assign("bitr", bitr_mock, envir = .GlobalEnv)
  on.exit(if (had_bitr) assign("bitr", old_bitr, envir = .GlobalEnv) else rm("bitr", envir = .GlobalEnv))

  fake_df <- data.frame(symbol = c("MYC", "TP53", "EGFR", "BRCA1", "GAPDH"),
                         score  = c(0.9, 0.1, 0.5, 0.7, 0.05))
  gl <- build_ranked_gene_list(fake_df, org_db = "fake_org_db")

  expect_length(gl, 5)
  expect_true(all(diff(gl) <= 0))
  expect_equal(unname(gl[1]), 0.9)
  expect_true(all(grepl("^EID", names(gl))))
})

test_that("plot_hallmark_heatmap filters to pathways significant in at least one group", {
  # testthat::test_file() doesn't carry attached-package state from the file's
  # top-level code into individual test_that() blocks - re-attach here (see the
  # note above the library-stubbing block near the top of this file).
  library(dplyr)
  library(ggplot2)

  out_dir <- tempfile()
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE))

  fake_hallmark <- list(
    clusterA = data.frame(ID = c("HALLMARK_MYC_TARGETS", "HALLMARK_APOPTOSIS"),
                           NES = c(2.1, -1.5), p.adjust = c(0.01, 0.2)),
    clusterB = data.frame(ID = c("HALLMARK_MYC_TARGETS", "HALLMARK_APOPTOSIS"),
                           NES = c(1.8, 1.2), p.adjust = c(0.03, 0.6)),
    clusterC = NULL
  )

  result <- plot_hallmark_heatmap(fake_hallmark, out_dir, "testalgo")

  expect_equal(nrow(result), 2)  # only MYC_TARGETS is significant anywhere -> 2 rows (one per cluster)
  expect_true(all(result$ID == "MYC_TARGETS"))
  expect_true(file.exists(file.path(out_dir, "testalgo_hallmark_heatmap.pdf")))
})

test_that("plot_hallmark_heatmap returns NULL when nothing is significant anywhere", {
  library(dplyr)
  library(ggplot2)

  out_dir <- tempfile()
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE))

  fake_hallmark <- list(clusterA = data.frame(ID = "HALLMARK_A", NES = 0.5, p.adjust = 0.9))
  expect_null(plot_hallmark_heatmap(fake_hallmark, out_dir, "testalgo"))
})

test_that("plot_hallmark_comparison merges two conditions and flags significance from either side", {
  library(dplyr)
  library(ggplot2)

  out_dir <- tempfile()
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE))

  hm_a <- data.frame(ID = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_APOPTOSIS"),
                      NES = c(2.5, 1.8, -1.2), p.adjust = c(0.001, 0.01, 0.2))
  hm_b <- data.frame(ID = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_MYC_TARGETS_V1"),
                      NES = c(1.9, 2.1), p.adjust = c(0.005, 0.02))

  res <- plot_hallmark_comparison(hm_a, hm_b, "wt_tnf", "ss_tnf", out_dir, "cyclum")

  expect_equal(nrow(res), 4)  # union of both ID sets
  expect_true(file.exists(file.path(out_dir, "hallmark_comparison_wt_tnf_vs_ss_tnf.pdf")))
  # a pathway present only in hm_b should have NES_a filled with 0, not NA
  myc_row <- res[res$ID == "MYC_TARGETS_V1", ]
  expect_equal(myc_row$NES_a, 0)
  expect_equal(myc_row$NES_b, 2.1)
  expect_true(myc_row$significant)
})

test_that("plot_hallmark_comparison returns NULL when both conditions are empty", {
  out_dir <- tempfile()
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE))

  expect_null(plot_hallmark_comparison(NULL, NULL, "x", "y", out_dir, "cyclum"))
})

test_that("run_pathway_analyses dispatches the right species string to each sub-function", {
  # run_pathway_analyses() was source()'d into .GlobalEnv, so its lexical scope
  # is .GlobalEnv - mocks defined inside this test_that() block (a child
  # environment) would NOT be found by it. Assign the mocks directly into
  # .GlobalEnv instead, with on.exit() removing them again so they don't leak
  # into other tests in this file.
  calls <- new.env()
  mocks <- list(
    resolve_org_db    = function(species) paste0("orgdb_", species),
    run_hallmark_gsea = function(geneList, species, figures_dir, fname, algorithm, hallmark = NULL) {
      calls$hallmark <- species; calls$hallmark_arg <- hallmark; "hallmark_result"
    },
    run_go_gsea = function(geneList, org_db, figures_dir, fname) {
      calls$go <- org_db; "go_result"
    },
    run_dorothea_gsea = function(geneList, species, org_db, figures_dir, fname, dorothea_term2gene = NULL) {
      calls$dorothea <- species; calls$dorothea_arg <- dorothea_term2gene; "dorothea_result"
    },
    run_tfea_chip = function(geneList, figures_dir, fname, species, tf_filter) {
      calls$tfea <- species; "tfea_result"
    }
  )
  originals <- mget(names(mocks), envir = .GlobalEnv)
  for (nm in names(mocks)) assign(nm, mocks[[nm]], envir = .GlobalEnv)
  on.exit(for (nm in names(originals)) assign(nm, originals[[nm]], envir = .GlobalEnv))

  out_dir <- tempfile()

  res <- run_pathway_analyses(c(A = 1), species = "human", figures_dir = out_dir,
                               fname = "f", algorithm = "cyclum")
  expect_equal(calls$hallmark, "Homo sapiens")
  expect_equal(calls$go, "orgdb_human")
  expect_equal(calls$dorothea, "human")
  expect_equal(calls$tfea, "human")
  expect_true(dir.exists(out_dir))
  expect_named(res, c("hallmark", "go", "dorothea", "tfea"))
  expect_false(is.null(res$go))
  expect_false(is.null(res$tfea))
  unlink(out_dir, recursive = TRUE)

  res2 <- run_pathway_analyses(c(A = 1), species = "mouse", figures_dir = tempfile(),
                                fname = "f2", algorithm = "scPrisma",
                                run_go = FALSE, run_tfea = FALSE)
  expect_equal(calls$hallmark, "Mus musculus")
  expect_null(res2$go)
  expect_null(res2$tfea)
  expect_null(calls$hallmark_arg)
  expect_null(calls$dorothea_arg)

  # Pre-built tables should be forwarded through untouched, not rebuilt.
  fake_hallmark <- data.frame(gs_name = "SET1", entrez_gene = "1")
  fake_dorothea <- data.frame(Geneset = "TF1", ENTREZID = "1")
  run_pathway_analyses(c(A = 1), species = "human", figures_dir = tempfile(),
                        fname = "f3", algorithm = "cyclum",
                        hallmark = fake_hallmark, dorothea_term2gene = fake_dorothea)
  expect_identical(calls$hallmark_arg, fake_hallmark)
  expect_identical(calls$dorothea_arg, fake_dorothea)
})
