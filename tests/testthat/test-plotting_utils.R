source(file.path(proj_root, "common/plotting_utils.R"))

test_that("theme_common returns a ggplot theme, with angle_x_labels toggling x-axis rotation", {
  t1 <- theme_common()
  t2 <- theme_common(base_size = 20, angle_x_labels = TRUE)

  expect_s3_class(t1, "theme")
  expect_s3_class(t2, "theme")
  expect_null(t1$axis.text.x$angle)
  expect_equal(t2$axis.text.x$angle, 45)
})

test_that("plot_score_comparison identifies the most score-divergent genes", {
  out_dir <- tempfile()
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE))

  set.seed(42)
  genes <- paste0("Gene", 1:200)
  scores_a <- data.frame(symbol = genes, score = runif(200))
  scores_b <- data.frame(symbol = genes, score = c(runif(190), scores_a$score[191:200] + 5))

  res <- plot_score_comparison(scores_a, scores_b, "wt_tnf", "ss_tnf", out_dir, "cyclum")

  expect_equal(nrow(res), 200)
  expect_true(file.exists(file.path(out_dir, "score_comparison_wt_tnf_vs_ss_tnf.pdf")))
  top20 <- res[order(-res$score_diff), ][1:20, ]
  expect_true("Gene200" %in% top20$symbol)
})

test_that("plot_score_comparison returns NULL gracefully when data is missing", {
  out_dir <- tempfile()
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE))

  scores_a <- data.frame(symbol = c("a", "b"), score = c(1, 2))
  scores_c <- data.frame(symbol = c("x", "y"), score = c(1, 2))

  expect_null(plot_score_comparison(scores_a, NULL, "p", "q", out_dir, "cyclum"))
  expect_null(plot_score_comparison(scores_a, scores_c, "p", "q", out_dir, "cyclum"))  # no shared genes
  expect_null(plot_score_comparison(scores_a, data.frame(symbol = character(0), score = numeric(0)),
                                     "p", "q", out_dir, "cyclum"))
})

test_that("plot_score_vs_expression highlights reference oscillating genes and saves a plot", {
  out_dir <- tempfile()
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE))

  cluster_df <- data.frame(
    score = c(0.9, 0.1, 0.5, 0.8),
    mean_expr = c(10, 100, 5, 50),
    is_oscillating = c(TRUE, FALSE, TRUE, FALSE)
  )

  p <- plot_score_vs_expression(cluster_df, "pPSM", "cyclum", "Cyclum Score", out_dir)

  expect_s3_class(p, "ggplot")
  expect_true(file.exists(file.path(out_dir, "pPSM_score_vs_expression_cyclum.pdf")))
})

test_that("plot_score_vs_expression returns NULL for empty input", {
  out_dir <- tempfile()
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE))

  expect_null(plot_score_vs_expression(data.frame(), "empty", "cyclum", "score", out_dir))
})
