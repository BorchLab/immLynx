# Tests for runGLIPH function

test_that("runGLIPH checks for turboGliph package", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if(requireNamespace("turboGliph", quietly = TRUE))

  data("immLynx_example", package = "immLynx")

  expect_error(runGLIPH(immLynx_example),
               "Package 'turboGliph' is required")
})

test_that("runGLIPH validates chain argument", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_not_installed("turboGliph")

  data("immLynx_example", package = "immLynx")

  expect_error(runGLIPH(immLynx_example, chains = "invalid"))
})

test_that("runGLIPH validates method arguments", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_not_installed("turboGliph")

  data("immLynx_example", package = "immLynx")

  expect_error(runGLIPH(immLynx_example, local_method = "invalid"))
  expect_error(runGLIPH(immLynx_example, global_method = "invalid"))
})

test_that("runGLIPH adds cluster column to Seurat object", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_not_installed("turboGliph")

  data("immLynx_example", package = "immLynx")

  result <- runGLIPH(immLynx_example, chains = "TRB")

  expect_s4_class(result, "Seurat")
  expect_true("gliph_cluster" %in% names(result@meta.data))
})

test_that("runGLIPH returns list when return_seurat=FALSE", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_not_installed("turboGliph")

  data("immLynx_example", package = "immLynx")

  result <- runGLIPH(immLynx_example, chains = "TRB", return_seurat = FALSE)

  expect_type(result, "list")
  expect_true("clusters" %in% names(result))
  expect_true("parameters" %in% names(result))
})

test_that("runGLIPH handles custom column prefix", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_not_installed("turboGliph")

  data("immLynx_example", package = "immLynx")

  result <- runGLIPH(immLynx_example, chains = "TRB",
                     column_prefix = "custom")

  expect_true("custom_cluster" %in% names(result@meta.data))
})

test_that("scoreGLIPH validates input", {
  expect_error(scoreGLIPH(list(clusters = NULL)))
})

test_that("scoreGLIPH returns scored clusters", {
  skip_if_not_installed("turboGliph")

  # Create mock GLIPH results
  mock_results <- list(
    clusters = data.frame(
      cluster_id = 1:5,
      size = c(10, 5, 20, 3, 15),
      motif = c("ASS", "SYG", "AGT", "LGQ", "RDN"),
      stringsAsFactors = FALSE
    )
  )

  result <- scoreGLIPH(mock_results, min_size = 2)

  expect_s3_class(result, "data.frame")
  expect_true("quality_score" %in% names(result))
  expect_true(all(result$quality_score >= 0 & result$quality_score <= 1))
})

test_that("scoreGLIPH filters by minimum size", {
  mock_results <- list(
    clusters = data.frame(
      cluster_id = 1:5,
      size = c(1, 2, 3, 4, 5),
      stringsAsFactors = FALSE
    )
  )

  result <- scoreGLIPH(mock_results, min_size = 3)

  expect_equal(nrow(result), 3)
  expect_true(all(result$size >= 3))
})

test_that("extractGLIPHmotifs handles missing data", {
  expect_null(extractGLIPHmotifs(list(motifs = NULL)))
})

test_that("extractGLIPHmotifs filters by FDR", {
  mock_results <- list(
    motifs = data.frame(
      motif = c("ASS", "SYG", "AGT"),
      count = c(10, 5, 20),
      pvalue = c(0.001, 0.1, 0.05),
      stringsAsFactors = FALSE
    )
  )

  result <- extractGLIPHmotifs(mock_results, fdr_threshold = 0.1)

  # After FDR adjustment, check filtering works
  expect_true(all(result$fdr < 0.1))
})
