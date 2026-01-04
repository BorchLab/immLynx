# Tests for runEmbeddings function

test_that("runEmbeddings validates chain argument", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  expect_error(runEmbeddings(immLynx_example, chains = "invalid"))
})

test_that("runEmbeddings adds dimensional reduction to Seurat object", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment and model download")

  data("immLynx_example", package = "immLynx")

  result <- runEmbeddings(immLynx_example, chains = "TRB")

  expect_s4_class(result, "Seurat")
  expect_true("tcr_esm" %in% Seurat::Reductions(result))
})

test_that("runEmbeddings returns list when return_seurat=FALSE", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment and model download")

  data("immLynx_example", package = "immLynx")

  result <- runEmbeddings(immLynx_example, chains = "TRB",
                          return_seurat = FALSE)

  expect_type(result, "list")
  expect_true("embeddings" %in% names(result))
  expect_true("metadata" %in% names(result))
  expect_s3_class(result$embeddings, "matrix")
})

test_that("runEmbeddings handles custom reduction name", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment and model download")

  data("immLynx_example", package = "immLynx")

  result <- runEmbeddings(immLynx_example, chains = "TRB",
                          reduction_name = "custom_embed")

  expect_true("custom_embed" %in% Seurat::Reductions(result))
})

test_that("runEmbeddings handles both chains", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment and model download")

  data("immLynx_example", package = "immLynx")

  result <- runEmbeddings(immLynx_example, chains = "both",
                          return_seurat = FALSE)

  expect_true(all(result$metadata$chain %in% c("TRA", "TRB", "TRA+TRB")))
})

test_that("runEmbeddings validates pool argument", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  expect_error(runEmbeddings(immLynx_example, pool = "invalid"))
})
