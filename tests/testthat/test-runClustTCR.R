# Tests for runClustTCR function

test_that("runClustTCR validates chain argument", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  expect_error(runClustTCR(immLynx_example, chains = "invalid"))
})

test_that("runClustTCR adds cluster column to Seurat object", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment")

  data("immLynx_example", package = "immLynx")

  result <- runClustTCR(immLynx_example, chains = "TRB", method = "mcl")

  expect_s4_class(result, "Seurat")
  expect_true("clustcr_TRB" %in% names(result@meta.data))
})

test_that("runClustTCR returns data.frame when return_seurat=FALSE", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment")

  data("immLynx_example", package = "immLynx")

  result <- runClustTCR(immLynx_example, chains = "TRB", return_seurat = FALSE)

  expect_s3_class(result, "data.frame")
  expect_true("barcode" %in% names(result))
  expect_true("cluster" %in% names(result))
})

test_that("runClustTCR handles custom column prefix", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment")

  data("immLynx_example", package = "immLynx")

  result <- runClustTCR(immLynx_example, chains = "TRB",
                        column_prefix = "custom")

  expect_true("custom_TRB" %in% names(result@meta.data))
})

test_that("runClustTCR handles both chains separately", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment")

  data("immLynx_example", package = "immLynx")

  result <- runClustTCR(immLynx_example, chains = "both",
                        combine_chains = FALSE)

  expect_true("clustcr_TRA" %in% names(result@meta.data) ||
              "clustcr_TRB" %in% names(result@meta.data))
})

test_that("runClustTCR handles combined chains", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment")

  data("immLynx_example", package = "immLynx")

  result <- runClustTCR(immLynx_example, chains = "both",
                        combine_chains = TRUE)

  expect_true("clustcr_combined" %in% names(result@meta.data))
})

test_that("runClustTCR stops when no valid sequences found", {
  skip_if_not_installed("Seurat")

  # Create a mock Seurat object with no TCR data
  mock_seurat <- Seurat::CreateSeuratObject(
    counts = matrix(1:100, nrow = 10),
    min.cells = 0,
    min.features = 0
  )

  expect_error(runClustTCR(mock_seurat, chains = "TRB"))
})
