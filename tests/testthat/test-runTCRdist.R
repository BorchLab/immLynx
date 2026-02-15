# Tests for runTCRdist function

test_that("runTCRdist validates chain argument", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  expect_error(runTCRdist(immLynx_example, chains = "invalid"))
})

test_that("runTCRdist returns correct structure", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment")

  data("immLynx_example", package = "immLynx")

  result <- runTCRdist(immLynx_example, chains = "beta")

  expect_type(result, "list")
  expect_true("distances" %in% names(result))
  expect_true("barcodes" %in% names(result))
  expect_true("tcr_data" %in% names(result))
})

test_that("runTCRdist returns distance matrix of correct dimensions", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment")

  data("immLynx_example", package = "immLynx")

  result <- runTCRdist(immLynx_example, chains = "beta")

  n_cells <- length(result$barcodes)
  expect_equal(nrow(result$distances$pw_beta), n_cells)
  expect_equal(ncol(result$distances$pw_beta), n_cells)
})

test_that("runTCRdist handles add_to_object parameter", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment")

  data("immLynx_example", package = "immLynx")

  result <- runTCRdist(immLynx_example, chains = "beta",
                       add_to_object = TRUE)

  expect_s4_class(result, "Seurat")
  expect_true(!is.null(result@misc$tcrdist))
})

test_that("runTCRdist handles both chains", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment")

  data("immLynx_example", package = "immLynx")

  result <- runTCRdist(immLynx_example, chains = c("alpha", "beta"))

  expect_true("pw_alpha" %in% names(result$distances) ||
              "pw_beta" %in% names(result$distances))
})

test_that("runTCRdist validates organism argument", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment")

  data("immLynx_example", package = "immLynx")

  # Should work with valid organism
  expect_no_error(
    runTCRdist(immLynx_example, chains = "beta", organism = "human")
  )
})
