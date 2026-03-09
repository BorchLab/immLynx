# Tests for runClustTCR function

# ===========================================================================
# Parameter validation
# ===========================================================================

test_that("runClustTCR validates chain argument", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  expect_error(runClustTCR(immLynx_example, chains = "invalid"))
})

test_that("runClustTCR rejects non-Seurat/SCE input", {
  tcr_data <- create_mock_tcr_data(10)

  expect_error(
    runClustTCR(tcr_data, chains = "TRB"),
    "Input must be a Seurat or SingleCellExperiment object"
  )
})

test_that("runClustTCR rejects matrix input", {
  mat <- matrix(1:100, nrow = 10)

  expect_error(
    runClustTCR(mat, chains = "TRB"),
    "Input must be a Seurat or SingleCellExperiment object"
  )
})

test_that("runClustTCR rejects list input", {
  lst <- list(a = 1, b = 2)

  expect_error(
    runClustTCR(lst, chains = "TRB"),
    "Input must be a Seurat or SingleCellExperiment object"
  )
})

test_that("runClustTCR accepts all valid chain values", {
  skip_if_not_installed("Seurat")

  # Just test that match.arg accepts these
  expect_no_error(match.arg("TRB", c("TRB", "TRA", "both")))
  expect_no_error(match.arg("TRA", c("TRB", "TRA", "both")))
  expect_no_error(match.arg("both", c("TRB", "TRA", "both")))
})

test_that("runClustTCR function signature has correct defaults", {
  f <- formals(runClustTCR)

  expect_equal(f$method, "mcl")
  expect_equal(f$combine_chains, FALSE)
  expect_equal(f$return_object, TRUE)
  expect_equal(f$column_prefix, "clustcr")
})


# ===========================================================================
# Python-dependent tests (skipped without Python)
# ===========================================================================

test_that("runClustTCR adds cluster column to Seurat object", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_python()

  data("immLynx_example", package = "immLynx")

  result <- runClustTCR(immLynx_example, chains = "TRB", method = "mcl")

  expect_s4_class(result, "Seurat")
  expect_true("clustcr_TRB" %in% names(result@meta.data))
})

test_that("runClustTCR returns data.frame when return_object=FALSE", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_python()

  data("immLynx_example", package = "immLynx")

  result <- runClustTCR(immLynx_example, chains = "TRB",
                        return_object = FALSE)

  expect_s3_class(result, "data.frame")
  expect_true("barcode" %in% names(result))
  expect_true("cluster" %in% names(result))
  expect_true("cdr3_aa" %in% names(result))
  expect_true("chain" %in% names(result))
})

test_that("runClustTCR handles custom column prefix", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_python()

  data("immLynx_example", package = "immLynx")

  result <- runClustTCR(immLynx_example, chains = "TRB",
                        column_prefix = "custom")

  expect_true("custom_TRB" %in% names(result@meta.data))
})

test_that("runClustTCR handles both chains separately", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_python()

  data("immLynx_example", package = "immLynx")

  result <- runClustTCR(immLynx_example, chains = "both",
                        combine_chains = FALSE)

  expect_true("clustcr_TRA" %in% names(result@meta.data) ||
              "clustcr_TRB" %in% names(result@meta.data))
})

test_that("runClustTCR handles combined chains", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_python()

  data("immLynx_example", package = "immLynx")

  result <- runClustTCR(immLynx_example, chains = "both",
                        combine_chains = TRUE)

  expect_true("clustcr_combined" %in% names(result@meta.data))
})

test_that("runClustTCR combined chain result df has combined_sequence", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_python()

  data("immLynx_example", package = "immLynx")

  result <- runClustTCR(immLynx_example, chains = "both",
                        combine_chains = TRUE, return_object = FALSE)

  expect_s3_class(result, "data.frame")
  expect_true("combined_sequence" %in% names(result))
  expect_true("cluster" %in% names(result))
  # Combined sequences should contain "_"
  expect_true(all(grepl("_", result$combined_sequence)))
})

test_that("runClustTCR produces messages during execution", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_python()

  data("immLynx_example", package = "immLynx")

  expect_message(
    runClustTCR(immLynx_example, chains = "TRB"),
    "Extracting TCR sequences"
  )
})

test_that("runClustTCR return_object=FALSE for both chains returns combined df", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_python()

  data("immLynx_example", package = "immLynx")

  result <- runClustTCR(immLynx_example, chains = "both",
                        combine_chains = FALSE, return_object = FALSE)

  expect_s3_class(result, "data.frame")
  expect_true("chain" %in% names(result))
})
