# Tests for runOLGA and generateOLGA functions

test_that("runOLGA validates chain argument", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  expect_error(runOLGA(immLynx_example, chains = "invalid"))
})

test_that("runOLGA infers correct model from organism and chain", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment")

  data("immLynx_example", package = "immLynx")

  # This should internally set model to "humanTRB"
  result <- runOLGA(immLynx_example, chains = "TRB", organism = "human")

  expect_s4_class(result, "Seurat")
})

test_that("runOLGA adds pgen columns to Seurat object", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment")

  data("immLynx_example", package = "immLynx")

  result <- runOLGA(immLynx_example, chains = "TRB")

  expect_true("olga_pgen_TRB" %in% names(result@meta.data))
  expect_true("olga_pgen_log10_TRB" %in% names(result@meta.data))
})

test_that("runOLGA returns data.frame when return_seurat=FALSE", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment")

  data("immLynx_example", package = "immLynx")

  result <- runOLGA(immLynx_example, chains = "TRB", return_seurat = FALSE)

  expect_s3_class(result, "data.frame")
  expect_true("barcode" %in% names(result))
  expect_true("pgen" %in% names(result))
  expect_true("log10_pgen" %in% names(result))
})

test_that("runOLGA handles custom column name", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment")

  data("immLynx_example", package = "immLynx")

  result <- runOLGA(immLynx_example, chains = "TRB",
                    column_name = "custom_pgen")

  expect_true("custom_pgen_TRB" %in% names(result@meta.data))
})

test_that("runOLGA stops for invalid model", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  expect_error(runOLGA(immLynx_example, chains = "TRB",
                       organism = "invalid"))
})

test_that("generateOLGA validates model argument", {
  expect_error(generateOLGA(n = 10, model = "invalid"))
})

test_that("generateOLGA returns correct number of sequences", {
  skip("Requires Python environment")

  result <- generateOLGA(n = 10, model = "humanTRB")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 10)
  expect_true("aa_seq" %in% names(result))
})
