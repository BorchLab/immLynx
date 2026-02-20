# Tests for runEmbeddings function

# ===========================================================================
# Parameter validation
# ===========================================================================

test_that("runEmbeddings validates chain argument", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  expect_error(runEmbeddings(immLynx_example, chains = "invalid"))
})

test_that("runEmbeddings rejects non-Seurat/SCE input", {
  tcr_data <- create_mock_tcr_data(10)

  expect_error(
    runEmbeddings(tcr_data, chains = "TRB"),
    "Input must be a Seurat or SingleCellExperiment object"
  )
})

test_that("runEmbeddings rejects data.frame input", {
  df <- data.frame(a = 1:5)

  expect_error(
    runEmbeddings(df, chains = "TRB"),
    "Input must be a Seurat or SingleCellExperiment object"
  )
})

test_that("runEmbeddings function signature has correct defaults", {
  f <- formals(runEmbeddings)

  expect_equal(f$model_name, "facebook/esm2_t12_35M_UR50D")
  expect_equal(f$pool, "mean")
  expect_equal(f$chunk_size, 32)
  expect_equal(f$reduction_name, "tcr_esm")
  expect_equal(f$reduction_key, "ESM_")
  expect_equal(f$return_object, TRUE)
})

test_that("runEmbeddings accepts valid chain values", {
  expect_no_error(match.arg("TRB", c("TRB", "TRA", "both")))
  expect_no_error(match.arg("TRA", c("TRB", "TRA", "both")))
  expect_no_error(match.arg("both", c("TRB", "TRA", "both")))
})


# ===========================================================================
# Python-dependent tests (skipped without Python)
# ===========================================================================

test_that("runEmbeddings adds dimensional reduction to Seurat object", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_transformers()

  data("immLynx_example", package = "immLynx")

  result <- runEmbeddings(immLynx_example, chains = "TRB")

  expect_s4_class(result, "Seurat")
  expect_true("tcr_esm" %in% Seurat::Reductions(result))
})

test_that("runEmbeddings returns list when return_object=FALSE", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_transformers()

  data("immLynx_example", package = "immLynx")

  result <- runEmbeddings(immLynx_example, chains = "TRB",
                          return_object = FALSE)

  expect_type(result, "list")
  expect_true("embeddings" %in% names(result))
  expect_true("metadata" %in% names(result))
  expect_true("model_name" %in% names(result))
  expect_true("pool" %in% names(result))
  expect_true(is.matrix(result$embeddings))
  expect_s3_class(result$metadata, "data.frame")
})

test_that("runEmbeddings handles custom reduction name", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_transformers()

  data("immLynx_example", package = "immLynx")

  result <- runEmbeddings(immLynx_example, chains = "TRB",
                          reduction_name = "custom_embed")

  expect_true("custom_embed" %in% Seurat::Reductions(result))
})

test_that("runEmbeddings handles both chains", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_transformers()

  data("immLynx_example", package = "immLynx")

  result <- runEmbeddings(immLynx_example, chains = "both",
                          return_object = FALSE)

  expect_true(all(result$metadata$chain %in% c("TRA", "TRB", "TRA+TRB")))
})

test_that("runEmbeddings adds chain metadata to Seurat object", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_transformers()

  data("immLynx_example", package = "immLynx")

  result <- runEmbeddings(immLynx_example, chains = "TRB")

  expect_true("tcr_esm_chain" %in% names(result@meta.data))
})

test_that("runEmbeddings embeddings have barcodes as rownames", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_transformers()

  data("immLynx_example", package = "immLynx")

  result <- runEmbeddings(immLynx_example, chains = "TRB",
                          return_object = FALSE)

  expect_true(!is.null(rownames(result$embeddings)))
  expect_equal(rownames(result$embeddings), result$metadata$barcode)
})

test_that("runEmbeddings produces informative messages", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_transformers()

  data("immLynx_example", package = "immLynx")

  expect_message(
    runEmbeddings(immLynx_example, chains = "TRB"),
    "Loading Hugging Face model"
  )
})
