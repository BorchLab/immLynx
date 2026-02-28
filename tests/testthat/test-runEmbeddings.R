# Tests for runEmbeddings function

# ===========================================================================
# Parameter validation
# ===========================================================================

test_that("runEmbeddings function exists and is exported", {
  expect_true(is.function(runEmbeddings))
})

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

test_that("runEmbeddings rejects matrix input", {
  mat <- matrix(1:10, nrow = 2)

  expect_error(
    runEmbeddings(mat, chains = "TRB"),
    "Input must be a Seurat or SingleCellExperiment object"
  )
})

test_that("runEmbeddings rejects list input", {
  lst <- list(a = 1, b = 2)

  expect_error(
    runEmbeddings(lst, chains = "TRB"),
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

test_that("runEmbeddings accepts valid chain values via match.arg", {
  expect_no_error(match.arg("TRB", c("TRB", "TRA", "both")))
  expect_no_error(match.arg("TRA", c("TRB", "TRA", "both")))
  expect_no_error(match.arg("both", c("TRB", "TRA", "both")))
})

test_that("runEmbeddings has all expected parameters", {
  f <- formals(runEmbeddings)

  expected_params <- c("input", "chains", "model_name", "pool",
                       "chunk_size", "reduction_name", "reduction_key",
                       "return_object")
  for (p in expected_params) {
    expect_true(p %in% names(f), info = paste("Missing parameter:", p))
  }
})


# ===========================================================================
# Python-dependent tests (require basilisk env with transformers + torch)
# ===========================================================================

test_that("runEmbeddings adds dimensional reduction to Seurat object", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_transformers()

  data("immLynx_example", package = "immLynx")

  result <- runEmbeddings(immLynx_example, chains = "TRB",
                          prefer_device = "cpu")

  expect_s4_class(result, "Seurat")
  expect_true("tcr_esm" %in% Seurat::Reductions(result))
})

test_that("runEmbeddings returns list when return_object=FALSE", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_transformers()

  data("immLynx_example", package = "immLynx")

  result <- runEmbeddings(immLynx_example, chains = "TRB",
                          return_object = FALSE,
                          prefer_device = "cpu")

  expect_type(result, "list")
  expect_true("embeddings" %in% names(result))
  expect_true("metadata" %in% names(result))
  expect_true("model_name" %in% names(result))
  expect_true("pool" %in% names(result))
  expect_true(is.matrix(result$embeddings))
  expect_s3_class(result$metadata, "data.frame")
})

test_that("runEmbeddings list output has correct model_name", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_transformers()

  data("immLynx_example", package = "immLynx")

  result <- runEmbeddings(immLynx_example, chains = "TRB",
                          return_object = FALSE,
                          prefer_device = "cpu")

  expect_equal(result$model_name, "facebook/esm2_t12_35M_UR50D")
  expect_equal(result$pool, "mean")
})

test_that("runEmbeddings handles custom reduction name", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_transformers()

  data("immLynx_example", package = "immLynx")

  result <- runEmbeddings(immLynx_example, chains = "TRB",
                          reduction_name = "custom_embed",
                          prefer_device = "cpu")

  expect_true("custom_embed" %in% Seurat::Reductions(result))
})

test_that("runEmbeddings handles both chains", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_transformers()

  data("immLynx_example", package = "immLynx")

  result <- runEmbeddings(immLynx_example, chains = "both",
                          return_object = FALSE,
                          prefer_device = "cpu")

  expect_true(all(result$metadata$chain %in% c("TRA", "TRB", "TRA+TRB")))
})

test_that("runEmbeddings adds chain metadata to Seurat object", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_transformers()

  data("immLynx_example", package = "immLynx")

  result <- runEmbeddings(immLynx_example, chains = "TRB",
                          prefer_device = "cpu")

  expect_true("tcr_esm_chain" %in% names(result@meta.data))
})

test_that("runEmbeddings embeddings have barcodes as rownames", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_transformers()

  data("immLynx_example", package = "immLynx")

  result <- runEmbeddings(immLynx_example, chains = "TRB",
                          return_object = FALSE,
                          prefer_device = "cpu")

  expect_true(!is.null(rownames(result$embeddings)))
  expect_equal(rownames(result$embeddings), result$metadata$barcode)
})

test_that("runEmbeddings produces loading model message", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_transformers()

  data("immLynx_example", package = "immLynx")

  expect_message(
    runEmbeddings(immLynx_example, chains = "TRB",
                  prefer_device = "cpu"),
    "Loading Hugging Face model"
  )
})

test_that("runEmbeddings produces sequence extraction message", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_transformers()

  data("immLynx_example", package = "immLynx")

  expect_message(
    runEmbeddings(immLynx_example, chains = "TRB",
                  prefer_device = "cpu"),
    "Extracting TCR sequences"
  )
})

test_that("runEmbeddings embedding matrix has correct dimensions", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_transformers()

  data("immLynx_example", package = "immLynx")

  result <- runEmbeddings(immLynx_example, chains = "TRB",
                          return_object = FALSE,
                          prefer_device = "cpu")

  # Number of rows should match metadata
  expect_equal(nrow(result$embeddings), nrow(result$metadata))
  # ESM-2 35M has 480 hidden dims
  expect_equal(ncol(result$embeddings), 480)
})

test_that("runEmbeddings metadata has expected columns", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_transformers()

  data("immLynx_example", package = "immLynx")

  result <- runEmbeddings(immLynx_example, chains = "TRB",
                          return_object = FALSE,
                          prefer_device = "cpu")

  expect_true("barcode" %in% names(result$metadata))
  expect_true("sequence" %in% names(result$metadata))
  expect_true("chain" %in% names(result$metadata))
})

test_that("runEmbeddings TRB-only chain info is correct", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_transformers()

  data("immLynx_example", package = "immLynx")

  result <- runEmbeddings(immLynx_example, chains = "TRB",
                          return_object = FALSE,
                          prefer_device = "cpu")

  expect_true(all(result$metadata$chain == "TRB"))
})

test_that("runEmbeddings with cls pooling works", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_transformers()

  data("immLynx_example", package = "immLynx")

  result <- runEmbeddings(immLynx_example, chains = "TRB",
                          pool = "cls",
                          return_object = FALSE,
                          prefer_device = "cpu")

  expect_true(is.matrix(result$embeddings))
  expect_equal(nrow(result$embeddings), nrow(result$metadata))
})

test_that("runEmbeddings with custom chunk_size works", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_transformers()

  data("immLynx_example", package = "immLynx")

  result <- runEmbeddings(immLynx_example, chains = "TRB",
                          chunk_size = 8,
                          return_object = FALSE,
                          prefer_device = "cpu")

  expect_true(is.matrix(result$embeddings))
  expect_true(nrow(result$embeddings) > 0)
})

test_that("runEmbeddings embedding values are finite", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_transformers()

  data("immLynx_example", package = "immLynx")

  result <- runEmbeddings(immLynx_example, chains = "TRB",
                          return_object = FALSE,
                          prefer_device = "cpu")

  expect_true(all(is.finite(result$embeddings)))
})

test_that("runEmbeddings Seurat object retains original cells", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_transformers()

  data("immLynx_example", package = "immLynx")
  original_ncells <- ncol(immLynx_example)

  result <- runEmbeddings(immLynx_example, chains = "TRB",
                          prefer_device = "cpu")

  # Should have same number of cells
  expect_equal(ncol(result), original_ncells)
})

test_that("runEmbeddings custom reduction key works", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_transformers()

  data("immLynx_example", package = "immLynx")

  result <- runEmbeddings(immLynx_example, chains = "TRB",
                          reduction_name = "my_esm",
                          reduction_key = "MYESM_",
                          prefer_device = "cpu")

  expect_true("my_esm" %in% Seurat::Reductions(result))
})
