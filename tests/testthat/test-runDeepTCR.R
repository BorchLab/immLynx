# Tests for runDeepTCR function

test_that("runDeepTCR validates chain argument", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  expect_error(runDeepTCR(immLynx_example, chains = "invalid"))
})

test_that("runDeepTCR creates output directory", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment")

  data("immLynx_example", package = "immLynx")

  temp_dir <- tempfile("deeptcr_test")
  on.exit(unlink(temp_dir, recursive = TRUE))

  result <- runDeepTCR(immLynx_example, chains = "TRB",
                       output_dir = temp_dir)

  expect_true(dir.exists(temp_dir))
  expect_true(dir.exists(file.path(temp_dir, "data")))
})

test_that("runDeepTCR adds dimensional reduction to Seurat object", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment")

  data("immLynx_example", package = "immLynx")

  temp_dir <- tempfile("deeptcr_test")
  on.exit(unlink(temp_dir, recursive = TRUE))

  result <- runDeepTCR(immLynx_example, chains = "TRB",
                       output_dir = temp_dir)

  expect_s4_class(result, "Seurat")
  expect_true("deeptcr" %in% Seurat::Reductions(result))
})

test_that("runDeepTCR returns matrix when return_seurat=FALSE", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment")

  data("immLynx_example", package = "immLynx")

  temp_dir <- tempfile("deeptcr_test")
  on.exit(unlink(temp_dir, recursive = TRUE))

  result <- runDeepTCR(immLynx_example, chains = "TRB",
                       output_dir = temp_dir,
                       return_seurat = FALSE)

  expect_true(is.matrix(result) || is.data.frame(result))
})

test_that("runDeepTCR handles custom reduction name", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment")

  data("immLynx_example", package = "immLynx")

  temp_dir <- tempfile("deeptcr_test")
  on.exit(unlink(temp_dir, recursive = TRUE))

  result <- runDeepTCR(immLynx_example, chains = "TRB",
                       output_dir = temp_dir,
                       reduction_name = "custom_deeptcr")

  expect_true("custom_deeptcr" %in% Seurat::Reductions(result))
})

test_that("runDeepTCR respects latent_dim parameter", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment")

  data("immLynx_example", package = "immLynx")

  temp_dir <- tempfile("deeptcr_test")
  on.exit(unlink(temp_dir, recursive = TRUE))

  result <- runDeepTCR(immLynx_example, chains = "TRB",
                       output_dir = temp_dir,
                       latent_dim = 50,
                       return_seurat = FALSE)

  expect_equal(ncol(result), 50)
})
