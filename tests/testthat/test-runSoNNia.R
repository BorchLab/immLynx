# Tests for runSoNNia function

test_that("runSoNNia validates chain argument", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  expect_error(runSoNNia(immLynx_example, chains = "invalid",
                         background_file = "fake.csv"))
})

test_that("runSoNNia checks for background file existence", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  expect_error(runSoNNia(immLynx_example, chains = "TRB",
                         background_file = "nonexistent_file.csv"),
               "Background file not found")
})

test_that("runSoNNia creates output directory", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment")

  data("immLynx_example", package = "immLynx")

  # Create a mock background file
  temp_bg <- tempfile(fileext = ".csv")
  bg_data <- data.frame(
    aa_seq = c("CASSLAPGATNEKLFF", "CASSLGQAYEQYF", "CASRLAGQETQYF"),
    v_index = c(1, 2, 3),
    j_index = c(1, 1, 2)
  )
  write.csv(bg_data, temp_bg, row.names = FALSE)

  temp_dir <- tempfile("sonia_test")
  on.exit({
    unlink(temp_bg)
    unlink(temp_dir, recursive = TRUE)
  })

  result <- runSoNNia(immLynx_example, chains = "TRB",
                      background_file = temp_bg,
                      save_folder = temp_dir)

  expect_true(dir.exists(temp_dir))
})

test_that("runSoNNia adds results to Seurat object", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment")

  data("immLynx_example", package = "immLynx")

  # Create a mock background file
  temp_bg <- tempfile(fileext = ".csv")
  bg_data <- data.frame(
    aa_seq = c("CASSLAPGATNEKLFF", "CASSLGQAYEQYF", "CASRLAGQETQYF"),
    v_index = c(1, 2, 3),
    j_index = c(1, 1, 2)
  )
  write.csv(bg_data, temp_bg, row.names = FALSE)
  on.exit(unlink(temp_bg))

  result <- runSoNNia(immLynx_example, chains = "TRB",
                      background_file = temp_bg,
                      return_seurat = TRUE)

  expect_s4_class(result, "Seurat")
})

test_that("runSoNNia returns results when return_seurat=FALSE", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment")

  data("immLynx_example", package = "immLynx")

  # Create a mock background file
  temp_bg <- tempfile(fileext = ".csv")
  bg_data <- data.frame(
    aa_seq = c("CASSLAPGATNEKLFF", "CASSLGQAYEQYF", "CASRLAGQETQYF"),
    v_index = c(1, 2, 3),
    j_index = c(1, 1, 2)
  )
  write.csv(bg_data, temp_bg, row.names = FALSE)
  on.exit(unlink(temp_bg))

  result <- runSoNNia(immLynx_example, chains = "TRB",
                      background_file = temp_bg,
                      return_seurat = FALSE)

  expect_type(result, "list")
})
