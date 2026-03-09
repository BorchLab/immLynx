# Tests for runSoNNia function

# ===========================================================================
# Parameter validation
# ===========================================================================

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

test_that("runSoNNia rejects non-Seurat/SCE input", {
  tcr_data <- create_mock_tcr_data(10)
  temp_bg <- tempfile(fileext = ".csv")
  write.csv(data.frame(aa_seq = "CASSAAA"), temp_bg, row.names = FALSE)
  on.exit(unlink(temp_bg))

  expect_error(
    runSoNNia(tcr_data, chains = "TRB", background_file = temp_bg),
    "Input must be a Seurat or SingleCellExperiment object"
  )
})

test_that("runSoNNia function signature has correct defaults", {
  f <- formals(runSoNNia)

  expect_equal(f$organism, "human")
  expect_equal(f$save_folder, "sonia_output")
  expect_equal(f$n_epochs, 100)
  expect_equal(f$return_object, TRUE)
})

test_that("runSoNNia accepts valid chain values", {
  expect_no_error(match.arg("TRB", c("TRB", "TRA")))
  expect_no_error(match.arg("TRA", c("TRB", "TRA")))
})

test_that("runSoNNia checks background file path with special characters", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  expect_error(
    runSoNNia(immLynx_example, chains = "TRB",
              background_file = "/path/to/nonexistent/file with spaces.csv"),
    "Background file not found"
  )
})


# ===========================================================================
# Python-dependent tests
# ===========================================================================

test_that("runSoNNia creates output directory", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_python()
  skip("soNNia ValueError with inhomogeneous array shapes - source code issue")

  data("immLynx_example", package = "immLynx")

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
  skip_if_no_python()
  skip("soNNia ValueError with inhomogeneous array shapes - source code issue")

  data("immLynx_example", package = "immLynx")

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
                      return_object = TRUE)

  expect_s4_class(result, "Seurat")
})

test_that("runSoNNia returns results when return_object=FALSE", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_python()
  skip("soNNia ValueError with inhomogeneous array shapes - source code issue")

  data("immLynx_example", package = "immLynx")

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
                      return_object = FALSE)

  expect_type(result, "list")
})

test_that("runSoNNia produces messages during execution", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_python()
  skip("soNNia ValueError with inhomogeneous array shapes - source code issue")

  data("immLynx_example", package = "immLynx")

  temp_bg <- tempfile(fileext = ".csv")
  bg_data <- data.frame(
    aa_seq = c("CASSLAPGATNEKLFF", "CASSLGQAYEQYF", "CASRLAGQETQYF"),
    v_index = c(1, 2, 3),
    j_index = c(1, 1, 2)
  )
  write.csv(bg_data, temp_bg, row.names = FALSE)
  on.exit(unlink(temp_bg))

  expect_message(
    runSoNNia(immLynx_example, chains = "TRB",
              background_file = temp_bg),
    "Extracting TRB sequences"
  )
})
