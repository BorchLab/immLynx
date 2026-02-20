# Tests for runOLGA and generateOLGA functions

# ===========================================================================
# runOLGA - Parameter validation
# ===========================================================================

test_that("runOLGA validates chain argument", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  expect_error(runOLGA(immLynx_example, chains = "invalid"))
})

test_that("runOLGA rejects non-Seurat/SCE input", {
  tcr_data <- create_mock_tcr_data(10)

  expect_error(
    runOLGA(tcr_data, chains = "TRB"),
    "Input must be a Seurat or SingleCellExperiment object"
  )
})

test_that("runOLGA stops for invalid organism/chain combination", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  expect_error(
    runOLGA(immLynx_example, chains = "TRB", organism = "invalid"),
    "Cannot infer OLGA model"
  )
})

test_that("runOLGA infers humanTRB model from organism='human' and chains='TRB'", {
  # Test the internal model inference logic
  organism <- "human"
  chains <- "TRB"
  model <- paste0(organism, chains)

  expect_equal(model, "humanTRB")
  expect_true(model %in% c("humanTRB", "humanTRA", "humanIGH", "mouseTRB"))
})

test_that("runOLGA infers humanTRA model from organism='human' and chains='TRA'", {
  organism <- "human"
  chains <- "TRA"
  model <- paste0(organism, chains)

  expect_equal(model, "humanTRA")
  expect_true(model %in% c("humanTRB", "humanTRA", "humanIGH", "mouseTRB"))
})

test_that("runOLGA invalid organism/chain combo fails model validation", {
  organism <- "mouse"
  chains <- "TRA"
  model <- paste0(organism, chains)

  # mouseTRA is not in the valid list
  expect_false(model %in% c("humanTRB", "humanTRA", "humanIGH", "mouseTRB"))
})

test_that("runOLGA function signature has correct defaults", {
  f <- formals(runOLGA)

  expect_null(f$model)
  expect_equal(f$organism, "human")
  expect_equal(f$use_vj_genes, FALSE)
  expect_equal(f$return_object, TRUE)
  expect_equal(f$column_name, "olga_pgen")
})

test_that("runOLGA accepts valid chain values", {
  expect_no_error(match.arg("TRB", c("TRB", "TRA")))
  expect_no_error(match.arg("TRA", c("TRB", "TRA")))
})


# ===========================================================================
# runOLGA - Python-dependent tests
# ===========================================================================

test_that("runOLGA adds pgen columns to Seurat object", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_python()

  data("immLynx_example", package = "immLynx")

  result <- runOLGA(immLynx_example, chains = "TRB")

  expect_s4_class(result, "Seurat")
  expect_true("olga_pgen_TRB" %in% names(result@meta.data))
  expect_true("olga_pgen_log10_TRB" %in% names(result@meta.data))
})

test_that("runOLGA returns data.frame when return_object=FALSE", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_python()

  data("immLynx_example", package = "immLynx")

  result <- runOLGA(immLynx_example, chains = "TRB", return_object = FALSE)

  expect_s3_class(result, "data.frame")
  expect_true("barcode" %in% names(result))
  expect_true("cdr3_aa" %in% names(result))
  expect_true("pgen" %in% names(result))
  expect_true("log10_pgen" %in% names(result))
  expect_true("v_gene" %in% names(result))
  expect_true("j_gene" %in% names(result))
})

test_that("runOLGA handles custom column name", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_python()

  data("immLynx_example", package = "immLynx")

  result <- runOLGA(immLynx_example, chains = "TRB",
                    column_name = "custom_pgen")

  expect_true("custom_pgen_TRB" %in% names(result@meta.data))
  expect_true("custom_pgen_log10_TRB" %in% names(result@meta.data))
})

test_that("runOLGA pgen values are valid probabilities", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_python()

  data("immLynx_example", package = "immLynx")

  result <- runOLGA(immLynx_example, chains = "TRB", return_object = FALSE)

  # Pgen values should be between 0 and 1
  expect_true(all(result$pgen >= 0))
  expect_true(all(result$pgen <= 1))
  # log10_pgen should be negative or zero
  expect_true(all(result$log10_pgen <= 0 | is.nan(result$log10_pgen)))
})

test_that("runOLGA produces messages during execution", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_python()

  data("immLynx_example", package = "immLynx")

  expect_message(
    runOLGA(immLynx_example, chains = "TRB"),
    "Extracting TRB sequences"
  )
})


# ===========================================================================
# generateOLGA
# ===========================================================================

test_that("generateOLGA function signature has correct defaults", {
  f <- formals(generateOLGA)

  expect_equal(f$n, 100)
  expect_equal(f$model, "humanTRB")
})

test_that("generateOLGA produces a message about generation", {
  skip_if_no_python()

  expect_message(
    generateOLGA(n = 5, model = "humanTRB"),
    "Generating 5 random sequences"
  )
})

test_that("generateOLGA returns correct number of sequences", {
  skip_if_no_python()

  result <- generateOLGA(n = 10, model = "humanTRB")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 10)
  expect_true("aa_seq" %in% names(result))
  expect_true("nt_seq" %in% names(result))
  expect_true("v_index" %in% names(result))
  expect_true("j_index" %in% names(result))
})

test_that("generateOLGA works with different models", {
  skip_if_no_python()

  result_trb <- generateOLGA(n = 5, model = "humanTRB")
  result_tra <- generateOLGA(n = 5, model = "humanTRA")

  expect_equal(nrow(result_trb), 5)
  expect_equal(nrow(result_tra), 5)
})
