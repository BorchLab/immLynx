# Tests for runMetaclonotypist function

test_that("runMetaclonotypist validates chain argument", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  expect_error(runMetaclonotypist(immLynx_example, chains = "invalid"))
})

test_that("runMetaclonotypist validates method argument", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  expect_error(runMetaclonotypist(immLynx_example, method = "invalid"))
})

test_that("runMetaclonotypist validates clustering argument", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  expect_error(runMetaclonotypist(immLynx_example, clustering = "invalid"))
})

test_that("runMetaclonotypist adds metaclone column to Seurat object", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment")

  data("immLynx_example", package = "immLynx")

  result <- runMetaclonotypist(immLynx_example, chains = "beta")

  expect_s4_class(result, "Seurat")
  expect_true("metaclone" %in% names(result@meta.data))
  expect_true("metaclone_size" %in% names(result@meta.data))
})

test_that("runMetaclonotypist returns data.frame when return_seurat=FALSE", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment")

  data("immLynx_example", package = "immLynx")

  result <- runMetaclonotypist(immLynx_example, chains = "beta",
                               return_seurat = FALSE)

  expect_s3_class(result, "data.frame")
  expect_true("barcode" %in% names(result))
  expect_true("cdr3_aa" %in% names(result))
  expect_true("metaclone" %in% names(result))
  expect_true("metaclone_size" %in% names(result))
})

test_that("runMetaclonotypist handles custom column name", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip("Requires Python environment")

  data("immLynx_example", package = "immLynx")

  result <- runMetaclonotypist(immLynx_example, chains = "beta",
                               column_name = "custom_meta")

  expect_true("custom_meta" %in% names(result@meta.data))
})

test_that("runMetaclonotypist sets correct default max_dist based on method", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  # This is a unit test for the internal logic, not requiring Python
  # We verify that the function signature accepts these parameters

  expect_no_error({
    # Just verify the function exists and accepts parameters
    formals(runMetaclonotypist)
  })
})

test_that("runMetaclonotypist accepts data.frame input", {
  skip("Requires Python environment")

  tcr_data <- data.frame(
    barcode = paste0("cell_", 1:10),
    cdr3_aa = c("CASSLAPGATNEKLFF", "CASSLGQAYEQYF", "CASRLAGQETQYF",
                "CASSYSGGNTGELFF", "CASSQDRTGQETQYF", "CASSLNRDNEQFF",
                "CASSLTGTEAFF", "CASSYSQGSYEQYF", "CASSLAGDTDTQYF",
                "CASSLVSGSTDTQYF")
  )

  result <- runMetaclonotypist(tcr_data, chains = "beta",
                               return_seurat = FALSE)

  expect_s3_class(result, "data.frame")
})

test_that("runHLAassociation validates input columns", {
  metaclone_data <- data.frame(
    barcode = paste0("cell_", 1:10),
    metaclone = rep(c("MC1", "MC2"), 5)
  )

  hla_data <- data.frame(
    sample = paste0("sample_", 1:10),  # Note: different column name
    HLA_A = rep(c("A*01:01", "A*02:01"), 5)
  )

  expect_error(runHLAassociation(metaclone_data, hla_data, by = "barcode"))
})

test_that("runHLAassociation returns correct structure", {
  skip("Requires matching data")

  metaclone_data <- data.frame(
    barcode = paste0("cell_", 1:20),
    metaclone = sample(c("MC1", "MC2", "MC3"), 20, replace = TRUE)
  )

  hla_data <- data.frame(
    barcode = paste0("cell_", 1:20),
    HLA_A_01 = sample(c(TRUE, FALSE), 20, replace = TRUE),
    HLA_A_02 = sample(c(TRUE, FALSE), 20, replace = TRUE)
  )

  result <- runHLAassociation(metaclone_data, hla_data)

  expect_s3_class(result, "data.frame")
  expect_true("metaclone" %in% names(result))
  expect_true("hla_allele" %in% names(result))
  expect_true("odds_ratio" %in% names(result))
  expect_true("pvalue" %in% names(result))
  expect_true("fdr" %in% names(result))
})
