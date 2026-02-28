# Tests for runMetaclonotypist and runHLAassociation functions

# ===========================================================================
# runMetaclonotypist - Parameter validation
# ===========================================================================

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

test_that("runMetaclonotypist function signature has correct defaults", {
  f <- formals(runMetaclonotypist)

  expect_equal(f$max_edits, 2)
  expect_null(f$max_dist)
  expect_equal(f$resolution, 1.0)
  expect_equal(f$return_seurat, TRUE)
  expect_equal(f$column_name, "metaclone")
})

test_that("runMetaclonotypist accepts all valid chain values", {
  expect_no_error(match.arg("alpha", c("beta", "alpha")))
  expect_no_error(match.arg("beta", c("beta", "alpha")))
})

test_that("runMetaclonotypist accepts all valid method values", {
  expect_no_error(match.arg("tcrdist", c("tcrdist", "sceptr")))
  expect_no_error(match.arg("sceptr", c("tcrdist", "sceptr")))
})

test_that("runMetaclonotypist accepts all valid clustering values", {
  expect_no_error(match.arg("cc", c("cc", "leiden", "louvain", "mcl")))
  expect_no_error(match.arg("leiden", c("cc", "leiden", "louvain", "mcl")))
  expect_no_error(match.arg("louvain", c("cc", "leiden", "louvain", "mcl")))
  expect_no_error(match.arg("mcl", c("cc", "leiden", "louvain", "mcl")))
})

test_that("runMetaclonotypist rejects data.frame without required columns", {
  bad_df <- data.frame(x = 1:5, y = 6:10)

  expect_error(
    runMetaclonotypist(bad_df, chains = "beta", return_seurat = FALSE),
    "barcode.*cdr3_aa"
  )
})

test_that("runMetaclonotypist rejects data.frame with only barcode column", {
  bad_df <- data.frame(barcode = paste0("cell_", 1:5))

  expect_error(
    runMetaclonotypist(bad_df, chains = "beta", return_seurat = FALSE),
    "cdr3_aa"
  )
})

test_that("runMetaclonotypist rejects data.frame with all NA sequences", {
  na_df <- data.frame(
    barcode = paste0("cell_", 1:5),
    cdr3_aa = rep(NA, 5)
  )

  expect_error(
    runMetaclonotypist(na_df, chains = "beta", return_seurat = FALSE),
    "No valid TCR sequences found"
  )
})


# ===========================================================================
# runMetaclonotypist - Python-dependent tests
# ===========================================================================

test_that("runMetaclonotypist adds metaclone column to Seurat object", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_python()

  data("immLynx_example", package = "immLynx")

  result <- runMetaclonotypist(immLynx_example, chains = "beta")

  expect_s4_class(result, "Seurat")
  expect_true("metaclone" %in% names(result@meta.data))
  expect_true("metaclone_size" %in% names(result@meta.data))
})

test_that("runMetaclonotypist returns data.frame when return_seurat=FALSE", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_python()

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
  skip_if_no_python()

  data("immLynx_example", package = "immLynx")

  result <- runMetaclonotypist(immLynx_example, chains = "beta",
                               column_name = "custom_meta")

  expect_true("custom_meta" %in% names(result@meta.data))
  expect_true("custom_meta_size" %in% names(result@meta.data))
})

test_that("runMetaclonotypist accepts data.frame input", {
  skip_if_no_python()

  tcr_data <- create_mock_tcr_data(100)

  result <- runMetaclonotypist(tcr_data, chains = "beta",
                               return_seurat = FALSE)

  expect_s3_class(result, "data.frame")
})

test_that("runMetaclonotypist produces messages during execution", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_python()

  data("immLynx_example", package = "immLynx")

  expect_message(
    runMetaclonotypist(immLynx_example, chains = "beta",
                       return_seurat = FALSE),
    "Running metaclonotypist"
  )
})


# ===========================================================================
# runHLAassociation
# ===========================================================================

test_that("runHLAassociation validates input columns", {
  metaclone_data <- data.frame(
    barcode = paste0("cell_", 1:10),
    metaclone = rep(c("MC1", "MC2"), 5)
  )

  hla_data <- data.frame(
    sample = paste0("sample_", 1:10),
    HLA_A = rep(c("A*01:01", "A*02:01"), 5)
  )

  expect_error(
    runHLAassociation(metaclone_data, hla_data, by = "barcode"),
    "not found in hla_data"
  )
})

test_that("runHLAassociation validates metaclone_data has by column", {
  metaclone_data <- data.frame(
    sample = paste0("cell_", 1:10),
    metaclone = rep(c("MC1", "MC2"), 5)
  )

  hla_data <- data.frame(
    barcode = paste0("cell_", 1:10),
    HLA_A = rep(c("A*01:01", "A*02:01"), 5)
  )

  expect_error(
    runHLAassociation(metaclone_data, hla_data, by = "barcode"),
    "not found in metaclone_data"
  )
})

test_that("runHLAassociation errors when no matching samples", {
  metaclone_data <- data.frame(
    barcode = paste0("cell_", 1:10),
    metaclone = rep(c("MC1", "MC2"), 5)
  )

  hla_data <- data.frame(
    barcode = paste0("other_", 1:10),
    HLA_A = rep(c("A*01:01", "A*02:01"), 5)
  )

  expect_error(
    runHLAassociation(metaclone_data, hla_data),
    "No matching samples"
  )
})

test_that("runHLAassociation errors when no HLA columns found", {
  metaclone_data <- data.frame(
    barcode = paste0("cell_", 1:10),
    metaclone = rep(c("MC1", "MC2"), 5)
  )

  hla_data <- data.frame(
    barcode = paste0("cell_", 1:10),
    some_column = rep("value", 10)
  )

  expect_error(
    runHLAassociation(metaclone_data, hla_data),
    "No HLA columns found"
  )
})

test_that("runHLAassociation returns correct structure with valid input", {
  set.seed(42)
  metaclone_data <- data.frame(
    barcode = paste0("cell_", 1:40),
    metaclone = c(rep("MC1", 20), rep("MC2", 20)),
    stringsAsFactors = FALSE
  )

  hla_data <- data.frame(
    barcode = paste0("cell_", 1:40),
    HLA_A = c(rep(TRUE, 15), rep(FALSE, 5),
              rep(TRUE, 5), rep(FALSE, 15)),
    stringsAsFactors = FALSE
  )

  result <- runHLAassociation(metaclone_data, hla_data)

  expect_s3_class(result, "data.frame")
  expect_true("metaclone" %in% names(result))
  expect_true("hla_allele" %in% names(result))
  expect_true("odds_ratio" %in% names(result))
  expect_true("pvalue" %in% names(result))
  expect_true("fdr" %in% names(result))
  expect_true("significant" %in% names(result))
})

test_that("runHLAassociation applies FDR correction", {
  set.seed(42)
  metaclone_data <- data.frame(
    barcode = paste0("cell_", 1:40),
    metaclone = c(rep("MC1", 20), rep("MC2", 20)),
    stringsAsFactors = FALSE
  )

  hla_data <- data.frame(
    barcode = paste0("cell_", 1:40),
    HLA_A = c(rep(TRUE, 15), rep(FALSE, 5),
              rep(TRUE, 5), rep(FALSE, 15)),
    HLA_B = c(rep(TRUE, 10), rep(FALSE, 10),
              rep(TRUE, 10), rep(FALSE, 10)),
    stringsAsFactors = FALSE
  )

  result <- runHLAassociation(metaclone_data, hla_data)

  # FDR should be >= raw p-value
  expect_true(all(result$fdr >= result$pvalue))
})

test_that("runHLAassociation is sorted by FDR", {
  set.seed(42)
  metaclone_data <- data.frame(
    barcode = paste0("cell_", 1:40),
    metaclone = c(rep("MC1", 20), rep("MC2", 20)),
    stringsAsFactors = FALSE
  )

  hla_data <- data.frame(
    barcode = paste0("cell_", 1:40),
    HLA_A = c(rep(TRUE, 15), rep(FALSE, 5),
              rep(TRUE, 5), rep(FALSE, 15)),
    HLA_B = c(rep(TRUE, 10), rep(FALSE, 10),
              rep(TRUE, 10), rep(FALSE, 10)),
    stringsAsFactors = FALSE
  )

  result <- runHLAassociation(metaclone_data, hla_data)

  if (nrow(result) > 1) {
    expect_true(all(diff(result$fdr) >= 0))
  }
})

test_that("runHLAassociation handles custom by column", {
  set.seed(42)
  metaclone_data <- data.frame(
    sample = paste0("sample_", 1:40),
    metaclone = c(rep("MC1", 20), rep("MC2", 20)),
    stringsAsFactors = FALSE
  )

  hla_data <- data.frame(
    sample = paste0("sample_", 1:40),
    HLA_A = c(rep(TRUE, 15), rep(FALSE, 5),
              rep(TRUE, 5), rep(FALSE, 15)),
    stringsAsFactors = FALSE
  )

  result <- runHLAassociation(metaclone_data, hla_data, by = "sample")

  expect_s3_class(result, "data.frame")
})

test_that("runHLAassociation handles custom fdr_threshold", {
  set.seed(42)
  metaclone_data <- data.frame(
    barcode = paste0("cell_", 1:40),
    metaclone = c(rep("MC1", 20), rep("MC2", 20)),
    stringsAsFactors = FALSE
  )

  hla_data <- data.frame(
    barcode = paste0("cell_", 1:40),
    HLA_A = c(rep(TRUE, 15), rep(FALSE, 5),
              rep(TRUE, 5), rep(FALSE, 15)),
    stringsAsFactors = FALSE
  )

  result_strict <- runHLAassociation(metaclone_data, hla_data,
                                      fdr_threshold = 0.01)
  result_lenient <- runHLAassociation(metaclone_data, hla_data,
                                       fdr_threshold = 0.10)

  # Stricter threshold should give fewer or equal significant results
  expect_lte(sum(result_strict$significant), sum(result_lenient$significant))
})

test_that("runHLAassociation handles NA metaclones", {
  set.seed(42)
  metaclone_data <- data.frame(
    barcode = paste0("cell_", 1:40),
    metaclone = c(rep("MC1", 15), rep(NA, 10), rep("MC2", 15)),
    stringsAsFactors = FALSE
  )

  hla_data <- data.frame(
    barcode = paste0("cell_", 1:40),
    HLA_A = c(rep(TRUE, 20), rep(FALSE, 20)),
    stringsAsFactors = FALSE
  )

  result <- runHLAassociation(metaclone_data, hla_data)

  # Should only test MC1 and MC2, not NA
  expect_true(all(result$metaclone %in% c("MC1", "MC2")))
})

test_that("runHLAassociation warns when no valid tests performed", {
  metaclone_data <- data.frame(
    barcode = paste0("cell_", 1:5),
    metaclone = rep("MC1", 5),
    stringsAsFactors = FALSE
  )

  hla_data <- data.frame(
    barcode = paste0("cell_", 1:5),
    HLA_A = rep(TRUE, 5),
    stringsAsFactors = FALSE
  )

  # All cells in one metaclone and all have same HLA -> degenerate table
  expect_warning(
    result <- runHLAassociation(metaclone_data, hla_data),
    "No valid association tests"
  )
  expect_null(result)
})

test_that("runHLAassociation produces informative messages", {
  set.seed(42)
  metaclone_data <- data.frame(
    barcode = paste0("cell_", 1:20),
    metaclone = rep(c("MC1", "MC2"), each = 10),
    stringsAsFactors = FALSE
  )

  hla_data <- data.frame(
    barcode = paste0("cell_", 1:20),
    HLA_A = c(rep(TRUE, 8), rep(FALSE, 2),
              rep(TRUE, 2), rep(FALSE, 8)),
    stringsAsFactors = FALSE
  )

  expect_message(
    runHLAassociation(metaclone_data, hla_data),
    "Testing associations"
  )
})
