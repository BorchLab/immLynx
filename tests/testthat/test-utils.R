# Tests for utility functions

test_that("extractTCRdata returns correct structure for single chain", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  result <- extractTCRdata(immLynx_example, chains = "TRB", format = "long")

  expect_s3_class(result, "data.frame")
  expect_true("barcode" %in% names(result))
  expect_true("cdr3_aa" %in% names(result))
  expect_true("chain" %in% names(result))
  expect_true(all(result$chain == "TRB"))
})

test_that("extractTCRdata returns wide format correctly", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  result <- extractTCRdata(immLynx_example, chains = "both", format = "wide")

  expect_s3_class(result, "data.frame")
  expect_true("barcode" %in% names(result))
  expect_true(any(grepl("_TRA$|_TRB$", names(result))))
})

test_that("extractTCRdata removes NA sequences when requested", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  result_no_na <- extractTCRdata(immLynx_example, chains = "TRB", remove_na = TRUE)
  result_with_na <- extractTCRdata(immLynx_example, chains = "TRB", remove_na = FALSE)

  expect_true(all(!is.na(result_no_na$cdr3_aa)))
  expect_gte(nrow(result_with_na), nrow(result_no_na))
})

test_that("validateTCRdata identifies missing columns", {
  bad_data <- data.frame(x = 1:5)
  result <- validateTCRdata(bad_data, strict = FALSE)

  expect_false(result$valid)
  expect_true(any(grepl("Missing required columns", result$errors)))
})

test_that("validateTCRdata passes for valid data", {
  good_data <- data.frame(
    barcode = paste0("cell_", 1:5),
    cdr3_aa = c("CASSLAPGATNEKLFF", "CASSLGQAYEQYF", "CASRLAGQETQYF",
                "CASSYSGGNTGELFF", "CASSQDRTGQETQYF"),
    v = paste0("TRBV", 1:5),
    j = paste0("TRBJ", 1:5),
    chain = "TRB"
  )

  result <- validateTCRdata(good_data, strict = FALSE)

  expect_true(result$valid)
  expect_length(result$errors, 0)
})

test_that("validateTCRdata detects invalid amino acids", {
  bad_seq_data <- data.frame(
    barcode = c("cell_1", "cell_2"),
    cdr3_aa = c("CASSXXX123", "CASSLGQAYEQYF")
  )

  result <- validateTCRdata(bad_seq_data, check_sequences = TRUE, strict = FALSE)

  expect_true(any(grepl("non-standard amino acids", result$warnings)))
})

test_that("validateTCRdata strict mode throws error", {
  bad_data <- data.frame(x = 1:5)

  expect_error(validateTCRdata(bad_data, strict = TRUE),
               "TCR data validation")
})

test_that("convertToTcrdist produces correct format for beta chain", {
  tcr_data <- data.frame(
    barcode = paste0("cell_", 1:3),
    cdr3_aa = c("CASSLAPGATNEKLFF", "CASSLGQAYEQYF", "CASRLAGQETQYF"),
    v = c("TRBV7-2", "TRBV6-1", "TRBV12-3"),
    j = c("TRBJ1-1", "TRBJ2-1", "TRBJ2-5"),
    chain = "TRB"
  )

  result <- convertToTcrdist(tcr_data, chains = "beta")

  expect_true("count" %in% names(result))
  expect_true("v_b_gene" %in% names(result))
  expect_true("j_b_gene" %in% names(result))
  expect_true("cdr3_b_aa" %in% names(result))
  expect_equal(result$count, rep(1L, 3))
})

test_that("convertToTcrdist handles wide format", {
  wide_data <- data.frame(
    barcode = c("cell_1", "cell_2"),
    cdr3_aa_TRA = c("CAVSEAPNQAGTALIF", "CAVRDSSYKLIF"),
    v_TRA = c("TRAV12-1", "TRAV8-1"),
    j_TRA = c("TRAJ15", "TRAJ12"),
    cdr3_aa_TRB = c("CASSLAPGATNEKLFF", "CASSLGQAYEQYF"),
    v_TRB = c("TRBV7-2", "TRBV6-1"),
    j_TRB = c("TRBJ1-1", "TRBJ2-1")
  )

  result <- convertToTcrdist(wide_data, chains = "both")

  expect_true("v_a_gene" %in% names(result))
  expect_true("cdr3_a_aa" %in% names(result))
  expect_true("v_b_gene" %in% names(result))
  expect_true("cdr3_b_aa" %in% names(result))
})

test_that("summarizeTCRrepertoire returns correct structure", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  result <- summarizeTCRrepertoire(immLynx_example, chains = "TRB")

  expect_s3_class(result, "TCR_summary")
  expect_true("total_cells" %in% names(result))
  expect_true("unique_clonotypes" %in% names(result))
  expect_true("diversity" %in% names(result))
  expect_true("top_clones" %in% names(result))
})

test_that("summarizeTCRrepertoire calculates diversity metrics", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  result <- summarizeTCRrepertoire(immLynx_example, chains = "TRB",
                                   calculate_diversity = TRUE)

  expect_true(!is.null(result$diversity))
  expect_true("shannon" %in% names(result$diversity))
  expect_true("simpson" %in% names(result$diversity))
  expect_true("clonality" %in% names(result$diversity))
  expect_gte(result$diversity$clonality, 0)
  expect_lte(result$diversity$clonality, 1)
})

test_that("print.TCR_summary works without error", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  result <- summarizeTCRrepertoire(immLynx_example, chains = "TRB")

  expect_output(print(result), "TCR Repertoire Summary")
})

test_that("summarizeTCRrepertoire works with data.frame input", {
  tcr_data <- data.frame(
    barcode = paste0("cell_", 1:100),
    cdr3_aa = sample(paste0("CASS", LETTERS[1:10]), 100, replace = TRUE),
    v = sample(paste0("TRBV", 1:5), 100, replace = TRUE),
    j = sample(paste0("TRBJ", 1:3), 100, replace = TRUE),
    chain = "TRB"
  )

  result <- summarizeTCRrepertoire(tcr_data, chains = "TRB")

  expect_s3_class(result, "TCR_summary")
  expect_equal(result$total_cells, 100)
})
