# Tests for utility functions

# ===========================================================================
# extractTCRdata
# ===========================================================================

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

test_that("extractTCRdata validates chains argument", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  expect_error(extractTCRdata(immLynx_example, chains = "invalid"))
})

test_that("extractTCRdata validates format argument", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  expect_error(extractTCRdata(immLynx_example, format = "invalid"))
})

test_that("extractTCRdata long format for both chains stacks data", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  result <- extractTCRdata(immLynx_example, chains = "both", format = "long")

  expect_s3_class(result, "data.frame")
  expect_true("chain" %in% names(result))
  chain_vals <- unique(result$chain)
  expect_true("TRA" %in% chain_vals || "TRB" %in% chain_vals)
})

test_that("extractTCRdata single chain adds chain column", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  result <- extractTCRdata(immLynx_example, chains = "TRA")

  expect_true("chain" %in% names(result))
  expect_true(all(result$chain == "TRA"))
})

test_that("extractTCRdata accepts default chains parameter", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  # Default should be TRB
  result <- extractTCRdata(immLynx_example)

  expect_s3_class(result, "data.frame")
  expect_true(all(result$chain == "TRB"))
})


# ===========================================================================
# validateTCRdata
# ===========================================================================

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

test_that("validateTCRdata strict mode returns invisible TRUE for valid data", {
  good_data <- data.frame(
    barcode = paste0("cell_", 1:5),
    cdr3_aa = c("CASSLAPGATNEKLFF", "CASSLGQAYEQYF", "CASRLAGQETQYF",
                "CASSYSGGNTGELFF", "CASSQDRTGQETQYF"),
    v = paste0("TRBV", 1:5),
    j = paste0("TRBJ", 1:5),
    chain = "TRB"
  )

  result <- validateTCRdata(good_data, strict = TRUE)
  expect_true(result)
})

test_that("validateTCRdata strict mode warns for missing recommended columns", {
  minimal_data <- data.frame(
    barcode = paste0("cell_", 1:3),
    cdr3_aa = c("CASSLAPGATNEKLFF", "CASSLGQAYEQYF", "CASRLAGQETQYF")
  )

  expect_warning(
    validateTCRdata(minimal_data, strict = TRUE),
    "Missing recommended columns"
  )
})

test_that("validateTCRdata warns about missing recommended columns in non-strict mode", {
  minimal_data <- data.frame(
    barcode = paste0("cell_", 1:3),
    cdr3_aa = c("CASSLAPGATNEKLFF", "CASSLGQAYEQYF", "CASRLAGQETQYF")
  )

  result <- validateTCRdata(minimal_data, strict = FALSE)

  expect_true(result$valid)
  expect_true(any(grepl("Missing recommended columns", result$warnings)))
})

test_that("validateTCRdata handles check_sequences=FALSE", {
  data_with_invalid <- data.frame(
    barcode = c("cell_1", "cell_2"),
    cdr3_aa = c("CASSXXX123", "INVALID!!!")
  )

  result <- validateTCRdata(data_with_invalid, check_sequences = FALSE, strict = FALSE)

  expect_false(any(grepl("non-standard amino acids", result$warnings)))
})

test_that("validateTCRdata check_genes validates V gene format", {
  data_bad_v <- data.frame(
    barcode = paste0("cell_", 1:3),
    cdr3_aa = c("CASSLAPGATNEKLFF", "CASSLGQAYEQYF", "CASRLAGQETQYF"),
    v = c("TRBV7-2", "bad_gene", "not_valid"),
    j = paste0("TRBJ", 1:3),
    chain = "TRB"
  )

  result <- validateTCRdata(data_bad_v, check_genes = TRUE, strict = FALSE)

  expect_true(any(grepl("V genes do not match IMGT", result$warnings)))
})

test_that("validateTCRdata check_genes validates J gene format", {
  data_bad_j <- data.frame(
    barcode = paste0("cell_", 1:3),
    cdr3_aa = c("CASSLAPGATNEKLFF", "CASSLGQAYEQYF", "CASRLAGQETQYF"),
    v = paste0("TRBV", 1:3),
    j = c("TRBJ1-1", "bad_gene", "not_valid"),
    chain = "TRB"
  )

  result <- validateTCRdata(data_bad_j, check_genes = TRUE, strict = FALSE)

  expect_true(any(grepl("J genes do not match IMGT", result$warnings)))
})

test_that("validateTCRdata check_genes handles NA genes gracefully", {
  data_na_genes <- data.frame(
    barcode = paste0("cell_", 1:3),
    cdr3_aa = c("CASSLAPGATNEKLFF", "CASSLGQAYEQYF", "CASRLAGQETQYF"),
    v = c("TRBV7-2", NA, "TRBV5-1"),
    j = c(NA, "TRBJ2-1", NA),
    chain = "TRB"
  )

  result <- validateTCRdata(data_na_genes, check_genes = TRUE, strict = FALSE)

  expect_true(result$valid)
})

test_that("validateTCRdata handles data with all NA cdr3_aa", {
  na_data <- data.frame(
    barcode = paste0("cell_", 1:3),
    cdr3_aa = c(NA, NA, NA)
  )

  result <- validateTCRdata(na_data, check_sequences = TRUE, strict = FALSE)

  expect_true(result$valid)
})

test_that("validateTCRdata summary stats are correct", {
  tcr_data <- data.frame(
    barcode = paste0("cell_", 1:5),
    cdr3_aa = c("CASSA", "CASSB", "CASSA", NA, "CASSC"),
    chain = "TRB"
  )

  result <- validateTCRdata(tcr_data, strict = FALSE)

  expect_equal(result$summary$n_rows, 5)
  expect_equal(result$summary$n_unique_barcodes, 5)
  expect_equal(result$summary$n_unique_cdr3, 3)
  expect_equal(result$summary$n_na_cdr3, 1)
  expect_equal(result$summary$chains, "TRB")
})

test_that("validateTCRdata reports missing barcode column as error", {
  no_barcode <- data.frame(
    cdr3_aa = c("CASSLAPGATNEKLFF", "CASSLGQAYEQYF")
  )

  result <- validateTCRdata(no_barcode, strict = FALSE)

  expect_false(result$valid)
  expect_true(any(grepl("barcode", result$errors)))
})

test_that("validateTCRdata reports missing cdr3_aa column as error", {
  no_cdr3 <- data.frame(
    barcode = paste0("cell_", 1:3)
  )

  result <- validateTCRdata(no_cdr3, strict = FALSE)

  expect_false(result$valid)
  expect_true(any(grepl("cdr3_aa", result$errors)))
})

test_that("validateTCRdata accepts valid amino acid sequences with stop codon", {
  data_star <- data.frame(
    barcode = paste0("cell_", 1:2),
    cdr3_aa = c("CASS*", "CASSLAPGATNEKLFF")
  )

  result <- validateTCRdata(data_star, check_sequences = TRUE, strict = FALSE)

  expect_false(any(grepl("non-standard amino acids", result$warnings)))
})


# ===========================================================================
# convertToTcrdist
# ===========================================================================

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

test_that("convertToTcrdist validates chains argument", {
  tcr_data <- create_mock_tcr_data(10)

  expect_error(convertToTcrdist(tcr_data, chains = "invalid"))
})

test_that("convertToTcrdist handles alpha chain from long format", {
  tcr_data <- data.frame(
    barcode = paste0("cell_", 1:3),
    cdr3_aa = c("CAVSEAPNQAGTALIF", "CAVRDSSYKLIF", "CAGQTGGFKTIF"),
    v = c("TRAV12-1", "TRAV8-1", "TRAV6"),
    j = c("TRAJ15", "TRAJ12", "TRAJ10"),
    chain = "TRA"
  )

  result <- convertToTcrdist(tcr_data, chains = "alpha")

  expect_true("v_a_gene" %in% names(result))
  expect_true("j_a_gene" %in% names(result))
  expect_true("cdr3_a_aa" %in% names(result))
  expect_false("v_b_gene" %in% names(result))
  expect_equal(nrow(result), 3)
})

test_that("convertToTcrdist handles both chains from long format", {
  tcr_data <- data.frame(
    barcode = c("cell_1", "cell_2", "cell_1", "cell_2"),
    cdr3_aa = c("CAVSEAPNQAGTALIF", "CAVRDSSYKLIF",
                "CASSLAPGATNEKLFF", "CASSLGQAYEQYF"),
    v = c("TRAV12-1", "TRAV8-1", "TRBV7-2", "TRBV6-1"),
    j = c("TRAJ15", "TRAJ12", "TRBJ1-1", "TRBJ2-1"),
    chain = c("TRA", "TRA", "TRB", "TRB")
  )

  result <- convertToTcrdist(tcr_data, chains = "both")

  expect_true("v_a_gene" %in% names(result))
  expect_true("v_b_gene" %in% names(result))
  expect_true("count" %in% names(result))
  expect_equal(nrow(result), 2)
})

test_that("convertToTcrdist errors on both chains without chain column", {
  tcr_data <- data.frame(
    barcode = paste0("cell_", 1:3),
    cdr3_aa = c("CASSLAPGATNEKLFF", "CASSLGQAYEQYF", "CASRLAGQETQYF"),
    v = c("TRBV7-2", "TRBV6-1", "TRBV12-3"),
    j = c("TRBJ1-1", "TRBJ2-1", "TRBJ2-5")
  )

  expect_error(convertToTcrdist(tcr_data, chains = "both"),
               "Cannot determine chain type")
})

test_that("convertToTcrdist without chain column assumes single chain", {
  tcr_data <- data.frame(
    barcode = paste0("cell_", 1:3),
    cdr3_aa = c("CASSLAPGATNEKLFF", "CASSLGQAYEQYF", "CASRLAGQETQYF"),
    v = c("TRBV7-2", "TRBV6-1", "TRBV12-3"),
    j = c("TRBJ1-1", "TRBJ2-1", "TRBJ2-5")
  )

  result <- convertToTcrdist(tcr_data, chains = "beta")

  expect_true("v_b_gene" %in% names(result))
  expect_equal(nrow(result), 3)
})

test_that("convertToTcrdist include_count=FALSE omits count column", {
  tcr_data <- create_mock_tcr_data(5)

  result <- convertToTcrdist(tcr_data, chains = "beta", include_count = FALSE)

  expect_false("count" %in% names(result))
})

test_that("convertToTcrdist count column is first position", {
  tcr_data <- create_mock_tcr_data(5)

  result <- convertToTcrdist(tcr_data, chains = "beta")

  expect_equal(names(result)[1], "count")
})

test_that("convertToTcrdist wide format alpha only", {
  wide_data <- data.frame(
    barcode = c("cell_1", "cell_2"),
    cdr3_aa_TRA = c("CAVSEAPNQAGTALIF", "CAVRDSSYKLIF"),
    v_TRA = c("TRAV12-1", "TRAV8-1"),
    j_TRA = c("TRAJ15", "TRAJ12"),
    cdr3_aa_TRB = c("CASSLAPGATNEKLFF", "CASSLGQAYEQYF"),
    v_TRB = c("TRBV7-2", "TRBV6-1"),
    j_TRB = c("TRBJ1-1", "TRBJ2-1")
  )

  result <- convertToTcrdist(wide_data, chains = "alpha")

  expect_true("v_a_gene" %in% names(result))
  expect_true("cdr3_a_aa" %in% names(result))
  expect_false("v_b_gene" %in% names(result))
})

test_that("convertToTcrdist wide format beta only", {
  wide_data <- data.frame(
    barcode = c("cell_1", "cell_2"),
    cdr3_aa_TRA = c("CAVSEAPNQAGTALIF", "CAVRDSSYKLIF"),
    v_TRA = c("TRAV12-1", "TRAV8-1"),
    j_TRA = c("TRAJ15", "TRAJ12"),
    cdr3_aa_TRB = c("CASSLAPGATNEKLFF", "CASSLGQAYEQYF"),
    v_TRB = c("TRBV7-2", "TRBV6-1"),
    j_TRB = c("TRBJ1-1", "TRBJ2-1")
  )

  result <- convertToTcrdist(wide_data, chains = "beta")

  expect_true("v_b_gene" %in% names(result))
  expect_true("cdr3_b_aa" %in% names(result))
  expect_false("v_a_gene" %in% names(result))
})

test_that("convertToTcrdist resets rownames", {
  tcr_data <- create_mock_tcr_data(5)
  rownames(tcr_data) <- paste0("custom_", 1:5)

  result <- convertToTcrdist(tcr_data, chains = "beta")

  # After rownames(result) <- NULL, R assigns default sequential rownames
  expect_equal(rownames(result), as.character(seq_len(nrow(result))))
})


# ===========================================================================
# summarizeTCRrepertoire
# ===========================================================================

test_that("summarizeTCRrepertoire returns correct structure from Seurat", {
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

test_that("summarizeTCRrepertoire validates chains argument", {
  tcr_data <- create_mock_tcr_data(10)

  expect_error(summarizeTCRrepertoire(tcr_data, chains = "invalid"))
})

test_that("summarizeTCRrepertoire calculate_diversity=FALSE returns NULL diversity", {
  tcr_data <- create_mock_tcr_data(20)

  result <- summarizeTCRrepertoire(tcr_data, chains = "TRB",
                                   calculate_diversity = FALSE)

  expect_null(result$diversity)
})

test_that("summarizeTCRrepertoire top_clones has correct structure", {
  tcr_data <- data.frame(
    barcode = paste0("cell_", 1:20),
    cdr3_aa = c(rep("CASSAAA", 10), rep("CASSBBB", 5),
                rep("CASSCCC", 3), rep("CASSDDD", 2)),
    chain = "TRB"
  )

  result <- summarizeTCRrepertoire(tcr_data, chains = "TRB")

  expect_s3_class(result$top_clones, "data.frame")
  expect_true("cdr3_aa" %in% names(result$top_clones))
  expect_true("count" %in% names(result$top_clones))
  expect_true("proportion" %in% names(result$top_clones))
  expect_equal(result$top_clones$cdr3_aa[1], "CASSAAA")
  expect_equal(result$top_clones$count[1], 10L)
})

test_that("summarizeTCRrepertoire computes clonotype_ratio correctly", {
  tcr_data <- data.frame(
    barcode = paste0("cell_", 1:10),
    cdr3_aa = c(rep("CASSAAA", 5), rep("CASSBBB", 5)),
    chain = "TRB"
  )

  result <- summarizeTCRrepertoire(tcr_data, chains = "TRB")

  expect_equal(result$clonotype_ratio, 2 / 10)
})

test_that("summarizeTCRrepertoire computes CDR3 length stats", {
  tcr_data <- data.frame(
    barcode = paste0("cell_", 1:3),
    cdr3_aa = c("CASS", "CASSLGT", "CASSLAPGATNEKLFF"),
    chain = "TRB"
  )

  result <- summarizeTCRrepertoire(tcr_data, chains = "TRB")

  expect_equal(result$cdr3_length$min, 4)
  expect_equal(result$cdr3_length$max, 16)
  expect_equal(result$cdr3_length$median, 7)
})

test_that("summarizeTCRrepertoire gene_usage has v_genes and j_genes", {
  tcr_data <- data.frame(
    barcode = paste0("cell_", 1:10),
    cdr3_aa = rep("CASSAAA", 10),
    v = sample(paste0("TRBV", 1:3), 10, replace = TRUE),
    j = sample(paste0("TRBJ", 1:2), 10, replace = TRUE),
    chain = "TRB"
  )

  result <- summarizeTCRrepertoire(tcr_data, chains = "TRB")

  expect_true("v_genes" %in% names(result$gene_usage))
  expect_true("j_genes" %in% names(result$gene_usage))
  expect_s3_class(result$gene_usage$v_genes, "data.frame")
  expect_true("gene" %in% names(result$gene_usage$v_genes))
  expect_true("count" %in% names(result$gene_usage$v_genes))
  expect_true("proportion" %in% names(result$gene_usage$v_genes))
})

test_that("summarizeTCRrepertoire handles data without v/j columns", {
  tcr_data <- data.frame(
    barcode = paste0("cell_", 1:10),
    cdr3_aa = rep("CASSAAA", 10),
    chain = "TRB"
  )

  result <- summarizeTCRrepertoire(tcr_data, chains = "TRB")

  expect_length(result$gene_usage, 0)
})

test_that("summarizeTCRrepertoire diversity metrics are mathematically correct", {
  single_clone <- data.frame(
    barcode = paste0("cell_", 1:5),
    cdr3_aa = rep("CASSAAA", 5),
    chain = "TRB"
  )

  result <- summarizeTCRrepertoire(single_clone, chains = "TRB")

  expect_equal(result$diversity$shannon, 0)
  expect_equal(result$diversity$shannon_normalized, 0)
  expect_equal(result$diversity$clonality, 1)
  expect_equal(result$diversity$simpson, 1)
})

test_that("summarizeTCRrepertoire diversity for even distribution", {
  even_data <- data.frame(
    barcode = paste0("cell_", 1:10),
    cdr3_aa = rep(c("CASSAAA", "CASSBBB"), each = 5),
    chain = "TRB"
  )

  result <- summarizeTCRrepertoire(even_data, chains = "TRB")

  expect_equal(result$diversity$shannon_normalized, 1, tolerance = 1e-10)
  expect_equal(result$diversity$clonality, 0, tolerance = 1e-10)
  expect_equal(result$diversity$simpson, 0.5, tolerance = 1e-10)
})

test_that("summarizeTCRrepertoire chains field is stored", {
  tcr_data <- create_mock_tcr_data(10)

  result <- summarizeTCRrepertoire(tcr_data, chains = "TRB")

  expect_equal(result$chains, "TRB")
})

test_that("summarizeTCRrepertoire top_clones limited to 10", {
  tcr_data <- data.frame(
    barcode = paste0("cell_", 1:20),
    cdr3_aa = paste0("CASS", LETTERS[1:20]),
    chain = "TRB"
  )

  result <- summarizeTCRrepertoire(tcr_data, chains = "TRB")

  expect_lte(nrow(result$top_clones), 10)
})

test_that("summarizeTCRrepertoire top_clones with fewer than 10 clones", {
  tcr_data <- data.frame(
    barcode = paste0("cell_", 1:6),
    cdr3_aa = c("CASSA", "CASSB", "CASSC", "CASSA", "CASSB", "CASSA"),
    chain = "TRB"
  )

  result <- summarizeTCRrepertoire(tcr_data, chains = "TRB")

  expect_equal(nrow(result$top_clones), 3)
})

test_that("summarizeTCRrepertoire diversity includes all expected metrics", {
  tcr_data <- create_mock_tcr_data(50)

  result <- summarizeTCRrepertoire(tcr_data, chains = "TRB")

  expect_true("shannon" %in% names(result$diversity))
  expect_true("shannon_normalized" %in% names(result$diversity))
  expect_true("simpson" %in% names(result$diversity))
  expect_true("inverse_simpson" %in% names(result$diversity))
  expect_true("gini_simpson" %in% names(result$diversity))
  expect_true("clonality" %in% names(result$diversity))
  expect_true("richness" %in% names(result$diversity))
})

test_that("summarizeTCRrepertoire gini_simpson equals 1 - simpson", {
  tcr_data <- create_mock_tcr_data(50)

  result <- summarizeTCRrepertoire(tcr_data, chains = "TRB")

  expect_equal(result$diversity$gini_simpson,
               1 - result$diversity$simpson,
               tolerance = 1e-10)
})

test_that("summarizeTCRrepertoire inverse_simpson equals 1/simpson", {
  tcr_data <- create_mock_tcr_data(50)

  result <- summarizeTCRrepertoire(tcr_data, chains = "TRB")

  expect_equal(result$diversity$inverse_simpson,
               1 / result$diversity$simpson,
               tolerance = 1e-10)
})


# ===========================================================================
# print.TCR_summary
# ===========================================================================

test_that("print.TCR_summary works without error", {
  tcr_data <- data.frame(
    barcode = paste0("cell_", 1:10),
    cdr3_aa = c("CASSLGTGELFF", "CASSIRSSYEQYF", "CASSLGTGELFF",
                "CASSYSTGELFF", "CASSIRSSYEQYF", "CASSLGTGELFF",
                "CASNQGLNEKLFF", "CASSYSTGELFF", "CASSLGTGELFF",
                "CASSIRSSYEQYF"),
    v = rep("TRBV7-2", 10),
    j = rep("TRBJ2-2", 10),
    chain = rep("TRB", 10),
    stringsAsFactors = FALSE
  )

  result <- summarizeTCRrepertoire(tcr_data)

  expect_output(print(result), "TCR Repertoire Summary")
  expect_output(print(result), "Chain\\(s\\):")
  expect_output(print(result), "Total cells with TCR:")
  expect_output(print(result), "Unique clonotypes:")
  expect_output(print(result), "Shannon entropy:")
  expect_output(print(result), "CDR3 Length Distribution")
})

test_that("print.TCR_summary returns object invisibly", {
  tcr_data <- create_mock_tcr_data(10)
  summary_obj <- summarizeTCRrepertoire(tcr_data)

  result <- withVisible(print(summary_obj))

  expect_false(result$visible)
  expect_s3_class(result$value, "TCR_summary")
})

test_that("print.TCR_summary works without diversity", {
  tcr_data <- create_mock_tcr_data(10)
  summary_obj <- summarizeTCRrepertoire(tcr_data, calculate_diversity = FALSE)

  expect_output(print(summary_obj), "TCR Repertoire Summary")
  output <- capture.output(print(summary_obj))
  expect_false(any(grepl("Shannon entropy", output)))
})
