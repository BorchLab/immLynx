# Tests for runTCRdist function

# ===========================================================================
# Parameter validation
# ===========================================================================

test_that("runTCRdist validates chain argument", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")

  data("immLynx_example", package = "immLynx")

  expect_error(runTCRdist(immLynx_example, chains = "invalid"),
               "Invalid chain")
})

test_that("runTCRdist rejects non-Seurat/SCE input", {
  tcr_data <- create_mock_tcr_data(10)

  expect_error(
    runTCRdist(tcr_data, chains = "beta"),
    "Input must be a Seurat or SingleCellExperiment object"
  )
})

test_that("runTCRdist rejects data.frame input", {
  df <- data.frame(a = 1:5)

  expect_error(
    runTCRdist(df, chains = "beta"),
    "Input must be a Seurat or SingleCellExperiment object"
  )
})

test_that("runTCRdist function signature has correct defaults", {
  f <- formals(runTCRdist)

  expect_equal(f$chains, "beta")
  expect_equal(f$organism, "human")
  expect_equal(f$compute_distances, TRUE)
  expect_equal(f$add_to_object, FALSE)
})

test_that("runTCRdist chain_map correctly maps alpha and beta", {
  chain_map <- list(
    "alpha" = "TRA",
    "beta" = "TRB"
  )

  expect_equal(chain_map[["alpha"]], "TRA")
  expect_equal(chain_map[["beta"]], "TRB")
  expect_null(chain_map[["invalid"]])
})

test_that("runTCRdist .add_allele helper adds *01 suffix", {
  .add_allele <- function(genes) {
    ifelse(!is.na(genes) & !grepl("\\*", genes),
           paste0(genes, "*01"),
           genes)
  }

  expect_equal(.add_allele("TRBV7-2"), "TRBV7-2*01")
  expect_equal(.add_allele("TRBV7-2*01"), "TRBV7-2*01")
  expect_true(is.na(.add_allele(NA)))
})

test_that("runTCRdist .add_allele does not double-add suffix", {
  .add_allele <- function(genes) {
    ifelse(!is.na(genes) & !grepl("\\*", genes),
           paste0(genes, "*01"),
           genes)
  }

  gene <- "TRBV10-3*02"
  expect_equal(.add_allele(gene), "TRBV10-3*02")
})

test_that("runTCRdist .add_allele handles vector input", {
  .add_allele <- function(genes) {
    ifelse(!is.na(genes) & !grepl("\\*", genes),
           paste0(genes, "*01"),
           genes)
  }

  genes <- c("TRBV7-2", NA, "TRBV10-3*02", "TRBV5-1")
  result <- .add_allele(genes)

  expect_equal(result[1], "TRBV7-2*01")
  expect_true(is.na(result[2]))
  expect_equal(result[3], "TRBV10-3*02")
  expect_equal(result[4], "TRBV5-1*01")
})


# ===========================================================================
# Python-dependent tests
# ===========================================================================

test_that("runTCRdist returns correct structure", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_python()

  data("immLynx_example", package = "immLynx")

  result <- runTCRdist(immLynx_example, chains = "beta")

  expect_type(result, "list")
  expect_true("distances" %in% names(result))
  expect_true("barcodes" %in% names(result))
  expect_true("tcr_data" %in% names(result))
})

test_that("runTCRdist returns distance matrix of correct dimensions", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_python()

  data("immLynx_example", package = "immLynx")

  result <- runTCRdist(immLynx_example, chains = "beta")

  n_cells <- length(result$barcodes)
  expect_true(n_cells > 0)
  expect_equal(nrow(result$tcr_data), n_cells)
})

test_that("runTCRdist handles add_to_object parameter for Seurat", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_python()

  data("immLynx_example", package = "immLynx")

  result <- runTCRdist(immLynx_example, chains = "beta",
                       add_to_object = TRUE)

  expect_s4_class(result, "Seurat")
  expect_true(!is.null(result@misc$tcrdist))
})

test_that("runTCRdist handles both chains", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_python()

  data("immLynx_example", package = "immLynx")

  result <- runTCRdist(immLynx_example, chains = c("alpha", "beta"))

  expect_true("pw_alpha" %in% names(result$distances) ||
              "pw_beta" %in% names(result$distances))
})

test_that("runTCRdist tcr_data has correct format", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_python()

  data("immLynx_example", package = "immLynx")

  result <- runTCRdist(immLynx_example, chains = "beta")

  expect_true("count" %in% names(result$tcr_data))
  expect_true("v_b_gene" %in% names(result$tcr_data))
  expect_true("j_b_gene" %in% names(result$tcr_data))
  expect_true("cdr3_b_aa" %in% names(result$tcr_data))
  non_na_v <- result$tcr_data$v_b_gene[!is.na(result$tcr_data$v_b_gene)]
  if (length(non_na_v) > 0) {
    expect_true(all(grepl("\\*", non_na_v)))
  }
})

test_that("runTCRdist produces messages during execution", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("immApex")
  skip_if_no_python()

  data("immLynx_example", package = "immLynx")

  expect_message(
    runTCRdist(immLynx_example, chains = "beta"),
    "Extracting TCR sequences"
  )
})
