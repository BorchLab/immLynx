# Tests for DeepTCR VAE featurization (issue #8).
#
# DeepTCR's Load_Data explicitly does not merge identical amino acid
# sequences, so the wrapper trains on unique sequences and expands the
# resulting feature rows back onto cells. That expansion is the pure-R logic
# worth testing without Python.

# ===========================================================================
# .deeptcrFeaturesToCells
# ===========================================================================

test_that(".deeptcrFeaturesToCells expands unique-sequence rows onto cells", {
  # Two unique sequences, three cells, the first sequence used twice.
  features <- matrix(c(1, 2,
                       3, 4), nrow = 2, byrow = TRUE)
  unique_seqs <- c("A", "B")
  cell_seqs <- c("A", "B", "A")

  out <- immLynx:::.deeptcrFeaturesToCells(features, unique_seqs, cell_seqs)

  expect_identical(dim(out), c(3L, 2L))
  expect_identical(out[1, ], c(1, 2))
  expect_identical(out[2, ], c(3, 4))
  # Cells sharing a sequence must share a feature vector.
  expect_identical(out[3, ], out[1, ])
})

test_that(".deeptcrFeaturesToCells preserves cell order", {
  features <- matrix(c(10, 20, 30), nrow = 3, ncol = 1)
  unique_seqs <- c("X", "Y", "Z")
  cell_seqs <- c("Z", "X", "Y")

  out <- immLynx:::.deeptcrFeaturesToCells(features, unique_seqs, cell_seqs)

  expect_identical(as.vector(out), c(30, 10, 20))
})

test_that(".deeptcrFeaturesToCells errors when a cell sequence is unmatched", {
  features <- matrix(1:2, nrow = 2, ncol = 1)

  expect_error(
    immLynx:::.deeptcrFeaturesToCells(features, c("A", "B"), c("A", "C")),
    "unique"
  )
})

test_that(".deeptcrFeaturesToCells errors on a length mismatch", {
  features <- matrix(1:2, nrow = 2, ncol = 1)

  expect_error(
    immLynx:::.deeptcrFeaturesToCells(features, c("A", "B", "C"), c("A")),
    "rows"
  )
})

# ===========================================================================
# Argument validation, no Python required
# ===========================================================================

test_that("runDeepTCR rejects non-SingleCellExperiment input", {
  expect_error(runDeepTCR(data.frame(x = 1)), "SingleCellExperiment")
})

test_that("runDeepTCR rejects an invalid latent_dim", {
  data("immLynx_example", package = "immLynx")

  expect_error(runDeepTCR(immLynx_example, latent_dim = 0), "latent_dim")
  expect_error(runDeepTCR(immLynx_example, latent_dim = -5), "latent_dim")
})

# ===========================================================================
# End-to-end, requires the DeepTCR environment
# ===========================================================================

test_that("runDeepTCR adds a reduction of the requested width", {
  skip_if_no_deeptcr()
  data("immLynx_example", package = "immLynx")

  sce <- runDeepTCR(immLynx_example, chains = "TRB",
                    latent_dim = 16, return_object = TRUE)

  expect_true("tcr_deeptcr" %in% SingleCellExperiment::reducedDimNames(sce))
  rd <- SingleCellExperiment::reducedDim(sce, "tcr_deeptcr")
  expect_identical(nrow(rd), ncol(immLynx_example))
  expect_lte(ncol(rd), 16L)
  # Cells carrying a TRB sequence must have finite features.
  expect_true(any(stats::complete.cases(rd)))
})

test_that("runDeepTCR returns raw features when return_object is FALSE", {
  skip_if_no_deeptcr()
  data("immLynx_example", package = "immLynx")

  res <- runDeepTCR(immLynx_example, chains = "TRB",
                    latent_dim = 16, return_object = FALSE)

  expect_type(res, "list")
  expect_true(all(c("features", "sequences", "explained_variance_ratio") %in%
                    names(res)))
  expect_true(is.matrix(res$features))
  # One feature row per unique sequence trained on.
  expect_identical(nrow(res$features), length(res$sequences))
})

test_that("runDeepTCR leaves no model directory behind", {
  skip_if_no_deeptcr()
  data("immLynx_example", package = "immLynx")

  before <- list.files(getwd(), all.files = TRUE)
  invisible(runDeepTCR(immLynx_example, chains = "TRB",
                       latent_dim = 16, return_object = FALSE))
  after <- list.files(getwd(), all.files = TRUE)

  expect_identical(setdiff(after, before), character(0))
})
