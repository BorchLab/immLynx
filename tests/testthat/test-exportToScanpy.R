# Tests for exportToScanpy() and its internal helpers

# ===========================================================================
# .pruneColData() — internal helper
# ===========================================================================

test_that(".pruneColData drops list-columns", {
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(0L, 2, 3)),
    colData = S4Vectors::DataFrame(
      keep = c("a", "b", "c"),
      drop_me = I(list(1:2, 1:3, 1:4))
    )
  )

  out <- immLynx:::.pruneColData(sce, obs_columns = NULL, drop_ct = FALSE)
  expect_true("keep" %in% colnames(SummarizedExperiment::colData(out)))
  expect_false("drop_me" %in% colnames(SummarizedExperiment::colData(out)))
})

test_that(".pruneColData drops CT* columns when drop_ct = TRUE", {
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(0L, 2, 2)),
    colData = S4Vectors::DataFrame(
      keep   = c("a", "b"),
      CTgene = c("TRBV1.TRBJ1", "TRBV2.TRBJ2"),
      CTaa   = c("CASS", "CASR")
    )
  )

  kept <- immLynx:::.pruneColData(sce, obs_columns = NULL, drop_ct = FALSE)
  expect_true(all(c("CTgene", "CTaa") %in%
                  colnames(SummarizedExperiment::colData(kept))))

  dropped <- immLynx:::.pruneColData(sce, obs_columns = NULL, drop_ct = TRUE)
  expect_false(any(c("CTgene", "CTaa") %in%
                   colnames(SummarizedExperiment::colData(dropped))))
})

test_that(".pruneColData honors obs_columns subset", {
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(0L, 2, 2)),
    colData = S4Vectors::DataFrame(a = 1:2, b = 3:4, c = 5:6)
  )

  out <- immLynx:::.pruneColData(sce, obs_columns = c("a", "c"),
                                 drop_ct = FALSE)
  expect_equal(colnames(SummarizedExperiment::colData(out)),
               c("a", "c"))
})

test_that(".pruneColData errors on unknown obs_columns", {
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(0L, 2, 2)),
    colData = S4Vectors::DataFrame(a = 1:2, b = 3:4)
  )

  expect_error(
    immLynx:::.pruneColData(sce, obs_columns = c("a", "missing"),
                            drop_ct = FALSE),
    "missing"
  )
})

test_that(".pruneColData preserves obs_columns ordering", {
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(0L, 2, 2)),
    colData = S4Vectors::DataFrame(a = 1:2, b = 3:4, c = 5:6)
  )

  out <- immLynx:::.pruneColData(sce, obs_columns = c("c", "a"),
                                 drop_ct = FALSE)
  expect_equal(colnames(SummarizedExperiment::colData(out)),
               c("c", "a"))
})

test_that(".pruneColData drops only exact CT columns, not lookalikes", {
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(0L, 2, 2)),
    colData = S4Vectors::DataFrame(
      CTgene = c("TRBV1.TRBJ1", "TRBV2.TRBJ2"),
      CTLA4  = c(0.5, 1.2),
      CTSL   = c(2.3, 0.1)
    )
  )

  out <- immLynx:::.pruneColData(sce, obs_columns = NULL, drop_ct = TRUE)
  cn <- colnames(SummarizedExperiment::colData(out))
  expect_false("CTgene" %in% cn)
  expect_true(all(c("CTLA4", "CTSL") %in% cn))
})

test_that(".pruneColData drops CT* before obs_columns subset", {
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(0L, 2, 2)),
    colData = S4Vectors::DataFrame(
      keep   = c("a", "b"),
      CTgene = c("TRBV1.TRBJ1", "TRBV2.TRBJ2")
    )
  )

  # When the user asks for a CT* column with drop_ct = TRUE, the column is
  # already gone by the time obs_columns is checked, so we error clearly.
  expect_error(
    immLynx:::.pruneColData(sce,
                            obs_columns = c("keep", "CTgene"),
                            drop_ct = TRUE),
    "CTgene"
  )
})

# ===========================================================================
# .buildAIRR() — internal helper
# ===========================================================================
# Mock SCE constructors live in tests/testthat/helper-immLynx.R:
#   - mock_tcr_sce(): paired/single TRA+TRB cells plus a chainless cell
#   - mock_bcr_sce(): IGH + IGK + IGL combinations

test_that(".buildAIRR returns AIRR-schema columns", {
  sce <- mock_tcr_sce()
  airr <- immLynx:::.buildAIRR(sce, chains = "both")

  required <- c("cell_id", "sequence_id", "locus", "productive",
                "v_call", "j_call", "junction_aa", "junction")
  expect_true(all(required %in% colnames(airr)))
})

test_that(".buildAIRR returns rows for cells with chains", {
  sce <- mock_tcr_sce()
  airr <- immLynx:::.buildAIRR(sce, chains = "both")

  # cell c1 has 2 chains, c2 has 1, c3 has 0, c4 has 1 -> 4 rows
  expect_equal(nrow(airr), 4L)
  expect_setequal(unique(airr$cell_id), c("c1", "c2", "c4"))
})

test_that(".buildAIRR sequence_ids are unique across loci", {
  sce <- mock_tcr_sce()
  airr <- immLynx:::.buildAIRR(sce, chains = "both")
  # c1 contributes both a TRA and a TRB row; sequence_ids must differ.
  expect_equal(anyDuplicated(airr$sequence_id), 0L)
})

test_that(".buildAIRR locus values are valid AIRR loci", {
  sce <- mock_tcr_sce()
  airr <- immLynx:::.buildAIRR(sce, chains = "both")
  expect_true(all(airr$locus %in%
                  c("TRA", "TRB", "TRG", "TRD",
                    "IGH", "IGK", "IGL")))
})

test_that(".buildAIRR returns NULL when no TCR data present", {
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(0L, 2, 2,
                                  dimnames = list(c("g1", "g2"),
                                                  c("c1", "c2"))))
  )
  expect_null(immLynx:::.buildAIRR(sce, chains = "both"))
})

test_that(".buildAIRR populates junction (nt) when CTnt is present", {
  sce <- mock_tcr_sce()
  airr <- immLynx:::.buildAIRR(sce, chains = "both")
  expect_true(any(!is.na(airr$junction) & nzchar(airr$junction)))
})

test_that(".buildAIRR derives productive from CDR3 sequence", {
  sce <- mock_tcr_sce()
  airr <- immLynx:::.buildAIRR(sce, chains = "both")
  # All mock CDR3s lack stop codons or frameshift markers, so they
  # should all be productive = TRUE.
  expect_true(all(airr$productive))
  expect_type(airr$productive, "logical")
})

test_that(".buildAIRR honors a single-locus chain request", {
  sce <- mock_tcr_sce()
  airr <- immLynx:::.buildAIRR(sce, chains = "TRB")
  expect_true(all(airr$locus == "TRB"))
  # c1 (TRB), c2 (TRB) — c4 has only TRA, c3 has nothing.
  expect_setequal(unique(airr$cell_id), c("c1", "c2"))
})

test_that(".buildAIRR returns NULL when all CDR3s are NA for a locus", {
  # TCR present in CTgene but no AA sequences anywhere
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(0L, 2, 2,
                                  dimnames = list(c("g1", "g2"),
                                                  c("c1", "c2")))),
    colData = S4Vectors::DataFrame(
      CTgene = c("TRBV1.TRBJ1", "TRBV2.TRBJ2"),
      CTaa   = c(NA_character_, NA_character_),
      CTnt   = c(NA_character_, NA_character_),
      row.names = c("c1", "c2")
    )
  )
  expect_null(immLynx:::.buildAIRR(sce, chains = "TRB"))
})

test_that(".buildAIRR maps IGH to immApex 'Heavy' and emits IGH locus", {
  sce <- mock_bcr_sce()
  airr <- immLynx:::.buildAIRR(sce, chains = "IGH")
  expect_true(all(airr$locus == "IGH"))
  # b1 (paired) and b2 (heavy only) both contribute IGH; b3 is lambda.
  expect_setequal(unique(airr$cell_id), c("b1", "b2"))
})

test_that(".buildAIRR disambiguates IGK vs IGL by V-gene prefix", {
  sce <- mock_bcr_sce()
  airr_k <- immLynx:::.buildAIRR(sce, chains = "IGK")
  airr_l <- immLynx:::.buildAIRR(sce, chains = "IGL")
  expect_true(all(airr_k$locus == "IGK"))
  expect_true(all(airr_l$locus == "IGL"))
  expect_true(all(startsWith(airr_k$v_call, "IGK")))
  expect_true(all(startsWith(airr_l$v_call, "IGL")))
})

test_that(".airr_locus_to_immapex translates BCR locus names", {
  expect_equal(immLynx:::.airr_locus_to_immapex("IGH"), "Heavy")
  expect_equal(immLynx:::.airr_locus_to_immapex("IGK"), "Light")
  expect_equal(immLynx:::.airr_locus_to_immapex("IGL"), "Light")
  expect_equal(immLynx:::.airr_locus_to_immapex("TRB"), "TRB")
})

# ===========================================================================
# exportToScanpy() — H5AD path (round-trip)
# ===========================================================================

# Tiny SCE constructor for round-trip integration tests.
.tiny_sce <- function() {
  m <- Matrix::rsparsematrix(50, 30, density = 0.1,
                             rand.x = function(n) rpois(n, lambda = 2))
  rownames(m) <- paste0("g", seq_len(50))
  colnames(m) <- paste0("c", seq_len(30))
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = m),
    colData = S4Vectors::DataFrame(
      sample  = rep(c("A", "B"), length.out = 30),
      n_count = colSums(as.matrix(m)),
      row.names = colnames(m)
    )
  )
  SingleCellExperiment::reducedDim(sce, "PCA") <-
    matrix(rnorm(30 * 5), 30, 5)
  SingleCellExperiment::reducedDim(sce, "UMAP") <-
    matrix(rnorm(30 * 2), 30, 2)
  sce
}

test_that("exportToScanpy writes a readable H5AD (round-trip)", {
  skip_if_not_installed("zellkonverter")
  skip_if_not_installed("Matrix")
  skip_on_bioc_build()  # writeH5AD spins basilisk on first call

  sce <- .tiny_sce()
  out_dir <- tempfile("expscan_")

  res <- exportToScanpy(sce, output_dir = out_dir,
                       format = "h5ad", write_airr = FALSE,
                       verbose = FALSE)

  expect_true(file.exists(res$h5ad))
  expect_null(res$airr)
  expect_null(res$h5mu)

  back <- zellkonverter::readH5AD(res$h5ad)
  expect_equal(dim(back), dim(sce))
  expect_true(all(c("sample", "n_count") %in%
                  colnames(SummarizedExperiment::colData(back))))
  # zellkonverter writes reductions to obsm with the same names.
  expect_true(any(c("PCA", "X_PCA") %in%
                  SingleCellExperiment::reducedDimNames(back)))
})

test_that("exportToScanpy errors when H5AD exists and overwrite = FALSE", {
  skip_if_not_installed("zellkonverter")
  skip_if_not_installed("Matrix")
  skip_on_bioc_build()

  sce <- .tiny_sce()
  out_dir <- tempfile("expscan_")
  dir.create(out_dir)
  file.create(file.path(out_dir, "adata.h5ad"))

  expect_error(
    exportToScanpy(sce, output_dir = out_dir,
                   format = "h5ad", write_airr = FALSE,
                   overwrite = FALSE, verbose = FALSE),
    "exists"
  )
})

test_that("exportToScanpy errors on unknown assay", {
  skip_if_not_installed("Matrix")
  sce <- .tiny_sce()
  out_dir <- tempfile("expscan_")

  expect_error(
    exportToScanpy(sce, output_dir = out_dir,
                   format = "h5ad", assay = "missing",
                   write_airr = FALSE, verbose = FALSE),
    "assay"
  )
})

test_that("exportToScanpy errors on unknown reductions", {
  skip_if_not_installed("Matrix")
  sce <- .tiny_sce()
  out_dir <- tempfile("expscan_")

  expect_error(
    exportToScanpy(sce, output_dir = out_dir,
                   format = "h5ad",
                   reductions = c("PCA", "MISSING"),
                   write_airr = FALSE, verbose = FALSE),
    "MISSING"
  )
})

test_that("exportToScanpy errors on h5mu request without TCR data", {
  skip_if_not_installed("Matrix")
  sce <- .tiny_sce()  # has no CT* columns
  out_dir <- tempfile("expscan_")

  expect_error(
    exportToScanpy(sce, output_dir = out_dir,
                   format = "h5mu", verbose = FALSE),
    "immune receptor"
  )
})

test_that("exportToScanpy rejects non-SCE / non-Seurat input", {
  expect_error(
    exportToScanpy(data.frame(x = 1), output_dir = tempfile(),
                   format = "h5ad", verbose = FALSE),
    "SingleCellExperiment"
  )
})

test_that("exportToScanpy rejects NULL or invalid output_dir", {
  skip_if_not_installed("Matrix")
  sce <- .tiny_sce()

  expect_error(
    exportToScanpy(sce, output_dir = NULL,
                   format = "h5ad", verbose = FALSE),
    "output_dir"
  )
  expect_error(
    exportToScanpy(sce, output_dir = NA_character_,
                   format = "h5ad", verbose = FALSE),
    "output_dir"
  )
  expect_error(
    exportToScanpy(sce, output_dir = c("a", "b"),
                   format = "h5ad", verbose = FALSE),
    "output_dir"
  )
})

test_that("exportToScanpy errors when output_dir is a file, not a dir", {
  skip_if_not_installed("Matrix")
  sce <- .tiny_sce()
  out_path <- tempfile("expscan_file_")
  file.create(out_path)
  on.exit(unlink(out_path), add = TRUE)

  expect_error(
    exportToScanpy(sce, output_dir = out_path,
                   format = "h5ad", verbose = FALSE),
    "not a directory"
  )
})

test_that("exportToScanpy rejects empty reductions", {
  skip_if_not_installed("Matrix")
  sce <- .tiny_sce()

  expect_error(
    exportToScanpy(sce, output_dir = tempfile("expscan_"),
                   format = "h5ad",
                   reductions = character(0),
                   write_airr = FALSE, verbose = FALSE),
    "reductions"
  )
})

test_that("exportToScanpy overwrite = TRUE replaces a stale H5AD", {
  skip_if_not_installed("zellkonverter")
  skip_if_not_installed("Matrix")
  skip_on_bioc_build()

  sce <- .tiny_sce()
  out_dir <- tempfile("expscan_overwrite_")
  dir.create(out_dir)
  stale <- file.path(out_dir, "adata.h5ad")
  writeLines("not a real h5ad", stale)
  stale_size <- file.size(stale)

  res <- exportToScanpy(sce, output_dir = out_dir,
                       format = "h5ad", write_airr = FALSE,
                       overwrite = TRUE, verbose = FALSE)

  expect_true(file.exists(res$h5ad))
  # The real H5AD must be much larger than the stub we planted, and must
  # actually round-trip through zellkonverter.
  expect_gt(file.size(res$h5ad), stale_size)
  back <- zellkonverter::readH5AD(res$h5ad)
  expect_equal(dim(back), dim(sce))
})

test_that("exportToScanpy round-trips factor columns in colData", {
  skip_if_not_installed("zellkonverter")
  skip_if_not_installed("Matrix")
  skip_on_bioc_build()

  sce <- .tiny_sce()
  SummarizedExperiment::colData(sce)$cluster <-
    factor(rep(c("Tcell", "Bcell", "NK"), length.out = ncol(sce)))
  out_dir <- tempfile("expscan_factor_")

  res <- exportToScanpy(sce, output_dir = out_dir,
                       format = "h5ad", write_airr = FALSE,
                       verbose = FALSE)

  back <- zellkonverter::readH5AD(res$h5ad)
  cd_back <- SummarizedExperiment::colData(back)
  expect_true("cluster" %in% colnames(cd_back))
  # zellkonverter encodes factors as anndata categoricals; the readback
  # is a factor with the same level set (order may differ).
  expect_setequal(levels(cd_back$cluster), c("Tcell", "Bcell", "NK"))
})

test_that("exportToScanpy accepts a Seurat input via .coerceToSCE", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("zellkonverter")
  skip_if_not_installed("Matrix")
  skip_on_bioc_build()

  m <- matrix(rpois(40 * 25, lambda = 2), 40, 25,
              dimnames = list(paste0("g", seq_len(40)),
                              paste0("c", seq_len(25))))
  sobj <- Seurat::CreateSeuratObject(counts = m)
  sobj$sample <- rep(c("A", "B"), length.out = 25)

  out_dir <- tempfile("expscan_seurat_")
  res <- exportToScanpy(sobj, output_dir = out_dir,
                       format = "h5ad", write_airr = FALSE,
                       verbose = FALSE)

  expect_true(file.exists(res$h5ad))
  back <- zellkonverter::readH5AD(res$h5ad)
  expect_equal(ncol(back), 25L)
  expect_true("sample" %in%
              colnames(SummarizedExperiment::colData(back)))
})

# ===========================================================================
# AIRR sidecar
# ===========================================================================

test_that("exportToScanpy writes a parseable AIRR TSV sidecar", {
  skip_if_not_installed("zellkonverter")
  skip_if_not_installed("data.table")
  skip_if_not_installed("Matrix")
  skip_on_bioc_build()

  data("immLynx_example", package = "immLynx")
  out_dir <- tempfile("expscan_airr_")

  res <- exportToScanpy(immLynx_example, output_dir = out_dir,
                       format = "h5ad", write_airr = TRUE,
                       verbose = FALSE)

  expect_true(file.exists(res$airr))
  airr <- data.table::as.data.table(
    read.delim(gzfile(res$airr), sep = "\t", check.names = FALSE,
               stringsAsFactors = FALSE)
  )
  required <- c("cell_id", "sequence_id", "locus", "productive",
                "v_call", "j_call", "junction_aa", "junction")
  expect_true(all(required %in% colnames(airr)))
  expect_gt(nrow(airr), 0L)
})

test_that("exportToScanpy warns and skips AIRR when no TCR present", {
  skip_if_not_installed("zellkonverter")
  skip_if_not_installed("Matrix")
  skip_on_bioc_build()

  sce <- .tiny_sce()  # no CT* fields
  out_dir <- tempfile("expscan_noairr_")

  expect_warning(
    res <- exportToScanpy(sce, output_dir = out_dir,
                         format = "h5ad", write_airr = TRUE,
                         verbose = FALSE),
    "no immune receptor"
  )
  expect_null(res$airr)
})

test_that("exportToScanpy AIRR row count matches .buildAIRR output", {
  skip_if_not_installed("zellkonverter")
  skip_if_not_installed("data.table")
  skip_if_not_installed("Matrix")
  skip_on_bioc_build()

  data("immLynx_example", package = "immLynx")
  out_dir <- tempfile("expscan_airrcount_")

  res <- exportToScanpy(immLynx_example, output_dir = out_dir,
                       format = "h5ad", write_airr = TRUE,
                       chains = "TRB", verbose = FALSE)

  written <- data.table::as.data.table(
    read.delim(gzfile(res$airr), sep = "\t", check.names = FALSE,
               stringsAsFactors = FALSE)
  )
  internal <- immLynx:::.buildAIRR(immLynx_example, chains = "TRB")
  expect_equal(nrow(written), nrow(internal))
  expect_equal(sort(unique(written$locus)),
               sort(unique(internal$locus)))
})

# ===========================================================================
# H5MU round-trip (gated on scanpyExportEnv)
# ===========================================================================

test_that("exportToScanpy writes a readable H5MU when TCR data is present", {
  skip_if_no_scanpy_env()
  skip_if_not_installed("zellkonverter")
  skip_if_not_installed("data.table")

  data("immLynx_example", package = "immLynx")
  out_dir <- tempfile("expscan_h5mu_")

  res <- exportToScanpy(immLynx_example, output_dir = out_dir,
                       format = "h5mu", write_airr = FALSE,
                       verbose = FALSE)

  expect_true(file.exists(res$h5mu))

  # Read back via the same env to confirm structure.
  proc <- basilisk::basiliskStart(immLynx:::scanpyExportEnv)
  on.exit(basilisk::basiliskStop(proc))
  ok <- basilisk::basiliskRun(proc, function(path) {
    mu <- reticulate::import("muon")
    mdata <- mu$read(path)
    list(
      modalities = names(mdata$mod),
      n_airr = mdata$mod$airr$n_obs
    )
  }, path = res$h5mu)

  expect_true(all(c("gex", "airr") %in% ok$modalities))
  expect_gt(ok$n_airr, 0L)
})

test_that("exportToScanpy h5mu + write_airr = TRUE writes both files", {
  skip_if_no_scanpy_env()
  skip_if_not_installed("zellkonverter")
  skip_if_not_installed("data.table")

  data("immLynx_example", package = "immLynx")
  out_dir <- tempfile("expscan_h5mu_both_")

  res <- exportToScanpy(immLynx_example, output_dir = out_dir,
                       format = "h5mu", write_airr = TRUE,
                       verbose = FALSE)

  expect_true(file.exists(res$h5ad))
  expect_true(file.exists(res$airr))
  expect_true(file.exists(res$h5mu))
})
