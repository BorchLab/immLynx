# Helper functions for immLynx tests

# Create mock TCR data for testing
create_mock_tcr_data <- function(n = 100) {
  # Common CDR3 sequences
  cdr3_pool <- c(
    "CASSLAPGATNEKLFF", "CASSLGQAYEQYF", "CASRLAGQETQYF",
    "CASSYSGGNTGELFF", "CASSQDRTGQETQYF", "CASSLNRDNEQFF",
    "CASSLTGTEAFF", "CASSYSQGSYEQYF", "CASSLAGDTDTQYF",
    "CASSLVSGSTDTQYF", "CASSQETQYF", "CASSLGANTGELFF"
  )

  v_pool <- paste0("TRBV", c("5-1", "6-1", "7-2", "12-3", "20-1", "28"))
  j_pool <- paste0("TRBJ", c("1-1", "1-2", "2-1", "2-3", "2-5", "2-7"))

  data.frame(
    barcode = paste0("cell_", seq_len(n)),
    cdr3_aa = sample(cdr3_pool, n, replace = TRUE),
    v = sample(v_pool, n, replace = TRUE),
    j = sample(j_pool, n, replace = TRUE),
    chain = "TRB",
    stringsAsFactors = FALSE
  )
}

# Create mock paired alpha-beta data
create_mock_paired_data <- function(n = 50) {
  alpha_cdr3 <- c(
    "CAVSEAPNQAGTALIF", "CAVRDSSYKLIF", "CAGQTGGFKTIF",
    "CALSDNNARLMF", "CAVSEGGSYIPTF", "CAVNGGSQGNLIF"
  )

  beta_cdr3 <- c(
    "CASSLAPGATNEKLFF", "CASSLGQAYEQYF", "CASRLAGQETQYF",
    "CASSYSGGNTGELFF", "CASSQDRTGQETQYF", "CASSLNRDNEQFF"
  )

  data.frame(
    barcode = paste0("cell_", seq_len(n)),
    cdr3_aa_TRA = sample(alpha_cdr3, n, replace = TRUE),
    v_TRA = sample(paste0("TRAV", 1:6), n, replace = TRUE),
    j_TRA = sample(paste0("TRAJ", 1:6), n, replace = TRUE),
    cdr3_aa_TRB = sample(beta_cdr3, n, replace = TRUE),
    v_TRB = sample(paste0("TRBV", 1:6), n, replace = TRUE),
    j_TRB = sample(paste0("TRBJ", 1:6), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

# Check if Python environment is available
python_available <- function() {
  tryCatch({
    proc <- basilisk::basiliskStart(immLynx:::immLynxEnv)
    on.exit(basilisk::basiliskStop(proc))
    TRUE
  }, error = function(e) FALSE)
}

# Skip on Bioconductor build machines (15-min timeout too tight for Python tests)
skip_on_bioc_build <- function() {
  if (nzchar(Sys.getenv("IS_BIOC_BUILD_MACHINE")) ||
      nzchar(Sys.getenv("BBS_HOME"))) {
    testthat::skip("Skipping Python tests on Bioconductor build machine")
  }
}

# Skip test if Python not available or on Bioc build machine
skip_if_no_python <- function() {
  skip_on_bioc_build()
  if (!python_available()) {
    testthat::skip("Python environment not available")
  }
}

# Alias: skip_if_no_transformers now delegates to skip_if_no_python
# since transformers and torch are included in the basilisk environment
skip_if_no_transformers <- function() {
  skip_if_no_python()
}

# Cache for module import checks (avoids repeated basiliskRun calls)
.module_check_cache <- new.env(parent = emptyenv())

# Check if a Python module can actually be imported (catches shared lib issues)
can_import_module <- function(module_name) {
  if (!is.null(.module_check_cache[[module_name]])) {
    return(.module_check_cache[[module_name]])
  }
  result <- tryCatch({
    basilisk::basiliskRun(env = immLynx:::immLynxEnv, fun = function(mod) {
      reticulate::import(mod)
      TRUE
    }, mod = module_name)
  }, error = function(e) FALSE)
  .module_check_cache[[module_name]] <- result
  result
}

# Skip if tcrdist (parasail shared library) can't be loaded
skip_if_no_tcrdist <- function() {
  skip_if_no_python()
  if (!can_import_module("tcrdist.repertoire")) {
    testthat::skip("tcrdist/parasail shared library not loadable")
  }
}

# Skip if metaclonotypist can't be loaded
skip_if_no_metaclonotypist <- function() {
  skip_if_no_python()
  if (!can_import_module("metaclonotypist")) {
    testthat::skip("metaclonotypist shared library not loadable")
  }
}

# Skip test if the scanpy export basilisk env is not available.
# This env is separate from immLynxEnv (anndata + scanpy + muon + scirpy
# stack), so it has its own first-use install cost.
skip_if_no_scanpy_env <- function() {
  skip_on_bioc_build()
  ok <- tryCatch({
    proc <- basilisk::basiliskStart(immLynx:::scanpyExportEnv)
    on.exit(basilisk::basiliskStop(proc))
    TRUE
  }, error = function(e) FALSE)
  if (!ok) testthat::skip("scanpyExportEnv not available")
}

# Skip test if the scXpand basilisk env is not available. This env carries
# python 3.11 + PyTorch and runs to several gigabytes, so its first-use
# install cost is much larger than the other two.
skip_if_no_scxpand_env <- function() {
  skip_on_bioc_build()
  ok <- tryCatch({
    proc <- basilisk::basiliskStart(immLynx:::scXpandEnv)
    on.exit(basilisk::basiliskStop(proc))
    TRUE
  }, error = function(e) FALSE)
  if (!ok) testthat::skip("scXpandEnv not available")
}

# Inference additionally downloads a pretrained model from figshare on
# first use. Gate that behind an explicit opt-in so it never fires
# unattended in CI.
skip_if_no_scxpand_model <- function() {
  skip_if_no_scxpand_env()
  if (!nzchar(Sys.getenv("IMMLYNX_TEST_SCXPAND"))) {
    testthat::skip("set IMMLYNX_TEST_SCXPAND=1 to run scXpand model tests")
  }
}

# Mock SCE indexed by Ensembl IDs with raw integer counts — the shape
# runScXpand expects. Clone calls follow the scRepertoire "_" contract.
mock_ensembl_sce <- function(n_genes = 50, n_cells = 20, seed = 1) {
  set.seed(seed)
  counts <- Matrix::Matrix(
    matrix(stats::rpois(n_genes * n_cells, lambda = 3), n_genes, n_cells),
    sparse = TRUE)
  rownames(counts) <- sprintf("ENSG%011d", seq_len(n_genes))
  colnames(counts) <- paste0("cell", seq_len(n_cells))

  # Skewed clone sizes: one dominant clone, a couple of medium ones, the
  # rest singletons. A flat distribution would make every cell fall on the
  # same side of the 1.5x median cutoff and leave nothing to score.
  n_big <- max(2L, n_cells %/% 4L)
  n_mid <- max(1L, n_cells %/% 10L)
  clone_idx <- c(rep(1L, n_big), rep(2L, n_mid), rep(3L, n_mid),
                 seq.int(4L, length.out = max(0L, n_cells - n_big -
                                                2L * n_mid)))
  clone_idx <- clone_idx[seq_len(n_cells)]
  clone <- paste0("CAS", clone_idx, "F_CAS", clone_idx, "F")

  SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(
      sample = rep(c("S1", "S2"), length.out = n_cells),
      CTstrict = clone,
      row.names = colnames(counts)
    )
  )
}

# Mock SCE whose clone layout makes the two median_basis readings diverge,
# so the tests can pin the arithmetic by hand.
#
#   sample A: c1 x 4, c2 x 1, c3 x 1
#     median over unique clones = median(4, 1, 1) = 1 -> cutoff 1.5
#       -> only c1 expanded
#     median over cells = median(4, 4, 4, 4, 1, 1) = 4 -> cutoff 6
#       -> nothing expanded
#   sample B: d1 x 2, d2 x 2, plus one cell with no clone call
mock_clonal_sce <- function() {
  clone <- c(rep("c1", 4), "c2", "c3", rep("d1", 2), rep("d2", 2), NA)
  sample <- c(rep("A", 6), rep("B", 5))
  n <- length(clone)

  # Dimnames set up front, and values varied: Matrix() collapses a square
  # constant matrix to a symmetric dsCMatrix, which ties rownames to
  # colnames and would drop the gene IDs.
  counts <- Matrix::Matrix(
    matrix(seq_len(4L * n), nrow = 4L,
           dimnames = list(sprintf("ENSG%011d", seq_len(4L)),
                           paste0("cell", seq_len(n)))),
    sparse = TRUE)

  SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(
      CTstrict = clone, sample = sample,
      row.names = colnames(counts)
    )
  )
}

# Mock SCE with hand-crafted scRepertoire CT* fields. Used by the
# exportToScanpy / .buildAIRR test suite. scRepertoire uses "_" as the
# chain-slot separator (TRA before, TRB after) with explicit "NA" tokens
# for missing chains; mocks must match that contract.
mock_tcr_sce <- function() {
  ct_gene <- c(
    "TRAV1.TRAJ1.TRAC_TRBV1.TRBD1.TRBJ1.TRBC1",  # paired
    "NA_TRBV2.TRBD2.TRBJ2.TRBC2",                # beta only
    NA,                                           # no TCR
    "TRAV3.TRAJ3.TRAC_NA"                         # alpha only
  )
  ct_aa <- c("CASA_CASB", "NA_CASB2", NA, "CASA3_NA")
  ct_nt <- c("TGTGCA_TGTGCB", "NA_TGTGCB2", NA, "TGTGCA3_NA")

  SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(0L, 2, 4,
                                  dimnames = list(c("g1", "g2"),
                                                  paste0("c", 1:4)))),
    colData = S4Vectors::DataFrame(
      CTgene = ct_gene, CTaa = ct_aa, CTnt = ct_nt,
      row.names = paste0("c", 1:4)
    )
  )
}

# Mock BCR SCE — heavy + light chain entries. Used by .buildAIRR BCR tests.
mock_bcr_sce <- function() {
  ct_gene <- c(
    "IGHV1-1.IGHJ1.IGHM_IGKV1-1.IGKJ1.IGKC",   # paired heavy + kappa
    "IGHV2-1.IGHJ2.IGHG1_NA",                  # heavy only
    "NA_IGLV1-1.IGLJ1.IGLC"                    # lambda only
  )
  ct_aa <- c("CARH1_CARK1", "CARH2_NA", "NA_CARL3")
  ct_nt <- c("TGTGCH1_TGTGCK1", "TGTGCH2_NA", "NA_TGTGCL3")

  SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(0L, 2, 3,
                                  dimnames = list(c("g1", "g2"),
                                                  paste0("b", 1:3)))),
    colData = S4Vectors::DataFrame(
      CTgene = ct_gene, CTaa = ct_aa, CTnt = ct_nt,
      row.names = paste0("b", 1:3)
    )
  )
}
