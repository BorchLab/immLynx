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

# ---------------------------------------------------------------------------
# Module availability checks (cached per test session)
# ---------------------------------------------------------------------------

# Cache to avoid redundant basilisk process starts during a test session
.module_cache <- new.env(parent = emptyenv())

# Check if Python environment is available
python_available <- function() {
  if (is.null(.module_cache$python)) {
    .module_cache$python <- tryCatch({
      proc <- basilisk::basiliskStart(immLynx:::immLynxEnv)
      on.exit(basilisk::basiliskStop(proc))
      TRUE
    }, error = function(e) FALSE)
  }
  .module_cache$python
}

# Check if transformers module is available in the basilisk environment
transformers_available <- function() {
  if (is.null(.module_cache$transformers)) {
    .module_cache$transformers <- tryCatch({
      proc <- basilisk::basiliskStart(immLynx:::immLynxEnv)
      on.exit(basilisk::basiliskStop(proc))
      basilisk::basiliskRun(proc, function() {
        reticulate::import("transformers")
        TRUE
      })
    }, error = function(e) FALSE)
  }
  .module_cache$transformers
}

# Check if tcrdist module (with submodules) is importable
tcrdist_available <- function() {
  if (is.null(.module_cache$tcrdist)) {
    .module_cache$tcrdist <- tryCatch({
      proc <- basilisk::basiliskStart(immLynx:::immLynxEnv)
      on.exit(basilisk::basiliskStop(proc))
      basilisk::basiliskRun(proc, function() {
        reticulate::import("tcrdist.repertoire")
        TRUE
      })
    }, error = function(e) FALSE)
  }
  .module_cache$tcrdist
}

# Check if metaclonotypist module is importable
metaclonotypist_available <- function() {
  if (is.null(.module_cache$metaclonotypist)) {
    .module_cache$metaclonotypist <- tryCatch({
      proc <- basilisk::basiliskStart(immLynx:::immLynxEnv)
      on.exit(basilisk::basiliskStop(proc))
      basilisk::basiliskRun(proc, function() {
        reticulate::import("metaclonotypist")
        TRUE
      })
    }, error = function(e) FALSE)
  }
  .module_cache$metaclonotypist
}

# Check if olga module is importable
olga_available <- function() {
  if (is.null(.module_cache$olga)) {
    .module_cache$olga <- tryCatch({
      proc <- basilisk::basiliskStart(immLynx:::immLynxEnv)
      on.exit(basilisk::basiliskStop(proc))
      basilisk::basiliskRun(proc, function() {
        reticulate::import("olga")
        TRUE
      })
    }, error = function(e) FALSE)
  }
  .module_cache$olga
}

# Check if clustcr module is importable
clustcr_available <- function() {
  if (is.null(.module_cache$clustcr)) {
    .module_cache$clustcr <- tryCatch({
      proc <- basilisk::basiliskStart(immLynx:::immLynxEnv)
      on.exit(basilisk::basiliskStop(proc))
      basilisk::basiliskRun(proc, function() {
        reticulate::import("clustcr")
        TRUE
      })
    }, error = function(e) FALSE)
  }
  .module_cache$clustcr
}

# Check if sonnia module is importable
sonnia_available <- function() {
  if (is.null(.module_cache$sonnia)) {
    .module_cache$sonnia <- tryCatch({
      proc <- basilisk::basiliskStart(immLynx:::immLynxEnv)
      on.exit(basilisk::basiliskStop(proc))
      basilisk::basiliskRun(proc, function() {
        reticulate::import("sonnia.sonnia")
        TRUE
      })
    }, error = function(e) FALSE)
  }
  .module_cache$sonnia
}

# ---------------------------------------------------------------------------
# Skip helpers
# ---------------------------------------------------------------------------

# Skip test if Python not available
skip_if_no_python <- function() {
  if (!python_available()) {
    testthat::skip("Python environment not available")
  }
}

# Skip test if transformers/torch not available in basilisk env
skip_if_no_transformers <- function() {
  skip_if_no_python()
  if (!transformers_available()) {
    testthat::skip("transformers module not available in Python environment")
  }
}

# Skip test if tcrdist is not importable (e.g., shared lib issues on macOS)
skip_if_no_tcrdist <- function() {
  skip_if_no_python()
  if (!tcrdist_available()) {
    testthat::skip("tcrdist module not available (possible shared library issue)")
  }
}

# Skip test if metaclonotypist is not importable (e.g., shared lib issues on macOS)
skip_if_no_metaclonotypist <- function() {
  skip_if_no_python()
  if (!metaclonotypist_available()) {
    testthat::skip("metaclonotypist module not available (possible shared library issue)")
  }
}

# Skip test if olga is not importable
skip_if_no_olga <- function() {
  skip_if_no_python()
  if (!olga_available()) {
    testthat::skip("olga module not available in Python environment")
  }
}

# Skip test if clustcr is not importable
skip_if_no_clustcr <- function() {
  skip_if_no_python()
  if (!clustcr_available()) {
    testthat::skip("clustcr module not available in Python environment")
  }
}

# Skip test if sonnia is not importable
skip_if_no_sonnia <- function() {
  skip_if_no_python()
  if (!sonnia_available()) {
    testthat::skip("sonnia module not available in Python environment")
  }
}
