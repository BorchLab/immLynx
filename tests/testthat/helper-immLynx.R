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

# Skip test if Python not available
skip_if_no_python <- function() {
  if (!python_available()) {
    testthat::skip("Python environment not available")
  }
}
