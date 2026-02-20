# Tests for internal calculate helper functions (calculate_helpers.R)

# ===========================================================================
# calculate.olga
# ===========================================================================

test_that("calculate.olga validates action argument", {
  expect_error(
    immLynx:::calculate.olga(action = "invalid"),
    "'arg' should be one of"
  )
})

test_that("calculate.olga accepts valid action values", {
  expect_no_error(match.arg("pgen", c("pgen", "generate")))
  expect_no_error(match.arg("generate", c("pgen", "generate")))
})

test_that("calculate.olga function has correct default parameters", {
  f <- formals(immLynx:::calculate.olga)

  expect_equal(f$model, "humanTRB")
  expect_equal(f$n, 100)
  expect_null(f$sequences)
  expect_null(f$v_genes)
  expect_null(f$j_genes)
})

test_that("calculate.olga pgen action generates probabilities", {
  skip_if_no_python()

  result <- immLynx:::calculate.olga(
    action = "pgen",
    model = "humanTRB",
    sequences = c("CASSLAPGATNEKLFF", "CASSLGQAYEQYF")
  )

  expect_type(result, "double")
  expect_length(result, 2)
  expect_true(all(result >= 0))
  expect_true(all(result <= 1))
})

test_that("calculate.olga generate action returns data.frame", {
  skip_if_no_python()

  result <- immLynx:::calculate.olga(
    action = "generate",
    model = "humanTRB",
    n = 5
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 5)
  expect_true("nt_seq" %in% names(result))
  expect_true("aa_seq" %in% names(result))
  expect_true("v_index" %in% names(result))
  expect_true("j_index" %in% names(result))
})

test_that("calculate.olga pgen with V/J genes", {
  skip_if_no_python()

  result <- immLynx:::calculate.olga(
    action = "pgen",
    model = "humanTRB",
    sequences = c("CASSLAPGATNEKLFF"),
    v_genes = c("TRBV7-2*01"),
    j_genes = c("TRBJ1-1*01")
  )

  expect_type(result, "double")
  expect_length(result, 1)
})

test_that("calculate.olga works with humanTRA model", {
  skip_if_no_python()

  result <- immLynx:::calculate.olga(
    action = "generate",
    model = "humanTRA",
    n = 3
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 3)
})


# ===========================================================================
# calculate.clustcr
# ===========================================================================

test_that("calculate.clustcr function has correct default parameters", {
  f <- formals(immLynx:::calculate.clustcr)

  expect_equal(f$method, "mcl")
  expect_equal(f$inflation, 2.0)
  expect_equal(f$eps, 0.5)
  expect_equal(f$min_samples, 2)
})

test_that("calculate.clustcr returns cluster assignments", {
  skip_if_no_python()

  sequences <- c(
    "CASSLAPGATNEKLFF", "CASSLGQAYEQYF", "CASRLAGQETQYF",
    "CASSLAPGATNEKLFF", "CASSYSGGNTGELFF", "CASSQDRTGQETQYF",
    "CASSLGQAYEQYF", "CASSLNRDNEQFF", "CASSLTGTEAFF",
    "CASSLAPGATNEKLFF"
  )

  result <- immLynx:::calculate.clustcr(sequences = sequences)

  expect_true(is.numeric(result))
  expect_length(result, length(sequences))
  # Some may be NA (singletons)
})


# ===========================================================================
# calculate.tcrDist
# ===========================================================================

test_that("calculate.tcrDist function has correct default parameters", {
  f <- formals(immLynx:::calculate.tcrDist)

  expect_equal(f$organism, "human")
  expect_equal(f$chains, "beta")
  expect_equal(f$compute_distances, TRUE)
})

test_that("calculate.tcrDist returns distance results for beta chain", {
  skip_if_no_python()

  df <- data.frame(
    count = rep(1L, 5),
    v_b_gene = paste0("TRBV", c("7-2", "6-1", "12-3", "5-1", "20-1"), "*01"),
    j_b_gene = paste0("TRBJ", c("1-1", "2-1", "2-5", "1-2", "2-3"), "*01"),
    cdr3_b_aa = c("CASSLAPGATNEKLFF", "CASSLGQAYEQYF", "CASRLAGQETQYF",
                  "CASSYSGGNTGELFF", "CASSQDRTGQETQYF"),
    stringsAsFactors = FALSE
  )

  result <- immLynx:::calculate.tcrDist(
    df = df,
    organism = "human",
    chains = "beta"
  )

  expect_type(result, "list")
  expect_true("pw_beta" %in% names(result) || "pw_cdr3_b_aa" %in% names(result))
})


# ===========================================================================
# calculate.sonia
# ===========================================================================

test_that("calculate.sonia function has correct default parameters", {
  f <- formals(immLynx:::calculate.sonia)

  expect_equal(f$organism, "human")
  expect_equal(f$dataset_type, "TCR")
  expect_equal(f$n_epochs, 100)
  expect_equal(f$save_folder, "sonia_output")
})


# ===========================================================================
# calculate.metaclonotypist
# ===========================================================================

test_that("calculate.metaclonotypist function has correct default parameters", {
  f <- formals(immLynx:::calculate.metaclonotypist)

  expect_equal(f$chain, "beta")
  expect_equal(f$method, "tcrdist")
  expect_equal(f$max_edits, 2)
  expect_equal(f$max_dist, 20)
  expect_equal(f$clustering, "leiden")
  expect_equal(f$resolution, 1.0)
})
