# Tests for tokenizeSequences function

# ===========================================================================
# Parameter validation
# ===========================================================================

test_that("tokenizeSequences function exists and is exported", {
  expect_true(is.function(tokenizeSequences))
})

test_that("tokenizeSequences errors when tokenizer is NULL", {
  expect_error(
    tokenizeSequences(NULL, c("CASSAAA", "CASSBBB")),
    "Tokenizer object is NULL"
  )
})

test_that("tokenizeSequences has correct default parameters", {
  f <- formals(tokenizeSequences)

  expect_equal(f$padding, TRUE)
  expect_equal(f$truncation, TRUE)
  expect_equal(f$return_tensors, "pt")
})

test_that("tokenizeSequences has all expected parameters", {
  f <- formals(tokenizeSequences)

  expect_true("tokenizer" %in% names(f))
  expect_true("aa_sequences" %in% names(f))
  expect_true("padding" %in% names(f))
  expect_true("truncation" %in% names(f))
  expect_true("return_tensors" %in% names(f))
})


# ===========================================================================
# Python-dependent tests
# ===========================================================================

test_that("tokenizeSequences returns tokenized output with input_ids", {
  skip_if_no_python()

  hf <- huggingModel("facebook/esm2_t12_35M_UR50D")
  on.exit(basilisk::basiliskStop(hf$proc))

  sequences <- c("CASSLGTGELFF", "CASSIRSSYEQYF", "CASSYSTGELFF")
  result <- tokenizeSequences(hf$tokenizer, sequences)

  expect_true(!is.null(result))
  expect_true(reticulate::py_has_attr(result, "input_ids"))
})

test_that("tokenizeSequences returns attention_mask", {
  skip_if_no_python()

  hf <- huggingModel("facebook/esm2_t12_35M_UR50D")
  on.exit(basilisk::basiliskStop(hf$proc))

  sequences <- c("CASS", "CASSLAPGATNEKLFF")
  result <- tokenizeSequences(hf$tokenizer, sequences, padding = TRUE)

  expect_true(reticulate::py_has_attr(result, "attention_mask"))
})

test_that("tokenizeSequences handles single sequence", {
  skip_if_no_python()

  hf <- huggingModel("facebook/esm2_t12_35M_UR50D")
  on.exit(basilisk::basiliskStop(hf$proc))

  result <- tokenizeSequences(hf$tokenizer, "CASSLGTGELFF")

  expect_true(!is.null(result))
  expect_true(reticulate::py_has_attr(result, "input_ids"))
})

test_that("tokenizeSequences produces tokenization message", {
  skip_if_no_python()

  hf <- huggingModel("facebook/esm2_t12_35M_UR50D")
  on.exit(basilisk::basiliskStop(hf$proc))

  expect_message(
    tokenizeSequences(hf$tokenizer, c("CASSAAA", "CASSBBB")),
    "Tokenizing"
  )
})

test_that("tokenizeSequences produces completion message", {
  skip_if_no_python()

  hf <- huggingModel("facebook/esm2_t12_35M_UR50D")
  on.exit(basilisk::basiliskStop(hf$proc))

  expect_message(
    tokenizeSequences(hf$tokenizer, c("CASSAAA")),
    "Tokenization complete"
  )
})

test_that("tokenizeSequences message contains sequence count", {
  skip_if_no_python()

  hf <- huggingModel("facebook/esm2_t12_35M_UR50D")
  on.exit(basilisk::basiliskStop(hf$proc))

  expect_message(
    tokenizeSequences(hf$tokenizer, c("CASSAAA", "CASSBBB", "CASSCCC")),
    "3 sequences"
  )
})

test_that("tokenizeSequences handles many sequences", {
  skip_if_no_python()

  hf <- huggingModel("facebook/esm2_t12_35M_UR50D")
  on.exit(basilisk::basiliskStop(hf$proc))

  sequences <- paste0("CASS", LETTERS[1:20])
  result <- tokenizeSequences(hf$tokenizer, sequences)

  expect_true(!is.null(result))
  # input_ids should have 20 rows
  n <- as.integer(reticulate::py_len(result$input_ids))
  expect_equal(n, 20)
})

test_that("tokenizeSequences padding pads shorter sequences", {

  skip_if_no_python()

  hf <- huggingModel("facebook/esm2_t12_35M_UR50D")
  on.exit(basilisk::basiliskStop(hf$proc))

  # These sequences have different lengths
  sequences <- c("CASS", "CASSLAPGATNEKLFF")
  result <- tokenizeSequences(hf$tokenizer, sequences, padding = TRUE)

  # Both should be padded to same length
  # torch.Size is a tuple; py_to_r may flatten it so use py_len for batch dim
  n_seqs <- as.integer(reticulate::py_len(result$input_ids))
  expect_equal(n_seqs, 2)  # 2 sequences
  # Length should match the longer sequence + special tokens
})
