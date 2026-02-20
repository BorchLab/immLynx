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


# ===========================================================================
# Python-dependent tests
# ===========================================================================

test_that("tokenizeSequences returns tokenized output", {
  skip_if_no_transformers()

  hf <- huggingModel("facebook/esm2_t12_35M_UR50D")
  sequences <- c("CASSLGTGELFF", "CASSIRSSYEQYF", "CASSYSTGELFF")

  result <- tokenizeSequences(hf$tokenizer, sequences)

  expect_true(!is.null(result))
  # Should have input_ids attribute
  expect_true(reticulate::py_has_attr(result, "input_ids"))
})

test_that("tokenizeSequences handles single sequence", {
  skip_if_no_transformers()

  hf <- huggingModel("facebook/esm2_t12_35M_UR50D")

  result <- tokenizeSequences(hf$tokenizer, "CASSLGTGELFF")

  expect_true(!is.null(result))
  expect_true(reticulate::py_has_attr(result, "input_ids"))
})

test_that("tokenizeSequences produces tokenization message", {
  skip_if_no_transformers()

  hf <- huggingModel("facebook/esm2_t12_35M_UR50D")

  expect_message(
    tokenizeSequences(hf$tokenizer, c("CASSAAA", "CASSBBB")),
    "Tokenizing"
  )
})

test_that("tokenizeSequences has attention_mask", {
  skip_if_no_transformers()

  hf <- huggingModel("facebook/esm2_t12_35M_UR50D")
  sequences <- c("CASS", "CASSLAPGATNEKLFF")

  result <- tokenizeSequences(hf$tokenizer, sequences, padding = TRUE)

  expect_true(reticulate::py_has_attr(result, "attention_mask"))
})
