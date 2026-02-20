# Tests for proteinEmbeddings function

# ===========================================================================
# Parameter validation
# ===========================================================================

test_that("proteinEmbeddings function exists and is exported", {
  expect_true(is.function(proteinEmbeddings))
})

test_that("proteinEmbeddings errors when model is NULL", {
  expect_error(
    proteinEmbeddings(NULL, list()),
    "!is.null\\(model\\)"
  )
})

test_that("proteinEmbeddings errors when tokenized.batch is NULL", {
  expect_error(
    proteinEmbeddings(list(), NULL),
    "!is.null\\(tokenized.batch\\)"
  )
})

test_that("proteinEmbeddings validates pool argument", {
  expect_no_error(match.arg("mean", c("none", "mean", "cls")))
  expect_no_error(match.arg("cls", c("none", "mean", "cls")))
  expect_no_error(match.arg("none", c("none", "mean", "cls")))
  expect_error(match.arg("invalid", c("none", "mean", "cls")))
})

test_that("proteinEmbeddings validates prefer_dtype argument", {
  expect_no_error(match.arg("float16", c("float16", "bfloat16", "float32")))
  expect_no_error(match.arg("bfloat16", c("float16", "bfloat16", "float32")))
  expect_no_error(match.arg("float32", c("float16", "bfloat16", "float32")))
  expect_error(match.arg("invalid", c("float16", "bfloat16", "float32")))
})

test_that("proteinEmbeddings validates prefer_device argument", {
  expect_no_error(match.arg("auto", c("auto", "cuda", "mps", "cpu")))
  expect_no_error(match.arg("cpu", c("auto", "cuda", "mps", "cpu")))
  expect_no_error(match.arg("cuda", c("auto", "cuda", "mps", "cpu")))
  expect_no_error(match.arg("mps", c("auto", "cuda", "mps", "cpu")))
  expect_error(match.arg("invalid", c("auto", "cuda", "mps", "cpu")))
})

test_that("proteinEmbeddings has correct default parameters", {
  f <- formals(proteinEmbeddings)

  expect_null(f$chunk_size)
})


# ===========================================================================
# Python-dependent tests
# ===========================================================================

test_that("proteinEmbeddings returns matrix with mean pooling", {
  skip_if_no_transformers()

  hf <- huggingModel("facebook/esm2_t12_35M_UR50D")
  sequences <- c("CASSLGTGELFF", "CASSIRSSYEQYF", "CASSYSTGELFF")
  tokenized <- tokenizeSequences(hf$tokenizer, sequences)

  result <- proteinEmbeddings(hf$model, tokenized,
                               pool = "mean", chunk_size = 32)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_true(ncol(result) > 0)
})

test_that("proteinEmbeddings returns matrix with cls pooling", {
  skip_if_no_transformers()

  hf <- huggingModel("facebook/esm2_t12_35M_UR50D")
  sequences <- c("CASSLGTGELFF", "CASSIRSSYEQYF")
  tokenized <- tokenizeSequences(hf$tokenizer, sequences)

  result <- proteinEmbeddings(hf$model, tokenized,
                               pool = "cls", chunk_size = 32)

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
})

test_that("proteinEmbeddings returns list with no pooling", {
  skip_if_no_transformers()

  hf <- huggingModel("facebook/esm2_t12_35M_UR50D")
  sequences <- c("CASSLGTGELFF", "CASSIRSSYEQYF")
  tokenized <- tokenizeSequences(hf$tokenizer, sequences)

  result <- proteinEmbeddings(hf$model, tokenized,
                               pool = "none", chunk_size = 32)

  expect_type(result, "list")
  expect_true(length(result) > 0)
})

test_that("proteinEmbeddings mean and cls produce same dimensions", {
  skip_if_no_transformers()

  hf <- huggingModel("facebook/esm2_t12_35M_UR50D")
  sequences <- c("CASSLGTGELFF", "CASSIRSSYEQYF")
  tokenized <- tokenizeSequences(hf$tokenizer, sequences)

  result_mean <- proteinEmbeddings(hf$model, tokenized,
                                    pool = "mean", chunk_size = 32)
  result_cls <- proteinEmbeddings(hf$model, tokenized,
                                   pool = "cls", chunk_size = 32)

  expect_equal(nrow(result_mean), nrow(result_cls))
  expect_equal(ncol(result_mean), ncol(result_cls))
})

test_that("proteinEmbeddings respects chunk_size", {
  skip_if_no_transformers()

  hf <- huggingModel("facebook/esm2_t12_35M_UR50D")
  sequences <- c("CASSLGTGELFF", "CASSIRSSYEQYF", "CASSYSTGELFF",
                  "CASSLNRDNEQFF", "CASSQDRTGQETQYF")
  tokenized <- tokenizeSequences(hf$tokenizer, sequences)

  # Small chunk size
  result_small <- proteinEmbeddings(hf$model, tokenized,
                                     pool = "mean", chunk_size = 2)
  # Large chunk size
  result_large <- proteinEmbeddings(hf$model, tokenized,
                                     pool = "mean", chunk_size = 10)

  # Results should be the same regardless of chunk size
  expect_equal(dim(result_small), dim(result_large))
  expect_equal(result_small, result_large, tolerance = 1e-4)
})
