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

test_that("proteinEmbeddings has all expected parameters", {
  f <- formals(proteinEmbeddings)

  expect_true("model" %in% names(f))
  expect_true("tokenized.batch" %in% names(f))
  expect_true("pool" %in% names(f))
  expect_true("chunk_size" %in% names(f))
  expect_true("prefer_dtype" %in% names(f))
  expect_true("prefer_device" %in% names(f))
})

test_that("proteinEmbeddings pool default includes all options", {
  f <- formals(proteinEmbeddings)

  pool_options <- eval(f$pool)
  expect_true("none" %in% pool_options)
  expect_true("mean" %in% pool_options)
  expect_true("cls" %in% pool_options)
})


# ===========================================================================
# Python-dependent tests (require basilisk env with torch + transformers)
# ===========================================================================

test_that("proteinEmbeddings returns matrix with mean pooling", {
  skip_if_no_python()

  hf <- huggingModel("facebook/esm2_t12_35M_UR50D")
  on.exit(basilisk::basiliskStop(hf$proc))

  sequences <- c("CASSLGTGELFF", "CASSIRSSYEQYF", "CASSYSTGELFF")
  tokenized <- tokenizeSequences(hf$tokenizer, sequences)

  result <- proteinEmbeddings(hf$model, tokenized,
                               pool = "mean", chunk_size = 32,
                               prefer_device = "cpu")

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_true(ncol(result) > 0)
})

test_that("proteinEmbeddings returns matrix with cls pooling", {
  skip_if_no_python()

  hf <- huggingModel("facebook/esm2_t12_35M_UR50D")
  on.exit(basilisk::basiliskStop(hf$proc))

  sequences <- c("CASSLGTGELFF", "CASSIRSSYEQYF")
  tokenized <- tokenizeSequences(hf$tokenizer, sequences)

  result <- proteinEmbeddings(hf$model, tokenized,
                               pool = "cls", chunk_size = 32,
                               prefer_device = "cpu")

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 2)
})

test_that("proteinEmbeddings returns list with no pooling", {
  skip_if_no_python()

  hf <- huggingModel("facebook/esm2_t12_35M_UR50D")
  on.exit(basilisk::basiliskStop(hf$proc))

  sequences <- c("CASSLGTGELFF", "CASSIRSSYEQYF")
  tokenized <- tokenizeSequences(hf$tokenizer, sequences)

  result <- proteinEmbeddings(hf$model, tokenized,
                               pool = "none", chunk_size = 32,
                               prefer_device = "cpu")

  expect_type(result, "list")
  expect_true(length(result) > 0)
})

test_that("proteinEmbeddings mean and cls produce same dimensions", {
  skip_if_no_python()

  hf <- huggingModel("facebook/esm2_t12_35M_UR50D")
  on.exit(basilisk::basiliskStop(hf$proc))

  sequences <- c("CASSLGTGELFF", "CASSIRSSYEQYF")
  tokenized <- tokenizeSequences(hf$tokenizer, sequences)

  result_mean <- proteinEmbeddings(hf$model, tokenized,
                                    pool = "mean", chunk_size = 32,
                                    prefer_device = "cpu")
  result_cls <- proteinEmbeddings(hf$model, tokenized,
                                   pool = "cls", chunk_size = 32,
                                   prefer_device = "cpu")

  expect_equal(nrow(result_mean), nrow(result_cls))
  expect_equal(ncol(result_mean), ncol(result_cls))
})

test_that("proteinEmbeddings respects chunk_size", {
  skip_if_no_python()

  hf <- huggingModel("facebook/esm2_t12_35M_UR50D")
  on.exit(basilisk::basiliskStop(hf$proc))

  sequences <- c("CASSLGTGELFF", "CASSIRSSYEQYF", "CASSYSTGELFF",
                  "CASSLNRDNEQFF", "CASSQDRTGQETQYF")
  tokenized <- tokenizeSequences(hf$tokenizer, sequences)

  # Small chunk size (force CPU for reproducibility)
  result_small <- proteinEmbeddings(hf$model, tokenized,
                                     pool = "mean", chunk_size = 2,
                                     prefer_device = "cpu")
  # Large chunk size
  result_large <- proteinEmbeddings(hf$model, tokenized,
                                     pool = "mean", chunk_size = 10,
                                     prefer_device = "cpu")

  # Results should be the same regardless of chunk size on CPU
  expect_equal(dim(result_small), dim(result_large))
  expect_equal(result_small, result_large, tolerance = 1e-5)
})

test_that("proteinEmbeddings embedding values are finite", {
  skip_if_no_python()

  hf <- huggingModel("facebook/esm2_t12_35M_UR50D")
  on.exit(basilisk::basiliskStop(hf$proc))

  sequences <- c("CASSLGTGELFF", "CASSIRSSYEQYF")
  tokenized <- tokenizeSequences(hf$tokenizer, sequences)

  result <- proteinEmbeddings(hf$model, tokenized,
                               pool = "mean", chunk_size = 32,
                               prefer_device = "cpu")

  expect_true(all(is.finite(result)))
})

test_that("proteinEmbeddings hidden dim matches ESM-2 35M", {
  skip_if_no_python()

  hf <- huggingModel("facebook/esm2_t12_35M_UR50D")
  on.exit(basilisk::basiliskStop(hf$proc))

  sequences <- c("CASSLGTGELFF")
  tokenized <- tokenizeSequences(hf$tokenizer, sequences)

  result <- proteinEmbeddings(hf$model, tokenized,
                               pool = "mean", chunk_size = 32,
                               prefer_device = "cpu")

  # ESM-2 35M has hidden dimension of 480
  expect_equal(ncol(result), 480)
})

test_that("proteinEmbeddings single sequence produces 1-row matrix", {
  skip_if_no_python()

  hf <- huggingModel("facebook/esm2_t12_35M_UR50D")
  on.exit(basilisk::basiliskStop(hf$proc))

  tokenized <- tokenizeSequences(hf$tokenizer, "CASSLGTGELFF")

  result <- proteinEmbeddings(hf$model, tokenized,
                               pool = "mean", chunk_size = 32,
                               prefer_device = "cpu")

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
})

test_that("proteinEmbeddings different sequences produce different embeddings", {
  skip_if_no_python()

  hf <- huggingModel("facebook/esm2_t12_35M_UR50D")
  on.exit(basilisk::basiliskStop(hf$proc))

  # Very different sequences should give different embeddings
  sequences <- c("CASSAAA", "CASSLAPGATNEKLFF")
  tokenized <- tokenizeSequences(hf$tokenizer, sequences)

  result <- proteinEmbeddings(hf$model, tokenized,
                               pool = "mean", chunk_size = 32,
                               prefer_device = "cpu")

  # The two rows should not be identical
  expect_false(isTRUE(all.equal(result[1, ], result[2, ])))
})

test_that("proteinEmbeddings cpu device works", {
  skip_if_no_python()

  hf <- huggingModel("facebook/esm2_t12_35M_UR50D")
  on.exit(basilisk::basiliskStop(hf$proc))

  sequences <- c("CASSLGTGELFF")
  tokenized <- tokenizeSequences(hf$tokenizer, sequences)

  result <- proteinEmbeddings(hf$model, tokenized,
                               pool = "mean", chunk_size = 32,
                               prefer_device = "cpu")

  expect_true(is.matrix(result))
  expect_equal(nrow(result), 1)
})

test_that("proteinEmbeddings float32 dtype works", {
  skip_if_no_python()

  hf <- huggingModel("facebook/esm2_t12_35M_UR50D")
  on.exit(basilisk::basiliskStop(hf$proc))

  sequences <- c("CASSLGTGELFF")
  tokenized <- tokenizeSequences(hf$tokenizer, sequences)

  result <- proteinEmbeddings(hf$model, tokenized,
                               pool = "mean", chunk_size = 32,
                               prefer_dtype = "float32",
                               prefer_device = "cpu")

  expect_true(is.matrix(result))
  expect_true(all(is.finite(result)))
})
