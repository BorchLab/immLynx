# Tests for huggingModel function

# ===========================================================================
# Parameter validation and function structure
# ===========================================================================

test_that("huggingModel function exists and is exported", {
  expect_true(is.function(huggingModel))
})

test_that("huggingModel has correct default parameter", {
  f <- formals(huggingModel)

  expect_equal(f$model_name, "facebook/esm2_t12_35M_UR50D")
})

test_that("huggingModel has single parameter", {
  f <- formals(huggingModel)
  expect_equal(length(f), 1)
  expect_true("model_name" %in% names(f))
})

# ===========================================================================
# Python-dependent tests (require basilisk env with transformers + torch)
# ===========================================================================

test_that("huggingModel returns list with model and tokenizer", {
  skip_if_no_python()

  result <- huggingModel("facebook/esm2_t12_35M_UR50D")
  on.exit(basilisk::basiliskStop(result$proc))

  expect_type(result, "list")
  expect_true("model" %in% names(result))
  expect_true("tokenizer" %in% names(result))
  expect_true("proc" %in% names(result))
  expect_true(!is.null(result$model))
  expect_true(!is.null(result$tokenizer))
})

test_that("huggingModel returns a basilisk process handle", {
  skip_if_no_python()

  result <- huggingModel("facebook/esm2_t12_35M_UR50D")
  on.exit(basilisk::basiliskStop(result$proc))

  # proc should be a basilisk process handle
  expect_true(!is.null(result$proc))
})

test_that("huggingModel produces initialization message", {
  skip_if_no_python()

  expect_message(
    {
      result <- huggingModel("facebook/esm2_t12_35M_UR50D")
      basilisk::basiliskStop(result$proc)
    },
    "Initializing Hugging Face"
  )
})

test_that("huggingModel produces tokenizer loading message", {
  skip_if_no_python()

  expect_message(
    {
      result <- huggingModel("facebook/esm2_t12_35M_UR50D")
      basilisk::basiliskStop(result$proc)
    },
    "Loading tokenizer"
  )
})

test_that("huggingModel produces model loading message", {
  skip_if_no_python()

  expect_message(
    {
      result <- huggingModel("facebook/esm2_t12_35M_UR50D")
      basilisk::basiliskStop(result$proc)
    },
    "Loading model"
  )
})

test_that("huggingModel produces success message", {
  skip_if_no_python()

  expect_message(
    {
      result <- huggingModel("facebook/esm2_t12_35M_UR50D")
      basilisk::basiliskStop(result$proc)
    },
    "Initialization successful"
  )
})

test_that("huggingModel errors with invalid model name", {
  skip_if_no_python()

  expect_error(
    huggingModel("completely/nonexistent-model-xyz-12345"),
    "Could not initialize Hugging Face"
  )
})

test_that("huggingModel tokenizer can tokenize sequences", {
  skip_if_no_python()

  result <- huggingModel("facebook/esm2_t12_35M_UR50D")
  on.exit(basilisk::basiliskStop(result$proc))

  sequences <- c("CASSLGTGELFF", "CASSIRSSYEQYF")
  tokenized <- result$tokenizer(sequences, padding = TRUE,
                                 truncation = TRUE, return_tensors = "pt")

  expect_true(reticulate::py_has_attr(tokenized, "input_ids"))
  expect_true(reticulate::py_has_attr(tokenized, "attention_mask"))
})

test_that("huggingModel model has eval method", {
  skip_if_no_python()

  result <- huggingModel("facebook/esm2_t12_35M_UR50D")
  on.exit(basilisk::basiliskStop(result$proc))

  expect_true(reticulate::py_has_attr(result$model, "eval"))
  expect_true(reticulate::py_has_attr(result$model, "to"))
})
