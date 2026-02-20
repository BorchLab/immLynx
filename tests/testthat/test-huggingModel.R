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


# ===========================================================================
# Python-dependent tests (require transformers)
# ===========================================================================

test_that("huggingModel returns list with model and tokenizer", {
  skip_if_no_transformers()

  result <- huggingModel("facebook/esm2_t12_35M_UR50D")

  expect_type(result, "list")
  expect_true("model" %in% names(result))
  expect_true("tokenizer" %in% names(result))
  expect_true(!is.null(result$model))
  expect_true(!is.null(result$tokenizer))
})

test_that("huggingModel produces informative messages", {
  skip_if_no_transformers()

  expect_message(
    huggingModel("facebook/esm2_t12_35M_UR50D"),
    "Initializing Hugging Face"
  )
})

test_that("huggingModel produces message about loading tokenizer", {
  skip_if_no_transformers()

  expect_message(
    huggingModel("facebook/esm2_t12_35M_UR50D"),
    "Loading tokenizer"
  )
})

test_that("huggingModel produces message about loading model", {
  skip_if_no_transformers()

  expect_message(
    huggingModel("facebook/esm2_t12_35M_UR50D"),
    "Loading model"
  )
})
