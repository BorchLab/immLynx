#' Initialize a Hugging Face Model and Tokenizer
#'
#' @description Downloads or loads a cached pre-trained model and its
#'   corresponding tokenizer from the Hugging Face Hub. Uses the basilisk-managed
#'   Python environment which includes \code{transformers} and \code{torch}.
#'
#' @param model_name A string specifying the model identifier from the
#'   Hugging Face Hub. For ESM-2 35M, this is "facebook/esm2_t12_35M_UR50D".
#'   Other options include "facebook/esm2_t33_650M_UR50D" for larger models.
#'
#' @return A list containing the R-wrapped Python objects for the 'model'
#'   and 'tokenizer', along with the basilisk process handle. The process
#'   must be stopped when no longer needed via \code{basilisk::basiliskStop()}.
#'
#' @export
#' @importFrom reticulate import
#'
#' @seealso \code{\link{tokenizeSequences}}, \code{\link{proteinEmbeddings}},
#'   \code{\link{runEmbeddings}}
#'
#' @examples
#' # Default model is ESM-2 35M
#' model_name <- "facebook/esm2_t12_35M_UR50D"
#' \donttest{
#'   # Load the default ESM-2 35M model
#'   hf_components <- huggingModel(model_name)
#'   names(hf_components)  # "model", "tokenizer", "proc"
#'
#'   # Use with tokenizeSequences and proteinEmbeddings
#'   sequences <- c("CASSLGTGELFF", "CASSIRSSYEQYF")
#'   tokenized <- tokenizeSequences(hf_components$tokenizer,
#'                                  sequences)
#'   embeddings <- proteinEmbeddings(hf_components$model,
#'                                   tokenized,
#'                                   pool = "mean",
#'                                   chunk_size = 32)
#'
#'   # Clean up the basilisk process when done
#'   basilisk::basiliskStop(hf_components$proc)
#' }
huggingModel <- function(model_name = "facebook/esm2_t12_35M_UR50D") {

  # Inform the user what's happening
  message("Initializing Hugging Face components for model: ", model_name)
  message("This may take a while if the model needs to be downloaded...")

  proc <- basilisk::basiliskStart(immLynxEnv)
  success <- FALSE
  on.exit(if (!success) basilisk::basiliskStop(proc))

  result <- basilisk::basiliskRun(proc, function(model_name) {
    transformers <- reticulate::import("transformers", convert = FALSE)

    message("Loading tokenizer...")
    tokenizer <- transformers$AutoTokenizer$from_pretrained(model_name)

    message("Loading model...")
    model <- transformers$AutoModel$from_pretrained(model_name)

    message("Initialization successful.")

    list(
      model = model,
      tokenizer = tokenizer
    )
  }, model_name = model_name)

  result$proc <- proc
  success <- TRUE

  return(result)
}
