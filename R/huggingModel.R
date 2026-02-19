#' Initialize a Hugging Face Model and Tokenizer
#'
#' @description Downloads or loads a cached pre-trained model and its
#'   corresponding tokenizer from the Hugging Face Hub.
#'
#' @param model_name A string specifying the model identifier from the
#'   Hugging Face Hub. For ESM-2 35M, this is "facebook/esm2_t12_35M_UR50D".
#'   Other options include "facebook/esm2_t33_650M_UR50D" for larger models.
#'
#' @return A list containing the R-wrapped Python objects for the 'model'
#'   and 'tokenizer'.
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
#' \dontrun{
#'   # Load the default ESM-2 35M model
#'   hf_components <- huggingModel(model_name)
#'   names(hf_components)  # "model" and "tokenizer"
#'
#'   # Load a larger ESM-2 model (650M parameters)
#'   hf_large <- huggingModel("facebook/esm2_t33_650M_UR50D")
#'
#'   # Use with tokenizeSequences and proteinEmbeddings
#'   sequences <- c("CASSLGTGELFF", "CASSIRSSYEQYF")
#'   tokenized <- tokenizeSequences(hf_components$tokenizer,
#'                                  sequences)
#'   embeddings <- proteinEmbeddings(hf_components$model,
#'                                   tokenized,
#'                                   pool = "mean",
#'                                   chunk_size = 32)
#' }
huggingModel <- function(model_name = "facebook/esm2_t12_35M_UR50D") {
  
  # Inform the user what's happening
  message("Initializing Hugging Face components for model:", model_name, "\n")
  message("This may take a while if the model needs to be downloaded...\n")
  
  tryCatch({
    # Import the 'transformers' Python library into the R session
    transformers <- import("transformers", convert = FALSE) 
    
    # Load the pre-trained tokenizer
    message("Loading tokenizer...\n")
    tokenizer <- transformers$AutoTokenizer$from_pretrained(model_name)
    
    # Load the pre-trained model
    message("Loading model...\n")
    model <- transformers$AutoModel$from_pretrained(model_name)
    
    message("Initialization successful.\n")
    
    # Return the components as a list
    return(list(
      model = model,
      tokenizer = tokenizer
    ))
    
  }, error = function(e) {
    stop(
      "Could not initialize Hugging Face components. ",
      "Please ensure that Python, 'transformers', and 'torch' ",
      "are installed and accessible by 'reticulate'. ",
      "Original error: ", e$message
    )
  })
}
