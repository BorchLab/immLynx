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
#' @seealso \code{\link{tokenizeSequences}}, \code{\link{proteinEmbeddings}}
#'
#' @examples
#' \dontrun{
#'   hf_components <- huggingModel()
#'   model <- hf_components$model
#'   tokenizer <- hf_components$tokenizer
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
    stop(paste(
      "Failed to initialize Hugging Face components.",
      "Please ensure that Python, 'transformers', and 'torch' are installed and accessible by 'reticulate'.",
      "Original error:", e$message
    ))
  })
}
