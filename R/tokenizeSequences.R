#' @title Tokenize Amino Acid Sequences
#' @description Takes a vector of amino acid sequences and uses a Hugging Face
#'   tokenizer to convert them into numerical input IDs.
#' @param tokenizer The tokenizer object returned by `initialize_protein_model_and_tokenizer`.
#' @param aa_sequences A character vector of amino acid sequences.
#' @param padding A logical or string. If TRUE, pads sequences to the length
#'   of the longest sequence in the batch. Defaults to TRUE.
#' @param truncation A logical. If TRUE, truncates sequences to the model's
#'   maximum input length. Defaults to TRUE.
#' @param return_tensors A string specifying the format for the returned tensors.
#'   "pt" for PyTorch tensors, "tf" for TensorFlow. Defaults to "pt".
#' @return The tokenized output, typically a dictionary-like object containing
#'   'input_ids' and 'attention_mask'.
#' @examples
#' \dontrun{
#'   hf_components <- initialize_protein_model_and_tokenizer()
#'   sequences <- c("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALP",
#'                  "MKKLFWRAVFLFLLAGLVACSP")
#'   tokenized_output <- tokenizeSequences(hf_components$tokenizer, sequences)
#'   print(tokenized_output)
#' }
tokenizeSequences <- function(tokenizer, 
                              aa_sequences, 
                              padding = TRUE, 
                              truncation = TRUE, 
                              return_tensors = "pt") {
  
  if (is.null(tokenizer)) {
    stop("Tokenizer object is NULL. Please initialize it first.")
  }
  
  message("Tokenizing", length(aa_sequences), "sequences...\n")
  
  # The tokenizer is a Python function, so we call it directly.
  tokenized_output <- tokenizer(
    aa_sequences,
    padding = padding,
    truncation = truncation,
    return_tensors = return_tensors
  )
  
  message("Tokenization complete.\n")
  
  return(tokenized_output)
}