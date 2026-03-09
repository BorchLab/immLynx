#' Tokenize Amino Acid Sequences
#'
#' @description Takes a vector of amino acid sequences and uses a Hugging Face
#'   tokenizer to convert them into numerical input IDs suitable for model input.
#'   The tokenizer should be obtained from \code{\link{huggingModel}}, which
#'   manages the Python environment via basilisk.
#'
#' @param tokenizer The tokenizer object returned by \code{\link{huggingModel}}.
#' @param aa_sequences A character vector of amino acid sequences (e.g., CDR3 sequences).
#' @param padding A logical or string. If TRUE, pads sequences to the length
#'   of the longest sequence in the batch. Defaults to TRUE.
#' @param truncation A logical. If TRUE, truncates sequences to the model's
#'   maximum input length. Defaults to TRUE.
#' @param return_tensors A string specifying the format for the returned tensors.
#'   "pt" for PyTorch tensors, "tf" for TensorFlow. Defaults to "pt".
#'
#' @return The tokenized output, typically a dictionary-like object containing
#'   'input_ids' and 'attention_mask'.
#'
#' @export
#'
#' @seealso \code{\link{huggingModel}}, \code{\link{proteinEmbeddings}},
#'   \code{\link{runEmbeddings}}
#'
#' @examples
#' sequences <- c("CASSLGTGELFF", "CASSIRSSYEQYF", "CASSYSTGELFF")
#' \dontrun{
#'   # Initialize model and tokenizer
#'   hf_components <- huggingModel()
#'
#'   # Tokenize CDR3 sequences
#'   tokenized <- tokenizeSequences(hf_components$tokenizer,
#'                                  sequences)
#'
#'   # Tokenize without padding (variable-length output)
#'   tokenized_nopad <- tokenizeSequences(
#'       hf_components$tokenizer,
#'       sequences,
#'       padding = FALSE)
#'
#'   # Pass tokenized output to proteinEmbeddings
#'   embeddings <- proteinEmbeddings(hf_components$model,
#'                                   tokenized,
#'                                   pool = "mean",
#'                                   chunk_size = 32)
#'
#'   # Clean up
#'   basilisk::basiliskStop(hf_components$proc)
#' }
tokenizeSequences <- function(tokenizer,
                              aa_sequences,
                              padding = TRUE,
                              truncation = TRUE,
                              return_tensors = "pt") {

  if (is.null(tokenizer)) {
    stop("Tokenizer object is NULL. Please initialize it first.")
  }

  message("Tokenizing ", length(aa_sequences), " sequences...")

  # The tokenizer is a Python function, so we call it directly.
  # It must be called within an active basilisk process (e.g., from
  # huggingModel() or runEmbeddings()).
  tokenized_output <- tokenizer(
    aa_sequences,
    padding = padding,
    truncation = truncation,
    return_tensors = return_tensors
  )

  message("Tokenization complete.")

  return(tokenized_output)
}
