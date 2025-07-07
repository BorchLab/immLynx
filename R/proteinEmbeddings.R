#' Get Protein Embeddings from a Model
#' @description Applies a pre-trained model to a batch of tokenized sequences
#'   to generate embeddings.
#' @param model The model object returned by `huggingModel`.
#' @param tokenized.batch The output from the `tokenizeSequences` function.
#' @param convert.to.matrix Return a matrix (TRUE) or a vector (FALSE)
#' @return The model's output, typically containing `last_hidden_state`.
#' @importFrom reticulate py_to_r array_reshape
#' @examples
#' \dontrun{
#'   hf_components <- initialize_protein_model_and_tokenizer()
#'   sequences <- c("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALP",
#'                  "MKKLFWRAVFLFLLAGLVACSP")
#'   tokenized_output <- tokenize_sequences(hf_components$tokenizer, sequences)
#'   model_output <- get_protein_embeddings(hf_components$model, tokenized_output)
#'   # Extract the per-residue embeddings
#'   embeddings <- model_output$last_hidden_state
#' }
proteinEmbeddings <- function(model, 
                              tokenized.batch, 
                              convert.to.matrix = TRUE) {
  if (is.null(model)) {
    stop("Model object is NULL. Please initialize it first.")
  }
  if (is.null(tokenized.batch)) {
    stop("Tokenized batch is NULL. Please tokenize sequences first.")
  }
  
  message("Applying model to generate embeddings...\n")
  
  # Import torch and tell it not to calculate gradients to save memory
  torch <- import("torch", convert = FALSE)
  with(torch$no_grad(), {
    model_output <- model(
      input_ids = tokenized.batch$input_ids,
      attention_mask = tokenized.batch$attention_mask
    )
  })
  
  last_hidden_states_array <- py_to_r(model_output$last_hidden_state$cpu()$numpy())
  if(length(dim(last_hidden_states_array)) == 3 & convert.to.matrix) { #Reshape 3D array into matrix
    last_hidden_states_array <- array_reshape(last_hidden_states_array, 
                                              c(dim(last_hidden_states_array)[1], prod(dim(last_hidden_states_array)[2:3])))
  } 
  
  return(last_hidden_states_array)
}