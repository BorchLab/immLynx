#' Generate Protein Language Model Embeddings for TCR Sequences
#'
#' @description Extracts TCR CDR3 sequences from a Seurat object and generates embeddings
#'   using a protein language model (e.g., ESM-2).
#'
#' @param seurat_obj A Seurat object containing scRepertoire TCR data.
#' @param chains Which chain(s) to embed: "TRB", "TRA", or "both". Default is "TRB".
#' @param model_name Hugging Face model name. Default is "facebook/esm2_t12_35M_UR50D".
#'   Other options: "facebook/esm2_t33_650M_UR50D", "facebook/esm2_t36_3B_UR50D"
#' @param pool Pooling method: "mean", "cls", or "none". Default is "mean".
#' @param chunk_size Number of sequences to process at once. Default is 32.
#' @param reduction_name Name for the dimensional reduction (stored in Seurat object).
#'   Default is "tcr_esm".
#' @param reduction_key Key prefix for embeddings in reduction. Default is "ESM_".
#' @param return_seurat Logical. If TRUE, adds embeddings as dimensional reduction.
#'   If FALSE, returns list with embeddings and metadata. Default is TRUE.
#' @param ... Additional arguments passed to proteinEmbeddings().
#'
#' @return If return_seurat=TRUE, returns Seurat object with embeddings added as a
#'   dimensional reduction (accessible via seurat_obj[[reduction_name]]).
#'   If FALSE, returns list with embeddings matrix and metadata.
#'
#' @details This function uses protein language models to generate dense vector
#'   representations of TCR CDR3 sequences. These embeddings can be used for:
#'   - Dimensionality reduction and visualization
#'   - Clustering TCRs by sequence similarity
#'   - Downstream machine learning tasks
#'
#' @export
#' @importFrom immApex getIR
#'
#' @examples
#' \dontrun{
#'   # Generate embeddings for TRB sequences
#'   seurat_obj <- runEmbeddings(seurat_obj, chains = "TRB")
#'
#'   # Use larger model for better embeddings
#'   seurat_obj <- runEmbeddings(seurat_obj, 
#'                               chains = "TRB",
#'                               model_name = "facebook/esm2_t33_650M_UR50D")
#'
#'   # Visualize embeddings
#'   seurat_obj <- RunUMAP(seurat_obj, reduction = "tcr_esm", dims = 1:30)
#'   DimPlot(seurat_obj, reduction = "umap")
#' }
runEmbeddings <- function(seurat_obj,
                          chains = c("TRB", "TRA", "both"),
                          model_name = "facebook/esm2_t12_35M_UR50D",
                          pool = "mean",
                          chunk_size = 32,
                          reduction_name = "tcr_esm",
                          reduction_key = "ESM_",
                          return_seurat = TRUE,
                          ...) {
  
  chains <- match.arg(chains)
  
  message("Loading Hugging Face model: ", model_name)
  hf_components <- huggingModel(model_name = model_name)
  
  message("Extracting TCR sequences from Seurat object...")
  
  # Extract sequences based on chains parameter
  if (chains == "both") {
    tra_data <- immApex::getIR(seurat_obj, chains = "TRA")
    trb_data <- immApex::getIR(seurat_obj, chains = "TRB")
    
    # Combine and filter
    combined_data <- merge(tra_data, trb_data, 
                          by = "barcode", 
                          suffixes = c("_TRA", "_TRB"),
                          all = TRUE)
    
    # Create combined sequences (using available chains)
    combined_data$combined_seq <- ifelse(
      !is.na(combined_data$cdr3_aa_TRA) & !is.na(combined_data$cdr3_aa_TRB),
      paste0(combined_data$cdr3_aa_TRA, combined_data$cdr3_aa_TRB),
      ifelse(!is.na(combined_data$cdr3_aa_TRA),
             combined_data$cdr3_aa_TRA,
             combined_data$cdr3_aa_TRB)
    )
    
    valid_data <- combined_data[!is.na(combined_data$combined_seq), ]
    sequences <- valid_data$combined_seq
    barcodes <- valid_data$barcode
    chain_info <- ifelse(!is.na(valid_data$cdr3_aa_TRA) & !is.na(valid_data$cdr3_aa_TRB),
                        "TRA+TRB",
                        ifelse(!is.na(valid_data$cdr3_aa_TRA), "TRA", "TRB"))
    
  } else {
    tcr_data <- immApex::getIR(seurat_obj, chains = chains)
    tcr_data <- tcr_data[!is.na(tcr_data$cdr3_aa), ]
    
    if (nrow(tcr_data) == 0) {
      stop("No valid ", chains, " sequences found.")
    }
    
    sequences <- tcr_data$cdr3_aa
    barcodes <- tcr_data$barcode
    chain_info <- rep(chains, length(sequences))
  }
  
  message("Found ", length(sequences), " valid sequences")
  message("Tokenizing sequences...")
  
  # Tokenize sequences
  tokenized <- tokenizeSequences(
    tokenizer = hf_components$tokenizer,
    aa_sequences = sequences,
    padding = TRUE,
    truncation = TRUE,
    return_tensors = "pt"
  )
  
  message("Generating embeddings...")
  
  # Generate embeddings
  embeddings <- proteinEmbeddings(
    model = hf_components$model,
    tokenized.batch = tokenized,
    pool = pool,
    chunk_size = chunk_size,
    ...
  )
  
  # Add rownames (barcodes)
  rownames(embeddings) <- barcodes
  
  # Create metadata data frame
  metadata_df <- data.frame(
    barcode = barcodes,
    sequence = sequences,
    chain = chain_info,
    stringsAsFactors = FALSE
  )
  
  if (return_seurat) {
    # Check if Seurat is available
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Package 'Seurat' is required to add embeddings to Seurat object.")
    }
    
    # Create dimensional reduction object
    # Only include cells that have embeddings
    cell_embeddings <- matrix(NA, 
                              nrow = ncol(seurat_obj), 
                              ncol = ncol(embeddings),
                              dimnames = list(colnames(seurat_obj), 
                                            paste0(reduction_key, 1:ncol(embeddings))))
    
    # Fill in embeddings for cells with TCR data
    cell_embeddings[barcodes, ] <- embeddings
    
    # Create DimReduc object
    tcr_reduction <- Seurat::CreateDimReducObject(
      embeddings = cell_embeddings,
      key = reduction_key,
      assay = Seurat::DefaultAssay(seurat_obj)
    )
    
    # Add to Seurat object
    seurat_obj[[reduction_name]] <- tcr_reduction
    
    # Also add chain info to metadata
    chain_meta <- rep(NA_character_, ncol(seurat_obj))
    names(chain_meta) <- colnames(seurat_obj)
    chain_meta[barcodes] <- chain_info
    seurat_obj[[paste0(reduction_name, "_chain")]] <- chain_meta
    
    message("Embeddings added as '", reduction_name, "' reduction")
    message("Use RunUMAP(obj, reduction='", reduction_name, "') to visualize")
    
    return(seurat_obj)
    
  } else {
    output <- list(
      embeddings = embeddings,
      metadata = metadata_df,
      model_name = model_name,
      pool = pool
    )
    return(output)
  }
}
