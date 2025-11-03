#' Run DeepTCR on scRepertoire Data
#'
#' @description Extract TCR sequences and run DeepTCR unsupervised feature extraction.
#'   Note: This function requires preparing data files in a specific format.
#'
#' @param seurat_obj A Seurat object containing scRepertoire TCR data
#' @param chains Which chain to use: "TRB" or "TRA". Default is "TRB".
#' @param output_dir Directory to save DeepTCR models and results
#' @param latent_dim Dimensionality of latent space. Default is 100.
#' @param epochs Number of training epochs. Default is 100.
#' @param reduction_name Name for dimensional reduction in Seurat. Default is "deeptcr".
#' @param return_seurat If TRUE, adds features to Seurat object. Default is TRUE.
#'
#' @details
#' DeepTCR requires creating intermediate CSV files with specific columns. This
#' function handles the data preparation automatically.
#'
#' @return If return_seurat=TRUE, returns Seurat object with DeepTCR features as
#'   a dimensional reduction. Otherwise returns the feature matrix.
#'
#' @export
#' @importFrom immApex getIR
#'
#' @examples
#' \dontrun{
#'   seurat_obj <- runDeepTCR(seurat_obj, 
#'                            chains = "TRB",
#'                            output_dir = "deeptcr_output")
#' }
runDeepTCR <- function(seurat_obj,
                       chains = c("TRB", "TRA"),
                       output_dir = "deeptcr_output",
                       latent_dim = 100,
                       epochs = 100,
                       reduction_name = "deeptcr",
                       return_seurat = TRUE) {
  
  chains <- match.arg(chains)
  
  message("Extracting ", chains, " sequences from Seurat object...")
  tcr_data <- immApex::getIR(seurat_obj, chains = chains)
  tcr_data <- tcr_data[!is.na(tcr_data$cdr3_aa), ]
  
  if (nrow(tcr_data) == 0) {
    stop("No valid ", chains, " sequences found.")
  }
  
  # Create data directory
  data_dir <- file.path(output_dir, "data")
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE)
  }
  
  message("Preparing data for DeepTCR...")
  
  # DeepTCR expects specific column names based on chain
  col_name <- ifelse(chains == "TRB", "beta", "alpha")
  
  # Create data file
  deeptcr_df <- data.frame(
    stringsAsFactors = FALSE
  )
  deeptcr_df[[col_name]] <- tcr_data$cdr3_aa
  
  # Add V and J genes if available
  if (!all(is.na(tcr_data$v))) {
    deeptcr_df[[paste0("v_", col_name)]] <- tcr_data$v
  }
  if (!all(is.na(tcr_data$j))) {
    deeptcr_df[[paste0("j_", col_name)]] <- tcr_data$j
  }
  
  # Save to file
  data_file <- file.path(data_dir, "tcr_data.csv")
  write.csv(deeptcr_df, data_file, row.names = FALSE)
  
  message("Running DeepTCR VAE...")
  
  # Run DeepTCR
  features <- calculate.deepTCR(
    output_dir = output_dir,
    data_dir = data_dir,
    files = c("tcr_data.csv"),
    latent_dim = latent_dim,
    epochs = epochs
  )
  
  # Add barcodes as rownames
  rownames(features) <- tcr_data$barcode
  
  if (return_seurat) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Package 'Seurat' is required.")
    }
    
    # Create dimensional reduction
    cell_features <- matrix(NA,
                           nrow = ncol(seurat_obj),
                           ncol = ncol(features),
                           dimnames = list(colnames(seurat_obj),
                                         paste0("DeepTCR_", 1:ncol(features))))
    
    cell_features[tcr_data$barcode, ] <- features
    
    deeptcr_reduction <- Seurat::CreateDimReducObject(
      embeddings = cell_features,
      key = "DeepTCR_",
      assay = Seurat::DefaultAssay(seurat_obj)
    )
    
    seurat_obj[[reduction_name]] <- deeptcr_reduction
    
    message("DeepTCR features added as '", reduction_name, "' reduction")
    return(seurat_obj)
    
  } else {
    return(features)
  }
}


