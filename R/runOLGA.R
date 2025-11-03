#' Calculate Generation Probability (Pgen) for TCRs in scRepertoire Data
#'
#' @description Extracts TCR sequences from a Seurat object and calculates their
#'   generation probability using OLGA.
#'
#' @param seurat_obj A Seurat object containing scRepertoire TCR data.
#' @param chains Which chain to analyze: "TRA" or "TRB". Default is "TRB".
#' @param model OLGA model to use. Options: "humanTRB", "humanTRA", "humanIGH", "mouseTRB".
#'   If NULL, will be inferred from organism and chains parameters.
#' @param organism Organism: "human" or "mouse". Used if model is NULL. Default is "human".
#' @param use_vj_genes Logical. If TRUE, includes V and J gene information in Pgen calculation.
#'   Default is FALSE (sequence-only Pgen).
#' @param return_seurat Logical. If TRUE, adds Pgen values to metadata. Default is TRUE.
#' @param column_name Name for the metadata column. Default is "olga_pgen".
#'
#' @return If return_seurat=TRUE, returns Seurat object with Pgen added to metadata.
#'   If FALSE, returns data.frame with barcodes, sequences, and Pgen values.
#'
#' @export
#' @importFrom immApex getIR
#'
#' @examples
#' \dontrun{
#'   # Calculate Pgen for TRB sequences
#'   seurat_obj <- runOLGA(seurat_obj, chains = "TRB", model = "humanTRB")
#'
#'   # Calculate Pgen with V/J gene info
#'   seurat_obj <- runOLGA(seurat_obj, chains = "TRB", use_vj_genes = TRUE)
#'
#'   # Just get results without adding to object
#'   pgen_results <- runOLGA(seurat_obj, chains = "TRB", return_seurat = FALSE)
#' }
runOLGA <- function(seurat_obj,
                    chains = c("TRB", "TRA"),
                    model = NULL,
                    organism = "human",
                    use_vj_genes = FALSE,
                    return_seurat = TRUE,
                    column_name = "olga_pgen") {
  
  chains <- match.arg(chains)
  
  # Infer model if not specified
  if (is.null(model)) {
    model <- paste0(organism, chains)
    if (!model %in% c("humanTRB", "humanTRA", "humanIGH", "mouseTRB")) {
      stop("Cannot infer OLGA model. Please specify explicitly using the 'model' parameter.")
    }
  }
  
  message("Extracting ", chains, " sequences from Seurat object...")
  tcr_data <- immApex::getIR(seurat_obj, chains = chains)
  tcr_data <- tcr_data[!is.na(tcr_data$cdr3_aa), ]
  
  if (nrow(tcr_data) == 0) {
    stop("No valid ", chains, " sequences found.")
  }
  
  message("Found ", nrow(tcr_data), " valid sequences")
  message("Calculating generation probabilities with OLGA (", model, ")...")
  
  # Prepare arguments for calculate.olga
  olga_args <- list(
    action = "pgen",
    model = model,
    sequences = tcr_data$cdr3_aa
  )
  
  # Add V and J gene info if requested
  if (use_vj_genes) {
    if (all(!is.na(tcr_data$v)) && all(!is.na(tcr_data$j))) {
      olga_args$v_genes <- tcr_data$v
      olga_args$j_genes <- tcr_data$j
      message("Including V and J gene information in Pgen calculation")
    } else {
      warning("use_vj_genes=TRUE but some V/J genes are NA. Proceeding without gene info.")
    }
  }
  
  # Calculate Pgen
  pgen_values <- do.call(calculate.olga, olga_args)
  
  # Create result data frame
  result_df <- data.frame(
    barcode = tcr_data$barcode,
    cdr3_aa = tcr_data$cdr3_aa,
    v_gene = tcr_data$v,
    j_gene = tcr_data$j,
    pgen = pgen_values,
    log10_pgen = log10(pgen_values),
    stringsAsFactors = FALSE
  )
  
  if (return_seurat) {
    # Add to metadata
    pgen_col <- rep(NA_real_, ncol(seurat_obj))
    names(pgen_col) <- colnames(seurat_obj)
    pgen_col[result_df$barcode] <- result_df$pgen
    
    log_pgen_col <- rep(NA_real_, ncol(seurat_obj))
    names(log_pgen_col) <- colnames(seurat_obj)
    log_pgen_col[result_df$barcode] <- result_df$log10_pgen
    
    seurat_obj[[paste0(column_name, "_", chains)]] <- pgen_col
    seurat_obj[[paste0(column_name, "_log10_", chains)]] <- log_pgen_col
    
    message("Pgen values added to metadata as '", 
            paste0(column_name, "_", chains), "' and '",
            paste0(column_name, "_log10_", chains), "'")
    return(seurat_obj)
  } else {
    return(result_df)
  }
}


#' Generate Random TCR Sequences using OLGA
#'
#' @description Generate random TCR sequences from OLGA's generative model.
#'
#' @param n Number of sequences to generate.
#' @param model OLGA model to use: "humanTRB", "humanTRA", "humanIGH", "mouseTRB".
#'
#' @return Data.frame with generated sequences (nt_seq, aa_seq, v_index, j_index).
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Generate 100 human TRB sequences
#'   random_seqs <- generateOLGA(n = 100, model = "humanTRB")
#' }
generateOLGA <- function(n = 100, model = "humanTRB") {
  message("Generating ", n, " random sequences with OLGA (", model, ")...")
  
  result <- calculate.olga(
    action = "generate",
    model = model,
    n = n
  )
  
  message("Generation complete.")
  return(result)
}
