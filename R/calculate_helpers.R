#' @title Internal Python Bridge Functions
#' @name calculate_helpers
#' @description Internal functions that use basilisk to call Python libraries.
#'   These functions are not exported and are used by the run* wrapper functions.
#' @keywords internal
NULL

#' Calculate TCR Distances using tcrdist3
#' @description Internal function that calls tcrdist3 via basilisk.
#' @param df A data.frame with TCR sequences in tcrdist3 format
#' @param organism Character: "human" or "mouse"
#' @param chains Character vector: "alpha", "beta", or c("alpha", "beta")
#' @param compute_distances Logical: whether to compute full distance matrix
#' @return A list containing distance matrices
#' @keywords internal
calculate.tcrDist <- function(df,
                              organism = "human",
                              chains = "beta",
                              compute_distances = TRUE) {

  proc <- basilisk::basiliskStart(immLynxEnv)
  on.exit(basilisk::basiliskStop(proc))

  results <- basilisk::basiliskRun(proc, function(df, organism, chains, compute_distances) {
    pd <- reticulate::import("pandas")
    tcrdist <- reticulate::import("tcrdist")
    tc <- tcrdist$repertoire

    # Convert R data.frame to pandas DataFrame
    df_py <- pd$DataFrame(df)

    # Build chains as R list so reticulate converts to Python list
    # (R list -> Python list, R character vector -> Python str)
    if (length(chains) == 2 && all(c("alpha", "beta") %in% chains)) {
      py_chains <- list("alpha", "beta")
    } else if ("alpha" %in% chains) {
      py_chains <- list("alpha")
    } else {
      py_chains <- list("beta")
    }

    tr <- tc$TCRrep(
      cell_df = df_py,
      organism = organism,
      chains = py_chains,
      compute_distances = compute_distances
    )

    # Extract distance matrices using py_has_attr to avoid AttributeError
    result <- list()

    .safe_extract <- function(obj, attr_name) {
      if (reticulate::py_has_attr(obj, attr_name)) {
        val <- obj[[attr_name]]
        if (!is.null(val)) {
          return(reticulate::py_to_r(val))
        }
      }
      NULL
    }

    if ("alpha" %in% chains) {
      result$pw_alpha <- .safe_extract(tr, "pw_alpha")
      result$pw_cdr3_a_aa <- .safe_extract(tr, "pw_cdr3_a_aa")
    }
    if ("beta" %in% chains) {
      result$pw_beta <- .safe_extract(tr, "pw_beta")
      result$pw_cdr3_b_aa <- .safe_extract(tr, "pw_cdr3_b_aa")
    }

    result
  }, df = df, organism = organism, chains = chains, compute_distances = compute_distances)

  return(results)
}


#' Cluster TCRs using clusTCR
#' @description Internal function that calls clusTCR via basilisk.
#' @param sequences Character vector of CDR3 amino acid sequences
#' @param method Clustering method: "mcl" or "dbscan"
#' @param inflation MCL inflation parameter (for mcl method)
#' @param eps DBSCAN epsilon parameter (for dbscan method)
#' @param min_samples DBSCAN minimum samples parameter (for dbscan method)
#' @return Integer vector of cluster assignments
#' @keywords internal
calculate.clustcr <- function(sequences,
                              method = "mcl",
                              inflation = 2.0,
                              eps = 0.5,
                              min_samples = 2) {

  proc <- basilisk::basiliskStart(immLynxEnv)
  on.exit(basilisk::basiliskStop(proc))

  cluster_df <- basilisk::basiliskRun(proc, function(sequences,
      method, inflation, eps, min_samples) {
    pd <- reticulate::import("pandas")
    clustcr <- reticulate::import("clustcr")

    # clusTCR fit() expects an iterable of CDR3 sequences
    cdr3_series <- pd$Series(sequences)

    # MCL clustering - mcl_params is [inflation, expansion]
    clustering <- clustcr$Clustering(
        method = "MCL",
        mcl_params = list(1.2, inflation)
    )

    result <- clustering$fit(cdr3_series)

    # clusters_df has columns: junction_aa, cluster
    # Only clustered sequences are returned (singletons dropped)
    reticulate::py_to_r(result$clusters_df)
  }, sequences = sequences, method = method,
     inflation = inflation, eps = eps,
     min_samples = min_samples)

  # Map cluster assignments back to input sequences
  # Sequences not in a cluster get NA
  cluster_map <- stats::setNames(
      cluster_df$cluster,
      cluster_df$junction_aa
  )
  labels <- unname(cluster_map[sequences])

  return(labels)
}


#' Calculate Generation Probability using OLGA
#' @description Internal function that calls OLGA via basilisk.
#' @param action Either "pgen" to calculate probability or "generate" to generate sequences
#' @param model OLGA model name
#' @param sequences Character vector of CDR3 sequences (for pgen)
#' @param v_genes Optional V gene annotations
#' @param j_genes Optional J gene annotations
#' @param n Number of sequences to generate (for generate action)
#' @return For pgen: numeric vector of probabilities. For generate: data.frame
#' @keywords internal
calculate.olga <- function(action = c("pgen", "generate"),
                           model = "humanTRB",
                           sequences = NULL,
                           v_genes = NULL,
                           j_genes = NULL,
                           n = 100) {

  action <- match.arg(action)

  proc <- basilisk::basiliskStart(immLynxEnv)
  on.exit(basilisk::basiliskStop(proc))

  result <- basilisk::basiliskRun(proc, function(action, model, sequences, v_genes, j_genes, n) {
    olga <- reticulate::import("olga")
    os <- reticulate::import("os")

    # Map model names to OLGA default model directories and chain types
    model_map <- list(
      humanTRB  = list(dir = "human_T_beta",  type = "VDJ"),
      humanTRA  = list(dir = "human_T_alpha", type = "VJ"),
      humanIGH  = list(dir = "human_B_heavy", type = "VDJ"),
      mouseTRB  = list(dir = "mouse_T_beta",  type = "VDJ")
    )

    if (!model %in% names(model_map)) {
      stop("Unknown model: ", model)
    }

    model_info <- model_map[[model]]
    olga_dir <- os$path$dirname(olga$`__file__`)
    model_path <- os$path$join(olga_dir, "default_models", model_info$dir)

    marginals_file <- os$path$join(model_path, "model_marginals.txt")
    params_file    <- os$path$join(model_path, "model_params.txt")
    v_anchor_file  <- os$path$join(model_path, "V_gene_CDR3_anchors.csv")
    j_anchor_file  <- os$path$join(model_path, "J_gene_CDR3_anchors.csv")

    # Load genomic data and generative model
    if (model_info$type == "VDJ") {
      genomic_data <- olga$load_model$GenomicDataVDJ()
      genomic_data$load_igor_genomic_data(params_file, v_anchor_file, j_anchor_file)
      gen_model_data <- olga$load_model$GenerativeModelVDJ()
      gen_model_data$load_and_process_igor_model(marginals_file)
      pgen_calc <- olga$generation_probability$GenerationProbabilityVDJ(gen_model_data, genomic_data)
      seq_gen   <- olga$sequence_generation$SequenceGenerationVDJ(gen_model_data, genomic_data)
    } else {
      genomic_data <- olga$load_model$GenomicDataVJ()
      genomic_data$load_igor_genomic_data(params_file, v_anchor_file, j_anchor_file)
      gen_model_data <- olga$load_model$GenerativeModelVJ()
      gen_model_data$load_and_process_igor_model(marginals_file)
      pgen_calc <- olga$generation_probability$GenerationProbabilityVJ(gen_model_data, genomic_data)
      seq_gen   <- olga$sequence_generation$SequenceGenerationVJ(gen_model_data, genomic_data)
    }

    if (action == "pgen") {
      # Calculate generation probabilities
      pgen_values <- numeric(length(sequences))

      for (i in seq_along(sequences)) {
        seq_aa <- sequences[i]

        # Calculate Pgen
        if (!is.null(v_genes) && !is.null(j_genes) &&
            !is.na(v_genes[i]) && !is.na(j_genes[i])) {
          # With V/J gene info
          pgen <- pgen_calc$compute_aa_CDR3_pgen(seq_aa, v_genes[i], j_genes[i])
        } else {
          # Without V/J gene info - marginalize over genes
          pgen <- pgen_calc$compute_aa_CDR3_pgen(seq_aa)
        }

        pgen_values[i] <- pgen
      }

      return(pgen_values)

    } else {
      # Generate sequences
      generated <- list(
        nt_seq = character(n),
        aa_seq = character(n),
        v_index = integer(n),
        j_index = integer(n)
      )

      for (i in seq_len(n)) {
        gen_seq <- seq_gen$gen_rnd_prod_CDR3()
        generated$nt_seq[i] <- gen_seq[[1]]
        generated$aa_seq[i] <- gen_seq[[2]]
        generated$v_index[i] <- gen_seq[[3]]
        generated$j_index[i] <- gen_seq[[4]]
      }

      return(as.data.frame(generated))
    }
  }, action = action, model = model, sequences = sequences,
     v_genes = v_genes, j_genes = j_genes, n = n)

  return(result)
}


#' Run soNNia Selection Analysis
#' @description Internal function that calls soNNia via basilisk.
#' @param data_folder Directory containing data files
#' @param data_filename Name of selected TCR data file
#' @param pgen_filename Name of background Pgen file
#' @param organism Organism: "human" or "mouse"
#' @param dataset_type Type of dataset: "TCR" or "BCR"
#' @param n_epochs Number of training epochs
#' @param save_folder Directory for saving outputs
#' @return List with selection model results
#' @keywords internal
calculate.sonia <- function(data_folder,
                            data_filename,
                            pgen_filename,
                            organism = "human",
                            dataset_type = "TCR",
                            n_epochs = 100,
                            save_folder = "sonia_output") {

  proc <- basilisk::basiliskStart(immLynxEnv)
  on.exit(basilisk::basiliskStop(proc))

  # Read and prepare data as lists of lists: [[cdr3, v_gene, j_gene], ...]
  data_path <- file.path(data_folder, data_filename)
  pgen_path <- file.path(data_folder, pgen_filename)

  selected_csv <- utils::read.csv(data_path, stringsAsFactors = FALSE)
  background_csv <- utils::read.csv(pgen_path, stringsAsFactors = FALSE)

  # Build data_seqs: list of character vectors [cdr3, v, j]
  sel_cols <- names(selected_csv)
  cdr3_col <- grep("cdr3", sel_cols, value = TRUE)[1]
  v_col <- grep("^v", sel_cols, value = TRUE)[1]
  j_col <- grep("^j", sel_cols, value = TRUE)[1]

  data_seqs <- lapply(seq_len(nrow(selected_csv)), function(i) {
    c(selected_csv[[cdr3_col]][i],
      if (!is.na(v_col)) selected_csv[[v_col]][i] else "",
      if (!is.na(j_col)) selected_csv[[j_col]][i] else "")
  })

  # Build gen_seqs from background (OLGA output has aa_seq column)
  bg_cols <- names(background_csv)
  aa_col <- grep("aa_seq|aa", bg_cols, value = TRUE)[1]

  gen_seqs <- lapply(seq_len(nrow(background_csv)), function(i) {
    c(background_csv[[aa_col]][i], "", "")
  })

  results <- basilisk::basiliskRun(proc, function(data_seqs, gen_seqs,
                                                   dataset_type, n_epochs, save_folder) {
    sonnia <- reticulate::import("sonnia.sonnia")

    # Determine model identifier for the chain type
    model_id <- if (dataset_type == "TCR") "humanTRB" else "humanIGH"

    # Initialize soNNia model — try pgen_model (>=0.4.0) then chain_type (<=0.1.x)
    qm <- tryCatch(
      sonnia$SoNNia(
        data_seqs = data_seqs,
        gen_seqs = gen_seqs,
        pgen_model = model_id
      ),
      error = function(e) {
        sonnia$SoNNia(
          data_seqs = data_seqs,
          gen_seqs = gen_seqs,
          chain_type = model_id
        )
      }
    )

    # Train the model
    qm$infer_selection(epochs = as.integer(n_epochs))

    # Get selection factors — try evaluate_selection_factors (>=0.4.0)
    # then compute_Q (<=0.1.x)
    selection_factors <- tryCatch(
      reticulate::py_to_r(qm$evaluate_selection_factors(qm$data_seqs)),
      error = function(e) {
        reticulate::py_to_r(qm$compute_Q(qm$data_seqs))
      }
    )

    # Save the model
    if (!is.null(save_folder)) {
      dir.create(save_folder, recursive = TRUE, showWarnings = FALSE)
      qm$save_model(save_folder)
    }

    list(
      selection_factors = selection_factors,
      model_path = save_folder
    )
  }, data_seqs = data_seqs, gen_seqs = gen_seqs,
     dataset_type = dataset_type, n_epochs = n_epochs, save_folder = save_folder)

  return(results)
}
