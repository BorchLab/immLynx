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

    # Determine which TCRrep class to use based on chains
    if (length(chains) == 2 && all(c("alpha", "beta") %in% chains)) {
      # Paired alpha-beta analysis
      tr <- tc$TCRrep(
        cell_df = df_py,
        organism = organism,
        chains = reticulate::py_eval("['alpha', 'beta']"),
        compute_distances = compute_distances
      )
    } else if ("alpha" %in% chains) {
      tr <- tc$TCRrep(
        cell_df = df_py,
        organism = organism,
        chains = reticulate::py_eval("['alpha']"),
        compute_distances = compute_distances
      )
    } else {
      tr <- tc$TCRrep(
        cell_df = df_py,
        organism = organism,
        chains = reticulate::py_eval("['beta']"),
        compute_distances = compute_distances
      )
    }

    # Extract distance matrices
    result <- list()

    if ("alpha" %in% chains && !is.null(tr$pw_alpha)) {
      result$pw_alpha <- reticulate::py_to_r(tr$pw_alpha)
    }
    if ("beta" %in% chains && !is.null(tr$pw_beta)) {
      result$pw_beta <- reticulate::py_to_r(tr$pw_beta)
    }
    if (!is.null(tr$pw_cdr3_a_aa)) {
      result$pw_cdr3_a_aa <- reticulate::py_to_r(tr$pw_cdr3_a_aa)
    }
    if (!is.null(tr$pw_cdr3_b_aa)) {
      result$pw_cdr3_b_aa <- reticulate::py_to_r(tr$pw_cdr3_b_aa)
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

  clusters <- basilisk::basiliskRun(proc, function(sequences, method, inflation, eps, min_samples) {
    pd <- reticulate::import("pandas")
    clustcr <- reticulate::import("clustcr")

    # Create a DataFrame with sequences
    df <- pd$DataFrame(list(CDR3 = sequences))

    if (method == "mcl") {
      # MCL clustering
      clustering <- clustcr$Clustering(method = "MCL", mcl_params = list(inflation = inflation))
    } else if (method == "dbscan") {
      # DBSCAN clustering
      clustering <- clustcr$Clustering(method = "DBSCAN", dbscan_params = list(eps = eps, min_samples = as.integer(min_samples)))
    } else {
      stop("Unknown clustering method: ", method)
    }

    # Fit the clustering
    result <- clustering$fit(df)

    # Get cluster labels
    labels <- reticulate::py_to_r(result$cluster_df$cluster)

    labels
  }, sequences = sequences, method = method, inflation = inflation, eps = eps, min_samples = min_samples)

  return(clusters)
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

    # Load the appropriate model
    if (model == "humanTRB") {
      load_model <- olga$load_model$load_model_parms("human_T_beta")
      gen_model <- olga$generation_probability$GenerationProbabilityVDJ(load_model[[1]], load_model[[2]])
      seq_gen <- olga$sequence_generation$SequenceGenerationVDJ(load_model[[1]], load_model[[2]])
    } else if (model == "humanTRA") {
      load_model <- olga$load_model$load_model_parms("human_T_alpha")
      gen_model <- olga$generation_probability$GenerationProbabilityVJ(load_model[[1]], load_model[[2]])
      seq_gen <- olga$sequence_generation$SequenceGenerationVJ(load_model[[1]], load_model[[2]])
    } else if (model == "humanIGH") {
      load_model <- olga$load_model$load_model_parms("human_B_heavy")
      gen_model <- olga$generation_probability$GenerationProbabilityVDJ(load_model[[1]], load_model[[2]])
      seq_gen <- olga$sequence_generation$SequenceGenerationVDJ(load_model[[1]], load_model[[2]])
    } else if (model == "mouseTRB") {
      load_model <- olga$load_model$load_model_parms("mouse_T_beta")
      gen_model <- olga$generation_probability$GenerationProbabilityVDJ(load_model[[1]], load_model[[2]])
      seq_gen <- olga$sequence_generation$SequenceGenerationVDJ(load_model[[1]], load_model[[2]])
    } else {
      stop("Unknown model: ", model)
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
          pgen <- gen_model$compute_aa_CDR3_pgen(seq_aa, v_genes[i], j_genes[i])
        } else {
          # Without V/J gene info - marginalize over genes
          pgen <- gen_model$compute_aa_CDR3_pgen(seq_aa)
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

  results <- basilisk::basiliskRun(proc, function(data_folder, data_filename, pgen_filename,
                                                   organism, dataset_type, n_epochs, save_folder) {
    sonnia <- reticulate::import("sonnia")
    pd <- reticulate::import("pandas")
    np <- reticulate::import("numpy")

    # Determine chain type for model loading
    if (dataset_type == "TCR") {
      chain_type <- "human_T_beta"
    } else {
      chain_type <- "human_B_heavy"
    }

    # Load data
    data_path <- file.path(data_folder, data_filename)
    pgen_path <- file.path(data_folder, pgen_filename)

    data_seqs <- pd$read_csv(data_path)
    pgen_seqs <- pd$read_csv(pgen_path)

    # Initialize and train soNNia model
    qm <- sonnia$SoNNia(
      data_seqs = reticulate::py_to_r(data_seqs$values$tolist()),
      gen_seqs = reticulate::py_to_r(pgen_seqs$values$tolist()),
      chain_type = chain_type
    )

    # Train the model
    qm$infer_selection(epochs = as.integer(n_epochs))

    # Get selection factors
    selection_factors <- qm$compute_Q(qm$data_seqs)

    # Save the model
    if (!is.null(save_folder)) {
      dir.create(save_folder, recursive = TRUE, showWarnings = FALSE)
      qm$save_model(save_folder)
    }

    list(
      selection_factors = reticulate::py_to_r(selection_factors),
      model_path = save_folder
    )
  }, data_folder = data_folder, data_filename = data_filename, pgen_filename = pgen_filename,
     organism = organism, dataset_type = dataset_type, n_epochs = n_epochs, save_folder = save_folder)

  return(results)
}
