# Predict T-cell clonal expansion from gene expression using scXpand.
# Public entry points: runScXpand(), listScXpandModels().
#
# Unlike every other tool wrapped by immLynx, scXpand consumes gene
# expression rather than receptor sequences: it infers whether a T cell
# belongs to an expanded clone without paired TCR sequencing.  When the
# object does carry scRepertoire clone calls, those are used to derive
# ground-truth labels so the prediction can be benchmarked.

#' @keywords internal
# Run a function inside the scXpandEnv basilisk env. Mirrors
# .run_in_scanpy_env in R/exportToScanpy.R but targets the torch/scxpand
# stack, which needs its own python (>= 3.11).
.run_in_scxpand_env <- function(FUN, ...) {
  proc <- basilisk::basiliskStart(scXpandEnv)
  on.exit(basilisk::basiliskStop(proc))
  basilisk::basiliskRun(proc, FUN, ...)
}

# Static mirror of scxpand.pretrained.model_registry.PRETRAINED_MODELS.
# Kept in R so listScXpandModels() can answer without building a
# multi-gigabyte conda environment; listScXpandModels(refresh = TRUE)
# queries the installed package and is authoritative.
.SCXPAND_MODELS <- data.frame(
  model_name = c("pan_cancer_autoencoder", "pan_cancer_mlp",
                 "pan_cancer_lightgbm", "pan_cancer_logistic",
                 "pan_cancer_svm"),
  model_type = c("autoencoder", "mlp", "lightgbm", "logistic", "svm"),
  version    = rep("1.0.0", 5L),
  description = c(
    "Pan-cancer autoencoder (scXpand default)",
    "Pan-cancer multi-layer perceptron",
    "Pan-cancer LightGBM gradient boosting",
    "Pan-cancer logistic regression",
    "Pan-cancer support vector machine"
  ),
  stringsAsFactors = FALSE
)

# ===========================================================================
# Gene identifier resolution
# ===========================================================================

#' @keywords internal
.isEnsembl <- function(x) {
  out <- grepl("^ENSG[0-9]{11}$", x)
  out[is.na(x)] <- FALSE
  out
}

#' @keywords internal
# Strip trailing version suffixes, but only from strings that already look
# like versioned Ensembl gene IDs. A blanket sub("\\.[0-9]+$", "", x) would
# mangle legitimate symbols such as MARCH1.2 or 7SK.2.
.stripEnsemblVersion <- function(x) {
  versioned <- grepl("^ENSG[0-9]{11}\\.[0-9]+$", x)
  versioned[is.na(x)] <- FALSE
  x[versioned] <- sub("\\.[0-9]+$", "", x[versioned])
  x
}

#' @keywords internal
# Map gene symbols to Ensembl gene IDs via org.Hs.eg.db, trying SYMBOL
# first and falling back to ALIAS for unmatched keys.
#
# One symbol routinely maps to several ENSG IDs (HLA-DRA -> 8, HSPA1A -> 5,
# alt haplotypes and patch scaffolds). scXpand silently zero-fills any gene
# it cannot find in its panel, so an unlucky pick deletes that gene rather
# than erroring. Every ambiguous case is therefore counted and reported.
.mapSymbolsToEnsembl <- function(symbols, multi_map = c("first", "expand",
                                                        "drop")) {
  multi_map <- match.arg(multi_map)

  empty <- list(map = data.frame(SYMBOL = character(0),
                                 ENSEMBL = character(0),
                                 stringsAsFactors = FALSE),
                n_query = 0L, n_symbol = 0L, n_alias = 0L,
                n_ambiguous = 0L, n_unmapped = 0L)

  uq <- unique(symbols[!is.na(symbols) & nzchar(symbols)])
  if (!length(uq)) return(empty)

  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE) ||
      !requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("Symbol to Ensembl mapping requires org.Hs.eg.db and AnnotationDbi. ",
         "Install with: BiocManager::install(c('org.Hs.eg.db', ",
         "'AnnotationDbi')), or pass map_symbols = 'never' and supply ",
         "Ensembl IDs via gene_ids.", call. = FALSE)
  }

  db <- org.Hs.eg.db::org.Hs.eg.db

  .lookup <- function(keys, keytype) {
    if (!length(keys)) return(NULL)
    res <- tryCatch(
      suppressWarnings(suppressMessages(
        AnnotationDbi::select(db, keys = keys, keytype = keytype,
                              columns = "ENSEMBL"))),
      error = function(e) NULL)
    if (is.null(res) || !nrow(res)) return(NULL)
    res <- res[!is.na(res$ENSEMBL), , drop = FALSE]
    if (!nrow(res)) return(NULL)
    data.frame(SYMBOL = as.character(res[[keytype]]),
               ENSEMBL = as.character(res$ENSEMBL),
               stringsAsFactors = FALSE)
  }

  m <- .lookup(uq, "SYMBOL")
  n_symbol <- if (is.null(m)) 0L else length(unique(m$SYMBOL))

  unmatched <- setdiff(uq, if (is.null(m)) character(0) else m$SYMBOL)
  a <- .lookup(unmatched, "ALIAS")
  n_alias <- if (is.null(a)) 0L else length(unique(a$SYMBOL))

  m <- rbind(m, a)
  if (is.null(m) || !nrow(m)) {
    empty$n_query <- length(uq)
    empty$n_unmapped <- length(uq)
    return(empty)
  }

  m <- m[!duplicated(m), , drop = FALSE]
  # Lexicographic order makes multi_map = "first" reproducible across
  # machines and org.Hs.eg.db releases; AnnotationDbi's own
  # multiVals = "first" follows database row order, which is not stable.
  m <- m[order(m$SYMBOL, m$ENSEMBL), , drop = FALSE]

  per_symbol <- table(m$SYMBOL)
  ambiguous <- names(per_symbol)[per_symbol > 1L]

  m <- switch(multi_map,
              first  = m[!duplicated(m$SYMBOL), , drop = FALSE],
              drop   = m[!(m$SYMBOL %in% ambiguous), , drop = FALSE],
              expand = m)

  list(map = m,
       n_query = length(uq),
       n_symbol = n_symbol,
       n_alias = n_alias,
       n_ambiguous = length(ambiguous),
       n_unmapped = length(setdiff(uq, m$SYMBOL)))
}

#' @keywords internal
# Turn whatever identifiers the object carries into Ensembl gene IDs, and
# work out which rows of the object survive.
#
# Returns:
#   ids             final rownames for the staged matrix
#   keep            row indices into the input object; may repeat when
#                   multi_map = "expand", and is longer than `ids` when
#                   collapse_groups is non-NULL
#   collapse_groups NULL, or a factor over `keep` whose levels are `ids`,
#                   meaning those rows must be summed together
#   report          named list summarising the mapping, surfaced to the
#                   user and stored in metadata()$scXpand$gene_mapping
.resolveGeneIDs <- function(sce, gene_ids = "rownames",
                            map_symbols = c("auto", "always", "never"),
                            ensembl_min_frac = 0.5,
                            multi_map = c("first", "expand", "drop"),
                            collapse = c("sum", "first", "drop"),
                            verbose = TRUE) {
  map_symbols <- match.arg(map_symbols)
  multi_map   <- match.arg(multi_map)
  collapse    <- match.arg(collapse)

  n_in <- nrow(sce)

  # --- Extract the raw identifiers -------------------------------------
  if (is.character(gene_ids) && length(gene_ids) == 1L &&
      identical(gene_ids, "rownames")) {
    raw <- rownames(sce)
    if (is.null(raw)) {
      stop("gene_ids = 'rownames' but the object has no rownames; pass a ",
           "rowData column name or a character vector instead.",
           call. = FALSE)
    }
    src <- "rownames"
  } else if (is.character(gene_ids) && length(gene_ids) == 1L) {
    rd <- SummarizedExperiment::rowData(sce)
    if (!(gene_ids %in% colnames(rd))) {
      stop("gene_ids must be \"rownames\", the name of a rowData column, ",
           "or a character vector of length nrow(input). ",
           "'", gene_ids, "' is not a rowData column",
           if (ncol(rd)) paste0(" (available: ",
                                paste(colnames(rd), collapse = ", "), ")")
           else " (the object has no rowData columns)", ".",
           call. = FALSE)
    }
    raw <- as.character(rd[[gene_ids]])
    src <- paste0("rowData$", gene_ids)
  } else if (is.character(gene_ids) && length(gene_ids) == n_in) {
    raw <- as.character(gene_ids)
    src <- "user-supplied vector"
  } else {
    stop("gene_ids must be \"rownames\", the name of a rowData column, or ",
         "a character vector of length nrow(input) (", n_in, "); got ",
         class(gene_ids)[1], " of length ", length(gene_ids), ".",
         call. = FALSE)
  }

  raw <- trimws(raw)
  raw[!nzchar(raw)] <- NA_character_
  raw <- .stripEnsemblVersion(raw)

  is_ens <- .isEnsembl(raw)
  frac <- mean(is_ens)

  # --- Decide whether to map -------------------------------------------
  do_map <- switch(map_symbols,
                   always = TRUE,
                   never  = FALSE,
                   auto   = frac < ensembl_min_frac)

  if (map_symbols == "never" && frac < ensembl_min_frac) {
    offenders <- utils::head(unique(raw[!is_ens & !is.na(raw)]), 5L)
    stop("Only ", sprintf("%.1f%%", 100 * frac), " of gene identifiers from ",
         src, " look like Ensembl gene IDs (need at least ",
         sprintf("%.1f%%", 100 * ensembl_min_frac), "). ",
         "Examples: ", paste(offenders, collapse = ", "), ". ",
         "scXpand's pretrained models are indexed by Ensembl ID. ",
         "Use map_symbols = 'auto' to map symbols via org.Hs.eg.db, or ",
         "supply Ensembl IDs through gene_ids.", call. = FALSE)
  }

  # --- Ensembl rows pass straight through -------------------------------
  keep <- which(is_ens)
  ids  <- raw[keep]

  mapinfo <- NULL
  if (do_map) {
    mapinfo <- .mapSymbolsToEnsembl(raw[!is_ens], multi_map = multi_map)
    m <- mapinfo$map
    if (nrow(m)) {
      by_sym <- split(m$ENSEMBL, m$SYMBOL)
      n_hits <- lengths(by_sym)[raw]
      n_hits[is.na(n_hits)] <- 0L
      n_hits[is_ens] <- 0L                 # already handled above
      hit_rows <- which(n_hits > 0L)
      keep <- c(keep, rep(hit_rows, times = n_hits[hit_rows]))
      ids  <- c(ids, unlist(by_sym[raw[hit_rows]], use.names = FALSE))
    }
    if (mapinfo$n_query > 0L &&
        (mapinfo$n_query - mapinfo$n_unmapped) / mapinfo$n_query < 0.25) {
      warning("Fewer than 25% of gene symbols mapped to Ensembl IDs. ",
              "Is this a human dataset, and are the identifiers really ",
              "gene symbols?", call. = FALSE)
    }
  }

  ord  <- order(keep)
  keep <- keep[ord]
  ids  <- ids[ord]

  if (!length(ids)) {
    stop("No gene identifiers could be resolved to Ensembl IDs from ", src,
         ". scXpand cannot run without them.", call. = FALSE)
  }

  # --- Post-mapping threshold re-check ----------------------------------
  post_frac <- length(unique(ids)) / max(1L, length(unique(raw[!is.na(raw)])))
  if (do_map && map_symbols == "auto" && post_frac < ensembl_min_frac) {
    stop("After symbol mapping only ", sprintf("%.1f%%", 100 * post_frac),
         " of input features resolved to Ensembl IDs (need at least ",
         sprintf("%.1f%%", 100 * ensembl_min_frac), "). ",
         "Lower ensembl_min_frac if this is expected (targeted panel), or ",
         "supply Ensembl IDs through gene_ids.", call. = FALSE)
  }

  # --- Collapse duplicate Ensembl IDs -----------------------------------
  n_precollapse <- length(ids)
  dup_ids <- unique(ids[duplicated(ids)])
  collapse_groups <- NULL

  if (length(dup_ids)) {
    if (collapse == "sum") {
      collapse_groups <- factor(ids)
      ids <- levels(collapse_groups)
    } else if (collapse == "first") {
      first <- !duplicated(ids)
      keep <- keep[first]
      ids  <- ids[first]
    } else {
      okay <- !(ids %in% dup_ids)
      keep <- keep[okay]
      ids  <- ids[okay]
    }
  }

  report <- list(
    source            = src,
    n_input           = n_in,
    n_already_ensembl = sum(is_ens),
    frac_ensembl_in   = frac,
    mapped            = do_map,
    n_mapped_symbol   = if (is.null(mapinfo)) 0L else mapinfo$n_symbol,
    n_mapped_alias    = if (is.null(mapinfo)) 0L else mapinfo$n_alias,
    n_unmapped        = if (is.null(mapinfo)) 0L else mapinfo$n_unmapped,
    n_ambiguous       = if (is.null(mapinfo)) 0L else mapinfo$n_ambiguous,
    multi_map         = if (do_map) multi_map else NA_character_,
    collapse          = if (length(dup_ids)) collapse else NA_character_,
    n_collapsed       = length(dup_ids),
    n_precollapse     = n_precollapse,
    n_output          = length(ids)
  )

  if (verbose) message(.formatGeneMapping(report))

  list(ids = ids, keep = keep, collapse_groups = collapse_groups,
       report = report)
}

#' @keywords internal
.formatGeneMapping <- function(r) {
  lines <- c(
    sprintf("scXpand gene ID resolution (%s):", r$source),
    sprintf("  input features:       %d", r$n_input),
    sprintf("  already Ensembl:      %d (%.1f%%)", r$n_already_ensembl,
            100 * r$frac_ensembl_in))
  if (isTRUE(r$mapped)) {
    lines <- c(lines,
      sprintf("  mapped via SYMBOL:    %d", r$n_mapped_symbol),
      sprintf("  mapped via ALIAS:     %d", r$n_mapped_alias),
      sprintf("  unmapped (dropped):   %d", r$n_unmapped),
      sprintf("  ambiguous 1:many:     %d symbols (multi_map = '%s')",
              r$n_ambiguous, r$multi_map))
  }
  if (r$n_collapsed > 0L) {
    lines <- c(lines,
      sprintf("  collapsed many:1:     %d Ensembl IDs received >1 row (collapse = '%s')",
              r$n_collapsed, r$collapse))
  }
  c(paste(c(lines, sprintf("  features written:     %d", r$n_output)),
          collapse = "\n"))
}

# ===========================================================================
# Ground-truth expansion labels from scRepertoire clone calls
# ===========================================================================

#' @keywords internal
# Derive scXpand's ground-truth obs fields from clone identity.
#
# scXpand defines a cell as expanded when its clone_id_size exceeds
# 1.5 x median_clone_size for that sample. Clone sizes are tabulated fresh
# here rather than read from scRepertoire's clonalFrequency column, because
# combineExpression() computes that under whatever group.by was in effect
# at combine time (or globally) and it is not guaranteed to be a per-sample
# count.
#
# Returns NULL when the object carries no usable clone data.
.deriveExpansionLabels <- function(sce, clone_col = NULL, sample_col = NULL,
                                   median_basis = c("clone", "cell"),
                                   verbose = TRUE) {
  median_basis <- match.arg(median_basis)
  cd <- SummarizedExperiment::colData(sce)

  if (is.null(clone_col)) {
    candidates <- intersect(c("CTstrict", "CTaa"), colnames(cd))
    if (!length(candidates)) return(NULL)
    clone_col <- candidates[1L]
  } else if (!(clone_col %in% colnames(cd))) {
    stop("clone_col '", clone_col, "' not found in colData.", call. = FALSE)
  }

  cid <- as.character(cd[[clone_col]])
  # scRepertoire writes literal "NA" tokens into the CT* columns for cells
  # with a missing chain; those are absent clone calls, not clone names.
  cid[!nzchar(cid) | cid %in% c("NA", "None", "NA_NA")] <- NA_character_
  has <- !is.na(cid)
  if (!any(has)) return(NULL)

  if (is.null(sample_col)) {
    warning("sample_col is NULL; treating all ", ncol(sce),
            " cells as a single sample. scXpand's expansion rule is defined ",
            "per sample -- pass sample_col = '<colData column>' (e.g. ",
            "'orig.ident' or 'Patient') for correct labels.", call. = FALSE)
    smp <- rep("__all__", ncol(sce))
  } else {
    if (!(sample_col %in% colnames(cd))) {
      stop("sample_col '", sample_col, "' not found in colData.",
           call. = FALSE)
    }
    smp <- as.character(cd[[sample_col]])
    if (anyNA(smp)) {
      warning(sum(is.na(smp)), " cells have a missing ", sample_col,
              " value; they are pooled into their own stratum.",
              call. = FALSE)
      smp[is.na(smp)] <- "__NA__"
    }
  }

  size <- rep(NA_real_, ncol(sce))
  med  <- rep(NA_real_, ncol(sce))
  samples <- unique(smp[has])
  median_by_sample <- stats::setNames(rep(NA_real_, length(samples)), samples)

  for (s in samples) {
    in_s <- has & smp == s
    tab  <- table(cid[in_s])
    sz   <- as.numeric(tab[cid[in_s]])
    size[in_s] <- sz
    m <- if (median_basis == "clone") {
      stats::median(as.numeric(tab))   # median over unique clones
    } else {
      stats::median(sz)                # median over cells
    }
    med[in_s] <- m
    median_by_sample[[s]] <- m
  }

  expansion <- ifelse(is.na(size), NA_character_,
                      ifelse(size > 1.5 * med, "expanded", "non-expanded"))

  if (verbose) {
    message("Derived expansion labels from '", clone_col, "'",
            if (is.null(sample_col)) "" else paste0(" within '", sample_col, "'"),
            " (median_basis = '", median_basis, "'): ",
            sum(has), "/", ncol(sce), " cells labelled, ",
            sum(expansion == "expanded", na.rm = TRUE), " expanded.")
  }

  list(clone_id_size = size,
       median_clone_size = med,
       expansion = expansion,
       clone_col = clone_col,
       sample_col = sample_col,
       median_basis = median_basis,
       median_by_sample = median_by_sample,
       n_labelled = sum(has),
       n_clones = length(unique(cid[has])),
       n_samples = length(samples))
}

#' @keywords internal
# Rank-based AUROC (Mann-Whitney U). Computed in R so the reported number
# never depends on how scXpand's evaluator handles NAs or ties.
.aurocR <- function(prob, label) {
  ok <- !is.na(prob) & !is.na(label)
  prob <- prob[ok]
  label <- as.logical(label[ok])
  n1 <- sum(label)
  n0 <- length(label) - n1
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(prob)
  (sum(r[label]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# ===========================================================================
# H5AD staging
# ===========================================================================

#' @keywords internal
.validateCountsAssay <- function(X, assay) {
  normalized <- c("logcounts", "lognorm", "lognormcounts", "data",
                  "normcounts", "scaledata", "scale.data")
  if (tolower(assay) %in% tolower(normalized)) {
    stop("assay '", assay, "' looks like normalized or log-transformed data. ",
         "scXpand requires raw UMI counts; pass assay = 'counts'.",
         call. = FALSE)
  }

  v <- if (methods::is(X, "dgCMatrix")) X@x else as.numeric(X)
  if (!length(v) || all(v == 0)) {
    stop("assay '", assay, "' contains no non-zero counts.", call. = FALSE)
  }
  if (min(v) < 0) {
    stop("assay '", assay, "' contains negative values; scXpand requires ",
         "raw UMI counts.", call. = FALSE)
  }
  # Integrality on a subsample: full checks on millions of cells are wasteful
  # and a single non-integer is enough to reject.
  if (length(v) > 1e5L) v <- v[seq_len(1e5L)]
  if (any(abs(v - round(v)) > 1e-8)) {
    stop("assay '", assay, "' contains non-integer values (it looks ",
         "normalized or log-transformed). scXpand requires raw UMI counts; ",
         "pass assay = 'counts'.", call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
# Sum rows sharing an Ensembl ID. Raw UMI counts are additive, so summing
# preserves the count semantics scXpand expects.
.collapseRowsSum <- function(X, groups) {
  groups <- droplevels(as.factor(groups))
  if (requireNamespace("Matrix", quietly = TRUE)) {
    g <- as.integer(groups)
    agg <- Matrix::sparseMatrix(i = g, j = seq_along(g), x = 1,
                                dims = c(nlevels(groups), length(g)))
    out <- agg %*% X
    rownames(out) <- levels(groups)
    colnames(out) <- colnames(X)
    return(out)
  }
  out <- rowsum(as.matrix(X), group = groups, reorder = TRUE)
  rownames(out) <- levels(groups)
  out
}

#' @keywords internal
# Build a minimal SCE holding only what scXpand reads (X and obs) and write
# it to H5AD. zellkonverter runs in its own basilisk env, so this happens on
# the R side and only the file path crosses into scXpandEnv.
.stageScXpandH5AD <- function(sce, res_ids, labels, obs_columns, assay,
                              h5ad_path, verbose = TRUE) {

  X <- SummarizedExperiment::assay(sce, assay)
  .validateCountsAssay(X, assay)

  X <- X[res_ids$keep, , drop = FALSE]
  if (!is.null(res_ids$collapse_groups)) {
    X <- .collapseRowsSum(X, res_ids$collapse_groups)
  } else {
    rownames(X) <- res_ids$ids
  }
  colnames(X) <- colnames(sce)

  stage <- SingleCellExperiment::SingleCellExperiment(
    assays  = list(counts = X),
    colData = SummarizedExperiment::colData(sce))

  # drop_ct = TRUE: the CT* columns are long strings scXpand never reads,
  # and their content is already distilled into clone_id_size.
  stage <- .pruneColData(stage, obs_columns = obs_columns, drop_ct = TRUE)

  # Attached after pruning so obs_columns can never remove them. Names come
  # from scxpand/data_util/constants.py. All three are written together or
  # not at all: scXpand reads them as a set during evaluation, and a partial
  # set is a plausible way to make its evaluator raise.
  if (!is.null(labels) && !is.null(labels$expansion)) {
    SummarizedExperiment::colData(stage)$clone_id_size <-
      labels$clone_id_size
    SummarizedExperiment::colData(stage)$median_clone_size <-
      labels$median_clone_size
    SummarizedExperiment::colData(stage)$expansion <- labels$expansion
  }

  if (verbose) message("Writing staged H5AD: ", h5ad_path)
  zellkonverter::writeH5AD(stage, file = h5ad_path, X_name = "counts")

  invisible(h5ad_path)
}

# ===========================================================================
# Result handling
# ===========================================================================

#' @keywords internal
# Match predictions back to R cells by barcode rather than trusting
# positional alignment through writeH5AD -> h5py -> AnnData -> DataLoader.
.alignPredictions <- function(preds, obs_names, expected_names) {
  if (length(preds) != length(obs_names)) {
    stop("scXpand returned ", length(preds), " predictions for ",
         length(obs_names), " cells in the staged H5AD.", call. = FALSE)
  }
  if (length(preds) != length(expected_names) ||
      !setequal(obs_names, expected_names)) {
    warning("H5AD obs_names do not match the staged cell names; falling ",
            "back to positional alignment.", call. = FALSE)
    names(preds) <- expected_names[seq_along(preds)]
    return(preds)
  }
  names(preds) <- obs_names
  preds[expected_names]
}

#' @keywords internal
# Flatten scXpand's nested metrics dict to a named numeric vector with
# dot-joined keys. Non-scalar leaves (arrays, confusion matrices) are
# dropped here; the unflattened list is retained as metrics_raw.
.flattenMetrics <- function(x, prefix = "") {
  if (is.null(x) || !length(x)) return(stats::setNames(numeric(0), character(0)))
  nms <- names(x)
  if (is.null(nms)) return(stats::setNames(numeric(0), character(0)))
  out <- stats::setNames(numeric(0), character(0))
  for (i in seq_along(x)) {
    nm <- nms[i]
    if (is.null(nm) || !nzchar(nm)) next
    v <- x[[i]]
    key <- if (nzchar(prefix)) paste(prefix, nm, sep = ".") else nm
    if (is.list(v)) {
      out <- c(out, .flattenMetrics(v, key))
    } else if (length(v) == 1L && (is.numeric(v) || is.logical(v)) &&
               !is.na(v) && is.finite(as.numeric(v))) {
      out <- c(out, stats::setNames(as.numeric(v), key))
    }
  }
  out
}

#' @keywords internal
# Fill a per-cell column, leaving NA for cells that were not scored.
.writeCellColumn <- function(obj, col_name, values, cell_names) {
  col_vec <- rep(NA, ncol(obj))
  names(col_vec) <- colnames(obj)
  col_vec[cell_names] <- values
  names(col_vec) <- NULL
  if (methods::is(obj, "SingleCellExperiment")) {
    SummarizedExperiment::colData(obj)[[col_name]] <- col_vec
  } else {
    obj[[col_name]] <- col_vec
  }
  obj
}

#' @keywords internal
# Stash the run summary. SingleCellExperiment has metadata(); Seurat keeps
# the equivalent in the misc slot.
.writeObjMetadata <- function(obj, key, value) {
  if (methods::is(obj, "SingleCellExperiment")) {
    md <- S4Vectors::metadata(obj)
    md[[key]] <- value
    S4Vectors::metadata(obj) <- md
  } else {
    misc <- methods::slot(obj, "misc")
    if (!is.list(misc)) misc <- list()
    misc[[key]] <- value
    methods::slot(obj, "misc") <- misc
  }
  obj
}

# ===========================================================================
# Exported functions
# ===========================================================================

#' List available scXpand pretrained models
#'
#' @description
#' Returns the pan-cancer models published with scXpand. By default this
#' reads a static table shipped with immLynx so it never triggers a
#' multi-gigabyte Python environment build. Set \code{refresh = TRUE} to
#' query the installed \code{scxpand} package, which is authoritative if
#' upstream has added models.
#'
#' @param refresh Logical. Query the installed \code{scxpand} package rather
#'   than the static table. Requires the scXpand basilisk environment and
#'   will build it on first use. Default \code{FALSE}.
#'
#' @return A \code{data.frame} with columns \code{model_name},
#'   \code{model_type}, \code{version} and \code{description}. When
#'   \code{refresh = TRUE}, \code{model_type} and \code{description} may be
#'   \code{NA} for models absent from the static table.
#'
#' @export
#'
#' @examples
#' listScXpandModels()
#'
#' \donttest{
#' # Authoritative, but builds the scXpand environment on first call.
#' listScXpandModels(refresh = TRUE)
#' }
listScXpandModels <- function(refresh = FALSE) {
  if (!isTRUE(refresh)) return(.SCXPAND_MODELS)

  live <- .run_in_scxpand_env(function() {
    sx <- reticulate::import("scxpand", convert = FALSE)
    bi <- reticulate::import_builtins(convert = FALSE)
    reg <- sx$PRETRAINED_MODELS
    nms <- as.character(reticulate::py_to_r(bi$list(reg$keys())))
    vers <- vapply(nms, function(n) tryCatch(
      as.character(reticulate::py_to_r(
        reticulate::py_get_attr(reg[[n]], "version"))),
      error = function(e) NA_character_), character(1))
    data.frame(model_name = nms, version = unname(vers),
               stringsAsFactors = FALSE)
  })

  idx <- match(live$model_name, .SCXPAND_MODELS$model_name)
  data.frame(model_name = live$model_name,
             model_type = .SCXPAND_MODELS$model_type[idx],
             version = live$version,
             description = .SCXPAND_MODELS$description[idx],
             stringsAsFactors = FALSE)
}

#' Predict T-cell clonal expansion from gene expression with scXpand
#'
#' @description
#' Runs inference with a scXpand pretrained pan-cancer model to estimate,
#' for each cell, the probability that it belongs to an expanded T-cell
#' clone. scXpand uses gene expression only, so this works on datasets with
#' no paired TCR sequencing. When the object does carry scRepertoire clone
#' calls, \code{runScXpand} also derives ground-truth expansion labels so
#' the prediction can be benchmarked against the observed repertoire.
#'
#' @details
#' \strong{Input requirements.} scXpand's pretrained models were trained on
#' raw UMI counts from T cells indexed by Ensembl gene ID. All three matter:
#' \itemize{
#'   \item \emph{Raw counts.} Normalized or log-transformed values are
#'     rejected. Use \code{assay = "counts"}.
#'   \item \emph{Ensembl IDs.} Most Seurat and scRepertoire objects carry
#'     gene symbols; see \code{gene_ids} and \code{map_symbols}.
#'   \item \emph{T cells only.} Filter to T cells before calling. There is
#'     no reliable programmatic check, and predictions on non-T cells are
#'     meaningless.
#' }
#' Genes the model expects but cannot find are zero-filled by scXpand, and
#' extra genes are ignored, so a partial overlap degrades quietly rather
#' than erroring. That is why the gene mapping summary is worth reading.
#'
#' \strong{Symbol mapping.} When identifiers are symbols, they are mapped
#' through \code{org.Hs.eg.db} (SYMBOL first, then ALIAS). One symbol often
#' maps to several Ensembl IDs because of alternate haplotypes:
#' \code{HLA-DRA} maps to eight. \code{multi_map = "first"} takes the
#' lexicographically first, which is reproducible but may pick a scaffold
#' the model does not know; \code{multi_map = "expand"} emits every
#' candidate, which maximizes overlap with the model's panel at the cost of
#' a larger staged matrix. Several input rows can also land on one Ensembl
#' ID, in which case \code{collapse} decides what happens (summing raw
#' counts is the default and preserves UMI semantics).
#'
#' \strong{Expansion labels.} scXpand calls a cell expanded when its clone
#' size exceeds 1.5 times the median clone size for its sample. Clone sizes
#' are tabulated fresh from \code{clone_col} within \code{sample_col};
#' scRepertoire's \code{clonalFrequency} column is deliberately not used,
#' because \code{combineExpression()} computes it under whatever
#' \code{group.by} was in effect and it is not necessarily a per-sample
#' count. "Median clone size" is genuinely ambiguous, so
#' \code{median_basis} selects the reading: \code{"clone"} takes the median
#' over unique clones, which in typical 10x data is 1 (most clones are
#' singletons) and therefore reduces the rule to \emph{clone size at least
#' 2}; \code{"cell"} takes the median over cells, which is dominated by
#' large clones and is much stricter.
#'
#' \strong{Cost.} The first call builds a dedicated basilisk environment
#' (Python 3.11, PyTorch CPU, scxpand) that occupies several gigabytes, and
#' downloads the selected model (roughly 35 MB) from figshare into
#' \code{cache_dir}. Both are reused afterwards. If a download does fail,
#' clear \code{cache_dir} before retrying: the downloader caches the failed
#' response and will otherwise keep reusing it.
#'
#' @param input A \code{SingleCellExperiment} or \code{Seurat} object.
#'   SingleCellExperiment is the native format and is used without
#'   conversion; Seurat objects are converted in and returned as Seurat.
#' @param model_name Pretrained model to use. See
#'   \code{\link{listScXpandModels}}. Default
#'   \code{"pan_cancer_autoencoder"}.
#' @param assay Assay holding raw UMI counts. Default \code{"counts"}.
#' @param gene_ids Where to find gene identifiers: \code{"rownames"}
#'   (default), the name of a \code{rowData} column, or a character vector
#'   of length \code{nrow(input)}.
#' @param map_symbols One of \code{"auto"} (map symbols to Ensembl IDs only
#'   if too few identifiers already look like Ensembl IDs), \code{"always"},
#'   or \code{"never"}. Mapping requires \code{org.Hs.eg.db} and
#'   \code{AnnotationDbi}.
#' @param ensembl_min_frac Minimum fraction of identifiers that must resolve
#'   to Ensembl IDs. Default \code{0.5}. Lower it for targeted panels.
#' @param multi_map How to resolve one symbol mapping to several Ensembl
#'   IDs: \code{"first"} (default, lexicographically first),
#'   \code{"expand"} (emit all), or \code{"drop"}.
#' @param collapse How to resolve several rows mapping to one Ensembl ID:
#'   \code{"sum"} (default), \code{"first"}, or \code{"drop"}.
#' @param derive_labels Logical. Derive ground-truth expansion labels from
#'   clone data when present. Default \code{TRUE}.
#' @param clone_col Column holding clone identity. Default \code{NULL}:
#'   use \code{"CTstrict"} if present, else \code{"CTaa"}.
#' @param sample_col Column identifying the sample or patient that clone
#'   sizes are tabulated within. Default \code{NULL}, which pools all cells
#'   into one sample and warns.
#' @param median_basis \code{"clone"} (default) or \code{"cell"}. See
#'   Details.
#' @param label_cells \code{"clonal"} (default) scores only cells with a
#'   clone call, which guarantees the ground-truth column has no missing
#'   values; \code{"all"} scores every cell.
#' @param threshold Probability cutoff for the binary call. Default
#'   \code{0.5}.
#' @param obs_columns Character vector of \code{colData} columns to carry
#'   into the staged H5AD. Default \code{NULL}: all of them, minus
#'   list-columns and the scRepertoire \code{CT*} columns.
#' @param batch_size Inference batch size. Default \code{1024}.
#' @param num_workers DataLoader workers. Default \code{0}, because worker
#'   processes forking out of embedded Python inside a basilisk child can
#'   hang on macOS. Raise it on Linux for speed.
#' @param work_dir Directory for the staged H5AD and scXpand's own output.
#'   Default \code{NULL}: a temporary directory, removed on exit.
#' @param cache_dir Directory for downloaded pretrained models. Default
#'   \code{NULL}: a \code{scxpand} subdirectory of
#'   \code{tools::R_user_dir("immLynx", "cache")}. Passing this explicitly
#'   matters because scXpand's own default would write a
#'   \code{.scxpand_cache} directory into the current working directory.
#' @param keep_files Logical. Keep the staged H5AD and scXpand outputs.
#'   Default \code{FALSE}.
#' @param overwrite Logical. Replace an existing \code{adata.h5ad} in an
#'   explicitly supplied \code{work_dir}. Default \code{FALSE}.
#' @param column_prefix Prefix for the columns written back. Default
#'   \code{"scXpand"}.
#' @param return_object Logical. Return the input object with new columns
#'   (default), or a \code{data.frame} of per-cell results.
#' @param verbose Logical. Emit progress messages. Default \code{TRUE}.
#'
#' @return If \code{return_object = TRUE}, the input object (same class)
#'   with \code{<prefix>_expansion_prob} and \code{<prefix>_expansion_pred}
#'   added to cell metadata, plus \code{<prefix>_clone_id_size},
#'   \code{<prefix>_median_clone_size} and \code{<prefix>_expansion_truth}
#'   when labels were derived. A run summary is stored under
#'   \code{metadata(x)$scXpand} for SingleCellExperiment or
#'   \code{x@misc$scXpand} for Seurat, holding the model info, the gene
#'   mapping report, scXpand's metrics and an independently computed AUROC.
#'   If \code{return_object = FALSE}, a \code{data.frame} with one row per
#'   scored cell and the run summary attached as
#'   \code{attr(x, "scXpand")}.
#'
#' @export
#' @importFrom methods is slot slot<-
#' @importFrom SummarizedExperiment colData colData<- rowData assay assayNames
#' @importFrom stats median setNames
#' @importFrom utils head
#'
#' @examples
#' data(immLynx_example)
#'
#' # The label derivation is pure R and runs without Python. Note that the
#' # derived clone sizes are tabulated per sample and so need not match the
#' # clonalFrequency column scRepertoire wrote.
#' labs <- immLynx:::.deriveExpansionLabels(
#'   immLynx_example,
#'   sample_col = "Patient"
#' )
#' table(labs$expansion, useNA = "ifany")
#'
#' \donttest{
#' # Full inference. The example object carries gene symbols and is not
#' # filtered to T cells, so this demonstrates the mechanics rather than a
#' # meaningful biological result.
#' sce <- runScXpand(
#'   immLynx_example,
#'   sample_col = "Patient",
#'   assay      = "counts"
#' )
#' summary(SummarizedExperiment::colData(sce)$scXpand_expansion_prob)
#' S4Vectors::metadata(sce)$scXpand$auroc_immLynx
#' }
runScXpand <- function(input,
                       model_name       = "pan_cancer_autoencoder",
                       assay            = "counts",
                       gene_ids         = "rownames",
                       map_symbols      = c("auto", "always", "never"),
                       ensembl_min_frac = 0.5,
                       multi_map        = c("first", "expand", "drop"),
                       collapse         = c("sum", "first", "drop"),
                       derive_labels    = TRUE,
                       clone_col        = NULL,
                       sample_col       = NULL,
                       median_basis     = c("clone", "cell"),
                       label_cells      = c("clonal", "all"),
                       threshold        = 0.5,
                       obs_columns      = NULL,
                       batch_size       = 1024L,
                       num_workers      = 0L,
                       work_dir         = NULL,
                       cache_dir        = NULL,
                       keep_files       = FALSE,
                       overwrite        = FALSE,
                       column_prefix    = "scXpand",
                       return_object    = TRUE,
                       verbose          = TRUE) {

  map_symbols  <- match.arg(map_symbols)
  multi_map    <- match.arg(multi_map)
  collapse     <- match.arg(collapse)
  median_basis <- match.arg(median_basis)
  label_cells  <- match.arg(label_cells)

  # --- Input dispatch ---------------------------------------------------
  # SingleCellExperiment passes through untouched; Seurat is converted in
  # here and the results are written back onto the original object below.
  sce <- .coerceToSCE(input)

  if (ncol(sce) == 0L) {
    stop("input has no cells.", call. = FALSE)
  }
  if (is.null(colnames(sce))) {
    stop("input must have cell names (colnames); predictions are matched ",
         "back by barcode.", call. = FALSE)
  }
  if (anyDuplicated(colnames(sce))) {
    stop("input has duplicated cell names; predictions cannot be matched ",
         "back by barcode. Make colnames unique first.", call. = FALSE)
  }

  # --- Scalar argument validation ---------------------------------------
  if (!is.character(model_name) || length(model_name) != 1L ||
      is.na(model_name) || !nzchar(model_name)) {
    stop("model_name must be a single non-empty string.", call. = FALSE)
  }
  if (!(model_name %in% .SCXPAND_MODELS$model_name)) {
    warning("Unknown model_name '", model_name, "'; known models: ",
            paste(.SCXPAND_MODELS$model_name, collapse = ", "),
            ". Passing it through to scXpand anyway.", call. = FALSE)
  }
  if (!is.character(assay) || length(assay) != 1L || is.na(assay)) {
    stop("assay must be a single string.", call. = FALSE)
  }
  if (!(assay %in% SummarizedExperiment::assayNames(sce))) {
    stop("assay '", assay, "' not found. Available: ",
         paste(SummarizedExperiment::assayNames(sce), collapse = ", "), ".",
         call. = FALSE)
  }
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold) ||
      threshold < 0 || threshold > 1) {
    stop("threshold must be a single number in [0, 1].", call. = FALSE)
  }
  if (!is.numeric(ensembl_min_frac) || length(ensembl_min_frac) != 1L ||
      is.na(ensembl_min_frac) || ensembl_min_frac < 0 ||
      ensembl_min_frac > 1) {
    stop("ensembl_min_frac must be a single number in [0, 1].", call. = FALSE)
  }
  .checkCount <- function(x, nm) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < 0 ||
        abs(x - round(x)) > 1e-8) {
      stop(nm, " must be a single non-negative whole number.", call. = FALSE)
    }
    as.integer(round(x))
  }
  batch_size  <- .checkCount(batch_size, "batch_size")
  num_workers <- .checkCount(num_workers, "num_workers")
  if (batch_size < 1L) {
    stop("batch_size must be at least 1.", call. = FALSE)
  }
  if (!is.character(column_prefix) || length(column_prefix) != 1L ||
      is.na(column_prefix)) {
    stop("column_prefix must be a single string.", call. = FALSE)
  }

  cd_names <- colnames(SummarizedExperiment::colData(sce))
  if (!is.null(clone_col) &&
      (!is.character(clone_col) || length(clone_col) != 1L ||
       !(clone_col %in% cd_names))) {
    stop("clone_col '", paste(clone_col, collapse = ", "),
         "' not found in colData.", call. = FALSE)
  }
  if (!is.null(sample_col) &&
      (!is.character(sample_col) || length(sample_col) != 1L ||
       !(sample_col %in% cd_names))) {
    stop("sample_col '", paste(sample_col, collapse = ", "),
         "' not found in colData.", call. = FALSE)
  }

  # --- Work directory ---------------------------------------------------
  transient <- is.null(work_dir)
  if (transient) {
    work_dir <- tempfile("scxpand_")
  } else if (!is.character(work_dir) || length(work_dir) != 1L ||
             is.na(work_dir) || !nzchar(work_dir)) {
    stop("work_dir must be a single non-empty path string.", call. = FALSE)
  }
  h5ad_path <- file.path(work_dir, "adata.h5ad")
  out_dir   <- file.path(work_dir, "scxpand_out")

  if (file.exists(h5ad_path)) {
    if (!isTRUE(overwrite)) {
      stop("Output exists (use overwrite = TRUE to replace): ", h5ad_path,
           call. = FALSE)
    }
    unlink(h5ad_path)
  }
  if (!dir.exists(work_dir)) {
    ok <- dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
    if (!ok || !dir.exists(work_dir)) {
      stop("Failed to create work_dir: ", work_dir, call. = FALSE)
    }
  }
  if (transient && !isTRUE(keep_files)) {
    on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)
  }

  # --- Model cache ------------------------------------------------------
  # scXpand would otherwise write ".scxpand_cache" into the current working
  # directory. Keep downloaded models in the standard per-user R cache so
  # they persist across sessions without touching the user's project.
  if (is.null(cache_dir)) {
    cache_dir <- file.path(tools::R_user_dir("immLynx", which = "cache"),
                           "scxpand")
  } else if (!is.character(cache_dir) || length(cache_dir) != 1L ||
             is.na(cache_dir) || !nzchar(cache_dir)) {
    stop("cache_dir must be a single non-empty path string.", call. = FALSE)
  }
  if (!dir.exists(cache_dir)) {
    ok <- dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    if (!ok || !dir.exists(cache_dir)) {
      stop("Failed to create cache_dir: ", cache_dir, call. = FALSE)
    }
  }
  cache_dir <- normalizePath(cache_dir, mustWork = TRUE)

  # --- Gene identifiers -------------------------------------------------
  res_ids <- .resolveGeneIDs(sce, gene_ids = gene_ids,
                             map_symbols = map_symbols,
                             ensembl_min_frac = ensembl_min_frac,
                             multi_map = multi_map, collapse = collapse,
                             verbose = verbose)

  # --- Ground-truth labels ----------------------------------------------
  labels <- NULL
  if (isTRUE(derive_labels)) {
    labels <- .deriveExpansionLabels(sce, clone_col = clone_col,
                                     sample_col = sample_col,
                                     median_basis = median_basis,
                                     verbose = verbose)
    if (is.null(labels)) {
      warning("No clone data found in colData (looked for CTstrict, CTaa). ",
              "Running inference without ground-truth labels; AUROC will be ",
              "unavailable.", call. = FALSE)
    }
  }

  # --- Cell subset ------------------------------------------------------
  scored <- colnames(sce)
  if (!is.null(labels) && label_cells == "clonal") {
    has <- !is.na(labels$expansion)
    if (!any(has)) {
      stop("label_cells = 'clonal' but no cell has a clone call.",
           call. = FALSE)
    }
    if (!all(has) && verbose) {
      message("Scoring ", sum(has), "/", ncol(sce),
              " cells with clone calls (label_cells = 'clonal').")
    }
    sce <- sce[, has, drop = FALSE]
    labels$clone_id_size     <- labels$clone_id_size[has]
    labels$median_clone_size <- labels$median_clone_size[has]
    labels$expansion         <- labels$expansion[has]
    scored <- colnames(sce)
  }

  # A single-class ground-truth column makes scXpand's evaluator raise
  # inside roc_auc_score, so drop it rather than lose the whole run.
  if (!is.null(labels)) {
    classes <- unique(labels$expansion[!is.na(labels$expansion)])
    if (length(classes) < 2L) {
      warning("Derived expansion labels have only one class ('",
              paste(classes, collapse = ", "),
              "'); omitting the ground-truth column so scXpand's evaluation ",
              "does not fail. Predictions are unaffected.", call. = FALSE)
      labels$expansion_dropped <- labels$expansion
      labels$expansion <- NULL
    }
  }

  # --- Stage the H5AD ---------------------------------------------------
  .stageScXpandH5AD(sce, res_ids = res_ids, labels = labels,
                    obs_columns = obs_columns, assay = assay,
                    h5ad_path = h5ad_path, verbose = verbose)

  # --- Inference --------------------------------------------------------
  if (verbose) {
    message("Running scXpand inference with '", model_name, "' on ",
            length(scored), " cells x ", length(res_ids$ids), " features...")
  }

  py_out <- tryCatch(
    .run_in_scxpand_env(function(h5ad, model, save_dir, cache, bsz, nw) {
      # convert = FALSE throughout: run_inference() returns an
      # InferenceResults dataclass whose fields are a numpy array, a nested
      # dict and a custom ModelInfo object. Auto-converting the container
      # yields a half-converted structure, so we hold a proxy and convert
      # each field explicitly, guarding them individually against upstream
      # field renames.
      sx <- reticulate::import("scxpand", convert = FALSE)
      ad <- reticulate::import("anndata", convert = FALSE)

      # Download as an explicit step rather than letting run_inference do
      # it, for two reasons.
      #
      # First, run_inference() calls download_pretrained_model() without a
      # cache_dir, and that defaults to ".scxpand_cache" in the *current
      # working directory* -- which would drop a few hundred megabytes into
      # whatever project the user happens to be sitting in.
      #
      # Second, scXpand's registry points at
      # https://figshare.com/ndownloader/articles/..., which answers HTTP
      # 202 with an empty body; pooch writes that empty file out and the
      # unzip then fails with "File is not a zip file". The same archive on
      # the ndownloader.figshare.com host returns 200 and a valid zip, so
      # we take the registry URL and rewrite the host. The substitution is
      # a no-op once upstream fixes its URLs.
      url <- as.character(reticulate::py_to_r(reticulate::py_get_attr(
        sx$get_pretrained_model_info(model), "url")))
      url <- sub("^https://figshare\\.com/ndownloader/",
                 "https://ndownloader.figshare.com/", url)

      model_path <- as.character(reticulate::py_to_r(
        sx$download_pretrained_model(model_url = url, cache_dir = cache)))

      res <- sx$run_inference(data_path = h5ad,
                              model_path = model_path,
                              save_path = save_dir,
                              batch_size = as.integer(bsz),
                              num_workers = as.integer(nw))

      preds <- as.numeric(reticulate::py_to_r(res$predictions))

      metrics <- tryCatch({
        if (isTRUE(reticulate::py_to_r(res$has_metrics))) {
          reticulate::py_to_r(res$metrics)
        } else NULL
      }, error = function(e) NULL)

      version <- tryCatch(
        as.character(reticulate::py_to_r(
          reticulate::py_get_attr(sx, "__version__"))),
        error = function(e) NA_character_)

      # Read obs_names back so R can match by barcode rather than trusting
      # positional alignment through the file round-trip.
      obs_names <- as.character(reticulate::py_to_r(
        ad$read_h5ad(h5ad, backed = "r")$obs_names$to_list()))

      list(predictions = preds, metrics = metrics, model_path = model_path,
           version = version, obs_names = obs_names)
    },
    h5ad = h5ad_path, model = model_name, save_dir = out_dir,
    cache = cache_dir, bsz = batch_size, nw = num_workers),
    error = function(e) e
  )

  if (inherits(py_out, "error")) {
    msg <- conditionMessage(py_out)
    hint <- if (grepl("not a zip file|Failed to download", msg)) paste0(
      "\nThis is a model download failure, not a problem with your data. ",
      "The failed response is cached, so clear it before retrying:\n  ",
      "unlink(\"", cache_dir, "\", recursive = TRUE)") else ""

    stop("scXpand inference failed: ", msg, "\n",
         "Staged H5AD: ",
         if (isTRUE(keep_files) || !transient) h5ad_path else
           "(deleted; re-run with keep_files = TRUE to inspect)", "\n",
         "Model cache: ", cache_dir, hint, call. = FALSE)
  }

  preds <- .alignPredictions(py_out$predictions, py_out$obs_names, scored)

  # --- Assemble the run summary ------------------------------------------
  truth <- if (!is.null(labels)) {
    if (is.null(labels$expansion)) labels$expansion_dropped else labels$expansion
  } else NULL

  auroc_immLynx <- if (!is.null(truth)) {
    .aurocR(as.numeric(preds), truth == "expanded")
  } else NA_real_

  metrics_flat <- .flattenMetrics(py_out$metrics)
  auroc_scxpand <- {
    hit <- grep("(^|\\.)AUROC$", names(metrics_flat), ignore.case = TRUE)
    if (length(hit)) unname(metrics_flat[[hit[1]]]) else NA_real_
  }

  # run_inference() only populates model_info for its registry branch, and
  # we deliberately take the local-path branch to control the cache. Rebuild
  # the same fields from the registry table plus the resolved path.
  reg <- .SCXPAND_MODELS[.SCXPAND_MODELS$model_name == model_name, ,
                         drop = FALSE]
  model_info <- list(
    model_name = model_name,
    model_type = if (nrow(reg)) reg$model_type[1] else NA_character_,
    version    = if (nrow(reg)) reg$version[1] else NA_character_,
    source     = "registry",
    path       = py_out$model_path)

  summary_list <- list(
    model = model_name,
    model_info = model_info,
    cache_dir = cache_dir,
    scxpand_version = py_out$version,
    n_cells = length(preds),
    n_features = length(res_ids$ids),
    threshold = threshold,
    metrics = metrics_flat,
    metrics_raw = py_out$metrics,
    auroc = auroc_scxpand,
    auroc_immLynx = auroc_immLynx,
    labels = if (is.null(labels)) NULL else list(
      clone_col = labels$clone_col,
      sample_col = labels$sample_col,
      median_basis = labels$median_basis,
      median_by_sample = labels$median_by_sample,
      n_labelled = labels$n_labelled,
      n_clones = labels$n_clones,
      n_samples = labels$n_samples),
    gene_mapping = res_ids$report,
    work_dir = if (isTRUE(keep_files) || !transient) work_dir else NULL,
    timestamp = Sys.time())

  if (verbose && !is.na(auroc_immLynx)) {
    message("AUROC against derived expansion labels: ",
            sprintf("%.3f", auroc_immLynx))
  }

  pred_call <- ifelse(is.na(preds), NA_character_,
                      ifelse(preds >= threshold, "expanded", "non-expanded"))

  # --- Return ------------------------------------------------------------
  if (isTRUE(return_object)) {
    p <- function(x) if (nzchar(column_prefix))
      paste0(column_prefix, "_", x) else x

    input <- .writeCellColumn(input, p("expansion_prob"),
                              as.numeric(preds), scored)
    input <- .writeCellColumn(input, p("expansion_pred"), pred_call, scored)
    if (!is.null(labels)) {
      input <- .writeCellColumn(input, p("clone_id_size"),
                                labels$clone_id_size, scored)
      input <- .writeCellColumn(input, p("median_clone_size"),
                                labels$median_clone_size, scored)
      input <- .writeCellColumn(input, p("expansion_truth"), truth, scored)
    }
    input <- .writeObjMetadata(input, "scXpand", summary_list)

    if (verbose) {
      message("Added ", p("expansion_prob"), " and ", p("expansion_pred"),
              " to cell metadata.")
    }
    return(input)
  }

  df <- data.frame(barcode = scored,
                   expansion_prob = as.numeric(preds),
                   expansion_pred = pred_call,
                   stringsAsFactors = FALSE)
  if (!is.null(labels)) {
    df$clone_id_size <- labels$clone_id_size
    df$median_clone_size <- labels$median_clone_size
    df$expansion_truth <- truth
  }
  rownames(df) <- NULL
  attr(df, "scXpand") <- summary_list
  df
}
