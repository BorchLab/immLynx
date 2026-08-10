library(testthat)

# ===========================================================================
# Signatures and the static model table
# ===========================================================================

test_that("runScXpand signature has correct defaults", {
  f <- formals(runScXpand)
  expect_equal(f$model_name, "pan_cancer_autoencoder")
  expect_equal(f$assay, "counts")
  expect_equal(f$gene_ids, "rownames")
  expect_equal(eval(f$map_symbols), c("auto", "always", "never"))
  expect_equal(f$ensembl_min_frac, 0.5)
  expect_equal(eval(f$multi_map), c("first", "expand", "drop"))
  expect_equal(eval(f$collapse), c("sum", "first", "drop"))
  expect_true(f$derive_labels)
  expect_null(eval(f$clone_col))
  expect_null(eval(f$sample_col))
  expect_equal(eval(f$median_basis), c("clone", "cell"))
  expect_equal(eval(f$label_cells), c("clonal", "all"))
  expect_equal(f$threshold, 0.5)
  expect_equal(f$batch_size, 1024L)
  expect_equal(f$num_workers, 0L)
  expect_null(eval(f$work_dir))
  expect_null(eval(f$cache_dir))
  expect_false(f$keep_files)
  expect_false(f$overwrite)
  expect_equal(f$column_prefix, "scXpand")
  expect_true(f$return_object)
  expect_true(f$verbose)
})

test_that("listScXpandModels returns the registry without touching Python", {
  m <- listScXpandModels()
  expect_s3_class(m, "data.frame")
  expect_equal(nrow(m), 5L)
  expect_equal(colnames(m),
               c("model_name", "model_type", "version", "description"))
  expect_true("pan_cancer_autoencoder" %in% m$model_name)
  expect_setequal(m$model_type,
                  c("autoencoder", "mlp", "lightgbm", "logistic", "svm"))
  expect_false(formals(listScXpandModels)$refresh)
})

# ===========================================================================
# Ensembl identifier helpers
# ===========================================================================

test_that(".isEnsembl matches only well-formed Ensembl gene IDs", {
  expect_true(immLynx:::.isEnsembl("ENSG00000198851"))
  expect_false(immLynx:::.isEnsembl("ENSG0000019885"))    # too short
  expect_false(immLynx:::.isEnsembl("ENST00000198851"))   # transcript
  expect_false(immLynx:::.isEnsembl("CD3E"))
  expect_false(immLynx:::.isEnsembl(NA_character_))
})

test_that(".stripEnsemblVersion strips only Ensembl versions", {
  expect_equal(immLynx:::.stripEnsemblVersion("ENSG00000198851.12"),
               "ENSG00000198851")
  # A blanket sub() would mangle these.
  expect_equal(immLynx:::.stripEnsemblVersion("MARCH1.2"), "MARCH1.2")
  expect_equal(immLynx:::.stripEnsemblVersion("7SK.2"), "7SK.2")
  expect_equal(immLynx:::.stripEnsemblVersion(NA_character_), NA_character_)
})

# ===========================================================================
# .resolveGeneIDs
# ===========================================================================

.gene_sce <- function(rn, id_col = NULL) {
  n <- length(rn)
  # Values must vary and dimnames must be set up front: a square matrix of
  # constant values makes Matrix() return a symmetric dsCMatrix, which ties
  # rownames to colnames and silently discards the gene identifiers.
  counts <- Matrix::Matrix(
    matrix(seq_len(n * 4L), nrow = n,
           dimnames = list(rn, paste0("c", 1:4))),
    sparse = TRUE)
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts))
  if (!is.null(id_col)) SummarizedExperiment::rowData(sce)$ID <- id_col
  sce
}

test_that(".resolveGeneIDs accepts all three gene_ids forms", {
  ens <- sprintf("ENSG%011d", 1:3)

  r1 <- immLynx:::.resolveGeneIDs(.gene_sce(ens), verbose = FALSE)
  expect_equal(r1$ids, ens)
  expect_equal(r1$keep, 1:3)
  expect_null(r1$collapse_groups)

  r2 <- immLynx:::.resolveGeneIDs(.gene_sce(c("a", "b", "c"), ens),
                                  gene_ids = "ID", verbose = FALSE)
  expect_equal(r2$ids, ens)

  r3 <- immLynx:::.resolveGeneIDs(.gene_sce(c("a", "b", "c")),
                                  gene_ids = ens, verbose = FALSE)
  expect_equal(r3$ids, ens)
})

test_that(".resolveGeneIDs strips Ensembl versions from rownames", {
  rn <- paste0(sprintf("ENSG%011d", 1:3), ".", 1:3)
  r <- immLynx:::.resolveGeneIDs(.gene_sce(rn), verbose = FALSE)
  expect_equal(r$ids, sprintf("ENSG%011d", 1:3))
})

test_that(".resolveGeneIDs rejects malformed gene_ids", {
  sce <- .gene_sce(c("a", "b", "c"), sprintf("ENSG%011d", 1:3))
  expect_error(immLynx:::.resolveGeneIDs(sce, gene_ids = "nope"),
               "not a rowData column")
  expect_error(immLynx:::.resolveGeneIDs(sce, gene_ids = c("x", "y")),
               "character vector of length nrow")
  expect_error(immLynx:::.resolveGeneIDs(sce, gene_ids = 1:3),
               "character vector of length nrow")
})

test_that("map_symbols = 'never' errors on symbol identifiers", {
  sce <- .gene_sce(c("CD3E", "CD8A", "GZMB"))
  expect_error(
    immLynx:::.resolveGeneIDs(sce, map_symbols = "never", verbose = FALSE),
    "0\\.0% of gene identifiers")
})

test_that(".resolveGeneIDs collapses duplicate Ensembl IDs", {
  ids <- c("ENSG00000000001", "ENSG00000000001", "ENSG00000000002")
  sce <- .gene_sce(c("g1", "g2", "g3"), ids)

  r_sum <- immLynx:::.resolveGeneIDs(sce, gene_ids = "ID", collapse = "sum",
                                     verbose = FALSE)
  expect_equal(r_sum$ids, c("ENSG00000000001", "ENSG00000000002"))
  expect_equal(r_sum$keep, 1:3)
  expect_equal(as.character(r_sum$collapse_groups), ids)

  r_first <- immLynx:::.resolveGeneIDs(sce, gene_ids = "ID",
                                       collapse = "first", verbose = FALSE)
  expect_equal(r_first$keep, c(1L, 3L))
  expect_null(r_first$collapse_groups)

  r_drop <- immLynx:::.resolveGeneIDs(sce, gene_ids = "ID", collapse = "drop",
                                      verbose = FALSE)
  expect_equal(r_drop$ids, "ENSG00000000002")
  expect_equal(r_drop$keep, 3L)
})

test_that(".resolveGeneIDs errors when nothing resolves", {
  sce <- .gene_sce(c("aa", "bb", "cc"))
  expect_error(
    immLynx:::.resolveGeneIDs(sce, map_symbols = "never",
                              ensembl_min_frac = 0, verbose = FALSE),
    "No gene identifiers could be resolved")
})

# ===========================================================================
# Symbol to Ensembl mapping
# ===========================================================================

test_that(".mapSymbolsToEnsembl maps symbols and counts what it dropped", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("AnnotationDbi")

  m <- immLynx:::.mapSymbolsToEnsembl(
    c("CD3E", "HLA-DRA", "NOTAGENEATALL"), multi_map = "first")

  expect_equal(m$map$ENSEMBL[m$map$SYMBOL == "CD3E"], "ENSG00000198851")
  expect_equal(m$n_query, 3L)
  expect_equal(m$n_unmapped, 1L)
  # HLA-DRA sits on several alt haplotypes.
  expect_gte(m$n_ambiguous, 1L)
})

test_that(".mapSymbolsToEnsembl multi_map modes behave differently", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("AnnotationDbi")

  first  <- immLynx:::.mapSymbolsToEnsembl("HLA-DRA", multi_map = "first")
  expand <- immLynx:::.mapSymbolsToEnsembl("HLA-DRA", multi_map = "expand")
  dropped <- immLynx:::.mapSymbolsToEnsembl("HLA-DRA", multi_map = "drop")

  expect_equal(nrow(first$map), 1L)
  expect_gt(nrow(expand$map), 1L)
  expect_equal(nrow(dropped$map), 0L)

  # Lexicographic ordering makes "first" reproducible run to run.
  again <- immLynx:::.mapSymbolsToEnsembl("HLA-DRA", multi_map = "first")
  expect_identical(first$map, again$map)
  expect_equal(first$map$ENSEMBL, min(expand$map$ENSEMBL))
})

test_that(".resolveGeneIDs maps symbols end to end", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("AnnotationDbi")

  sce <- .gene_sce(c("CD3E", "CD8A", "GZMB", "NOTAGENEATALL"))
  r <- immLynx:::.resolveGeneIDs(sce, map_symbols = "auto", verbose = FALSE)

  expect_true(all(immLynx:::.isEnsembl(r$ids)))
  expect_equal(length(r$ids), 3L)          # the bogus symbol is dropped
  expect_true(r$report$mapped)
  expect_equal(r$report$n_unmapped, 1L)

  # "expand" duplicates rows rather than dropping candidates.
  sce2 <- .gene_sce(c("CD3E", "HLA-DRA"))
  re <- immLynx:::.resolveGeneIDs(sce2, map_symbols = "auto",
                                  multi_map = "expand", verbose = FALSE)
  expect_gt(length(re$ids), 2L)
  expect_true(any(duplicated(re$keep)))
})

# ===========================================================================
# Expansion label derivation
# ===========================================================================

test_that(".deriveExpansionLabels reproduces scXpand's 1.5x median rule", {
  sce <- mock_clonal_sce()

  l <- immLynx:::.deriveExpansionLabels(sce, sample_col = "sample",
                                        median_basis = "clone",
                                        verbose = FALSE)

  # Sample A: sizes 4,4,4,4,1,1; median over unique clones = 1; cutoff 1.5.
  expect_equal(l$clone_id_size[1:6], c(4, 4, 4, 4, 1, 1))
  expect_equal(unname(l$median_clone_size[1:6]), rep(1, 6))
  expect_equal(l$expansion[1:6],
               c(rep("expanded", 4), "non-expanded", "non-expanded"))

  # Sample B: two clones of size 2; median 2; cutoff 3; nothing expanded.
  expect_equal(l$clone_id_size[7:10], rep(2, 4))
  expect_equal(l$expansion[7:10], rep("non-expanded", 4))

  # The cell with no clone call stays NA in all three vectors.
  expect_true(is.na(l$clone_id_size[11]))
  expect_true(is.na(l$median_clone_size[11]))
  expect_true(is.na(l$expansion[11]))

  expect_equal(l$n_labelled, 10L)
  expect_equal(l$n_clones, 5L)
  expect_equal(l$n_samples, 2L)
  expect_equal(unname(l$median_by_sample[["A"]]), 1)
})

test_that("median_basis = 'cell' gives a stricter cutoff", {
  sce <- mock_clonal_sce()
  l <- immLynx:::.deriveExpansionLabels(sce, sample_col = "sample",
                                        median_basis = "cell",
                                        verbose = FALSE)
  # Sample A median over cells = median(4,4,4,4,1,1) = 4; cutoff 6.
  expect_equal(unname(l$median_by_sample[["A"]]), 4)
  expect_equal(l$expansion[1:6], rep("non-expanded", 6))
  expect_false(any(l$expansion == "expanded", na.rm = TRUE))
})

test_that(".deriveExpansionLabels warns when sample_col is NULL", {
  sce <- mock_clonal_sce()
  expect_warning(
    l <- immLynx:::.deriveExpansionLabels(sce, verbose = FALSE),
    "single sample")
  # Pooled: sizes 4,1,1,2,2 -> median over unique clones = 2 -> cutoff 3.
  expect_equal(unname(l$median_clone_size[1]), 2)
  expect_equal(l$expansion[1], "expanded")
  expect_equal(l$expansion[7], "non-expanded")
})

test_that(".deriveExpansionLabels falls back from CTstrict to CTaa", {
  sce <- mock_clonal_sce()
  SummarizedExperiment::colData(sce)$CTaa <-
    SummarizedExperiment::colData(sce)$CTstrict
  SummarizedExperiment::colData(sce)$CTstrict <- NULL

  l <- immLynx:::.deriveExpansionLabels(sce, sample_col = "sample",
                                        verbose = FALSE)
  expect_equal(l$clone_col, "CTaa")
  expect_equal(l$n_labelled, 10L)
})

test_that(".deriveExpansionLabels returns NULL without usable clone data", {
  sce <- mock_clonal_sce()
  SummarizedExperiment::colData(sce)$CTstrict <- NULL
  expect_null(immLynx:::.deriveExpansionLabels(sce, sample_col = "sample",
                                               verbose = FALSE))

  # scRepertoire writes literal "NA" tokens for cells with no clone call.
  sce2 <- mock_clonal_sce()
  SummarizedExperiment::colData(sce2)$CTstrict <- "NA"
  expect_null(immLynx:::.deriveExpansionLabels(sce2, sample_col = "sample",
                                               verbose = FALSE))
})

test_that(".deriveExpansionLabels rejects unknown columns", {
  sce <- mock_clonal_sce()
  expect_error(immLynx:::.deriveExpansionLabels(sce, clone_col = "nope"),
               "clone_col")
  expect_error(immLynx:::.deriveExpansionLabels(sce, sample_col = "nope"),
               "sample_col")
})

test_that("clone sizes are tabulated fresh from the clone column", {
  data(immLynx_example, envir = environment())

  l <- immLynx:::.deriveExpansionLabels(immLynx_example,
                                        sample_col = "Patient",
                                        verbose = FALSE)
  cd <- SummarizedExperiment::colData(immLynx_example)
  ok <- !is.na(l$clone_id_size)
  expect_true(any(ok))

  # Pin the arithmetic against a hand-rolled per-patient tabulation.
  key <- paste(cd$Patient, cd$CTstrict, sep = "||")
  expected <- as.numeric(table(key[ok])[key[ok]])
  expect_equal(l$clone_id_size[ok], expected)
})

test_that("per-sample tabulation diverges from a pooled count", {
  # The shipped example has no clone shared across patients, so per-patient
  # and pooled counts happen to coincide there. Construct a case where a
  # clone does span samples: that is the situation in which reading
  # scRepertoire's clonalFrequency (computed under whatever group.by was in
  # effect) would give the wrong per-sample size.
  clone <- c("cX", "cX", "cX", "cY", "cZ", "cX", "cX", "cW", "cV", "cU")
  sample <- c(rep("A", 5), rep("B", 5))
  n <- length(clone)
  counts <- Matrix::Matrix(
    matrix(seq_len(3L * n), nrow = 3L,
           dimnames = list(sprintf("ENSG%011d", 1:3),
                           paste0("cell", seq_len(n)))),
    sparse = TRUE)
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(CTstrict = clone, sample = sample,
                                   row.names = colnames(counts)))

  per_sample <- immLynx:::.deriveExpansionLabels(sce, sample_col = "sample",
                                                 verbose = FALSE)
  pooled <- suppressWarnings(
    immLynx:::.deriveExpansionLabels(sce, verbose = FALSE))

  # cX has 3 cells in A and 2 in B, but 5 overall.
  expect_equal(per_sample$clone_id_size[1], 3)
  expect_equal(per_sample$clone_id_size[6], 2)
  expect_equal(pooled$clone_id_size[1], 5)
  expect_false(isTRUE(all.equal(per_sample$clone_id_size,
                                pooled$clone_id_size)))

  # Within B the unique clone sizes are 2, 1, 1, 1, so the median is 1 and
  # the cutoff 1.5: cX's two cells there clear it. Pooled, cX is counted as
  # a single clone of 5 and B's cells inherit A's count.
  expect_equal(per_sample$median_clone_size[6], 1)
  expect_equal(per_sample$expansion[6], "expanded")
  expect_equal(per_sample$median_clone_size[1], 1)
})

# ===========================================================================
# Counts validation and row collapsing
# ===========================================================================

test_that(".validateCountsAssay accepts raw counts and rejects everything else", {
  X <- Matrix::Matrix(matrix(as.numeric(rpois(40, 3)), 10, 4), sparse = TRUE)
  expect_true(immLynx:::.validateCountsAssay(X, "counts"))

  expect_error(immLynx:::.validateCountsAssay(X, "logcounts"),
               "normalized or log-transformed")

  Xn <- X; Xn@x <- Xn@x + 0.5
  expect_error(immLynx:::.validateCountsAssay(Xn, "counts"),
               "non-integer values")

  Xneg <- X; Xneg@x[1] <- -1
  expect_error(immLynx:::.validateCountsAssay(Xneg, "counts"),
               "negative values")

  expect_error(
    immLynx:::.validateCountsAssay(Matrix::Matrix(0, 4, 4, sparse = TRUE),
                                   "counts"),
    "no non-zero counts")
})

test_that(".collapseRowsSum adds colliding rows together", {
  Y <- Matrix::Matrix(matrix(c(1, 2, 3, 4, 5, 6), nrow = 3), sparse = TRUE)
  colnames(Y) <- c("a", "b")
  out <- immLynx:::.collapseRowsSum(Y, factor(c("E1", "E1", "E2")))

  expect_equal(rownames(out), c("E1", "E2"))
  expect_equal(colnames(out), c("a", "b"))
  expect_equal(as.matrix(unname(as.matrix(out))),
               matrix(c(3, 3, 9, 6), nrow = 2))
})

# ===========================================================================
# Result handling
# ===========================================================================

test_that(".aurocR matches the Mann-Whitney definition", {
  expect_equal(immLynx:::.aurocR(c(0.1, 0.2, 0.8, 0.9), c(0, 0, 1, 1)), 1)
  expect_equal(immLynx:::.aurocR(c(0.9, 0.8, 0.2, 0.1), c(0, 0, 1, 1)), 0)
  expect_equal(immLynx:::.aurocR(rep(0.5, 4), c(0, 0, 1, 1)), 0.5)
  expect_true(is.na(immLynx:::.aurocR(c(0.1, 0.2), c(0, 0))))
  # NAs are dropped pairwise rather than poisoning the result.
  expect_equal(immLynx:::.aurocR(c(0.1, NA, 0.8, 0.9), c(0, 1, 1, 1)), 1)
})

test_that(".flattenMetrics flattens scalars and drops the rest", {
  f <- immLynx:::.flattenMetrics(
    list(overall = list(AUROC = 0.81, n = 10L),
         arr = c(1, 2, 3), txt = "x", empty = list()))
  expect_equal(names(f), c("overall.AUROC", "overall.n"))
  expect_equal(unname(f), c(0.81, 10))
  expect_length(immLynx:::.flattenMetrics(NULL), 0L)
  expect_length(immLynx:::.flattenMetrics(list(1, 2)), 0L)  # unnamed
})

test_that(".alignPredictions matches by barcode", {
  p <- c(0.1, 0.9, 0.5)
  a <- immLynx:::.alignPredictions(p, c("c", "a", "b"), c("a", "b", "c"))
  expect_equal(names(a), c("a", "b", "c"))
  expect_equal(unname(a), c(0.9, 0.5, 0.1))

  expect_error(immLynx:::.alignPredictions(p, c("a", "b"), c("a", "b")),
               "3 predictions for 2 cells")
  expect_warning(
    immLynx:::.alignPredictions(p, c("x", "y", "z"), c("a", "b", "c")),
    "positional alignment")
})

test_that(".writeCellColumn leaves NA for unscored cells", {
  sce <- mock_clonal_sce()
  out <- immLynx:::.writeCellColumn(sce, "prob", c(0.2, 0.8),
                                    c("cell3", "cell1"))
  v <- SummarizedExperiment::colData(out)$prob
  expect_equal(v[1], 0.8)
  expect_equal(v[3], 0.2)
  expect_true(is.na(v[2]))
  expect_length(v, ncol(sce))
})

# ===========================================================================
# runScXpand argument validation (no Python)
# ===========================================================================

test_that("runScXpand rejects malformed input objects", {
  expect_error(runScXpand(list(a = 1)),
               "SingleCellExperiment or Seurat")

  sce <- mock_ensembl_sce(n_genes = 5, n_cells = 4)
  colnames(sce) <- c("a", "a", "b", "c")
  expect_error(runScXpand(sce), "duplicated cell names")
})

test_that("runScXpand validates its scalar arguments", {
  sce <- mock_ensembl_sce(n_genes = 5, n_cells = 6)

  expect_error(runScXpand(sce, assay = "logcounts"), "not found")
  expect_error(runScXpand(sce, threshold = 1.5), "threshold")
  expect_error(runScXpand(sce, ensembl_min_frac = 2), "ensembl_min_frac")
  expect_error(runScXpand(sce, batch_size = -1), "batch_size")
  expect_error(runScXpand(sce, num_workers = 1.5), "num_workers")
  expect_error(runScXpand(sce, model_name = c("a", "b")), "model_name")
  expect_error(runScXpand(sce, clone_col = "nope"), "clone_col")
  expect_error(runScXpand(sce, sample_col = "nope"), "sample_col")
})

test_that("runScXpand warns rather than errors on an unknown model", {
  sce <- mock_ensembl_sce(n_genes = 5, n_cells = 6)
  # An unrecognised name must not hard-stop: upstream may add models. Pair
  # it with a bad assay so validation aborts before any Python work.
  expect_warning(
    expect_error(runScXpand(sce, model_name = "future_model",
                            assay = "nope", verbose = FALSE),
                 "not found"),
    "Unknown model_name")
})

test_that("runScXpand refuses to clobber an existing staged H5AD", {
  sce <- mock_ensembl_sce(n_genes = 5, n_cells = 6)
  wd <- tempfile("scxpand_test_")
  dir.create(wd)
  writeLines("placeholder", file.path(wd, "adata.h5ad"))

  expect_error(runScXpand(sce, work_dir = wd, verbose = FALSE),
               "overwrite = TRUE")
  unlink(wd, recursive = TRUE)
})

# ===========================================================================
# H5AD staging (zellkonverter only, no scXpand env)
# ===========================================================================

test_that(".stageScXpandH5AD writes what scXpand expects", {
  skip_if_not_installed("zellkonverter")

  sce <- mock_clonal_sce()
  SummarizedExperiment::colData(sce)$junk <-
    I(replicate(ncol(sce), 1:2, simplify = FALSE))

  res <- immLynx:::.resolveGeneIDs(sce, verbose = FALSE)
  labs <- immLynx:::.deriveExpansionLabels(sce, sample_col = "sample",
                                           verbose = FALSE)
  f <- tempfile(fileext = ".h5ad")
  immLynx:::.stageScXpandH5AD(sce, res, labs, NULL, "counts", f,
                              verbose = FALSE)
  expect_true(file.exists(f))

  back <- zellkonverter::readH5AD(f)
  expect_true(all(grepl("^ENSG", rownames(back))))
  expect_equal(colnames(back), colnames(sce))
  expect_true(all(c("clone_id_size", "median_clone_size", "expansion") %in%
                    colnames(SummarizedExperiment::colData(back))))
  # CT* and list-columns must not survive into the H5AD.
  expect_false(any(c("CTstrict", "CTaa", "junk") %in%
                     colnames(SummarizedExperiment::colData(back))))
  unlink(f)
})

test_that(".stageScXpandH5AD sums collapsed rows", {
  skip_if_not_installed("zellkonverter")

  sce <- mock_ensembl_sce(n_genes = 4, n_cells = 5)
  SummarizedExperiment::rowData(sce)$ID <- c("ENSG00000000001",
                                             "ENSG00000000001",
                                             "ENSG00000000002",
                                             "ENSG00000000003")
  res <- immLynx:::.resolveGeneIDs(sce, gene_ids = "ID", collapse = "sum",
                                   verbose = FALSE)
  f <- tempfile(fileext = ".h5ad")
  immLynx:::.stageScXpandH5AD(sce, res, NULL, NULL, "counts", f,
                              verbose = FALSE)

  back <- zellkonverter::readH5AD(f)
  expect_equal(nrow(back), 3L)
  raw <- SummarizedExperiment::assay(sce, "counts")
  expect_equal(
    as.numeric(SummarizedExperiment::assay(back)["ENSG00000000001", ]),
    as.numeric(raw[1, ] + raw[2, ]))
  unlink(f)
})

# ===========================================================================
# Python-dependent tests
#
# These build the scXpand basilisk env (several GB) and download a
# pretrained model from figshare. Opt in with IMMLYNX_TEST_SCXPAND=1.
# ===========================================================================

test_that("runScXpand adds predictions to a SingleCellExperiment", {
  skip_if_no_scxpand_model()

  sce <- mock_ensembl_sce(n_genes = 200, n_cells = 40)
  out <- runScXpand(sce, sample_col = "sample", verbose = FALSE)

  expect_s4_class(out, "SingleCellExperiment")
  prob <- SummarizedExperiment::colData(out)$scXpand_expansion_prob
  expect_length(prob, ncol(sce))
  expect_true(all(prob >= 0 & prob <= 1, na.rm = TRUE))

  pred <- SummarizedExperiment::colData(out)$scXpand_expansion_pred
  expect_true(all(pred %in% c("expanded", "non-expanded", NA)))

  md <- S4Vectors::metadata(out)$scXpand
  expect_equal(md$model, "pan_cancer_autoencoder")
  expect_equal(md$threshold, 0.5)
  expect_true(!is.null(md$gene_mapping))
})

test_that("runScXpand returns a data.frame when asked", {
  skip_if_no_scxpand_model()

  sce <- mock_ensembl_sce(n_genes = 200, n_cells = 40)
  df <- runScXpand(sce, sample_col = "sample", return_object = FALSE,
                   verbose = FALSE)

  expect_s3_class(df, "data.frame")
  expect_true(all(c("barcode", "expansion_prob", "expansion_pred") %in%
                    colnames(df)))
  expect_true(all(df$barcode %in% colnames(sce)))
  expect_false(is.null(attr(df, "scXpand")))
})

test_that("runScXpand computes AUROC against derived labels", {
  skip_if_no_scxpand_model()

  sce <- mock_ensembl_sce(n_genes = 200, n_cells = 60)
  out <- runScXpand(sce, sample_col = "sample", verbose = FALSE)
  md <- S4Vectors::metadata(out)$scXpand

  expect_true(!is.null(md$labels))
  expect_true(is.na(md$auroc_immLynx) ||
                (md$auroc_immLynx >= 0 && md$auroc_immLynx <= 1))
})

test_that("runScXpand returns a Seurat object for Seurat input", {
  skip_if_no_scxpand_model()
  skip_if_not_installed("Seurat")

  sce <- mock_ensembl_sce(n_genes = 200, n_cells = 40)
  # Seurat chatters about coercing the matrix and about empty layers while
  # building the fixture; none of it comes from runScXpand.
  seu <- suppressWarnings(Seurat::CreateSeuratObject(
    counts = as.matrix(SummarizedExperiment::assay(sce, "counts")),
    meta.data = as.data.frame(SummarizedExperiment::colData(sce))))

  out <- suppressWarnings(
    runScXpand(seu, sample_col = "sample", verbose = FALSE))
  expect_s4_class(out, "Seurat")
  expect_true("scXpand_expansion_prob" %in% colnames(out[[]]))
  expect_false(is.null(out@misc$scXpand))
})

test_that("listScXpandModels(refresh = TRUE) matches the static table", {
  skip_if_no_scxpand_env()

  live <- listScXpandModels(refresh = TRUE)
  expect_true(all(immLynx:::.SCXPAND_MODELS$model_name %in% live$model_name))
})

test_that("keep_files preserves the staged H5AD", {
  skip_if_no_scxpand_model()

  sce <- mock_ensembl_sce(n_genes = 200, n_cells = 30)
  wd <- tempfile("scxpand_keep_")
  out <- runScXpand(sce, sample_col = "sample", work_dir = wd,
                    keep_files = TRUE, verbose = FALSE)
  expect_true(file.exists(file.path(wd, "adata.h5ad")))
  expect_equal(S4Vectors::metadata(out)$scXpand$work_dir, wd)
  unlink(wd, recursive = TRUE)
})
