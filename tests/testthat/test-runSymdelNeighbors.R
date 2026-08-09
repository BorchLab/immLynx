# Tests for symdel neighbor search (issue #6).
#
# The pure-R helpers are tested without Python. symdel returns a list of
# (i, j, distance) tuples using 0-based indices, and it emits both directions
# of every pair, so the conversion helpers carry the real logic.

# ===========================================================================
# .symdelTripletsToEdges
# ===========================================================================

test_that(".symdelTripletsToEdges converts 0-based indices to sequences", {
  seqs <- c("CASSA", "CASSB", "CASSC")
  # Both directions of a single pair (0, 1), as symdel emits them.
  triplets <- list(c(0, 1, 1), c(1, 0, 1))

  edges <- immLynx:::.symdelTripletsToEdges(triplets, seqs)

  expect_s3_class(edges, "data.frame")
  expect_identical(names(edges), c("from_seq", "to_seq", "distance"))
  expect_identical(nrow(edges), 1L)
  expect_identical(edges$from_seq, "CASSA")
  expect_identical(edges$to_seq, "CASSB")
  expect_identical(edges$distance, 1L)
})

test_that(".symdelTripletsToEdges collapses both directions into one edge", {
  seqs <- c("A", "B", "C")
  triplets <- list(c(0, 1, 1), c(1, 0, 1), c(1, 2, 1), c(2, 1, 1))

  edges <- immLynx:::.symdelTripletsToEdges(triplets, seqs)

  expect_identical(nrow(edges), 2L)
  # Every edge is stored in a canonical orientation, so no pair repeats.
  keys <- paste(edges$from_seq, edges$to_seq)
  expect_identical(anyDuplicated(keys), 0L)
})

test_that(".symdelTripletsToEdges drops self-pairs", {
  seqs <- c("A", "B")
  triplets <- list(c(0, 0, 0), c(0, 1, 1), c(1, 0, 1))

  edges <- immLynx:::.symdelTripletsToEdges(triplets, seqs)

  expect_identical(nrow(edges), 1L)
  expect_false(any(edges$from_seq == edges$to_seq))
})

test_that(".symdelTripletsToEdges returns an empty frame for no neighbors", {
  edges <- immLynx:::.symdelTripletsToEdges(list(), c("A", "B"))

  expect_s3_class(edges, "data.frame")
  expect_identical(nrow(edges), 0L)
  expect_identical(names(edges), c("from_seq", "to_seq", "distance"))
})

# ===========================================================================
# .symdelDegree
# ===========================================================================

test_that(".symdelDegree counts neighbors for each unique sequence", {
  seqs <- c("A", "B", "C", "D")
  edges <- data.frame(
    from_seq = c("A", "B"),
    to_seq   = c("B", "C"),
    distance = c(1L, 1L),
    stringsAsFactors = FALSE
  )

  deg <- immLynx:::.symdelDegree(edges, seqs)

  # B touches both edges, A and C touch one each, D is isolated.
  expect_identical(deg[["A"]], 1L)
  expect_identical(deg[["B"]], 2L)
  expect_identical(deg[["C"]], 1L)
  expect_identical(deg[["D"]], 0L)
})

test_that(".symdelDegree returns all zeros when there are no edges", {
  seqs <- c("A", "B")
  edges <- data.frame(from_seq = character(0), to_seq = character(0),
                      distance = integer(0), stringsAsFactors = FALSE)

  deg <- immLynx:::.symdelDegree(edges, seqs)

  expect_identical(unname(deg), c(0L, 0L))
  expect_identical(names(deg), seqs)
})

test_that(".symdelDegree names every unique sequence in input order", {
  seqs <- c("Z", "Y", "X")
  edges <- data.frame(from_seq = "Z", to_seq = "X", distance = 1L,
                      stringsAsFactors = FALSE)

  deg <- immLynx:::.symdelDegree(edges, seqs)

  expect_identical(names(deg), seqs)
})

# ===========================================================================
# Argument validation
# ===========================================================================

test_that("runSymdelNeighbors rejects non-SingleCellExperiment input", {
  expect_error(runSymdelNeighbors(data.frame(x = 1)),
               "SingleCellExperiment")
})

test_that("runSymdelNeighbors rejects an invalid max_edits", {
  data("immLynx_example", package = "immLynx")

  expect_error(runSymdelNeighbors(immLynx_example, max_edits = 0),
               "max_edits")
  expect_error(runSymdelNeighbors(immLynx_example, max_edits = -1),
               "max_edits")
})

# ===========================================================================
# End-to-end, requires the Python environment
# ===========================================================================

test_that("runSymdelNeighbors returns an edge list", {
  skip_if_no_python()
  data("immLynx_example", package = "immLynx")

  edges <- runSymdelNeighbors(immLynx_example, chains = "TRB",
                              max_edits = 1, return_object = FALSE)

  expect_s3_class(edges, "data.frame")
  expect_identical(names(edges), c("from_seq", "to_seq", "distance"))
  expect_true(all(edges$distance <= 1))
  expect_false(any(edges$from_seq == edges$to_seq))
})

test_that("runSymdelNeighbors adds a degree column to colData", {
  skip_if_no_python()
  data("immLynx_example", package = "immLynx")

  sce <- runSymdelNeighbors(immLynx_example, chains = "TRB",
                            max_edits = 1, return_object = TRUE)

  expect_true("symdel_degree" %in% colnames(SummarizedExperiment::colData(sce)))
  deg <- SummarizedExperiment::colData(sce)$symdel_degree
  expect_identical(length(deg), ncol(immLynx_example))
  # Cells without a TRB chain get NA, everything else is a non-negative count.
  expect_true(all(deg[!is.na(deg)] >= 0))
})

test_that("runSymdelNeighbors honors column_prefix", {
  skip_if_no_python()
  data("immLynx_example", package = "immLynx")

  sce <- runSymdelNeighbors(immLynx_example, chains = "TRB",
                            column_prefix = "nn", return_object = TRUE)

  expect_true("nn_degree" %in% colnames(SummarizedExperiment::colData(sce)))
})

test_that("runSymdelNeighbors finds more neighbors at a larger max_edits", {
  skip_if_no_python()
  data("immLynx_example", package = "immLynx")

  e1 <- runSymdelNeighbors(immLynx_example, chains = "TRB",
                           max_edits = 1, return_object = FALSE)
  e2 <- runSymdelNeighbors(immLynx_example, chains = "TRB",
                           max_edits = 2, return_object = FALSE)

  expect_gte(nrow(e2), nrow(e1))
  expect_true(all(e2$distance <= 2))
})
