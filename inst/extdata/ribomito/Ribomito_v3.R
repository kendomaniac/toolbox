#!/usr/bin/env Rscript

  library(data.table)
  library(Matrix)


#' ribomito_qc
#'
#' Compute per-cell QC metrics: N.genes (counts > 2), ribo.score (%), mito.score (%)
#' Supports:
#'   1) Dense table input (first column = gene identifiers; remaining columns = cells)
#'   2) 10x-style sparse MatrixMarket input (matrix.mtx + barcodes.tsv + features.tsv)
#'
#' @param infile Path to dense counts file OR to a MatrixMarket .mtx file OR to a directory containing matrix.mtx / barcodes.tsv / features.tsv
#' @param out_metrics Output CSV for per-cell metrics
#' @param out_freq Output CSV for N.genes frequency bins
#' @param rownames_format One of c("ensemblID:symbol","symbol")
#' @param mtx_file Optional explicit path to matrix.mtx (overrides inference)
#' @param barcodes_file Optional explicit path to barcodes.tsv (overrides inference)
#' @param features_file Optional explicit path to features.tsv (overrides inference)
#' @param counts_threshold Threshold for "gene present" (default 2, matching original function)
#' @return Invisibly returns a list(metrics=..., freq=...)
ribomito_qc <- function(infile,
                        out_metrics = "ribomito_metrics.csv",
                        out_freq    = "ribomito_Ngenes_freq.csv",
                        rownames_format = c("ensemblID:symbol", "symbol"),
                        mtx_file = NULL,
                        barcodes_file = NULL,
                        features_file = NULL,
                        counts_threshold = 2) {

  rownames_format <- match.arg(rownames_format)

  if (is.null(infile) || !nzchar(infile)) {
    stop("infile must be provided.")
  }
  if (!file.exists(infile)) {
    stop("Input path does not exist: ", infile)
  }

  # -----------------------------
  # Helpers
  # -----------------------------
  is_mtx_path <- function(x) {
    grepl("\\.mtx(\\.gz)?$", x, ignore.case = TRUE)
  }

  # Determine paths for MatrixMarket trio if applicable
  infer_10x_paths <- function(path) {
    if (dir.exists(path)) {
      mtx  <- file.path(path, "matrix.mtx")
      mtxg <- file.path(path, "matrix.mtx.gz")
      bc   <- file.path(path, "barcodes.tsv")
      bcg  <- file.path(path, "barcodes.tsv.gz")
      feat <- file.path(path, "features.tsv")
      featg<- file.path(path, "features.tsv.gz")

      mtx_file0 <- if (file.exists(mtx)) mtx else if (file.exists(mtxg)) mtxg else NA_character_
      bc_file0  <- if (file.exists(bc)) bc else if (file.exists(bcg)) bcg else NA_character_
      ft_file0  <- if (file.exists(feat)) feat else if (file.exists(featg)) featg else NA_character_

      return(list(mtx_file = mtx_file0, barcodes_file = bc_file0, features_file = ft_file0))
    }

    if (is_mtx_path(path)) {
      d <- dirname(path)
      # If only .mtx provided, try infer barcodes/features alongside it
      bc   <- file.path(d, "barcodes.tsv")
      bcg  <- file.path(d, "barcodes.tsv.gz")
      feat <- file.path(d, "features.tsv")
      featg<- file.path(d, "features.tsv.gz")

      bc_file0  <- if (file.exists(bc)) bc else if (file.exists(bcg)) bcg else NA_character_
      ft_file0  <- if (file.exists(feat)) feat else if (file.exists(featg)) featg else NA_character_

      return(list(mtx_file = path, barcodes_file = bc_file0, features_file = ft_file0))
    }

    return(list(mtx_file = NA_character_, barcodes_file = NA_character_, features_file = NA_character_))
  }

  # Efficient N.genes for dgCMatrix without creating a full logical sparse matrix
  # Counts number of nonzero entries > threshold per column.
  n_genes_sparse_gt <- function(mat_dgC, threshold) {
    if (!inherits(mat_dgC, "dgCMatrix")) {
      mat_dgC <- as(mat_dgC, "dgCMatrix")
    }
    nnz <- length(mat_dgC@x)
    if (nnz == 0L) return(integer(ncol(mat_dgC)))

    gt <- mat_dgC@x > threshold
    p <- mat_dgC@p # 0-based, length = ncol+1
    # cumulative sum over nonzeros
    cgt <- cumsum(as.integer(gt))

    # per-column counts using p
    res <- integer(ncol(mat_dgC))
    for (j in seq_len(ncol(mat_dgC))) {
      start <- p[j]     # 0..nnz
      end   <- p[j + 1] # 0..nnz
      if (end == start) {
        res[j] <- 0L
      } else if (start == 0L) {
        res[j] <- cgt[end]
      } else {
        res[j] <- cgt[end] - cgt[start]
      }
    }
    res
  }

  # -----------------------------
  # Load counts as matrix (dense or sparse)
  # -----------------------------
  counts_mat <- NULL
  genes <- NULL
  cell_ids <- NULL

  use_sparse <- FALSE

  # If user provided explicit mtx paths or infile points to .mtx/dir containing matrix.mtx -> use sparse
  if (!is.null(mtx_file) || !is.null(barcodes_file) || !is.null(features_file) ||
      dir.exists(infile) || is_mtx_path(infile)) {

    paths <- infer_10x_paths(infile)
    mtx_file <- if (!is.null(mtx_file)) mtx_file else paths$mtx_file
    barcodes_file <- if (!is.null(barcodes_file)) barcodes_file else paths$barcodes_file
    features_file <- if (!is.null(features_file)) features_file else paths$features_file

    if (is.na(mtx_file) || !file.exists(mtx_file)) {
      stop("MatrixMarket input detected/selected, but matrix.mtx not found. Provide mtx_file or a directory containing matrix.mtx.")
    }
    if (is.na(barcodes_file) || !file.exists(barcodes_file)) {
      stop("MatrixMarket input detected/selected, but barcodes.tsv not found. Provide barcodes_file or a directory containing barcodes.tsv.")
    }
    if (is.na(features_file) || !file.exists(features_file)) {
      stop("MatrixMarket input detected/selected, but features.tsv not found. Provide features_file or a directory containing features.tsv.")
    }

    use_sparse <- TRUE

    # Read sparse matrix
    counts_mat <- readMM(mtx_file)
    # Coerce to dgCMatrix for fast column operations
    counts_mat <- as(counts_mat, "dgCMatrix")

    # Read barcodes (cell IDs)
    bc <- fread(barcodes_file, header = FALSE)
    if (ncol(bc) < 1) stop("barcodes.tsv must have at least 1 column.")
    cell_ids <- as.character(bc[[1]])

    # Read features (gene annotations)
    feat <- fread(features_file, header = FALSE)
    if (ncol(feat) < 1) stop("features.tsv must have at least 1 column (ensembl IDs).")

    # Typical 10x features.tsv has: V1=ensembl, V2=symbol, V3=feature_type
    if (ncol(feat) >= 2) {
      ensembl <- as.character(feat[[1]])
      symbol  <- as.character(feat[[2]])
    } else {
      ensembl <- as.character(feat[[1]])
      symbol  <- ensembl
    }

    if (rownames_format == "symbol") {
      genes <- symbol
    } else { # "ensemblID:symbol"
      genes <- paste0(ensembl, ":", symbol)
    }

    # Attach dimnames
    if (length(genes) != nrow(counts_mat)) {
      stop("Row count mismatch: features.tsv has ", length(genes),
           " rows but matrix has ", nrow(counts_mat), " rows.")
    }
    if (length(cell_ids) != ncol(counts_mat)) {
      stop("Column count mismatch: barcodes.tsv has ", length(cell_ids),
           " rows but matrix has ", ncol(counts_mat), " columns.")
    }

    rownames(counts_mat) <- genes
    colnames(counts_mat) <- cell_ids

  } else {
    # Dense input mode (original behavior)
    dt <- fread(infile, check.names = TRUE)

    if (ncol(dt) < 2) {
      stop("Dense input must have at least 2 columns (gene column + >=1 cell): ", infile)
    }

    gene_col <- names(dt)[1]
    genes <- dt[[gene_col]]

    count_dt <- dt[, -1, with = FALSE]
    counts_mat <- as.matrix(count_dt)
    storage.mode(counts_mat) <- "numeric"
    rownames(counts_mat) <- genes
    cell_ids <- colnames(counts_mat)
  }

  # -----------------------------
  # QC computations
  # -----------------------------
  # Total counts per cell
  total_counts <- if (use_sparse) {
    Matrix::colSums(counts_mat)
  } else {
    colSums(counts_mat, na.rm = TRUE)
  }

  total_counts_safe <- ifelse(total_counts == 0, NA_real_, as.numeric(total_counts))

  # Identify mito and ribo genes by prefix rules
  if (rownames_format == "symbol") {
    mito_idx <- grepl("^MT-", genes)
    ribo_idx <- grepl("^RP", genes)
  } else { # ensemblID:symbol
    mito_idx <- grepl(":MT-", genes)
    ribo_idx <- grepl(":RP", genes)
  }

  mito_counts <- if (any(mito_idx)) {
    if (use_sparse) Matrix::colSums(counts_mat[mito_idx, , drop = FALSE]) else colSums(counts_mat[mito_idx, , drop = FALSE], na.rm = TRUE)
  } else {
    rep(0, ncol(counts_mat))
  }

  ribo_counts <- if (any(ribo_idx)) {
    if (use_sparse) Matrix::colSums(counts_mat[ribo_idx, , drop = FALSE]) else colSums(counts_mat[ribo_idx, , drop = FALSE], na.rm = TRUE)
  } else {
    rep(0, ncol(counts_mat))
  }

  mito_score <- 100 * (as.numeric(mito_counts) / total_counts_safe)
  ribo_score <- 100 * (as.numeric(ribo_counts) / total_counts_safe)

  # Number of genes called present: counts > counts_threshold
  n_genes <- if (use_sparse) {
    n_genes_sparse_gt(counts_mat, counts_threshold)
  } else {
    colSums(counts_mat > counts_threshold, na.rm = TRUE)
  }

  out_dt <- data.table(
    cell.ID    = colnames(counts_mat),
    N.genes    = as.integer(n_genes),
    ribo.score = as.numeric(ribo_score),
    mito.score = as.numeric(mito_score)
  )

  fwrite(out_dt, out_metrics)

  # Frequency bins for N.genes
  bin <- ifelse(out_dt$N.genes == 0, "0",
                as.character(cut(out_dt$N.genes,
                                 breaks = c(0, 100, 200, 500, 1000, Inf),
                                 right = TRUE,
                                 include.lowest = TRUE,
                                 labels = c("1-100", "101-200", "201-500", "501-1000", ">1000"))))

  freq_dt <- data.table(
    N.genes.range = c("0", "1-100", "101-200", "201-500", "501-1000", ">1000")
  )
  freq_dt[, Frequency := as.integer(table(factor(bin, levels = N.genes.range)))]

  fwrite(freq_dt, out_freq)

  invisible(list(metrics = out_dt, freq = freq_dt))
}

# -----------------------------
# Example usage:
# Dense:
# ribomito_qc("counts_dense.csv", out_metrics="metrics.csv", out_freq="freq.csv", rownames_format="ensemblID:symbol")
#
# Sparse 10x directory (contains matrix.mtx, barcodes.tsv, features.tsv):
# ribomito_qc("filtered_feature_bc_matrix/", out_metrics="metrics.csv", out_freq="freq.csv", rownames_format="ensemblID:symbol")
#
# Sparse explicit paths:
# ribomito_qc("matrix.mtx", barcodes_file="barcodes.tsv", features_file="features.tsv")
# -----------------------------

#test
ribomito_qc("real_dataset.csv", out_metrics="dense_metrics.csv", out_freq="dense_freq.csv", rownames_format="ensemblID:symbol")
ribomito_qc("filtered_feature_bc_matrix/", out_metrics="sparse1_metrics.csv", out_freq="sparse_freq.csv", rownames_format="symbol")
