#!/usr/bin/env Rscript

  library(data.table)
  library(Matrix)
  library(ggplot2)

# versione 3.1 mat_dgC <- as(mat_dgC, "dgCMatrix") was changed in mat_dgC <- as(mat_dgC, "CsparseMatrix"), since dgCMatrix is deprecated
# expression binning was inserted in the in the output file
# Version 3.2 a scatter plot for mito.score vs ribo.score was built coloring cells on the basis of the number of detected genes.
# Version 3.3 the function was inserted in a docker ubuntu and modified to be run directly from the file
# Versione 3.4 handle 10Xgenomics scMultiomics data, provide a QC referring to the GeneExpression subset. 
# NOTE for scMultiomics10X only sparse matrix are supported
# Version 3.4 has a MALAT1 versus % MT protein genes data and plot
#Version 3.5 support Takara seeker spatial transcriptomics output: a comma separate file with as rownames gene symbols and as column names the cellID, 
# and a comma spearated file with the cellID and X,Y spatial coordinates

setwd("/scratch")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript Ribomito_v3.3.R <infile> <rownames_format> <kit> [beads_location] [counts_threshold] [out_metrics] [out_freq] [mtx_file] [barcodes_file] [features_file]")
}

infile <- ifelse(length(args) >= 1, args[1])
rownames_format <- ifelse(length(args) >= 2, args[2], "symbol")
kit <- ifelse(length(args) >= 3, args[3], "scRNAseq")
beads_location <- if (length(args) >= 4) args[4] else NULL
counts_threshold <- ifelse(length(args) >= 5, as.numeric(args[5]), 2)
out_metrics <- ifelse(length(args) >= 6, args[6], "ribomito_metrics.csv")
out_freq <- ifelse(length(args) >= 7, args[7], "ribomito_freq.csv")
mtx_file <- if (length(args) >= 8) args[8] else NULL
barcodes_file <- if (length(args) >= 9) args[9] else NULL
features_file <- if (length(args) >= 10) args[10] else NULL


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
#' @param rownames_format One of c("ensemblID:symbol","symbol"), symbol is selected for sparse matrix and ensemblID:symbol is the default for dense matrix from rCASC output
#' @param kit One of c("scRNAseq", "scMultiomics10X"), scRNAseq wirks for Illumina, 10Xgenomics sparse and rCASC dense matrix
#' @param counts_threshold Threshold for "gene present" (default 2, as used in rCASC ribomito function)
#' @param beads_location The comma separated file which contains the beads locations required for Seeker specific visualisation
#' @param mtx_file Optional explicit path to matrix.mtx (overrides inference)
#' @param barcodes_file Optional explicit path to barcodes.tsv (overrides inference)
#' @param features_file Optional explicit path to features.tsv (overrides inference)
#' @return Invisibly returns a list(metrics=..., freq=...)
ribomito_qc <- function(infile,
						    rownames_format = c("ensemblID:symbol", "symbol"),
							  kit = c("scRNAseq", "scMultiomics10X", "Seeker"),
							  beads_location = NULL,
						    counts_threshold = 2,
                out_metrics = "ribomito_metrics.csv",
                out_freq = "ribomito_Ngenes_freq.csv",
                mtx_file = NULL,
                barcodes_file = NULL,
                features_file = NULL) {

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
      mat_dgC <- as(mat_dgC, "CsparseMatrix")
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
  if(kit == "Seeker" & is.null(beads_location)){
    stop("You need to add the name of the beads location file in the beads_location parameter")
  }

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
    if(use_sparse == TRUE & kit == "scMultiomics10X"){
    	# Read sparse matrix
    	counts_mat <- readMM(mtx_file)
    	# Coerce to dgCMatrix for fast column operations
    	counts_mat <- as(counts_mat, "CsparseMatrix")

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
			    feature_type <- as.character(feat[[3]]) #needed to depict the subset of GeneExpression
			    n_gene_expression = length(which(feature_type == "Gene Expression"))
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
		  counts_mat = counts_mat[1:n_gene_expression,]
		  genes = genes[1:n_gene_expression]
		} 
    if (use_sparse == TRUE & kit == "scRNAseq"){
    	# Read sparse matrix
    	counts_mat <- readMM(mtx_file)
    	# Coerce to dgCMatrix for fast column operations
    	counts_mat <- as(counts_mat, "CsparseMatrix")

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
    }
  }
  else {
    if(use_sparse == FALSE & kit == "scMultiomics10X"){
      stop("Dense matrix is not supported for scMultiomics10X. Provide a directory containing matrix.mtx.")
    }
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
    malat1_idx <- grepl("^malat1$", genes, ignore.case = TRUE) #generalization for Hs and Mm
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

  # Frequency bins for N.genes
  bin <- ifelse(as.integer(n_genes) == 0, "0",
                as.character(cut(as.integer(n_genes),
                                 breaks = c(0, 100, 200, 500, 1000, Inf),
                                 right = TRUE,
                                 include.lowest = TRUE,
                                 labels = c("1-100", "101-200", "201-500", "501-1000", ">1000"))))
  
  
  if(kit == "scMultiomics10X"){
    out_dt <- data.table(
      cell.ID    = colnames(counts_mat),
      N.genes    = as.integer(n_genes),
      ribo.score = as.numeric(ribo_score),
      mito.score = as.numeric(mito_score),
      expression.groups = bin,
      MALAT1.counts = Matrix::colSums(counts_mat[malat1_idx, , drop = FALSE])
    )
  } else if(kit == "Seeker"){
    out_dt <- data.table(
      cell.ID    = colnames(counts_mat),
      N.genes    = as.integer(n_genes),
      ribo.score = as.numeric(ribo_score),
      mito.score = as.numeric(mito_score),
      expression.groups = bin
    )
    beads = read.csv(beads_location, row.names=1)
    beads = beads[order(rownames(beads)),]
    setorder(out_dt, cell.ID)
    if(identical(rownames(beads), out_dt$cell.ID)){
      out_dt <- cbind(out_dt, as.data.table(beads[,1:2]))
    }
  } else{
    out_dt <- data.table(
      cell.ID    = colnames(counts_mat),
      N.genes    = as.integer(n_genes),
      ribo.score = as.numeric(ribo_score),
      mito.score = as.numeric(mito_score),
      expression.groups = bin
    )
  }

  fwrite(out_dt, out_metrics)
  
  
  #plots
  
  cols <- c(
    "0" = "black",
    "1-100" = "blue",
    "101-200" = "green",
    "201-500" = "gold",
    "501-1000" = "red",
    ">1000" = "darkred"
  )
  
  if(kit == "scMultiomics10X"){
    out_dt$expression.groups = as.factor(out_dt$expression.groups)
    p = ggplot(out_dt, aes(x=mito.score, y=MALAT1.counts, color=expression.groups)) + 
      geom_point(size=0.8, alpha=0.5, shape = 16) + 
      scale_color_manual(values = cols)
    
    ggsave(
      filename = "ribomito_scMultiomics.jpg",
      plot = p,
      device = "jpeg",
      width = 10, height = 5, units = "cm",
      dpi = 300
    )
  } else if(kit == "Seeker"){
    out_dt$expression.groups = as.factor(out_dt$expression.groups)
    p = ggplot(out_dt, aes(x=SPATIAL_1, y=SPATIAL_2, color=expression.groups)) + 
      geom_point(size=0.8, alpha=0.3, shape = 16) + 
      scale_color_manual(values = cols)
    
    ggsave(
      filename = "ribomito_seeker.jpg",
      plot = p,
      device = "jpeg",
      width = 10, height = 5, units = "cm",
      dpi = 300
    )
  }
  
  
  
  
  out_dt$expression.groups = as.factor(out_dt$expression.groups)
  p = ggplot(out_dt, aes(x=ribo.score, y=mito.score, color=expression.groups)) + 
    geom_point(size=0.8, alpha=0.5, shape = 16) + 
    scale_color_manual(values = cols)
  
  ggsave(
    filename = "ribomito.jpg",
    plot = p,
    device = "jpeg",
    width = 10, height = 5, units = "cm",
    dpi = 300
  )

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
ribomito_qc(infile=infile, rownames_format=rownames_format, kit=kit, beads_location=beads_location, out_metrics=out_metrics, out_freq=out_freq, 
	mtx_file=mtx_file, barcodes_file=barcodes_file, features_file=features_file, counts_threshold=counts_threshold)
#debug(ribomito_qc)
#ribomito_qc("filtered_feature_bc_matrix/", out_metrics="sparse2_metrics.csv", out_freq="sparse2_freq.csv", rownames_format="symbol", kit = "scMultiomics10X")
#ribomito_qc("filtered_feature_bc_matrix/", rownames_format="symbol", kit = "scMultiomics10X")
#ribomito_qc(infile="GSM9428579_dv90pc9.csv", rownames_format="symbol", kit = "scRNAseq")
#ribomito_qc(infile="dv90pc9.csv", rownames_format="symbol", kit = "Seeker", beads_location="dv90pc9_MatchedBeadLocation.csv")




