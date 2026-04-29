#' ribomitoQC
#'
#' @description This function perform a set of QC on single cell RNAseq, single cell nuclear RNAseq, Takara Seeker dense matrix, visium HD. It returns statistics data in csv format as well as plots QC plots.
#' @param input_dir_path, a character string indicating the path of the directory containing data to be used
#' @param infile Path to dense counts file OR to a MatrixMarket .mtx file OR to a directory containing matrix.mtx / barcodes.tsv / features.tsv
#' @param out_metrics Output CSV for per-cell metrics
#' @param out_freq Output CSV for N.genes frequency bins
#' @param rownames_format One of c("ensemblID:symbol","symbol"), symbol is selected for sparse matrix and ensemblID:symbol is the default for dense matrix from rCASC output
#' @param kit One of c("scRNAseq", "scMultiomics10X", "Seeker", "visiumHD"). Note: scRNAseq works for Illumina, 10Xgenomics sparse and rCASC dense matrix
#' @param organism One of c("hs", "mm"). Note: organism is required to the correct selection of ribosomal and mitochondrial annotation
#' @param counts_threshold Threshold for "gene present" (default 2, as used in rCASC ribomito function)
#' @param beads_location The comma separated file which contains the beads locations required for Seeker specific visualisation. For visiumHD this is  the binary visium parquet that requires parquet-tools to be installed in the system. Running parquet-tools generates a comma separate file of visium_parquet, which contains the association among barcodes as different resolution 0.02, 0.08 and 0.16 um and cellID
#' @param mtx_file Optional explicit path to matrix.mtx (overrides inference)
#' @param barcodes_file Optional explicit path to barcodes.tsv (overrides inference)
#' @param features_file Optional explicit path to features.tsv (overrides inference)
#' @author Raffaele A Calogero
#'
#' @examples
#' \dontrun{
#' ribomitoQC(
#'   infile = "filtered_feature_bc_matrix/",
#'   rownames_format="symbol",
#'   kit = "scRNAseq"
#'	 organism = "hs"
#' )
#' }
#' @export
ribomitoQC <- function(input_dir_path, infile, rownames_format, kit, organism, beads_location=NULL, counts_threshold=NULL, 
	out_metrics=NULL, out_freq=NULL, mtx_file=NULL, barcodes_file=NULL, features_file=NULL) {
  # Type checking.
  if (typeof(input_dir_path) != "character") {
    stop(paste("input_dir_path type is", paste0(typeof(input_dir_path), "."), "It should be \"character\""))
  }
    # Type checking.
  if (typeof(infile) != "character") {
    stop(paste("infile type is", paste0(typeof(infile), "."), "It should be \"character\""))
  }
  if (typeof(rownames_format) != "character") {
    stop(paste("rownames_format type is", paste0(typeof(rownames_format), "."), "It should be \"character\""))
  }
  if (typeof(kit) != "character") {
    stop(paste("kit type is", paste0(typeof(kit), "."), "It should be \"character\""))
  }
  if (typeof(organism) != "character") {
    stop(paste("kit type is", paste0(typeof(organism), "."), "It should be \"character\""))
  }
  # Executing the docker job
  run_in_docker(
    image_name = paste0("repbioinfo/ribomito:latest"),
    volumes = list(
      c(input_dir_path, "/scratch")
    ),
    additional_arguments = c(
      "Rscript /home/Ribomito_v3.6.3.R",
      infile,
      rownames_format,
      kit,
	  organism,
      beads_location,
	  counts_threshold,
	  out_metrics,
	  out_freq,
	  mtx_file,
	  barcodes_file,
	  features_file
    )
  )
}


