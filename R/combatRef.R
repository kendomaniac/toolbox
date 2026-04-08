#' combatRef
#'
#' @description This function is the front-end to the Combat-ref tool, Zhang Comput Struct Biotechnol J. 2024 Dec 16:27:58-64. This tool seems to perform sligly better tha Combat-seq from Bioconductor sva packae.
#' @param input_dir_path, a character string indicating the path of the directory containing the data to be analysed.
#' @param countmatrix_name, name of the count matrix file obtained fastq mapping. This must be a tab separated file, with genes on rows and samples on columns
#' @param metadata_name, name of the metadata file, which contains three columns: i. the names of the samples (empty header, which correspond to the column names of countmatrix file); ii. group_sub (header), which contains the groups of sample replicates; iii. batch_sub (header), which contains the batches to be removed.
#' @author Raffaele A. Calogero
#'
#' @examples
#' \dontrun{
#' combatRef(
#'   input_dir_path = "/the/input/dir",
#'   countmatrix_name = "_counts.txt",
#'   metadata_name = "covar.csv",
#' )
#' }
#' @export
combatRef <- function(input_dir_path, countmatrix_name, metadata_name) {
  # Type checking.
  if (typeof(input_dir_path) != "character") {
    stop(paste("input_dir_path type is", paste0(typeof(input_dir_path), "."), "It should be \"character\""))
  }
  if (typeof(countmatrix_name) != "character") {
    stop(paste("countmatrix_name type is", paste0(typeof(countmatrix_name), "."), "It should be \"character\""))
  }
  if (typeof(metadata_name) != "character") {
    stop(paste("metadata_name type is", paste0(typeof(metadata_name), "."), "It should be \"character\""))
  }

  # Check if input_dir_path exists
  if (is_running_in_docker()) {
    if (!dir.exists(input_dir_path)) {
      stop(paste("input_dir_path:", input_dir_path, "does not exist."))
    }
  }

  # Executing the docker job HERE IS WHERE TO CHANGE THE DOCKER NAME AND THE NAME OF THE SCRIPT IN THE DOCKER
  run_in_docker(
    image_name = paste0("repbioinfo/batch_effect_correction:latest"),
    volumes = list(
      c(input_dir_path, "/scratch")
    ),
    additional_arguments = c(
      "Rscript /Combat-ref/batch_correct.R",
      countmatrix_name,
      metadata_name
    )
  )
}


