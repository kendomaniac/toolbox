#' skeleton
#'
#' @description This function is the basic for designing a local front end for function based on Rscript embedded in docker containers.
#' @param input_dir_path, a character string indicating the path of the directory containing the fastq files and the csv files obtained from the indexing.Important this must be always present, the other parameters can be modified
#' @param countmatrix_name, name of the count matrix file obtained after the genome indexing
#' @param metadata_name, name of the metadata file obtained after the genome indexing
#' @param reference_group, name of the reference group inside the meetadata (wt, cpes,...)
#' @param organism, name of the organism subject of the analysis. Supported organisms are: 'Homo sapiens', 'Mus musculus' or 'Drosophila melanogaster'
#' @author Luca Alessandri, Agata D'Onofrio, Eliseo Martelli
#'
#' @examples
#' \dontrun{
#' skeleton(
#'   input_dir_path = "/the/input/dir",
#'   countmatrix_name = "gene_count_matrix.csv",
#'   metadata_name = "Covariatesstat.csv",
#'   reference_group = "wt",
#'   organism = "Drosophilamelanogaster"
#' )
#' }
#' @export
skeleton <- function(input_dir_path, countmatrix_name, metadata_name, reference_group, organism) {
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
  if (typeof(reference_group) != "character") {
    stop(paste("reference_group type is", paste0(typeof(reference_group), "."), "It should be \"character\""))
  }
  if (typeof(organism) != "character") {
    stop(paste("organism type is", paste0(typeof(organism), "."), "It should be \"character\""))
  }

  # Check if input_dir_path exists
  if (is_running_in_docker()) {
    if (!dir.exists(input_dir_path)) {
      stop(paste("input_dir_path:", input_dir_path, "does not exist."))
    }
  }

  # Executing the docker job HERE IS WHERE TO CHANGE THE DOCKER NAME AND THE NAME OF THE SCRIPT IN THE DOCKER
  run_in_docker(
    image_name = paste0("repbioinfo/rnaseqstar_v2:latest"),
    volumes = list(
      c(input_dir_path, "/scratch")
    ),
    additional_arguments = c(
      "Rscript /home/Skeleton.R",
      countmatrix_name,
      metadata_name,
      reference_group,
      organism
    )
  )
}


