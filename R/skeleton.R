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



#DO NOT TOUCH ANYTHING OVER HERE
#' Check and Identify the Container Runtime
#'
#' The `detect_container_runtime` function checks if the current R process
#' is running  inside a Docker or Singularity container. It returns the name
#' of the container runtime (either "Docker", "Singularity", or "None") along
#' with a boolean indicating whether it is running inside any of these
#' containers.
#'
#' @return
#' A list with two elements:
#' \item{name}{A character string indicating the container runtime: "Docker",
#' "Singularity", or "None".}
#' \item{is_running}{A logical value: `TRUE` if the process is running inside
#' a recognized container runtime, `FALSE` otherwise.}
#'
#' @examples
#' result <- detect_container_runtime()
#' cat("Container runtime:", result$name, "\n")
#' cat("Is running in container:", result$is_running, "\n")
#'
#' @export
detect_container_runtime <- function() {
  detection_map <- list(
    Docker = is_running_in_docker,
    Singularity = is_running_in_singularity
  )
  for (runtime in names(detection_map)) {
    if (detection_map[[runtime]]()) {
      return(list(name = runtime, is_running = TRUE))
    }
  }
  return(list(name = "None", is_running = FALSE))
}

#' Run a container.
#'
#' @param image_name The image you want to run.
#' @param volumes The list of volumes to mount to the container.
#' @param additional_arguments Vector of arguments to pass to the container.
#'
#' @export
run_in_container <- function(image_name,
                             volumes = list(),
                             additional_arguments = c()) {
  detection_map <- list(
    Docker = has_docker(),
    Singularity = has_singularity()
  )

  filtered_providers <- Filter(function(x) x$found, detection_map)
  if (length(filtered_providers) < 1) {
    return("in run_in_container function filtered_providers < 1")
  }
  first_found_name <- names(filtered_providers)[1]
  first_found_fn <- filtered_providers[[first_found_name]]$fn

  first_found_fn(image_name, volumes, additional_arguments)
}

#' Check if Docker is Available and Return Its Path
#'
#' The `has_docker` function checks if the Docker executable is available
#' in the system's `PATH`. It returns a list with a boolean indicating if Docker
#' is found, and the path to the Docker executable if it is available.
#'
#' @return
#' A list with two elements:
#' \item{found}{A logical value: `TRUE` if Docker is available in the system's
#' `PATH`, `FALSE` otherwise.}
#' \item{path}{A character string containing the full path to the Docker
#' executable if found,  or an empty string if not found.}
#'
#' @examples
#' result <- has_docker()
#' if (result$found) {
#'   cat("Docker is available at:", result$path, "\n")
#' } else {
#'   cat("Docker is not available.\n")
#' }
#'
#' @export
has_docker <- function() {
  path <- Sys.which("docker")
  return(list(found = nzchar(path), path = path, fn = run_in_docker))
}

#' Check if the script is running in Docker.
#'
#' @returns A truthy value indicating the state.
#' @export
is_running_in_docker <- function() {
  dockerenv_exists <- file.exists("/.dockerenv")
  cgroup_exists <- file.exists("/proc/1/cgroup")
  in_container_runtime <- FALSE
  if (cgroup_exists) {
    in_container_runtime <- any(
      grepl("docker", readLines("/proc/1/cgroup", warn = FALSE))
    )
  }
  return(dockerenv_exists || in_container_runtime)
}


#' Run a docker container.
#'
#' @param image_name The docker image you want to run.
#' @param volumes The list of volumes to mount to the container.
#' @param additional_arguments Vector of arguments to pass to the container.
#'
#' @export
run_in_docker <- function(image_name,
                          volumes = list(),
                          additional_arguments = c()) {
  base_command <- "run --privileged=true --platform linux/amd64 --rm"
  for (volume in volumes) {
    volume[1] <- normalize_path(volume[1],
      path_mappers = c(docker_mount_mapper)
    )
    base_command <- paste(base_command, "-v", paste(
      volume[1],
      volume[2],
      sep = ":"
    ))
  }
  base_command <- paste(base_command, image_name)
  for (argument in additional_arguments) {
    base_command <- paste(base_command, argument)
  }
  # Using an empty string to output to R's console according to documentation.
  system2("docker", args = base_command, stdout = "", stderr = "")
}

#' Check if Singularity is Available and Return Its Path
#'
#' The `has_singularity` function checks if the Singularity executable is
#' available  in the system's `PATH`. It returns a list with a boolean
#' indicating if Singularity is found, and the path to the Singularity
#' executable if it is available.
#'
#' @return
#' A list with two elements:
#' \item{found}{A logical value: `TRUE` if Singularity is available in the
#' system's `PATH`, `FALSE` otherwise.}
#' \item{path}{A character string containing the full path to the Singularity
#' executable if found, or an empty string if not found.}
#'
#' @examples
#' result <- has_singularity()
#' if (result$found) {
#'   cat("Singularity is available at:", result$path, "\n")
#' } else {
#'   cat("Singularity is not available.\n")
#' }
#'
#' @export
has_singularity <- function() {
  path <- Sys.which("singularity")
  return(list(found = nzchar(path), path = path, fn = run_in_singularity))
}

#' Check if the script is running in a container.
#'
#' @returns A truthy value indicating the state.
#' @export
is_running_in_singularity <- function() {
  cgroup_exists <- file.exists("/proc/1/cgroup")
  in_container_runtime <- FALSE
  if (cgroup_exists) {
    in_container_runtime <- any(
      grepl("singularity", readLines("/proc/1/cgroup", warn = FALSE))
    )
  }
  return(in_container_runtime)
}

#' Run a singularity container.
#'
#' @param image_name The singularity image you want to run.
#' @param volumes The list of volumes to mount to the container.
#' @param additional_arguments Vector of arguments to pass to the container.
#'
#' @export
run_in_singularity <- function(image_name,
                               volumes = list(),
                               additional_arguments = c()) {
  build_mapping <- function(volume) {
    volume[1] <- normalize_path(volume[1])
    return(paste(
      volume[1],
      volume[2],
      sep = ":"
    ))
  }
  base_command <- "run"
  if (length(volumes) > 0) {
    base_command <- paste(base_command, "--bind", build_mapping(volumes[[1]]))
  }
  for (volume in volumes[-1]) {
    base_command <- paste(base_command, build_mapping(volume), sep = ",")
  }
  base_command <- paste(base_command, image_name)
  for (argument in additional_arguments) {
    base_command <- paste(base_command, argument)
  }
  system2("singularity", args = base_command, stdout = "", stderr = "")
}

#' Check if the script is running in a container.
#'
#' @returns A truthy value indicating the state.
is_running_in_docker <- function() {
  dockerenv_exists <- file.exists("/.dockerenv")
  cgroup_exists <- file.exists("/proc/1/cgroup")
  in_container_runtime <- FALSE
  if (cgroup_exists) {
    in_container_runtime <- any(
      grepl("docker", readLines("/proc/1/cgroup", warn = FALSE))
    )
  }
  return(dockerenv_exists || in_container_runtime)
}

#' Gets the absolute path of a file.
#'
#' @param path a normalized, absolute path.
#' @return The absolute host path, if it exists.
#' @export
absolute_path_mapper <- function(path) {
  return(normalizePath(path, mustWork = FALSE))
}

#' Normalizes a path
#'
#' @param path The path to normalize.
#' @param path_mappers The mappers to be utilized to normalize the path.
#' @returns The normalized path.
#' @export
normalize_path <- function(path, path_mappers = c()) {
  path_mappers <- c(
    absolute_path_mapper, # Adds absolute_path_mapper at the beginning.
    path_mappers
  )
  for (mapper in path_mappers) {
    path <- mapper(path)
  }
  return(path)
}

#' Maps a path to host volumes.
#'
#' @param path a normalized, absolute path.
#' @return The absolute host path, if it exists.
#' @export
docker_mount_mapper <- function(path) {
  if (!is_running_in_docker()) {
    return(path)
  }

  # Since we are running in docker, we can use our hostname as an heuristic
  # to obtain our id.
  # TODO: Check if this assumption is right for other container engines and
  #       generalize this implementation.
  hostname <- Sys.info()["nodename"]
  output <- system2("docker",
    args = paste("inspect -f '{{ json .Mounts }}'", hostname),
    stdout = TRUE,
  )
  parsed_output <- jsonlite::fromJSON(output)

  # Iterate over mounts, return first match.
  for (i in seq_len(nrow(parsed_output))) {
    destination <- parsed_output[i, ]$Destination
    if (startsWith(path, destination)) {
      source <- parsed_output[i, ]$Source
      path <- sub(destination, source, path)
      return(path)
    }
  }
  return(path)
}


#' Execute code inside a scratch directory context.
#'
#' @param source_directory The source directory.
#' @param target_directory The target directory.
#' @param modifier A target path modifier.
#' @param callback_function The function to run inside the scratch directory
#' context.
#' @param cleanup_after Clean up after execution.
#' @param copy_pattern Pattern used to copy files from the target_directory back
#' to the source directory.
#' @return Timed directory with the format %Y%m%d_%H%M%S appended to
#' target_directory path.
#' @export
with_scratch <- function(
    source_directory,
    target_directory,
    modifier = modifier_timed_directory,
    callback_function,
    cleanup_after = FALSE,
    copy_pattern = NULL) {
  if (!dir.exists(source_directory)) {
    stop(paste("source_directory", source_directory, "doesn't exist."))
  }
  if (typeof(callback_function) != "closure") {
    stop("callback_function is not a closure.")
  }
  if (modifier != NULL) {
    if (typeof(modifier) != "closure") {
      stop("modifier is not a closure.")
    }
    target_directory <- modifier(target_directory)
  }
  # Create directory if doesn't exist.
  dir.create(target_directory, recursive = TRUE)

  # Calls the callback function.
  callback_function(target_directory)

  if (!cleanup_after) {
    return(cat("in with_scratch !cleanup_after"))
  }
  # Copy back from target_directory to source_directory.
  file.copy(
    list.files(
      pattern = copy_pattern,
      path = target_directory,
      full.names = TRUE,
      all.files = TRUE,
      recursive = TRUE,
    ),
    source_directory
  )
  unlink(target_directory, recursive = TRUE)
}

#' Append a timed directory with the format %Y%m%d_%H%M%S to the
#' target_directory path.
#'
#' @param target_directory The target directory.
#' @return Timed directory with the format %Y%m%d_%H%M%S appended to
#' target_directory path.
#' @export
modifier_timed_directory <- function(target_directory) {
  return(
    file.path(target_directory, strftime(Sys.time(), "%Y%m%d_%H%M%S"))
  )
}


