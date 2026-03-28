# Check and Identify the Container Runtime

The `detect_container_runtime` function checks if the current R process
is running inside a Docker or Singularity container. It returns the name
of the container runtime (either "Docker", "Singularity", or "None")
along with a boolean indicating whether it is running inside any of
these containers.

## Usage

``` r
detect_container_runtime()
```

## Value

A list with two elements:

- name:

  A character string indicating the container runtime: "Docker",
  "Singularity", or "None".

- is_running:

  A logical value: `TRUE` if the process is running inside a recognized
  container runtime, `FALSE` otherwise.

## Examples

``` r
result <- detect_container_runtime()
cat("Container runtime:", result$name, "\n")
#> Container runtime: None 
cat("Is running in container:", result$is_running, "\n")
#> Is running in container: FALSE 
```
