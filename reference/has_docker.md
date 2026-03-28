# Check if Docker is Available and Return Its Path

The `has_docker` function checks if the Docker executable is available
in the system's `PATH`. It returns a list with a boolean indicating if
Docker is found, and the path to the Docker executable if it is
available.

## Usage

``` r
has_docker()
```

## Value

A list with two elements:

- found:

  A logical value: `TRUE` if Docker is available in the system's `PATH`,
  `FALSE` otherwise.

- path:

  A character string containing the full path to the Docker executable
  if found, or an empty string if not found.

## Examples

``` r
result <- has_docker()
if (result$found) {
  cat("Docker is available at:", result$path, "\n")
} else {
  cat("Docker is not available.\n")
}
#> Docker is available at: /usr/bin/docker 
```
