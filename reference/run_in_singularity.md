# Run a singularity container.

Run a singularity container.

## Usage

``` r
run_in_singularity(image_name, volumes = list(), additional_arguments = c())
```

## Arguments

- image_name:

  The singularity image you want to run.

- volumes:

  The list of volumes to mount to the container.

- additional_arguments:

  Vector of arguments to pass to the container.
