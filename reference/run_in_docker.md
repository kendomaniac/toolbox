# Run a docker container.

Run a docker container.

## Usage

``` r
run_in_docker(image_name, volumes = list(), additional_arguments = c())
```

## Arguments

- image_name:

  The docker image you want to run.

- volumes:

  The list of volumes to mount to the container.

- additional_arguments:

  Vector of arguments to pass to the container.
