# Execute code inside a scratch directory context.

Execute code inside a scratch directory context.

## Usage

``` r
with_scratch(
  source_directory,
  target_directory,
  modifier = modifier_timed_directory,
  callback_function,
  cleanup_after = FALSE,
  copy_pattern = NULL
)
```

## Arguments

- source_directory:

  The source directory.

- target_directory:

  The target directory.

- modifier:

  A target path modifier.

- callback_function:

  The function to run inside the scratch directory context.

- cleanup_after:

  Clean up after execution.

- copy_pattern:

  Pattern used to copy files from the target_directory back to the
  source directory.

## Value

Timed directory with the format %Y%m%d\_%H%M%S appended to
target_directory path.
