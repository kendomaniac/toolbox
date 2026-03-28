# skeleton

This function is the basic for designing a local front end for function
based on Rscript embedded in docker containers.

## Usage

``` r
skeleton(
  input_dir_path,
  countmatrix_name,
  metadata_name,
  reference_group,
  organism
)
```

## Arguments

- input_dir_path, :

  a character string indicating the path of the directory containing the
  fastq files and the csv files obtained from the indexing.Important
  this must be always present, the other parameters can be modified

- countmatrix_name, :

  name of the count matrix file obtained after the genome indexing

- metadata_name, :

  name of the metadata file obtained after the genome indexing

- reference_group, :

  name of the reference group inside the meetadata (wt, cpes,...)

- organism, :

  name of the organism subject of the analysis. Supported organisms are:
  'Homo sapiens', 'Mus musculus' or 'Drosophila melanogaster'

## Author

Luca Alessandri, Agata D'Onofrio, Eliseo Martelli

## Examples

``` r
if (FALSE) { # \dontrun{
skeleton(
  input_dir_path = "/the/input/dir",
  countmatrix_name = "gene_count_matrix.csv",
  metadata_name = "Covariatesstat.csv",
  reference_group = "wt",
  organism = "Drosophilamelanogaster"
)
} # }
```
