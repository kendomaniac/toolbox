# combatRef

This function is the front-end to the Combat-ref tool, Zhang Comput
Struct Biotechnol J. 2024 Dec 16:27:58-64. This tool seems to perform
sligly better tha Combat-seq from Bioconductor sva packae.

## Usage

``` r
combatRef(input_dir_path, countmatrix_name, metadata_name)
```

## Arguments

- input_dir_path, :

  a character string indicating the path of the directory containing the
  data to be analysed.

- countmatrix_name, :

  name of the count matrix file obtained fastq mapping. This must be a
  tab separated file, with genes on rows and samples on columns

- metadata_name, :

  name of the metadata file, which contains three columns: i. the names
  of the samples (empty header, which correspond to the column names of
  countmatrix file); ii. group_sub (header), which contains the groups
  of sample replicates; iii. batch_sub (header), which contains the
  batches to be removed.

## Author

Raffaele A. Calogero

## Examples

``` r
if (FALSE) { # \dontrun{
combatRef(
  input_dir_path = "/the/input/dir",
  countmatrix_name = "_counts.txt",
  metadata_name = "covar.csv",
)
} # }
```
