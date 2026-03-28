# ribomitoQC

This function perform a set of QC on single cell RNAseq, single cell
nuclear RNAseq, Takara Seeker dense matrix, visium HD. It returns
statistics data in csv format as well as plots QC plots.

## Usage

``` r
ribomitoQC(
  input_dir_path,
  infile,
  rownames_format,
  kit,
  beads_location = NULL,
  counts_threshold = NULL,
  out_metrics = NULL,
  out_freq = NULL,
  mtx_file = NULL,
  barcodes_file = NULL,
  features_file = NULL
)
```

## Arguments

- input_dir_path, :

  a character string indicating the path of the directory containing
  data to be used

- infile:

  Path to dense counts file OR to a MatrixMarket .mtx file OR to a
  directory containing matrix.mtx / barcodes.tsv / features.tsv

- rownames_format:

  One of c("ensemblID:symbol","symbol"), symbol is selected for sparse
  matrix and ensemblID:symbol is the default for dense matrix from rCASC
  output

- kit:

  One of c("scRNAseq", "scMultiomics10X", "Seeker", "visiumHD"). Note:
  scRNAseq works for Illumina, 10Xgenomics sparse and rCASC dense matrix

- beads_location:

  The comma separated file which contains the beads locations required
  for Seeker specific visualisation. For visiumHD this is the binary
  visium parquet that requires parquet-tools to be installed in the
  system. Running parquet-tools generates a comma separate file of
  visium_parquet, which contains the association among barcodes as
  different resolution 0.02, 0.08 and 0.16 um and cellID

- counts_threshold:

  Threshold for "gene present" (default 2, as used in rCASC ribomito
  function)

- out_metrics:

  Output CSV for per-cell metrics

- out_freq:

  Output CSV for N.genes frequency bins

- mtx_file:

  Optional explicit path to matrix.mtx (overrides inference)

- barcodes_file:

  Optional explicit path to barcodes.tsv (overrides inference)

- features_file:

  Optional explicit path to features.tsv (overrides inference)

## Author

Raffaele A Calogero

## Examples

``` r
if (FALSE) { # \dontrun{
ribomitoQC(
  infile = "filtered_feature_bc_matrix/",
  rownames_format="symbol",
  kit = "scRNAseq"
)
} # }
```
