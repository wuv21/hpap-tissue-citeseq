### Immune perturbations in human pancreas lymphatic tissues prior to and after type 1 diabetes onset

------------------------------------------------------------------------

This repo covers the scripts necessary to reproduce the downstream analysis of HPAP tissue CITEseq and is not a standalone tool for other datasets as the scripts are tailored to this specific analysis and dataset.

### Preprint

The preprint can be found [here](https://www.biorxiv.org/content/10.1101/2024.04.23.590798v2) on biorxiv.

### Requirements

1.  The `renv` package/environment manager is used to keep track of packages. Inspect the `renv.lock` file to see the necessary packages with version numbers that are used for this analysis. Afterwards, do the following to install the packages needed for this analysis:

``` r
renv::init()
renv::restore(lockfile = "renv.lock")
```

2.  A compute environment with at least 128GB of RAM due to the size of this project (and potentially more is needed with hdWGCNA analysis depending on number of cells used). Alternatively, one can set up the Seurat objects with the hdf5 file format to reduce memory usage.

### Downstream analysis

1.  Download the processed RNA, ADT, and HTO files from GEO (or alternatively pre-process from the raw files available on GEO/SRA see below) using the series accession number: [GSE221787](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221787). These files are placed in `../hpap-citeseq-data/rna`, `../hpap-citeseq-data/adt_v2`, `../hpap-citeseq-data/hto_v2` respectively.

2.  Make the following directories: `./outs/rds`, `./rds`, and `matrixForDoubletScoring` which are not loaded in this repo due to size constraints.

3.  Run `main_merge.R` which will generate the necessary `.RData` image files and `.RDS` object files for figure generation.

4.  Run all of the figure R scripts in `./figures` to generate the figures used for the manuscript.

### Expected outputs

-   PDF and PNG files of all figures in the manuscript will be outputted in `./figures/outs/pdf` and `./figures/outs/png` respectively.

### Preprocessing

*(optional if running from raw FASTQ files)*

All fastq files are also found in the GEO series [GSE221787](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221787).

#### For RNA

1.  Raw BCL files were converted to fastq and demultiplexed using cellranger mkfastq (v7.0., 10X Genomics).
2.  Fastq files were aligned to hg38 and processed into count matrices using cellranger count.

#### For ADT and HTO

1.  Raw BCL files were converted to fastq and demultiplexed using bcl2fastq2 (Illumina).
2.  Using the [scc-proc](https://github.com/betts-lab/scc-proc) Snakemake pipeline that is based off the [KITE pipeline](https://github.com/pachterlab/kite) which uses [kallisto and bustools](https://www.kallistobus.tools/getting_started.html). Briefly, Fastq files were aligned to a reference index that either contains the barcodes found in the TotalSeqA Universal Human Cocktail (v1) or the TotalSeaA hashtag barcodes used in the experiment. Resultant alignments were collated with corrected cell barcodes and UMIs to produce a cell x feature matrix. Please refer to [scc-proc](https://github.com/betts-lab/scc-proc) for more detailed information.
