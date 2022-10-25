# hpap-tissue-citeseq

This repo covers the scripts necessary to reproduce the downstream analysis of HPAP tissue CITEseq. More info to be added...

## Preprint
Lorem ipsum

## Requirements
1. The `renv` package/environment manager is used to keep track of packages. Inspect the `renv.lock` file to see the necessary packages that are used for this analysis. Afterwards, do the following: 
  ```r
  renv::init()
  renv::restore(lockfile = "renv.lock")
  ```

2. A compute environment with at least 64GB of RAM due to the size of this project. Alternatively, one can set up the Seurat objects with the hdf5 file format to reduce memory usage.

## Upstream analysis/preprocessing
### For RNA
1. Raw BCL files were converted to fastq and demultiplexed using cellranger mkfastq (v7.0., 10X Genomics).
2. Fastq files were aligned to hg38 and processed into count matrices using cellranger count.

### For ADT and HTO
1. Raw BCL files were converted to fastq and demultiplexed using bcl2fastq2 (Illumina).
2. Using the [scc-proc](https://github.com/betts-lab/scc-proc) Snakemake pipeline that is based off the [KITE pipeline](https://github.com/pachterlab/kite) which uses [kallisto and bustools](https://www.kallistobus.tools/getting_started.html). Briefly, Fastq files were aligned to a reference index that either contains the barcodes found in the TotalSeqA Universal Human Cocktail (v1) or the TotalSeaA hashtag barcodes used in the experiment. Resultant alignments were collated with corrected cell barcodes and UMIs to produce a cell x feature matrix. Please refer to the `scc-proc` for more detailed information.

## Downstream analysis
Lorem ipsum