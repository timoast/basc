# BASC: barcoded single-cell data processing pipeline.

BASC is a snakemake workflow for processing single-cell
chromatin data containing multiple barcodes per sequenced
fragment. Most single-cell chromatin datasets contain a single
barcode (the cell barcode), and can be handled by other
processing pipelines such as cellranger-atac, or simpler
methods. This workflow is designed for experiments where
there are 1 or 2 additional barcodes sequenced for each 
DNA fragment, which may encode different information such 
as the data modality measured by the DNA fragment, droplet
combinatorial indexing well of origin, or both.

## Installation

The workflow is implemented using [Snakemake](https://snakemake.readthedocs.io),
and contains several dependencies. These dependencies can be
easily installed using conda (or [mamba](https://github.com/mamba-org/mamba)).

To get started, clone the git repository:

```
git clone git@github.com:timoast/basc.git
```

To install dependencies in a new conda environment:

```
# using conda
conda env create -f env.yaml
```

```
# using mamba
mamba env create -f env.yaml
```

## Configuring the BASC pipeline

To configure BASC, create a config file containing the barcode combinations for
each sample. See `config.example` for the list of required parameters in the config
file. Alternatively, individual configuration parameters can be passed to snakemake
on the command line, for example:

```
snakemake --config reference=/path/to/reference name=NTT
```

Create a tab-delimited `splitcode` configuration file describing the barcode
sequences that were used for each condition, their location, and the
Hamming distance to tolerate. Example `splitcode.tsv` file:

<!-- | sample_name | mark | i5_index | i7_index | sample_index |
| ----------- | ---- | -------- | -------- | ------------ |
| whole_cell  | H3K27me3 | GGATTGCT | AACAACAC | GGACTCCT,TAGGCATG |
| whole_cell  | H3K27ac | GTGTGACC | ACGTATGG | GGACTCCT,TAGGCATG |
| nuclei      | S2S5P | CGTCTATG | AACATTCC | CACATCGG,GGTTGGCA | -->

| groups  | ids      | tags      | distances |    locations
| ------- | -------- | --------- | --------- | ----------
| bc1     | H3K27me3 | GTGTTGTT  |  2        | 0:0:8
| bc1     | H3K27ac  | AGCATGTG  |  2        | 0:0:8
| bc2     | H3K27ac  | ACGCTCAC  |  2        | 1:0:8
| bc2     | H3K27me3 | AGCAATCC  |  2        | 1:0:8

## Running the workflow

To execute the snakemake workflow first activate the conda environment
containing all of the dependencies: `conda activate basc`. Next, run
the snakemake workflow specifying the `configfile` parameter.

Setting the `-n` option will perform a dry-run, and allow you to see
what steps will be executed in the workflow without computing anyting.
Setting the `-j` parameter will start the workflow with the desired
number of cores.

```
snakemake --configfile /path/to/config -j 24
```

See the [snakemake documentation](https://snakemake.readthedocs.io/en/stable/index.html)
for more information about running snakemake and all the available options.

## Workflow steps

The BASC pipeline executes the following steps:

1. Barcode extraction using [`splitcode`](https://github.com/Yenaled/splitcode)
2. Map reads to the genome using [`bwa-mem2`](https://github.com/bwa-mem2/bwa-mem2)
3. BAM file sorted and indexed using [`samtools`](https://github.com/samtools/samtools)
4. Fragment file created using [`sinto`](https://github.com/timoast/sinto)
