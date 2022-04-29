# basc

Barcoded single-cell data processing pipeline.

To install dependencies in a new conda environment:

```
mamba env create -f env.yaml
```

# Configuring the `basc` pipeline

To configure `basc`, create a config file containing the barcode combinations for
each sample. See `config.example` for the list of required parameters in the config
file. Alternatively, individual configuration parameters can be passed to snakemake
on the command line, for example:

```
snakemake --config reference=/path/to/reference name=NTT
```

Create a tab-delimited file describing the barcode sequences that were used for
each sample or assay. Example `samples.tsv` file:

```
sample_name	well	mark  i5_index	i7_index	sample_index
whole_cell	1	H3K27me3  GGATTGCT	AACAACAC	GGACTCCT,TAGGCATG
whole_cell	1	H3K27ac GTGTGACC	ACGTATGG	GGACTCCT,TAGGCATG
whole_cell	1	S2S5P cCGTCTATG	gAACATTCC	GGACTCCT,TAGGCATG
```

# Running the workflow

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

# Workflow steps

The `basc` pipeline executes the following steps:

1. Barcode fasta and sample sheet are generated from the config file
2. BCL converstion to fastq files using `bcl2fastq` 
3. Combinatorial barcode demultiplexing using a custom python script 
4. Map reads to the genome using `bwa-mem2`
5. BAM file sorted and indexed using `samtools`
6. Fragment file created using `sinto`
