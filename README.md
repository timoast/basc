# basc

Barcoded single-cell data processing pipeline.

To install dependencies in a new conda environment:

```
mamba env create -f env.yaml
```

# Running `basc` pipeline

To run `basc`, create a config file containing the barcode combinations for
each sample.

Example `samples.tsv` file:

```
sample_name	i5_index	i7_index	sample_index
H3K27me3	AAAAAAAA	AAAAAAAA	AAAAAAAA
H3K27ac	AAAAAAAA	AAAAAAAA	AAAAAAAA
S2S5P	AAAAAAAA	AAAAAAAA	AAAAAAAA
```

Example `config.yaml` file:

```
reference: path/to/reference
samples: config/samples.tsv
reads: path/to/runfolder
```

Execute the snakemake workflow:

```
snakemake -j 12
```

# Steps

1. Barcode fasta and sample sheet are generated from the config file
2. BCL converstion to fastq files using `bcl2fastq`. Demultiplex based on sample index.
3. Combinatorial barcode demultiplexing using `cutadapt`
4. Extraction of cell barcode using `sinto`
5. Map reads to the genome using `bwa-mem2`
6. BAM file sorted and indexed using `samtools`
7. Fragment file created using `sinto`
