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
sample_name	well	mark  i5_index	i7_index	sample_index
whole_cell	1	H3K27me3  GGATTGCT	AACAACAC	GGACTCCT,TAGGCATG
whole_cell	1	H3K27ac GTGTGACC	ACGTATGG	GGACTCCT,TAGGCATG
whole_cell	1	S2S5P cCGTCTATG	gAACATTCC	GGACTCCT,TAGGCATG
```

Example `config.yaml` file:

```
reference: path/to/bwa-mem2/reference.fa
samples: config/samples.tsv
reads: path/to/runfolder
name: <experiment_name>
reverse_complement: true
outdir: path/to/output
```

Execute the snakemake workflow:

```
conda activate basc
snakemake -j 12
```

# Steps

1. Barcode fasta and sample sheet are generated from the config file
2. BCL converstion to fastq files using `bcl2fastq`  
3. Combinatorial barcode demultiplexing using a custom python script  
4. Map reads to the genome using `bwa-mem2`
5. BAM file sorted and indexed using `samtools`
6. Fragment file created using `sinto`
