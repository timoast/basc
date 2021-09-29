import pandas as pd

configfile: "config.yaml"
samples = pd.read_csv(config["samples"], sep="\t")
unique_groups = list(set(samples['group_name'].to_list()))

rule create_samplesheet:
    output:
        "{}_barcodes/samplesheet.csv".format(config['name']),
        expand("{x}_barcodes/{y}_i5.fa", x=config['name'], y=unique_groups),
        expand("{x}_barcodes/{y}_i7.fa", x=config['name'], y=unique_groups)
    params:
        samplesheet=config['samples'],
        samplename=config['name'],
        revcomp=config['reverse_complement']
    message: "Generate samplesheet"
    shell:
        """
        Rscript scripts/generate_samplesheet.R \
            {params.samplesheet} \
            {params.samplename} \
            {params.revcomp}
        """

rule bcl2fastq:
    input: "{}_barcodes/samplesheet.csv".format(config['name'])
    message: "BCL conversion"
    output: directory("{}/fastq".format(config['reads']))
    params:
        bases_mask=config['bases_mask'],
        rundir=config['reads']
    threads: 12
    shell:
        """
        bcl2fastq \
          --use-bases-mask={params.bases_mask} \
          --create-fastq-for-index-reads \
          --minimum-trimmed-read-length=8 \
          --mask-short-adapter-reads=8 \
          --ignore-missing-positions \
          --ignore-missing-controls \
          --ignore-missing-filter \
          --ignore-missing-bcls \
          -r {threads} -w {threads} \
          -R {params.rundir} \
          --output-dir={output} \
          --sample-sheet={input}
        """

## After BCL conversion:
# I1 = i7 sample index
# R1 = genomic read 1
# R2 = Tn5 barcode 1, intervening sequence (15)
# R3 (novaseq) = Cell barcode (16), intervening sequence (14), Tn5 barcode 2
# R3 (nextseq) = Tn5 barcode 2, intervening sequence (14), cell barcode (16)
# R4 = genomic read 2

# i7 intervening seq: CCGAGCCCACGAGAC
# i5 intervening seq: AGCAGCCGTCGCAG

rule condense_lanes:
    input:
        fastqs=directory("{}/fastq".format(config['reads'])),
    message: "Combine output from all lanes"
    threads: 1
    output:
    shell:
        """
        # concatenate fastq files from each lane
        # rename with simplified naming: group_sample_R1.fastq.gz, etc
        """

# sinto barcodes just takes the whole read and puts it in the read name
# need to extract cell barcode using cutadapt first
rule extract_cell_barcodes:
    input:
    output:
    message: "Extract cell barcode sequences"
    threads: 4
    shell:
        """
        # need to remove 14 + 8 from start of read 3
        # for samples with 8 bp index, the cell barcode will have an extra T base
        # TODO need to deal with extra base downstream so cell barcode matches
        # across marks
        cutadapt -u 
        """

rule attach_barcodes:
    input:
    output:
    message: "Attach cell barcodes"
    threads: 1
    shell:
        """
        sinto barcode ...
        """

rule demux:
    input:
        fastqs=directory("{}/fastq".format(config['reads'])),
        bc1="{group}_barcodes/{sample}_i5.fa",
        bc2="{group}_barcodes/{sample}_i7.fa"
    output: directory(config['outdir'])
    message: "Demultiplex FASTQ files based on barcode combinations"
    threads: 12
    shell:
        """
        cutadapt \
          -e 0.15 \
          -g file:{input.bc1} \
          -G file:{input.bc2} \
          -o {output}/{{name1}}-{{name2}}.1.fastq.gz \
          -p {output}/{{name1}}-{{name2}}.2.fastq.gz \
          {input.fastqs}/{wildcards.group}_S1_L001_I1_001.fastq.gz \ # TODO
          {input.fastqs}/{wildcards.group}_S1_L001_I1_001.fastq.gz \
        """

# after demux we should have R1, R2 with cell barcode in read name
# split by Tn5 index combination

rule map:
    input:
        genome=config['reference'],
        read1=,
        read2=
    output: 
    message: "Map reads to genome"
    threads: 24
    shell:
        """
        bwa-mem2 mem {input.genome}  \
          | samtools sort - > aln.bam
        
        samtools index ...
        """

rule fragments:
    input:
    output:
    message: "Create fragment file"
    threads: 12
    shell:
        """
        sinto fragment -b {input}
        sort -k1,1 -k2,2n 
        bgzip -@ {threads}
        tabix index -p bed 
        """
