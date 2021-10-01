import pandas as pd

configfile: "config.yaml"
samples = pd.read_csv(config["samples"], sep="\t")
unique_groups = list(set(samples['group_name'].to_list()))

rule create_samplesheet:
    output:
        "{}_barcodes/samplesheet.csv".format(config['name']),
        "{}_barcodes/tn5_i5.fa".format(config['name']),
        "{}_barcodes/tn5_i7.fa".format(config['name']),
        "{}_barcodes/sampleindex.fa".format(config['name'])
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
    output:
        read1="{}/fastq/unassigned_R1.fastq.gz".format(config['reads']),
        read2="{}/fastq/unassigned_R2.fastq.gz".format(config['reads']),
        read3="{}/fastq/unassigned_R3.fastq.gz".format(config['reads']),
        read4="{}/fastq/unassigned_R4.fastq.gz".format(config['reads'])
    params:
        rundir=config['reads']
    threads: 12
    shell:
        """
        bcl2fastq \
          --use-bases-mask=Y*,Y*,Y*,Y* \
          --create-fastq-for-index-reads \
          --minimum-trimmed-read-length=8 \
          --mask-short-adapter-reads=8 \
          --ignore-missing-positions \
          --ignore-missing-controls \
          --ignore-missing-filter \
          --ignore-missing-bcls \
          -r {threads} -w {threads} \
          -R {params.rundir} \
          --output-dir={params.rundir}/fastq \
          --sample-sheet={input}
         
        # collapse lanes
        # read1
        cat {params.rundir}/fastq/unassigned_S1_L001_R1_001.fastq.gz {params.rundir}/fastq/unassigned_S1_L002_R1_001.fastq.gz {params.rundir}/fastq/unassigned_S1_L003_R1_001.fastq.gz {params.rundir}/fastq/unassigned_S1_L004_R1_001.fastq.gz > {output.read1}

        rm {params.rundir}/fastq/unassigned_S1_L001_R1_001.fastq.gz \
           {params.rundir}/fastq/unassigned_S1_L002_R1_001.fastq.gz \
           {params.rundir}/fastq/unassigned_S1_L003_R1_001.fastq.gz \
           {params.rundir}/fastq/unassigned_S1_L004_R1_001.fastq.gz

        # read2
        cat {params.rundir}/fastq/unassigned_S1_L001_R2_001.fastq.gz {params.rundir}/fastq/unassigned_S1_L002_R2_001.fastq.gz {params.rundir}/fastq/unassigned_S1_L003_R2_001.fastq.gz {params.rundir}/fastq/unassigned_S1_L004_R2_001.fastq.gz > {output.read2}

        rm {params.rundir}/fastq/unassigned_S1_L001_R2_001.fastq.gz \
           {params.rundir}/fastq/unassigned_S1_L002_R2_001.fastq.gz \
           {params.rundir}/fastq/unassigned_S1_L003_R2_001.fastq.gz \
           {params.rundir}/fastq/unassigned_S1_L004_R2_001.fastq.gz

        # read3
        cat {params.rundir}/fastq/unassigned_S1_L001_R3_001.fastq.gz {params.rundir}/fastq/unassigned_S1_L002_R3_001.fastq.gz {params.rundir}/fastq/unassigned_S1_L003_R3_001.fastq.gz {params.rundir}/fastq/unassigned_S1_L004_R3_001.fastq.gz > {output.read3}

        rm {params.rundir}/fastq/unassigned_S1_L001_R3_001.fastq.gz \
           {params.rundir}/fastq/unassigned_S1_L002_R3_001.fastq.gz \
           {params.rundir}/fastq/unassigned_S1_L003_R3_001.fastq.gz \
           {params.rundir}/fastq/unassigned_S1_L004_R3_001.fastq.gz

        # read4
        cat {params.rundir}/fastq/unassigned_S1_L001_R4_001.fastq.gz {params.rundir}/fastq/unassigned_S1_L002_R4_001.fastq.gz {params.rundir}/fastq/unassigned_S1_L003_R4_001.fastq.gz {params.rundir}/fastq/unassigned_S1_L004_R4_001.fastq.gz > {output.read4}

        rm {params.rundir}/fastq/unassigned_S1_L001_R4_001.fastq.gz \
           {params.rundir}/fastq/unassigned_S1_L002_R4_001.fastq.gz \
           {params.rundir}/fastq/unassigned_S1_L003_R4_001.fastq.gz \
           {params.rundir}/fastq/unassigned_S1_L004_R4_001.fastq.gz
        """

## After BCL conversion:
# R1 = genomic read 1
# R2 = Tn5 barcode 1, intervening sequence (15), sample index
# R3 (novaseq) = Cell barcode (16), intervening sequence (14), Tn5 barcode 2
# R3 (nextseq) = Tn5 barcode 2, intervening sequence (14), cell barcode (16)
# R4 = genomic read 2

# i7 intervening seq: CCGAGCCCACGAGAC
# i5 intervening seq: AGCAGCCGTCGCAG

rule demux:
    input:
        read1="{}/fastq/unassigned_R1.fastq.gz".format(config['reads']),
        read2="{}/fastq/unassigned_R2.fastq.gz".format(config['reads']),
        read3="{}/fastq/unassigned_R3.fastq.gz".format(config['reads']),
        read4="{}/fastq/unassigned_R4.fastq.gz".format(config['reads']),
        tn5_i5="{}_barcodes/tn5_i5.fa".format(config['name']),
        tn5_i7="{}_barcodes/tn5_i7.fa".format(config['name']),
        sampleindex="{}_barcodes/sampleindex.fa".format(config['name'])
    output: directory("{}/demux".format(config['outdir']))
    message: "Extract Tn5 barcode sequences"
    threads: 1
    shell:
        """
        python scripts/demux.py \
          --read1 {input.read1} \
          --read2 {input.read2} \
          --read3 {input.read3} \
          --read4 {input.read4} \
          --sampleindex {input.sampleindex} \
          --tn5_i5 {input.tn5_i5} \
          --tn5_i7 {input.tn5_i7} \
          --output {output}
        """

rule map:
    input:
        genome=config['reference'],
        dir=directory("{}/demux".format(config['outdir']))
    output: directory("{}/mapped".format(config['outdir']))
    message: "Map reads to genome"
    threads: 24
    shell:
        """
        
        # get all conditions
        
        # iterate over fastq file pairs and map
        
        bwa-mem2 mem {input.genome}  \
          | samtools sort - > aln.bam

        samtools index ...
        """

# rule fragments:
#     input:
#     output:
#     message: "Create fragment file"
#     threads: 12
#     shell:
#         """
#         sinto fragment -b {input}
#         sort -k1,1 -k2,2n 
#         bgzip -@ {threads}
#         tabix index -p bed 
#         """
