import sys


if len(config) == 0:
    print("Please supply a configuration using --configfile", file=sys.stderr)
    sys.exit()

rule all:
    input: "{}/mapped/aln.bam".format(config['outdir'])
    # input: directory("{}/fragments".format(config['outdir']))


rule index_genome:
    message: "Building bwa-mem2 genome index"
    output: "{}.0123".format(config['reference'])
    params:
        genome=config['reference']
    threads: 1
    shell:
        """
        bwa-mem2 index {params.genome}
        """

rule splitcode_extract:
    input:
        read1="{}/{}.fastq.gz".format(config['reads'], config['read1']),
        index1="{}/{}.fastq.gz".format(config['reads'], config['index1']),
        index2="{}/{}.fastq.gz".format(config['reads'], config['index2']),
        read2="{}/{}.fastq.gz".format(config['reads'], config['read2'])
    output:
        r1="{}/read1.fastq.gz".format(config['outdir']),
        r2="{}/read2.fastq.gz".format(config['outdir'])
    params:
        cfg=config['barcodes'],
        cb=config['cb'],
        keep=config['keep_group']
    message: "Extract barcodes from index reads"
    threads: 12
    shell:
        """
        splitcode \
          -t {threads} \
          --config={params.cfg} \
          --seq-names \
          --nFastqs=4 \
          --trim-only \
          --extract={params.cb} \
          --x-names \
          --no-x-out \
          --output={output.r1},{output.r2} \
          --gzip \
          --mod-names \
          --keep-grp=<(echo {params.keep}) \
          --select=2,3 \
          --sam-tags=TB:Z:,CB:Z: \
          {input.index1} {input.index2} {input.read1} {input.read2}
        """

rule map:
    input:
        r1="{}/read1.fastq.gz".format(config['outdir']),
        r2="{}/read2.fastq.gz".format(config['outdir']),
        idx="{}.0123".format(config['reference'])
    output: "{}/mapped/aln.bam".format(config['outdir'])
    message: "Map reads to genome"
    params:
        genome=config['reference']
    threads: 24
    shell:
        """
        bwa-mem2 mem \
          -t {threads} \
          -C \
          {params.genome} \
          {input.r1} \
          {input.r2} \
            | samtools sort -@ {threads} -O bam - \
            > {output}
        samtools index -@ {threads} {output}
        """

rule fragments:
    input:
        bam="{}/mapped/aln.bam".format(config['outdir']),
        samples=config['samples']
    output: directory("{}/fragments".format(config['outdir']))
    message: "Create fragment files"
    threads: 12
    shell:
        """
        mkdir {output}
        sinto fragments \
          -p {threads} \
          -b {input.bam} \
          -f {output} \
          --header \
          --split \
          --collapse_within
        
        cd {output}
        for fragfile in *.tsv; do
          fname=(${fragfile//.tsv/ })
          sinto sort -i $fragfile -o $fname.bed
          bgzip -@ {threads} $fname.bed
          tabix -p bed $fname.bed.gz
        done
        """
