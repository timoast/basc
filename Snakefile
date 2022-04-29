import pandas as pd

configfile: "config.example"
samples = pd.read_csv(config["samples"], sep="\t")
unique_groups = list(set(samples['sample_name'].to_list()))

rule all:
    input: directory("{}/fragments".format(config['outdir']))

rule create_fasta:
    output:
        "{}_barcodes/tn5_i5.fa".format(config['name']),
        "{}_barcodes/tn5_i7.fa".format(config['name']),
        "{}_barcodes/sampleindex.fa".format(config['name'])
    params:
        samples=config['samples'],
        revcomp_i5=config['reverse_complement_i5'],
        revcomp_i7=config['reverse_complement_i7'],
        samplename=config['name']
    message: "Generate barcode FASTA files"
    shell:
        """
        Rscript scripts/generate_fasta.R \
            {params.samples} \
            {params.samplename} \
            {params.revcomp_i5} \
            {params.revcomp_i7}
        """

rule demux:
    input:
        read1="{}/{}.fastq.gz".format(config['reads'], config['read1']),
        read2="{}/{}.fastq.gz".format(config['reads'], config['read2']),
        read3="{}/{}.fastq.gz".format(config['reads'], config['read3']),
        read4="{}/{}.fastq.gz".format(config['reads'], config['read4']),
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
          --min_barcode_len 14 \
          --output {output}
        """

rule map:
    input: dir=directory("{}/demux".format(config['outdir']))
    output: directory("{}/mapped".format(config['outdir']))
    message: "Map reads to genome"
    params:
        group=unique_groups,
        genome=config['reference']
    threads: 24
    shell:
        """
        mkdir {output}
        for i in {params.group}; do
            cd {input}/$i
            mkdir {output}/$i
            for R1 in $(ls -d *.R1.fastq); do
                fname=${{R1%.R1.fastq}}
                R2=$fname.R2.fastq
                bwa-mem2 mem -t {threads} -C {params.genome} $R1 $R2 \
                | samtools sort -@ {threads} -O bam - \
                > "{output}/"$i"/"$fname".bam"
                samtools index -@ {threads} "{output}/"$i"/"$fname".bam"
            done
            cd -
        done
        """

rule fragments:
    input: directory("{}/mapped".format(config['outdir']))
    output: directory("{}/fragments".format(config['outdir']))
    message: "Create fragment file"
    params: unique_groups
    threads: 12
    shell:
        """
        mkdir {output}
        for i in {params[0]}; do
            cd {input}/$i
            mkdir {output}/$i
            for bam in $(ls -d *.bam); do
                fname=${{bam%.bam}}
                sinto fragments -p {threads} -b $bam -f {output}/$i/$fname.tmp
                sort -k1,1 -k2,2n {output}/$i/$fname.tmp > {output}/$i/$fname.tsv
                bgzip -@ {threads} {output}/$i/$fname.tsv
                tabix -p bed {output}/$i/$fname.tsv.gz
                rm {output}/$i/$fname.tmp
            done
            cd -
        done
        """
