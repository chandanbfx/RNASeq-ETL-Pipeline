# Process RNA-seq data (align → count)
rule trim_fastq:
    input:
        r1="data/raw/{sample}_1.fastq.gz",
        r2="data/raw/{sample}_2.fastq.gz"
        
    output:
        trimR1="results/trimmed/{sample}_1.fastq",
        trimR2="results/trimmed/{sample}_2.fastq"
    conda:
        "../envs/rnaseq.yaml"
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.trimR1} -O {output.trimR2}"

rule align_bam:
    input:
        r1="results/trimmed/{sample}_1.fastq",
        r2="results/trimmed/{sample}_2.fastq",
        index=expand("resources/genomes/grch38/genome.{n}.ht2", n=range(1, 9))  # Updated path
    output:
        bam="results/aligned/{sample}.bam",
        log="results/aligned/{sample}.hisat2.log"
    params:
        genome="resources/genomes/grch38/genome"  # HISAT2 index prefix (no .ht2)
    threads: 4  # Recommended for HISAT2
    conda:
        "../envs/rnaseq.yaml"
    shell:
        "hisat2 -x {params.genome} -1 {input.r1} -2 {input.r2} 2> {output.log} | samtools view -Sb > {output.bam}"

rule samtools_stats:
    input: "results/aligned/{sample}.bam"
    output: "results/aligned/{sample}.bam.stats"
    shell: "samtools stats {input} > {output}"

rule organize_bams:
    input:
        expand("results/aligned/{sample}.bam", sample=config["samples"])
    output:
        expand("results/aligned/final/{sample}.sorted.bam", sample=config["samples"]),  # Individual files
        expand("results/aligned/final/{sample}.sorted.bam.bai", sample=config["samples"]),  # Index files
        dir_done = touch("results/aligned/final/.dir_complete")
    params:
        dir="results/aligned/final"
    threads: 4
    shell:
        """
        mkdir -p {params.dir}
        for bam in {input}; do
            sample=$(basename ${{bam}} .bam)
            samtools sort -@ {threads} -o {params.dir}/${{sample}}.sorted.bam ${{bam}}
            samtools index {params.dir}/${{sample}}.sorted.bam
        done
        """

rule generate_counts:
    input:
        bam="results/aligned/{sample}.bam",
        gtf=config["gtf"]
    output:
        "results/counts/{sample}_counts.txt"
    conda:
        "../envs/rnaseq.yaml"
    shell:
        "featureCounts -a {input.gtf} -o {output} {input.bam} -p"  # -p for paired-end

rule merge_counts:
    input:
        expand("results/counts/{sample}_counts.txt", sample=config["samples"])
    output:
        "results/counts/gene_counts_matrix.csv"
    conda:
        "../envs/rnaseq.yaml"
    script:
        "../scripts/merge_counts.py"  # Uses pandas to merge counts
