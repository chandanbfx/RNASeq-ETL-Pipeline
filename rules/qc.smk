# Aggregate QC reports
rule multiqc:
    input:
        expand("results/trimmed/{sample}_1.fastq", sample=config["samples"]),  # From trim_fastq
        expand("results/aligned/{sample}.bam", sample=config["samples"])       # From align_bam
    output:
        "results/multiqc_report.html"
    conda:
        "../envs/rnaseq.yaml"
    shell:
        "multiqc results/ -o results/"
