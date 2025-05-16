rule geo_fetch:
    output:
        "results/geo/geo_metadata.csv"
    conda:
        "../envs/geo_etl.yaml"
    script:
        "../scripts/geo_fetch.py"

rule load_metadata:
    input:
        "results/geo/geo_metadata.csv"
    output:
        "results/metadata.db"
    conda:
        "../envs/geo_etl.yaml"
    shell:
        """
        sqlite3 {output} <<EOF
        .mode csv
        .import {input} samples
        EOF
        """
rule download_sra:
    input:
        "results/geo/geo_metadata.csv"
    output:
        # For paired-end data
        temp("data/raw/{sample}_1.fastq.gz"),
        temp("data/raw/{sample}_2.fastq.gz")
    params:
        sra_id=lambda wildcards: config["sra_ids"][wildcards.sample]
    conda:
        "../envs/geo_etl.yaml"
    shell:
        """
        fasterq-dump {params.sra_id} --outdir data/raw/ --split-files
        gzip data/raw/{params.sra_id}_1.fastq
        gzip data/raw/{params.sra_id}_2.fastq
        """
