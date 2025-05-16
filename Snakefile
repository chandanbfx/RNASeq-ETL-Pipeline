# Master workflow to chain GEO ETL + RNA-seq processing
configfile: "config/config.yaml"
configfile: "config/geo_query.yaml"

include: "rules/geo_etl.smk"
include: "rules/alignment.smk"
include: "rules/qc.smk"
include: "rules/genome_index.smk"
include: "rules/upload _s3.smk"

rule all:
    input:
        "results/counts/gene_counts_matrix.csv",  # Merged counts
        "results/metadata.db",                   # SQLite DB
        "results/aligned/final/.dir_complete",	# Sorted BAMs
        "results/multiqc_report.html",            # QC
        ".s3_upload_complete"
