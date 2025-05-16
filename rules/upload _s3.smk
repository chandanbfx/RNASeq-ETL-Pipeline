from datetime import datetime

def get_s3_path(wildcards):
    """Generate S3 path with current date"""
    return f"s3://project.etl.rnaseq/project_{datetime.now().strftime('%Y%m%d')}"

rule upload_to_s3:
    input:
        "results/counts/gene_counts_matrix.csv",
        "results/metadata.db",
        "results/multiqc_report.html",
        "results/aligned/final/.dir_complete"
    output:
        touch(".s3_upload_complete")
    params:
        s3_path=get_s3_path,
        exclude="*.tmp *.log"
    conda:
        "../envs/aws_cli.yaml"
    log:
        "logs/s3_upload.log"
    shell:
        """
        # Sync with checksum verification and user metadata
        aws s3 sync results/ {params.s3_path} --profile rnaseq-s3-user \
            --exclude "{params.exclude}" \
            --checksum \
            --metadata "project={config[project_id]}" \
            --no-progress \
            >> {log} 2>&1

        # Verify upload completeness
        aws s3 ls {params.s3_path}/ > /dev/null || (echo "S3 upload failed!"; exit 1)
        """
