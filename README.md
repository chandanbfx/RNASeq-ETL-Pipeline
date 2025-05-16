# ETL-RNASeq

A modular and reproducible ETL (Extract-Transform-Load) pipeline for RNA-Seq data using Snakemake, AWS CLI, and Python.

---

## Overview

**ETL-RNASeq** automates the full RNA-seq data lifecycle:

- **Extract**: Download public RNA-seq datasets from NCBI GEO/SRA
- **Transform**: Preprocess via QC, trimming, alignment, quantification
- **Load**: Upload final outputs (counts, BAM files) to an S3 bucket

Built using Snakemake for reproducibility and scalability, this pipeline is suitable for batch processing on both local and cloud environments.

---

## Workflow Stages

1. **geo_etl.smk** – Download and extract FASTQ files using `fasterq-dump`
2. **qc.smk** – Perform quality control using FastQC and optional trimming
3. **genome_index.smk** – Build reference genome index (HISAT2/STAR)
4. **alignment.smk** – Align reads to reference genome
5. **upload_s3.smk** – Upload processed data to AWS S3 using `aws-cli`

---

## Directory Structure

```

ETL-RNASeq/
├── config      # Configuration files (YAML)
│   ├── config.yaml
│   └── geo_query.yaml
├── data        # Raw and intermediate data files
│   └── raw
├── envs        # Conda environment YAMLs
│   ├── aws_cli.yaml
│   ├── genome_index.yaml
│   ├── geo_etl.yaml
│   └── rnaseq.yaml
├── License
├── resources   # Reference genome and annotation files
│   ├── annotations 
│   └── genomes
├── results     # Final outputs (BAM, counts)
├── rules       # Modular Snakemake rule files
│   ├── alignment.smk
│   ├── genome_index.smk
│   ├── geo_etl.smk
│   ├── qc.smk
│   └── upload _s3.smk
├── scripts     # Custom Python scripts
│   ├── geo_fetch.py
│   ├── merge_counts.py
│   └── sra_download.py
└── Snakefile

```

---

## Requirements

- [Snakemake](https://snakemake.readthedocs.io)
- [AWS CLI](https://aws.amazon.com/cli/)
- [SRA Toolkit](https://github.com/ncbi/sra-tools)
- MultiQC
- HISAT2
- featureCounts
- Python 3.8+
- Conda (recommended)

---

## Getting Started

### Clone and Setup

```bash
git clone https://github.com/your-username/ETL-RNASeq.git
cd ETL-RNASeq
```

### Configure AWS

```bash
aws configure
```

### Modify the config

Update `config/config.yaml` with:

```yaml
samples: "config/samples.csv"
reference_genome: "resources/genome.fa"
annotation_gtf: "resources/genes.gtf"
s3_bucket: "your-s3-bucket-name"
```

### Run Pipeline

```bash
snakemake --cores 8 --use-conda
```

---

## Upload to S3

The rule `upload_s3.smk` uploads final results to your specified S3 bucket.

```bash
aws s3 sync results/ s3://your-s3-bucket/results/
```
Note: Ensure you have appropriate access and permissions to the specified S3 bucket (read/write) before running the upload step.
---

## Modular Snakemake Rules

Each pipeline step is in `rules/`:

- `geo_etl.smk`: Download from SRA
- `qc.smk`: MultiQC and trimming
- `genome_index.smk`: Build index
- `alignment.smk`: Map reads to genome
- `upload_s3.smk`: Sync output to S3

---

## Custom Scripts

All custom processing logic (e.g., metadata parsing, file renaming, post-processing) is in `scripts/`.

---

## Conda Environments

Each module uses its own environment defined in `envs/`. Snakemake ensures reproducibility by using:

```bash
--use-conda
```

---

## To Do

- Add support for paired-end/single-end auto-detection
- Integrate MultiQC reports
- Add Nextflow support (optional)

---

## License

This project is licensed under the MIT License - see the `License` file for details.

---

## Acknowledgements

- NCBI GEO / SRA
- Bioconda community
- Snakemake developers
- MultiQC, HISAT2, featurecounts, AWS
