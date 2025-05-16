rule download_genome:
    output:
        "resources/genomes/GRCh38.fa"
    params:
        url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
    conda:
        "../envs/genome_index.yaml"
    shell:
        "wget -O - {params.url} | gunzip > {output}"

rule download_gtf:
    output:
        "resources/annotations/genes.gtf"
    params:
	    url="https://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz"
    shell:
	    "wget -O - {params.url} | gunzip > {output}"

rule download_hisat2_index:
    """
    Downloads and extracts a pre-built HISAT2 index for GRCh38.
    Source: https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
    Note: The tarball extracts into a 'grch38/' subdirectory.
    """
    output:
        expand("resources/genomes/grch38/genome.{n}.ht2", n=range(1, 9))
    params:
        url="https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz",
        dir="resources/genomes",
        extracted_dir="resources/genomes/grch38"# Where files are actually extracted
    shell:
        """
        mkdir -p {params.dir}
        wget -O {params.dir}/grch38_genome.tar.gz {params.url}
        tar -xzf {params.dir}/grch38_genome.tar.gz -C {params.dir}
        rm {params.dir}/grch38_genome.tar.gz  # Cleanup
        """
