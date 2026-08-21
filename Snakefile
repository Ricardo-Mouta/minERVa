import os
import glob
from pathlib import Path

configfile: "config.yaml"

GENOMES = []
if os.path.exists(config["genomes_list"]):
    with open(config["genomes_list"]) as f:
        GENOMES = [line.strip() for line in f if line.strip()]

GBK_DIR = config["gbk_dir"]
GBKS = [Path(f).stem for f in glob.glob(os.path.join(GBK_DIR, "*.gbk"))]
TARGETS = config["targets"]

rule all:
    input:
        "results/summary_hits.tsv"

rule download_genome:
    output:
        fasta = temp("resources/genomes/{genome}.fna")
    log:
        "logs/download_genome/{genome}.log"
    conda:
        "envs/datasets.yaml"
    retries: 3
    shell:
        """
        datasets download genome accession {wildcards.genome} --include genome --filename resources/genomes/{wildcards.genome}.zip > {log} 2>&1
        unzip -p resources/genomes/{wildcards.genome}.zip ncbi_dataset/data/{wildcards.genome}/*.fna > {output.fasta} 2>> {log}
        rm resources/genomes/{wildcards.genome}.zip
        """

rule makeblastdb:
    input:
        fasta = "resources/genomes/{genome}.fna"
    output:
        db_dir = temp(directory("resources/blast_dbs/{genome}"))
    log:
        "logs/makeblastdb/{genome}.log"
    conda:
        "envs/blast.yaml"
    shell:
        """
        mkdir -p {output.db_dir}
        makeblastdb -in {input.fasta} -dbtype nucl -out {output.db_dir}/{wildcards.genome} -title {wildcards.genome} > {log} 2>&1
        """

rule extract_queries:
    input:
        gbk = os.path.join(GBK_DIR, "{gbk}.gbk")
    output:
        fasta = "results/queries/{gbk}_{target}.fasta"
    log:
        "logs/extract_queries/{gbk}_{target}.log"
    conda:
        "envs/python.yaml"
    shell:
        "python scripts/extract_queries.py --gbk {input.gbk} --target {wildcards.target} --out {output.fasta} > {log} 2>&1"

rule local_blastn:
    input:
        fasta = "results/queries/{gbk}_{target}.fasta",
        db_dir = "resources/blast_dbs/{genome}"
    output:
        tsv = "results/blast/{genome}_{gbk}_{target}.tsv"
    log:
        "logs/blast/{genome}_{gbk}_{target}.log"
    conda:
        "envs/blast.yaml"
    params:
        db = "resources/blast_dbs/{genome}/{genome}",
        hspcov = config["blast_params"]["hspcoverage"],
        ident = config["blast_params"]["identity"],
        evalue = config["blast_params"]["evalue"]
    shell:
        """
        # Se o FASTA de query estiver vazio, criamos um output vazio para não quebrar o pipeline
        if [ ! -s {input.fasta} ]; then
            echo "Aviso: {input.fasta} vazio, ignorando BLAST." > {log} 2>&1
            touch {output.tsv}
        else
            echo "Rodando BLAST local para {wildcards.target} contra {wildcards.genome}..." > {log} 2>&1
            blastn -db {params.db} \\
                   -query {input.fasta} \\
                   -out {output.tsv} \\
                   -outfmt "6 qseqid sseqid sstart send sstrand pident qcovs" \\
                   -qcov_hsp_perc {params.hspcov} \\
                   -evalue {params.evalue} >> {log} 2>&1 || touch {output.tsv}
        fi
        """

rule summarize_hits:
    input:
        blasts = expand("results/blast/{genome}_{gbk}_{target}.tsv", genome=GENOMES, gbk=GBKS, target=TARGETS)
    output:
        summary = "results/summary_hits.tsv"
    log:
        "logs/summarize_hits.log"
    conda:
        "envs/python.yaml"
    params:
        ltr_range = config["blast_params"]["LTR_range"]
    shell:
        "python scripts/parse_blast.py --blast_files {input.blasts} --out {output.summary} --ltr_range {params.ltr_range} > {log} 2>&1"
