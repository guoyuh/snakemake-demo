config:"config.yaml"

bwa="/mnt/db/software/bwa/bwa"
samtools="/mnt/db/software/samtools-1.9/samtools"
bcftools="/mnt/db/software/bcftools/bcftools-v1.9/bcftools"

rule all:
    input:
        "calls/all.vcf"


def get_bwa_map_input_fastqs(wildcards):
    return config["samples"][wildcards.sample]

# rule bwa:
#     input:
#         "data/genome.fa",
#         "data/samples/{sample}.fastq"
#     output:
#         temp("mapped/{sample}.bam")
#     conda:
#         "envs/mapping.yaml"
#     threads: 8
#     shell:
#         "{bwa} mem -t {threads} {input} | {samtools} view -Sb - > {output}"



rule bwa:
    input:
        "data/genome.fa",
        get_bwa_map_input_fastqs
    output:
        temp("mapped/{sample}.bam")
    conda:
        "envs/mapping.yaml"
    threads: 8
    shell:
        "{bwa} mem -t {threads} {input} | {samtools} view -Sb - > {output}"


rule sort:
    input:
        "mapped/{sample}.bam"
    output:
        "mapped/{sample}.sorted.bam"
    conda:
        "envs/mapping.yaml"
    shell:
        "{samtools} sort -o {output} {input}"



rule call:
    input:
        fa="data/genome.fa",
        bam=expand("mapped/{sample}.sorted.bam", sample=config["samples"])
    output:
        "calls/all.vcf"
    conda:
        "envs/calling.yaml"
    shell:
        "{samtools} mpileup -g -f {input.fa} {input.bam} | "
        "{bcftools} call -mv - > {output}"

rule stats:
    input:
        "calls/all.vcf"
    output:
        report("plots/quals.svg", caption="report/calling.rst")
    conda:
        "envs/stats.yaml"
    script:
        "scripts/plot-quals.py"