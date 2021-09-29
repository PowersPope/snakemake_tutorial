configfile: "config.yaml"

rule all:
    input:
        "plots/quals.svg"

def get_bwa_map_input_fastqs(wildcards):
    return config["samples"][wildcards.sample]


rule bwa_map:
    input:
        "data/genome.fa",
        get_bwa_map_input_fastqs
    output:
        "mapped_reads/{sample}.bam"
    # Add in paramter features to the snakefile
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}"

    #Create a log of the output
    log:
        "logs/bwa_mem/{sample}.log"
    # Threads can be specified, so if I give 10 cores, only 8 will be used for this.
    # The other two will be used on the other rules to run stuff in parallel.
    threads: 8
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {input} | samtools view -Sb - > {output}) 2> {log}t"

rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"

rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"

rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
    output:
        "calls/all.vcf"
    params:
        newg=config['prior_mutation_rate']
    log:
        "logs/bcftools_call/{input.bam}.log"
    shell:
        "(samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -P '{params.newg}' -mv - > {output}) 2> {log}"
        
rule plot_quals:
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    script:
        "scripts/plot-quals.py"