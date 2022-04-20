from snakemake.utils import validate
import pandas as pd

configfile: "config/config_downsampling.yaml"
#validate(config, schema = "../schemas/config.schema.downsampling.yaml"")
samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")
      
        
rule calculate_coverage:
    input: 
        BAMFILE=lambda wildcards: samples["path"][wildcards.SAMPLE],
    output: 
        coverageStats="results/downsampling/{SAMPLE}.{GENOME}.mosdepth.summary.txt"
    conda:
        "../envs/downsampling.yml"
    shell: 
        """
        mosdepth -n {wildcards.SAMPLE} {input.BAMFILE}
        mv {wildcards.SAMPLE}.mosdepth.summary.txt results/downsampling/{wildcards.SAMPLE}.{wildcards.GENOME}.mosdepth.summary.txt
        mv {wildcards.SAMPLE}.mosdepth.global.dist.txt  results/downsampling/{wildcards.SAMPLE}.{wildcards.GENOME}.mosdepth.global.dist.txt
        """


rule downsampling:
    input: 
        coverageStats="results/downsampling/{SAMPLE}.{GENOME}.mosdepth.summary.txt",
        BAMFILE=lambda wildcards: samples["path"][wildcards.SAMPLE],
    output: 
        normalizedBAM="results/downsampling/{GENOME}/{SAMPLE}.{COV}x.bam"
    conda:
        "../envs/downsampling.yml"
    shell:
        """
        factor=$(tail -n 1 results/downsampling/{wildcards.SAMPLE}.{wildcards.GENOME}.mosdepth.summary.txt | awk '{{print {wildcards.COV}/$4}}')
        echo "${{factor}}"
        samtools view -s $factor -b {input.BAMFILE} > {output.normalizedBAM}
        """

rule build_index:
    input:
        normalizedBAM="results/downsampling/{GENOME}/{SAMPLE}.{COV}x.bam"
    output: 
        indexFile="results/downsampling/{GENOME}/{SAMPLE}.{COV}x.bam.bai"
    conda:
        "../envs/downsampling.yml"
    shell:
        """ 
        samtools index -b {input.normalizedBAM}
        """