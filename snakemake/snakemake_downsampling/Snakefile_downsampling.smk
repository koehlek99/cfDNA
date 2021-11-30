SAMPLES = ["BH01","IH01","IH02","IC15","IC17","IC35","IC37"]

rule all: 
    input: 
        expand("BAMs/37/{SAMPLE}_10x.bam.bai", 
        SAMPLE = SAMPLES)
        
        
rule calculate_coverage:
    input: 
        BAMFILE="/fast/groups/ag_kircher/cfDNA-data/Cell_2016/BAMs/{SAMPLE}.bam",
    output: 
        coverageStats="{SAMPLE}.mosdepth.summary.txt"
    shell: 
        """
        mosdepth -n {wildcards.SAMPLE} {input.BAMFILE}
        """


rule downsampling:
    input: 
        coverageStats="{SAMPLE}.mosdepth.summary.txt",
        BAMFILE="/fast/groups/ag_kircher/cfDNA-data/Cell_2016/BAMs/{SAMPLE}.bam",
    output: 
        normalizedBAM="BAMs/37/{SAMPLE}_10x.bam"
    run:
        import pandas as pd 
        stats = pd.read_csv(input.coverageStats, sep="\t", header=0)
        total_cov = stats[stats["chrom"]=="total"]["mean"].item()
        print(input.coverageStats, ": ",total_cov)
        factor = 10/total_cov
        print(factor)
        shell("samtools view -s {factor} -b {input.BAMFILE} > {output.normalizedBAM}")

rule build_index:
    input:
        normalizedBAM="BAMs/37/{SAMPLE}_10x.bam"
    output: 
        indexFile="BAMs/37/{SAMPLE}_10x.bam.bai"
    shell:
        """ 
        samtools index -b {input.normalizedBAM}
        """