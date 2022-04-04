from snakemake.utils import validate
import pandas as pd

configfile: "config/config_GE.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv("config/samples_downsampling.tsv", sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="../schemas/samples.schema.yaml")

def get_refsample(refsample):
    return "results/intermediate/{ID}/FFT_table/transcriptanno-{SAMPLE}-FFT_table.{GENOME}.tsv".format(ID = samples["ID"].loc[samples["sample"] == refsample].values[0], SAMPLE = refsample, GENOME = samples["genome_build"].loc[samples["sample"] == refsample].values[0])


rule join: 
    input:
        transcriptAnno= lambda wc: config[wc.GENOME]["transcriptAnno"],
        proteinAtlas= expand("resources/protein_atlas/RNAtable{SOURCE}.tsv.gz",
                            SOURCE=config["proteinAtlas"]),
    output:
        filteredTranscriptAnno="results/intermediate/transcriptAnno/transcriptAnno-{GENOME}.103.filtered.tsv.gz"
    conda: "../envs/cfDNA.yml"
    shell:
        """
        set +o pipefail;
        (zcat {input.transcriptAnno} | \
        head -n 1; join -t"$(echo -e "\\t")" \
        <(\zcat {input.transcriptAnno}| tail -n +2 | sort -k1,1 ) \
        <(zcat {input.proteinAtlas} | tail -n +2 |cut -f 1 | sort )) | \
        gzip -c > {output.filteredTranscriptAnno}
        """

rule prep:
    input:
        transcriptAnno="results/intermediate/transcriptAnno/transcriptAnno-{GENOME}.103.filtered.tsv.gz"
    output:
        body="results/intermediate/transcriptAnno/transcriptAnno-{GENOME}.103.body.bed.gz"
    conda: "../envs/cfDNA.yml"
    shell:
        """
        zcat {input.transcriptAnno} | tail -n +2 | \
        awk 'BEGIN{{ FS="\\t"; OFS="\\t" }}{{ if ($5 == "+") {{ print $2,$3-1-1000,$3-1+10000,$1,0,$5 }} else {{ print $2,$4-1-10000,$4-1+1000,$1,0,$5 }} }}'| \
        gzip -c > {output.body}
        """

rule generate_random_background:
    input:
        region=(
            "results/intermediate/transcriptAnno/transcriptAnno-{GENOME}.103.body.bed.gz"
        ),
        genome=lambda wildcards: config[wildcards.GENOME]["genome_autosomes"], 
        gap=lambda wildcards: config[wildcards.GENOME]["UCSC_gap"],
    output:
        "results/intermediate/transcriptAnno/transcriptAnno_background-{GENOME}.103.body.bed.gz"
    params:
        length=10000, #lambda wildcards, input: get_length(input.region)
    conda:
        "../envs/background.yml"
    shell:
        """
        bedtools random -n 1000 -l {params.length} -g {input.genome} | \
        bedtools shuffle -i stdin -g {input.genome} -excl {input.gap} -noOverlapping | \
        gzip -c > {output}
        """

rule extract_counts:
    input:
        target=(
            "results/intermediate/transcriptAnno/transcriptAnno-{GENOME}.103.body.bed.gz"
        ),
        BAMFILE="results/downsampling/{GENOME}/{SAMPLE}.{COV}x.bam",
    output:
        WPS="results/intermediate/{ID}/table/transcriptanno_{SAMPLE}_WPS.{COV}x.{GENOME}.csv.gz",
        COV="results/intermediate/{ID}/table/transcriptanno_{SAMPLE}_COV.{COV}x.{GENOME}.csv.gz",
        STARTS="results/intermediate/{ID}/table/transcriptanno_{SAMPLE}_STARTS.{COV}x.{GENOME}.csv.gz",
    params:
        minRL=config["minRL"],
        maxRL=config["maxRL"],
        out_pre="results/intermediate/{ID}/table/transcriptanno_{SAMPLE}_%s.{COV}x.{GENOME}.csv.gz",
    conda:
        "../envs/cfDNA.yml"
    shell:
        """
        workflow/scripts/expression_analysis/extractFromBAM_RegionBed_WPS_Cov.py \
        --minInsert={params.minRL} \
        --maxInsert={params.maxRL} \
        -i {input.target} \
        -o {params.out_pre} {input.BAMFILE}
        """

rule extract_counts_background:
    input:
        background="results/intermediate/transcriptAnno/transcriptAnno_background-{GENOME}.103.body.bed.gz",
        BAMFILE="results/downsampling/{GENOME}/{SAMPLE}.{COV}x.bam",
    output:
        WPS="results/intermediate/{ID}/background_region/table/transcriptanno_{SAMPLE}_WPS_background.{COV}x.{GENOME}.csv.gz",
        COV="results/intermediate/{ID}/background_region/table/transcriptanno_{SAMPLE}_COV_background.{COV}x.{GENOME}.csv.gz",
        STARTS="results/intermediate/{ID}/background_region/table/transcriptanno_{SAMPLE}_STARTS_background.{COV}x.{GENOME}.csv.gz",
    params:
        minRL=config["minRL"],
        maxRL=config["maxRL"],
        out_pre="results/intermediate/{ID}/background_region/table/transcriptanno_{SAMPLE}_%s_background.{COV}x.{GENOME}.csv.gz",
    conda:
        "../envs/cfDNA.yml"
    shell:
        """
        workflow/scripts/expression_analysis/extractFromBAM_RegionBed_WPS_Cov.py \
        --minInsert={params.minRL} \
        --maxInsert={params.maxRL} \
        -i {input.background} \
        -o {params.out_pre} {input.BAMFILE}
        """


rule normalize_WPS:
    input:
        target_WPS="results/intermediate/{ID}/table/transcriptanno_{SAMPLE}_WPS.{COV}x.{GENOME}.csv.gz",
        background_WPS="results/intermediate/{ID}/background_region/table/transcriptanno_{SAMPLE}_WPS_background.{COV}x.{GENOME}.csv.gz",
        target_COV="results/intermediate/{ID}/table/transcriptanno_{SAMPLE}_COV.{COV}x.{GENOME}.csv.gz",
        background_COV="results/intermediate/{ID}/background_region/table/transcriptanno_{SAMPLE}_COV_background.{COV}x.{GENOME}.csv.gz",
    output:
        output_WPS="results/intermediate/{ID}/table/transcriptanno_{SAMPLE}_WPS_normalized.{COV}x.{GENOME}.tsv.gz",
        output_COV="results/intermediate/{ID}/table/transcriptanno_{SAMPLE}_COV_normalized.{COV}x.{GENOME}.tsv.gz",
    conda:
        "../envs/overlays.yml"
    script:
        """../scripts/expression_analysis/normalize.py"""


rule FFT_table:
    input:
        normalized_WPS="results/intermediate/{ID}/table/transcriptanno_{SAMPLE}_WPS_normalized.{COV}x.{GENOME}.tsv.gz",
    output:
        FFT_table="results/intermediate/{ID}/FFT_table/transcriptanno-{SAMPLE}-FFT_table.{COV}x.{GENOME}.tsv.gz",
    conda:
        "../envs/overlays.yml"
    script:
        """../scripts/expression_analysis/fft_table.py"""

