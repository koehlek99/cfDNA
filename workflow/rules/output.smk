configfile: "config/config_downsampling.yaml"

def get_final_output():
    final_output = []
    final_output.extend(
        expand("results/downsampling/{SAMPLE}.{GENOME}.mosdepth.summary.txt",
            SAMPLE=samples["sample"],
            COV=config["coverage"],
            GENOME=samples["genome_build"],
        ),
    ),
    final_output.extend(
        expand("results/downsampling/{GENOME}/{SAMPLE}.{COV}x.bam",
            SAMPLE=samples["sample"],
            COV=config["coverage"],
            GENOME=samples["genome_build"],
        ), 
    ),
    final_output.extend(
        expand("results/downsampling/{GENOME}/{SAMPLE}.{COV}x.bam.bai",
            SAMPLE=samples["sample"],
            COV=config["coverage"],
            GENOME=samples["genome_build"],
        ),
    ),
    final_output.extend(
        expand(
            "results/intermediate/table/transcriptanno_{SAMPLE}_WPS_normalized.{COV}x.{GENOME}.tsv.gz",
            zip,
            SAMPLE=samples["sample"],
            COV=config["coverage"],
            GENOME=samples["genome_build"],
        ),
    ),
    final_output.extend(
        expand("results/intermediate/FFT_table/transcriptanno-{SAMPLE}-FFT_table.{COV}x.{GENOME}.tsv.gz",
            zip,
            SAMPLE=samples["sample"],
            COV=config["coverage"],
            GENOME=samples["genome_build"],
        ),
    ),
    final_output.extend(
        expand("results/intermediate/table/transcriptanno_{SAMPLE}_COV_normalized.{COV}x.{GENOME}.tsv.gz",
            zip,
            SAMPLE=samples["sample"],
            COV=config["coverage"],
            GENOME=samples["genome_build"],
        ),
    ), 
    final_output.extend(
        expand("results/rf/{GENOME}/features_{SAMPLE}.{COV}x.csv.gz",
            SAMPLE=samples["sample"],
            COV=config["coverage"],
            GENOME=samples["genome_build"],
        ),
    ),
    final_output.extend(
        expand("results/rf/{GENOME}/predictions_{SAMPLE}.{COV}x.csv.gz",
            SAMPLE=samples["sample"],
            COV=config["coverage"],
            GENOME=samples["genome_build"],
        ),
    ),



    return final_output