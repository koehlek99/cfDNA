rule extract_features:
    input:
        WPS = "results/intermediate/table/transcriptanno_{SAMPLE}_WPS_normalized.{COV}x.{GENOME}.tsv.gz",
        COVERAGE = "results/intermediate/table/transcriptanno_{SAMPLE}_COV_normalized.{COV}x.{GENOME}.tsv.gz",
        FFT = "results/intermediate/FFT_table/transcriptanno-{SAMPLE}-FFT_table.{COV}x.{GENOME}.tsv.gz",
        GENES = "resources/rf/gene_ids.csv",
        MEAN_EXP = "resources/rf/mean_expression.csv",
        GC = "resources/rf/GC_content.csv",
        MONOCYTES = "resources/rf/RNAtableExtended.tsv.gz",
    output:
        FEATURES = "results/rf/{GENOME}/features_{SAMPLE}.{COV}x.csv.gz"
    conda: 
        "../envs/rf.yml"
    script:
        """../scripts/rf/extract_features.py"""


rule predict_expression:
    input:
        FEATURES = "results/rf/{GENOME}/features_{SAMPLE}.{COV}x.csv.gz",
        MODEL = "resources/rf/rf_model.joblib",
    output:
        PREDICTIONS = "results/rf/{GENOME}/predictions_{SAMPLE}.{COV}x.csv.gz"
    conda: 
        "../envs/rf.yml"
    script:
        """../scripts/rf/predict_expression.py"""
        