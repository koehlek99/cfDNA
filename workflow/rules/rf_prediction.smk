rule extract_features:
    input:
        WPS = "results/intermediate/table/transcriptanno_{SAMPLE}_WPS_normalized.{COV}x.{GENOME}.tsv.gz",
        COVERAGE = "results/intermediate/table/transcriptanno_{SAMPLE}_COV_normalized.{COV}x.{GENOME}.tsv.gz",
        FFT = "results/intermediate/FFT_table/transcriptanno-{SAMPLE}-FFT_table.{COV}x.{GENOME}.tsv.gz",
        GC = "resources/rf/GC_content.csv",
        EXPRESSIONATLAS = "resources/rf/RNAtableExtended.tsv.gz",
    output:
        FEATURES = "results/rf/{GENOME}/features_{SAMPLE}.{COV}x.csv.gz"
    conda: 
        "../envs/rf.yml"
    script:
        """../scripts/rf/extract_features.py"""


rule predict_expression:
    input:
        FEATURES = "results/rf/{GENOME}/features_{SAMPLE}.{COV}x.csv.gz",
        MODEL = "resources/rf/BH01IH01_originalCoverage_{MODEL}.joblib",
    output:
        PREDICTIONS = "results/rf/{GENOME}/predictions_{SAMPLE}_{MODEL}.{COV}x.csv.gz"
    params: 
        model = "{MODEL}"
    conda: 
        "../envs/rf.yml"
    script:
        """../scripts/rf/predict_expression.py"""