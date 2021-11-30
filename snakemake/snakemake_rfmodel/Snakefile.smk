#SAMPLES = ["BH01","IH01","IH02","IC15","IC17","IC20","IC35","IC37"]
SAMPLES = ["NPH001","NPH002","NPH003","NPH004","NPH005","P147_1","P147_3","P148_1","P148_3","P40_1","P40_2","IC20","C2_6","C2_7"]


rule all: 
    input: 
        expand("results/predictions/37/predictions_{SAMPLE}.csv", 
        SAMPLE = SAMPLES),
        expand("results/features/37/features_{SAMPLE}.csv", 
        SAMPLE = SAMPLES)


rule extract_features: 
    input: 
        WPS = "resources/WPS/transcriptanno_{SAMPLE}_WPS_normalized.GRCh37.tsv",
        COVERAGE = "resources/COV/transcriptanno_{SAMPLE}_COV_normalized.GRCh37.tsv",
        FFT = "resources/FFT/transcriptanno-{SAMPLE}-FFT_table.GRCh37.tsv",
        GENES = "resources/gene_ids.csv",
        MEAN_EXP = "resources/mean_expression.csv",
        GC = "resources/GC_content.csv",
        MONOCYTES = "resources/RNAtableExtended.tsv.gz",
    output:
        FEATURES = "results/features/37/features_{SAMPLE}.csv"
    script:
        """scripts/extract_features.py"""


rule predict_expression: 
    input: 
        FEATURES = "results/features/37/features_{SAMPLE}.csv",
        MODEL = "resources/rf_model_37_BH01IH01.joblib",
    output:
        PREDICTIONS = "results/predictions/37/predictions_{SAMPLE}.csv"
    script:
        """scripts/predict_expression.py"""