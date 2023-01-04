import joblib
import pandas as pd
import numpy as np
import sys

rf = joblib.load(snakemake.input['MODEL'])
features = pd.read_csv(snakemake.input['FEATURES'], sep = "\t", index_col = 0, header = 0).fillna(0)

features = features[rf.feature_names_in_]

predictions = pd.DataFrame(index = features.index)
predictions['predictedExpression'] = rf.predict(features)

predictions.to_csv(snakemake.output['PREDICTIONS'])