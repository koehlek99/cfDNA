import joblib
import pandas as pd
import numpy as np

rf = joblib.load(snakemake.input['MODEL'])
features = pd.read_csv(snakemake.input['FEATURES'], sep = "\t", index_col = 0, header = 0).fillna(0)

features.loc[:,'predicted'] = rf.predict(np.array(features.drop(["monocytes"], axis = 1)))

features.to_csv(snakemake.output['PREDICTIONS'])