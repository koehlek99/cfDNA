import joblib
import pandas as pd
import numpy as np

rf = joblib.load(snakemake.input['MODEL'])

feature_names_37 = ['mean_expression', 'mean_cov_5k', 'GC_promotor', 'mean_cov','mean_cov_2k', 'mean_wps_5k', 'mean_wps_2k', 'mean_wps', '197','195', 'ndr_mean_cov', '200', '202', 'mean_cov_body1kb','mean_wps_body1kb']

feature_names_38 = ['mean_expression', 'GC_promotor', 'mean_cov_5k', 'mean_wps_5k','198', '201', 'amp2', 'amp1', '195', '205', '166', '190','ndr_mean_cov', '192', 'ndr_width2']

feature_names = feature_names_37

features = pd.read_csv(snakemake.input['FEATURES'], sep = "\t", index_col = 0, header = 0) 

features_reduced = features.loc[:, feature_names]
features_reduced['monocytes'] = features['monocytes']
features_reduced = features_reduced.fillna(0)

features_reduced.loc[:,'predicted'] = rf.predict(np.array(features_reduced.drop(["monocytes"], axis = 1)))

features_reduced.to_csv(snakemake.output['PREDICTIONS'])
