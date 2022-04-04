import gzip
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
import scipy.signal as sps
#import pysam


##read datasets
fft = pd.read_csv(snakemake.input['FFT'], header = 0, sep = "\t", index_col = 0)
#print(fft.columns)
coverage = pd.read_csv(snakemake.input['COVERAGE'], sep = '\t', header = None, index_col = 0)
wps = pd.read_csv(snakemake.input['WPS'], sep = '\t', header = None, index_col = 0)
genes = pd.read_csv(snakemake.input['GENES'], sep = '\t', header = None, index_col = 0)
mean_expression = pd.read_csv(snakemake.input['MEAN_EXP'], header = 0, index_col = 0)
gc_promotor = pd.read_csv(snakemake.input['GC'], header = 0, sep = "\t", index_col = 0)
monocytes = pd.read_csv(snakemake.input['MONOCYTES'], header = 0, sep = "\t", index_col = 0, compression = "gzip")

##intersect gene id's
union = genes.index.intersection(fft.index)

monocytes = monocytes["monocytes"]
monocytes = monocytes.loc[union]
monocytes = np.log2(monocytes + 1)

fft = fft.loc[union]
coverage = coverage.loc[union]
wps = wps.loc[union]

##calculate most important features
features = pd.DataFrame(index = fft.index)

features["mean_expression"] = mean_expression
features['mean_cov_5k'] = coverage.loc[:,0:5000].mean(axis = 1)
features['mean_cov_2k'] = coverage.loc[:,0:2000].mean(axis = 1)
features = features.merge(gc_promotor, left_index = True, right_index = True)
features['mean_wps_5k'] = wps.loc[:,0:5000].mean(axis = 1)
features['mean_cov_body1kb'] = coverage.loc[:,1000:2000].mean(axis = 1)
features['mean_wps_2k'] = wps.loc[:,0:2000].mean(axis = 1)
features["177"] = fft["177"]
features["198"] = fft["198"]
features["195"] = fft["195"]
features['ndr_mean_cov'] = 0
features["208"] = fft["208"]
features["192"] = fft["192"]
features["201"] = fft["201"]
features["205"] = fft["205"]


##smooth wps
wps_2k = wps.loc[:,0:2000]
wps_2k_smoothed = wps_2k.copy()

for gene in wps_2k_smoothed.index:
    wps_2k_smoothed.loc[gene] = sps.savgol_filter(wps_2k.loc[gene], window_length = 65, polyorder = 3)

import numpy as np
import math

for gene in wps_2k_smoothed.index:

    ##calculate peak distances
    row = wps_2k_smoothed.loc[gene,0:5000]
    max_peaks = np.array(sps.argrelextrema(np.array(row), comparator = np.greater, order = 50)).flatten()
    
    try:
        peak1_upstream = max(list(filter(lambda x: x < 1000, max_peaks)))
        peak1_body = min(list(filter(lambda x: x > 1000, max_peaks)))

        features.loc[gene, 'ndr_mean_cov'] = coverage.loc[gene,int(peak1_upstream):int(peak1_body)].mean()

    except:
        continue

##add monocyte expression
features['monocytes'] = monocytes

features.to_csv(snakemake.output['FEATURES'], sep = '\t', header = True)



##unused features

#features['median_cov_body1kb'] = coverage.loc[:,1000:2000].median(axis = 1)
#features['var_cov_body1kb'] = coverage.loc[:,1000:2000].var(axis = 1)
#features['mean_cov_upstream1kb'] = coverage.loc[:,0:1000].mean(axis = 1)
#features['median_cov_upstream1kb'] = coverage.loc[:,0:1000].median(axis = 1)
#features['var_cov_upstream1kb'] = coverage.loc[:,0:1000].var(axis = 1)
#features['variance_cov_2k'] = coverage.loc[:,0:2000].var(axis = 1)
#features['median_cov_2k'] = coverage.loc[:,0:2000].median(axis = 1)
#features['mean_wps_body1kb'] = wps.loc[:,1000:2000].mean(axis = 1)
#features['median_wps_body1kb'] = wps.loc[:,1000:2000].median(axis = 1)
#features['var_wps_body1kb'] = wps.loc[:,1000:2000].var(axis = 1)
#features['mean_wps_upstream1kb'] = wps.loc[:,0:1000].mean(axis = 1)
#features['median_wps_upstream1kb'] = wps.loc[:,0:1000].median(axis = 1)
#features['var_wps_upstream1kb'] = wps.loc[:,0:1000].var(axis = 1)
#features['variance_wps_2k'] = wps.loc[:,0:2000].var(axis = 1)
#features['median_wps_2k'] = wps.loc[:,0:2000].median(axis = 1)
#features['mean_cov'] = coverage.loc[:,1000:12700].mean(axis = 1)
#features['mean_wps'] = wps.loc[:,1000:12700].mean(axis = 1)
#features['median_cov'] = coverage.loc[:,1000:12700].median(axis = 1)
#features['median_wps'] = wps.loc[:,1000:12700].median(axis = 1)
#features['amp1'] = 0
#features['amp2'] = 0
#features['ndr_width1'] = 0
#features['ndr_width2'] = 0