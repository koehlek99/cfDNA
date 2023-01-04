import gzip
import pandas as pd
import numpy as np
import scipy.signal as sps


##read datasets
fft = pd.read_csv(snakemake.input['FFT'], header = 0, sep = "\t", index_col = 0)
#print(fft.columns)
coverage = pd.read_csv(snakemake.input['COVERAGE'], sep = '\t', header = None, index_col = 0)
wps = pd.read_csv(snakemake.input['WPS'], sep = '\t', header = None, index_col = 0)
gc_promotor = pd.read_csv(snakemake.input['GC'], header = 0, sep = "\t", index_col = 0)
expression = pd.read_csv(snakemake.input['EXPRESSIONATLAS'], header = 0, sep = "\t", index_col = 0, compression = "gzip")

##intersect gene id's
union = expression.index.intersection(fft.index)

##filter and log-transform RNAtable 
expression = expression.loc[union]
expression = np.log2(expression.iloc[:,69:]+1)

fft = fft.loc[union]
coverage = coverage.loc[union]
wps = wps.loc[union]


##calculate all features
features = fft.copy()

features['mean_cov_body1kb'] = coverage.loc[:,1000:2000].mean(axis = 1)
features['median_cov_body1kb'] = coverage.loc[:,1000:2000].median(axis = 1)
features['var_cov_body1kb'] = coverage.loc[:,1000:2000].var(axis = 1)
features['mean_cov_upstream1kb'] = coverage.loc[:,0:1000].mean(axis = 1)
features['median_cov_upstream1kb'] = coverage.loc[:,0:1000].median(axis = 1)
features['var_cov_upstream1kb'] = coverage.loc[:,0:1000].var(axis = 1)

features['variance_cov_2k'] = coverage.loc[:,0:2000].var(axis = 1)
features['mean_cov_2k'] = coverage.loc[:,0:2000].mean(axis = 1)
features['median_cov_2k'] = coverage.loc[:,0:2000].median(axis = 1)

features['mean_wps_body1kb'] = wps.loc[:,1000:2000].mean(axis = 1)
features['median_wps_body1kb'] = wps.loc[:,1000:2000].median(axis = 1)
features['var_wps_body1kb'] = wps.loc[:,1000:2000].var(axis = 1)

features['mean_wps_upstream1kb'] = wps.loc[:,0:1000].mean(axis = 1)
features['median_wps_upstream1kb'] = wps.loc[:,0:1000].median(axis = 1)
features['var_wps_upstream1kb'] = wps.loc[:,0:1000].var(axis = 1)

features['variance_wps_2k'] = wps.loc[:,0:2000].var(axis = 1)
features['mean_wps_2k'] = wps.loc[:,0:2000].mean(axis = 1)
features['median_wps_2k'] = wps.loc[:,0:2000].median(axis = 1)

features['mean_cov_5k'] = coverage.loc[:,0:5000].mean(axis = 1)
features['mean_wps_5k'] = wps.loc[:,0:5000].mean(axis = 1)

features['mean_cov'] = coverage.loc[:,0:11000].mean(axis = 1)
features['mean_wps'] = wps.loc[:,0:11000].mean(axis = 1)

features['median_cov'] = coverage.loc[:,0:11000].median(axis = 1)
features['median_wps'] = wps.loc[:,0:11000].median(axis = 1)

features['amp1'] = 0
features['amp2'] = 0 
features['ndr_width1'] = 0
features['ndr_width2'] = 0 
features['ndr_mean_cov'] = 0 


##smooth wps
wps_2k = wps.loc[:,0:2000]
wps_2k_smoothed = wps_2k.copy()

for gene in wps_2k_smoothed.index:
    wps_2k_smoothed.loc[gene] = sps.savgol_filter(wps_2k.loc[gene], window_length = 65, polyorder = 3)

import numpy as np
import math

for gene in wps_2k_smoothed.index:
    
    ##calculate amplitudes
    row = wps_2k_smoothed.loc[gene,1000:2000]
    maxima = np.array(sps.argrelextrema(np.array(row), comparator = np.greater, order = 50)).flatten()
    minima = np.array(sps.argrelextrema(np.array(row), comparator = np.less, order = 50)).flatten()
    
    if maxima.size > 1 and minima.size > 1: 
        amp1 = max(row.iloc[maxima[0]] - row.iloc[minima[0]], row.iloc[minima[0]] - row.iloc[maxima[0]])
        amp2 = max(row.iloc[maxima[1]] - row.iloc[minima[1]], row.iloc[minima[1]] - row.iloc[maxima[1]])
    elif maxima.size == 1 and minima.size == 1:
        amp1 = max(row.iloc[maxima[0]] - row.iloc[minima[0]], row.iloc[minima[0]] - row.iloc[maxima[0]])
        amp2 = 0
    else: 
        amp1 = 0
        amp2 = 0
    
    features.loc[gene,'amp1'] = amp1
    features.loc[gene,'amp2'] = amp2
    
    ##calculate peak distances 
    row = wps_2k_smoothed.loc[gene,0:5000]
    
    max_peaks = np.array(sps.argrelextrema(np.array(row), comparator = np.greater, order = 50)).flatten()
    try: 
        peak1_upstream = max(list(filter(lambda x: x < 1000, max_peaks)))
        peak1_body = min(list(filter(lambda x: x > 1000, max_peaks))) 
        peak2_body = list(filter(lambda x: x > 1000, max_peaks))[1]
        features.loc[gene, 'ndr_width1'] = int(peak1_body) - int(peak1_upstream)
        features.loc[gene, 'ndr_width2'] = int(peak2_body) - int(peak1_upstream)
        features.loc[gene, 'ndr_mean_cov'] = coverage.loc[gene,int(peak1_upstream):int(peak1_body)].mean()
        
    except: 
        continue


features["mean_expression"] = expression.mean(axis=1)

features = features.merge(gc_promotor, left_index = True, right_index = True)

features.to_csv(snakemake.output['FEATURES'], sep = '\t', header = True)
