#!/usr/bin/env python

import numpy as np
#import import_ipynb
#import numpy.linalg as LA
#import LMM as lmm
import random
import time
import os
import pandas as pd
from scipy.io import mmread
import scipy.sparse as sp
import matplotlib.pyplot as plt
from IPython.display import Image
import scanpy as sc
import matplotlib.pyplot as plt
import statsmodels.api as sm

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline

np.random.seed(10)

#### import the immune subpopulation of the rat samples
data = sc.read('/home/delaram/scLMM/input_data_designMat/inputdata_rat_set1_countData_2.h5ad') ## attributes removed
data.var_names_make_unique()
# a.obs['orig.ident'].head()
### renaming the meta info column names: https://github.com/theislab/scvelo/issues/255
data.__dict__['_raw'].__dict__['_var'] = data.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})

data_numpy = data.X.toarray()
cell_sums = np.sum(data_numpy,axis=1) # row sums - library size
gene_sums = np.sum(data_numpy,axis=0) # col sums - sum reads in a gene
gene_vars = np.var(data_numpy, axis=0)
data_numpy = data_numpy[:,gene_sums != 0]

data_sub = data_numpy
strain = data.obs.strain
#### sample metadata
y_cluster =data.obs.cluster 
y_strain = np.zeros((strain.shape[0]))

for i in range(len(data.obs.strain)):
    if data.obs.strain[i] != 'DA':
        y_strain[i] = 1

## working with the rat data
num_cells = data_sub.shape[0]
num_genes = data_sub.shape[1]
num_vars = 1 # number of variables in the design matrix - strain

### sample 100 genes and subset the data
# num_genes = 100
# gene_idx = random.sample(range(0, data_numpy.shape[1]), num_genes)
# data_numpy = data_numpy[:, gene_idx]


y = data_numpy
x = y_strain

### fit a poisson regression model to each gene and save the results

### make an empty array to store the p-values and coefficients
pvalue = []
coefficient = []
fittedvalues = []
deviance = []
null_deviance = []

#yhat = []
#tvalues = []
#resid_pearson = []
#resid_deviance = []
#resid_response = []
#resid_working = []
#nobs = []
#models = []
#pearson_chi2 = []

### time the fitting process
start_time = time.time()

for i in range(len(y[0])):
    y_a_gene = y[:, i]
    model = sm.GLM(y_a_gene, x, family=sm.families.Poisson())
    result = model.fit()
    #print(result.summary())
    #models.append([result])
    coefficient.append([result.params])
    pvalue.append([result.pvalues]) ## yhat == fittedvalue == mu
    fittedvalues.append([result.fittedvalues])
    deviance.append([result.deviance])
    null_deviance.append([result.null_deviance])

    '''
    yhat.append([result.predict()])
    nobs.append([result.nobs])
    tvalues.append([result.tvalues])
    resid_pearson.append([result.resid_pearson])
    resid_deviance.append([result.resid_deviance])
    resid_response.append([result.resid_response])
    resid_working.append([result.resid_working])
    pearson_chi2.append([result.pearson_chi2])
    '''

end_time = time.time()
print('time to fit the model: ', end_time - start_time)

pvalue = np.asarray(pvalue).reshape(num_genes, num_vars)
coefficient = np.asarray(coefficient).reshape(num_genes, num_vars)
fittedvalues = np.asarray(fittedvalues).reshape(num_genes, num_cells)
deviance = np.asarray(deviance).reshape(num_genes, 1)
null_deviance = np.asarray(null_deviance).reshape(num_genes, 1)

'''
yhat = np.asarray(yhat).reshape(num_genes, num_cells)
tvalues = np.asarray(tvalues).reshape(num_genes, num_vars)
resid_pearson = np.asarray(resid_pearson).reshape(num_genes, num_cells)
resid_deviance = np.asarray(resid_deviance).reshape(num_genes, num_cells)
resid_response = np.asarray(resid_response).reshape(num_genes, num_cells)
resid_working = np.asarray(resid_working).reshape(num_genes, num_cells)
#nobs = np.asarray(nobs).reshape(num_genes, 1)

pearson_chi2 = np.asarray(pearson_chi2).reshape(num_genes, 1)
'''


#### save the results to csv files
variable_lists = [coefficient, pvalue, fittedvalues, deviance, null_deviance]
variable_names = ['coefficient', 'pvalue', 'fittedvalues', 'deviance', 'null_deviance']

for i in range(len(variable_lists)):
    np.savetxt('GLM_FA_res/' + variable_names[i] + ".csv", variable_lists[i], delimiter=",")
    print(variable_names[i] + ' saved')
