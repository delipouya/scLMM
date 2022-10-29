#!/usr/bin/env python3
import numpy as np
#import import_ipynb
import numpy.linalg as LA
import LMM as lmm
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
from time import time as unix

np.random.seed(10)

##### Simulation of data
num_random_effects = 20
num_samples = 3000
num_fixed_effects = 0
num_features = 1000
sig_1 = 0.2
sig_2 = 0.8

# Generate data
X = np.ones((num_samples, 1))  ### only intercepts
Z = np.random.normal(size=(num_samples, num_random_effects))

# Generate kernel
Z = (Z - Z.mean(axis=0)) / Z.std(axis=0, ddof=1)
print(Z.shape)

K = Z @ Z.T / Z.shape[1]
I = np.eye(K.shape[0])
Kernel = [K]

# Generate phenotype assuming the generative process being y ~ MVN(Xb, V)
b = np.random.normal(size=(X.shape[1], 1))
sigma = [sig_1, sig_2]
V = sigma[0] * K + sigma[1] * I
L = LA.cholesky(V)

# Sample from spherical normal dist
R = np.random.normal(size=(L.shape[1], num_features))
# r = MVN(0, I)

# reshape it to have the V covariance structure
Y = X @ b + L @ R
# L @ r ~ MVN(0, V)

# Y shape now is 1000 by 1000 meaning that we have 10000 simulated reponse with covariance of V

print('X shape is: ', X.shape)
print('Z shape is: ', Z.shape)
print('K shape is: ', K.shape)
print('V shape is: ', V.shape)
print('L shape is: ', L.shape)
print('R shape is: ', R.shape)
print('Y shape is: ', Y.shape)

a_model_result = []
results = []

start_time = unix()
for i in range(Y.shape[1]):
    # print('------------------', i, '------------------')
    a_model_result = ([0], [0], [0], [0])
    if Y.sum(axis=0)[i] != 0:
        a_model_result = lmm.minqe_unbiased(Kernel, X, Y[:, i][:, np.newaxis])

    results.append(a_model_result)

execution_time = unix() - start_time

#### making a dict with gene names as keys and sigmas as the values
sig_1_list = [a_model[0][0] for a_model in results]
sig_2_list = [a_model[0][1] for a_model in results]
convergence = [a_model[3] for a_model in results]
features_name = ['gene_' + str(i + 1) for i in range(len(results))]

sig_df = pd.DataFrame(list(zip(features_name, sig_1_list, sig_2_list, convergence)),
                      columns=['genes', 'sig_1', 'sig_2', 'convergence'])
sig_df = sig_df.sort_values(by='sig_1', axis=0, ascending=False)
print('execution time: ', execution_time)
sig_df.to_csv('sharedsig_sumSamp' + str(num_samples) + '_numFeat' + str(num_features) + '.csv',
              encoding='utf-8', index=False)
