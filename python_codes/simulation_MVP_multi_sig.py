#!/usr/bin/env python3

import numpy as np
import import_ipynb
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


def get_standard(Z):
    return (Z - Z.mean(axis=0)) / Z.std(axis=0, ddof=1)


def sum_sigma_K(sigma, Kernel):
    sum_val = 0
    for i in range(len(Kernel)):
        sum_val += sigma[i] * Kernel[i]
    return sum_val


##### Simulation of data

num_simulation = 50
num_random_effects = 2
num_samples = 30000
num_fixed_effects = 0
num_features = 10

sig_1 = 0.2
sig_2 = 0.5
sig_res = 0.3

simulation_results = []
simulation_runtime = []
total_start = unix()

for sim in range(num_simulation):
    print('------------------------------------------------------------')
    print('simulation #', str(sim), ' started...')

    # Generate data
    X = np.ones((num_samples, 1))  ### only intercepts

    # Generate kernel
    Z_1 = get_standard(np.random.normal(size=(num_samples, num_random_effects)))
    Z_2 = get_standard(np.random.normal(size=(num_samples, num_random_effects)))

    K_1 = Z_1 @ Z_1.T / Z_1.shape[1]
    K_2 = Z_2 @ Z_2.T / Z_2.shape[1]

    I = np.eye(K_1.shape[0])

    Kernel = [K_1, K_2]
    sigma = [sig_1, sig_2, sig_res]

    # Generate phenotype assuming the generative process being y ~ MVN(Xb, V)
    ## b = np.random.normal(size=(X.shape[1],1)) ## equal b for all genes
    b = np.random.normal(size=(X.shape[1], num_features))  ### different b for each gene

    V = sum_sigma_K(sigma, Kernel) + sigma[len(sigma) - 1] * I
    L = LA.cholesky(V)

    # Sample from spherical normal dist
    R = np.random.normal(size=(L.shape[1], num_features))
    # r = MVN(0, I)

    # reshape it to have the V covariance structure
    Y = X @ b + L @ R
    # L @ r ~ MVN(0, V)

    a_model_result = []
    results = []
    time_log = []

    start_time = unix()
    for i in range(Y.shape[1]):  #
        a_model_result = ([-1], [-1], [-1], [-1])
        if Y.sum(axis=0)[i] != 0:
            a_model_result = lmm.reml_ei(Kernel, X, Y[:, i][:, np.newaxis],verbose=False)

        results.append(a_model_result)

    execution_time = unix() - start_time
    print('------------------------------------------------------------')
    print('simulation #', str(sim), ' execution time is ', str(execution_time), ' seconds')
    print('------------------------------------------------------------')
    simulation_results.append(results)
    simulation_runtime.append(execution_time)

total_runtime = unix() - total_start


gene_sim_list = [np.empty((num_simulation, num_random_effects + 2), dtype=object) for gene in range(num_features)]
print(gene_sim_list[0].shape)

sim_count = 0
for sim in simulation_results:
    for a_gene_index in range(len(sim)):
        gene_sim_list[a_gene_index][sim_count, 0] = sim[a_gene_index][0][0]  ## sig 1
        gene_sim_list[a_gene_index][sim_count, 1] = sim[a_gene_index][0][1]  ## sig 2
        gene_sim_list[a_gene_index][sim_count, 2] = sim[a_gene_index][0][2]  ## sig res
        gene_sim_list[a_gene_index][sim_count, 3] = sim[a_gene_index][len(sim[0]) - 1]  # convergence

    sim_count += 1


print('total_runtime is: ', str(total_runtime), ' seconds')
print('total_runtime is: ', str(total_runtime/60), ' mins')

#### saving the results
folder_name = 'sim_' + str(num_simulation)+ '_numFeat_'+ str(num_features) + '_numSamp_' + str(num_samples)
os.mkdir(folder_name)

# print(gene_sim_list[0]) ### View the results for the first feature - list contains one matrix per gene
runtime_df = pd.DataFrame(simulation_runtime, columns=['run_time'])
runtime_df.to_csv(folder_name + '/' +'runtime_info_' + str(num_simulation) +'_simulation_' + str(num_samples) + '_numFeat' + str(num_features) + '.csv',encoding='utf-8', index=False)

for i in range(num_features):
    a_gene_df = pd.DataFrame(gene_sim_list[i], columns=['sig1', 'sig2', 'sig_res', 'convergence'])
    a_gene_df.to_csv(folder_name + '/' + 'gene_' + str(i) + '_sim_' + str(num_simulation) + '_numFeat_'+ str(num_features) + '_numSamp_' + str(num_samples) + '.csv',
                     encoding = 'utf-8', index = False)



