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

np.random.seed(14)


##### Simulation of data
'''
num_random_effects = 2000
num_samples = 1000
num_fixed_effects = 10
num_features = 10000
sig_1 = 0.2
sig_2 = 0.8
'''
num_random_effects = 2
num_samples = 100
num_fixed_effects = 6
num_features = 25
sig_1 = 0.2
sig_2 = 0.8



# Generate data
X = np.random.normal(size=(num_samples,num_fixed_effects))
X = np.c_[np.ones((num_samples,1)), X]

#### simulating Z as a continious normal varirable
Z = np.random.normal(size=(num_samples,num_random_effects))
Z = (Z -  Z.mean(axis=0)) / Z.std(axis=0, ddof=1)
print(Z.shape)

#### simulating Z as a categorical random varirable
Z = np.zeros((num_samples,num_random_effects))
for i in range(num_samples):
    if i < round(num_samples/2):
        Z[i,0] = 1
    else:
        Z[i,1] = 1

print(Z.shape)

# Generate kernel
K = Z @ Z.T / Z.shape[1]
I = np.eye(K.shape[0])
Kernel = [K]

# Generate phenotype assuming the generative process being y ~ MVN(Xb, V)

b = np.random.normal(size=(X.shape[1],1))
sigma = [sig_1, sig_2]
V = sigma[0] * K + sigma[1] * I

L = LA.cholesky(V)

# Sample from spherical normal dist
R = np.random.normal(size=(L.shape[1],num_features))
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
print(type(Y))

num_features_list = [10, 50, 100, 200, 500, 1000, 5000]  # , 2000, 5000
execution_time = []

for num_features in num_features_list:

    selected_features = random.sample(range(immune_data.shape[1]), num_features)
    selected_features_name = immune_data.var.index[selected_features]

    Y = immune_data.X[:, selected_features].toarray()  #### subsetting the features
    print(Y.shape)  ## data matrix

    a_model_result = []
    results = []

    start_time = time.clock()
    for i in range(Y.shape[1]):
        # print('------------------', i, '------------------')
        a_model_result = ([0], [0], [0], [0])
        if Y.sum(axis=0)[i] != 0:
            a_model_result = lmm.minqe_unbiased(Kernel, X, Y[:, i][:, np.newaxis])

        results.append(a_model_result)

    execution_time.append(time.clock() - start_time)

plt.scatter(execution_time, num_features_list)
plt.xlabel("Execution time (sec)")
plt.ylabel("Number of features")
plt.show()

#### number of slipped genes
sum(Y.sum(axis=0) == 0)

print(len(results))
for a_model in results[1:10]:
    print(a_model[0])


#### making a dict with gene names as keys and sigmas as the values
sig_1_list = [a_model[0][0] for a_model in results]
{selected_features_name[i]: sig_1_list[i] for i in range(len(selected_features_name))}

#### making a dataframe with gene names and sigmas

sig_df = pd.DataFrame(list(zip(selected_features_name, sig_1_list)),
               columns =['genes', 'sig'])
sig_df.sort_values(by='sig',axis=0,ascending=False)