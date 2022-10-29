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


#### import the immune subpopulation of the rat samples
immune_data = sc.read('/home/delaram/RatLiver/ImmuneSub_files/immuneSub_raw.h5ad') ## attributes removed
immune_data.var_names_make_unique()
print(immune_data)
print(type(immune_data.obs))
print(immune_data.obs.head())
# a.obs['orig.ident'].head()
### renaming the meta info column names: https://github.com/theislab/scvelo/issues/255
immune_data.__dict__['_raw'].__dict__['_var'] = immune_data.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})


#### adding strain information to the meta data of the dataframe
immune_data.obs['strain'] = [UMI.split('_')[1] for UMI in list(immune_data.obs.index)]
immune_data.obs.head()

### Generating a Z matrix based on the strain information
num_random_effects = 2
num_samples = immune_data.obs.shape[0]
Z = np.zeros((num_samples,num_random_effects))


for i in range(num_samples):
    if immune_data.obs['strain'][i] == 'DA':
        Z[i,0] = 1
    else:
        Z[i,1] = 1

K = Z @ Z.T / Z.shape[1]
I = np.eye(K.shape[0])
Kernel = [K]

### Generating a X matrix - includes only the intercept
X = np.ones((num_samples,1))

print(Z.shape) ### random effect
print(K.shape) ### kernel
print(immune_data.shape)
print(X.shape) ## fixed effects


#### evaluating the efficiency of the model by running it with different number of features
num_features_list = [10, 50, 100, 200, 500, 1000, 5000]
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


print(len(results)) ### number of final models
print(results[0]) ### how one of the model results looks
convergence_status = [result[3] for result in results]
sum(Y.sum(axis=0) == 0) #### number of slipped genes

#### making a dataframe with gene names and sigmas
sig_1_list = [a_model[0][0] for a_model in results]
sig_df = pd.DataFrame(list(zip(selected_features_name, sig_1_list)),
               columns =['genes', 'sig'])
sig_df = sig_df.sort_values(by='sig',axis=0,ascending=False)
sig_df.to_csv('sig_df.csv', encoding='utf-8', index=False)

#### visualizing the top genes using a heatmap
num_2vis = 50
df = sig_df['sig'].iloc[1:num_2vis]
df = df.rename(index=sig_df['genes'][1:num_2vis])

df.to_frame().style.background_gradient(cmap ='viridis').set_properties(**{'font-size': '20px'})

