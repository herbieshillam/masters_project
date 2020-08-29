# samplecov.py

import pandas as pd
import numpy as np
    
def samplecov(sample_std):    
    # Creates sample covariance matrix from standardised p x n dataset
    
    # sample_std = standardised p x n dataset
    
    # Initialise sample covariance:    
    S = pd.DataFrame(data=np.zeros([len(sample_std),len(sample_std)]))
    
    n = len(sample_std[0]) # number of observations
    
    for i in range(len(S)):
        for j in range(len(S[0])):
            counter = 0
            for k in range(n):
                counter += sample_std[i][k]*sample_std[j][k] 
                
            S.iloc[i][j] = counter / (n-1)    
            
    return S