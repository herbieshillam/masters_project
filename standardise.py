# standardise.py 

import pandas as pd
import numpy as np
import statistics

def standardise(sample):
    # Normalises each row so that each gene has zero mean and unit variance
    
    # sample = p x n dataset
    
    if isinstance(sample,pd.DataFrame):
        sample = sample.to_numpy()
        
    mean = np.zeros(len(sample))
    std = np.zeros(len(sample))
    
    for j in range(len(sample)):
        mean[j] = np.mean(sample[j])
        std[j] = np.sqrt(statistics.variance(sample[j]))
        
    sample_std = sample.copy()
    
    for i in range(len(sample)):
        for j in range(len(sample[0])):
            sample_std[i][j] = (sample[i][j]-mean[i])/std[i]
            
    return sample_std