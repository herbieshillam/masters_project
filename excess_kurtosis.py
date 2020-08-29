# excess_kurtosis.py

import numpy as np
import pandas as pd

def excess_kurtosis(x):
    # Calculates excess kurtosis for each gene
    
    # x = all observations for single gene
    
    if isinstance(x,pd.DataFrame):
        x = x.to_numpy()
        
    mu = np.mean(x) # mean
    std = np.std(x) # overall standard deviation
    mu4 = np.mean((x - mu)**4) # 4th standardised moment
    KURT = mu4 / (std**4) # kurtosis value
    return KURT - 3