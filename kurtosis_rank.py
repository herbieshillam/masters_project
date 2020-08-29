# kurtosis_rank.py

import numpy as np

def kurtosis_rank(gene_data_ttest):
    # Ranks genes by the absolute value of their excess kurtosis, largest first
    
    # gene_data_ttest = data with significantly small p-values, smallest first
    
    dataset = gene_data_ttest.iloc[:,1:-1].to_numpy() # removes non-numeric columns
    xs_kurt = np.zeros(len(dataset)) # initialise excess kurtosis
    
    for i in range(len(dataset)):
        xs_kurt[i] = abs(excess_kurtosis(dataset[i]))
    
    gene_data_ttest['Abs of Excess Kurtosis'] = xs_kurt # make abs of excess kurtosis a DataFrame column
        
    gene_data_ttest = gene_data_ttest.sort_values(by=['Abs of Excess Kurtosis'],ascending=False) # ranking
    
    glasso_label = gene_data_ttest.iloc[0:100] # selects 100 genes with top excess kurtosis values, contains gene ID's
    data_for_glasso = glasso_label.drop(glasso_label.columns[[0,-1]],axis=1) # drops non-numeric columns
    
    return data_for_glasso, glasso_label, gene_data_ttest

