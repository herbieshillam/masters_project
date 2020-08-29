# StARS.py

import random
import numpy as np
import pandas as pd

def StARS(glasso_data,rho_values,b,num):
    # Identifies optimal rho to use in glasso
    
    # glasso_data = table of standardised gene expression data for p interesting genes
    # rho_values = values of tuning parameter that are to be considered e.g. [0.4,0.6,0.8]
    # NOTE: rho_values needs to be an increasing list
    # b = number between 1 and n, size of each subsample (in terms of number of
    #                                                         samples selected)
    # num = number of subsamples 
    
    if isinstance(glasso_data,pd.DataFrame):
        glasso_data = glasso_data.to_numpy()
    n = len(glasso_data[0]) # number of samples
    p = len(glasso_data) # number of genes
    
    # Generate the samples used for each subsample:
    N = []
    for i in range(num): 
        N.append(random.sample(range(0,n),b))

    # Initiate empty lists:
    glasso_data_sample = [ [[] for i in range(len(N))] for i in range(len(rho_values))]
    glasso_data_std_sample = [ [[] for i in range(len(N))] for i in range(len(rho_values))]
    S_sample = [ [[] for i in range(len(N))] for i in range(len(rho_values))]
    Theta_sample = [ [[] for i in range(len(N))] for i in range(len(rho_values))]
    theta_hat = [ np.zeros([p,p]) for i in range(len(rho_values))]
    xi_hat = [ np.zeros([p,p]) for i in range(len(rho_values))]
    D_hat = [[] for i in range(len(rho_values))]
    
    for l in range(len(rho_values)): # for each value of rho
        for i in range(len(N)):  # for each subsample
            print('rho number:',l,', subsample number:',i) # to identify codes progress since can take time
            glasso_data_sample[l][i] = glasso_data[:,N[i]]
            glasso_data_std_sample[l][i] = standardise(glasso_data_sample[l][i])
            S_sample[l][i] = samplecov(glasso_data_std_sample[l][i])
            Theta_sample[l][i], W, it = glasso(S_sample[l][i],rho_values[l],0.1,1000)
            
            # Turning Theta into adjacency matrix and turning into fraction:
            for j in range(len(Theta_sample[l][i])):
                for k in range(len(Theta_sample[l][i][0])):
                    if Theta_sample[l][i][j][k] != 0:
                        Theta_sample[l][i][j][k] = 1
                        theta_hat[l][j][k] +=1
        theta_hat[l] = theta_hat[l] / len(N)
        
        # Measure of instability across subsamples:
        xi_hat[l] = 2*theta_hat[l]*(np.ones([p,p])-theta_hat[l])

        # Averaging over all edges:
        counter = 0 
        for s in range(len(xi_hat[l])):
            for t in range(len(xi_hat[l])):
                if s < t:
                    counter += xi_hat[l][s][t]
        D_hat[l] = counter / (nCr(p,2))
    
    print('D_hat:',D_hat)
    # Monotonise D_hat:
    D_bar = D_hat.copy()
    for i in range(len(rho_values)-1):
        D_bar[-1-i] = max(D_hat[-1-i:])
    print('D_bar:',D_bar)
    
    # Select optimum rho:
    valid_rho = []
    for i in range(len(D_bar)):
        if D_bar[i] <= 0.05:
            valid_rho.append(D_bar[i])
    D_bar_star = valid_rho[0]
    for i in range(len(D_bar)):
        if D_bar[i] == D_bar_star:
            print('Optimal value of rho:',rho_values[i])
            return rho_values[i]