# glasso.py

import numpy as np
import pandas as pd

def glasso(S,rho,tol,max_it): 
    # Min_{Theta > 0}: tr(S * Theta) - log(det(Theta)) + rho * ||Theta||_{1}
    
    # S = sample covariance matrix
    # rho = tuning parameter
    # tol = tolerance: to converge, requires the difference between consecutive iterations < tol
    # max_it = maximum number of iterations 
    
    if isinstance(S,pd.DataFrame):
        S = S.to_numpy()

    W = S + rho*np.eye(len(S)) # initialise estimated covariance matrix W
    W_old = np.zeros([len(S),len(S)])
    it = 0 
    B = np.ones([len(S)-1,len(S)]) # each col represents betahat for respective target col
    Theta = np.zeros([len(S),len(S)]) # initialise estimated precision matrix
    
    while np.sum(abs(W - W_old)) > tol: # loops until W converges
        if it > max_it:
            print('W did not converge within maximum number of iterations')
            Theta = 'error - glasso failed to converge'
            return W,Theta,it
        
        for k in range(len(S)): # iterates so each column will be target once
            W_old = W.copy()
            start = list(range(len(W))) # list of each column number
            other = start.copy()
            target = other[-1-k] # gives target column
            del other[-1-k] # gives remaining columns, make up W_11 block matrix
            
            W_11 = W[other][:,other] # (p-1) x (p-1) matrix, p = no. features
            s_12 = S[other][:,target] # (p-1) x 1 vector
            
            betahat = lasso(W_11,S,B[:,len(S)-1-k],s_12,rho,tol,max_it) # lasso regression
            
            # Updating w_12 and w_12_T:
            W[other,target] = W_11.dot(betahat) # w_12
            W[target,other] = W[other,target].T # w_12_T
            
            # Updating theta_22 and theta_12:
            Theta[target,target] = 1 / ( W[target,target] - W[target,other].dot(betahat) ) # theta_22
            Theta[other,target] = - ( betahat*Theta[target,target] ) # theta_12
        
        it += 1

    return Theta, W ,it