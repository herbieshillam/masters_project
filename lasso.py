# lasso.py

import numpy as np
from numpy.linalg import norm
from scipy.linalg import sqrtm

def lasso(W_11,S,b,s_12,rho,tol,max_it):
    # Lasso regression using coordinate descent and soft threshold operator
    
    # W_11 = upper block of W
    # S = sample covariance matrix
    # b = column of B, updates with betahat, used in calculation of Theta
    # s_12 = target column of S
    # rho = tuning parameter
    # tol = tolerance
    # max_it = maximum number of iterations
    
    it = 0 
    W_11_half = sqrtm(W_11) # square root of W_11
    b_old = np.zeros([len(S)-1,1]).T
    
    while np.sum(abs(b - b_old)) > tol:  # loops until b converges
        b_old = b.copy()
        if it > max_it:
            print('Betahat did not converge within maximum number of iterations')
            betahat = b
            return betahat 
        
        for j in range(len(S)-1): # iterates over each entry in b vector
            sum_ = 0 # initialise counter used in lasso regression
            for k in range(len(S)-1):
                if k != j:
                    sum_ += W_11[k][j]*b[k]         
            b[j] = STO(s_12[j] - sum_,rho)/(norm(W_11_half[:,j])**2) # L2 norm
            # STO = soft threshold operator 
        it += 1   
    betahat = b.copy()
    return betahat