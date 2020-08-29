# sparse_bayes.py

import statistics
import numpy as np
from numpy.linalg import norm
import pandas as pd
import random
from random import randint

def bayes_regression(data,rho):    
    # Applies Bayesian model to gene expression data to find maximum
    # likelihood estimates for weights and hence parent set for each gene
    
    # data = standardised gene expression data (p x n)
    # rho = value that multiplies var to get sigma: larger rho increases sparsity
    # of resulting graph
    
    if isinstance(data,pd.DataFrame):
        data = data.to_numpy()
        
    p = len(data) # number of genes
    n = len(data[0]) # number of observations 
    W = np.zeros([p,p]) # empty weight matrix  

    for i in range(p): # finds parents for 1 gene at a time
        
        e_i = data[i,:] # target gene
        
        # Remove gene i from data:
        start = list(range(p)) 
        other = start.copy()
        del other[i]                     
        data1 = data[other,:]          
        
        # Initialise inverse variance:
        alpha = np.ones(p-1)*np.inf # initially all set to +infinity
        alpha_diff = alpha.copy()
        
        # Initialise variance, sparsity factors, quality factors and theta:
        var = statistics.variance(e_i)
        sigma = var * rho
        s = np.zeros(p-1) # sparsity factor, note p-1 used since removed 1 gene / row
        S = np.zeros(p-1)
        q = np.zeros(p-1) # quality factor
        Q = np.zeros(p-1)
        theta = np.zeros(p-1)
             
        it = 0 # initialise iteration counter
        A = np.zeros(p-1) # binary parent vector: 1 = in basis, 0 = not in basis
        
        check = np.zeros(p-1)
        left = [] # needed for generating a gene that is different from previous few
        
        # Choose first potential parent j:
        j = randint(0,p-2)
        check[j] = 1
        e_j = data1[j,:] # parent j
        e_j = np.reshape(e_j,[1,len(e_j)]) # so dotting below works - row vector 1xn      
        
        alpha[j] = norm(e_j)**2/ ( ((norm(e_j.dot(e_i))**2)/norm(e_j)**2) - sigma) 
        A[j] = 1 # j in basis
        basis = e_j # gene expression for genes currently parents
        list_ = [j] # list of genes in basis, in order of addition to basis
              
        while max(alpha_diff) > 0.5: # until alpha converges         
                                      
            alpha_old = alpha.copy()
            
            # Initiate / update Cminusi_inv:  
            mid = np.diag(np.ones(len(list_))/alpha[list_])            
            C = sigma*np.eye(n) + (basis.T.dot(mid)).dot(basis)  
            Cminusi_inv = np.linalg.inv(C)  
                  
            # Check a new parent: 
            for k in range(len(check)):
                if check[k] != 1:
                    left.append(k)
            if left == []: # if all have been checked it resets checklist 
                check = np.zeros(p-1)
                j = randint(0,p-2)
                check[j] = 1
            else:
                j = random.choice(left)  
                check[j] = 1
            left = []         
        
            e_j = data1[j,:] 
            e_j = np.reshape(e_j,[1,len(e_j)]) # so dotting works
            
            S[j] = (e_j.dot(Cminusi_inv)).dot(e_j.T)
            Q[j] = (e_j.dot(Cminusi_inv)).dot(e_i)   
            
            s[j] = S[j].copy()
            q[j] = Q[j].copy()
                
            if A[j] == 1: # if j is currently in basis set s and j differently
                s[j] = (alpha[j]*S[j])/(alpha[j]-S[j])
                q[j] = (alpha[j]*Q[j])/(alpha[j]-S[j])

            theta[j] = q[j]**2 - s[j] 
        
            if theta[j] > 0 and alpha[j] < np.inf: # j is already a parent
                alpha[j] = (s[j]**2)/((q[j]**2)-s[j]) 
                alpha_diff[j] = abs(alpha[j]-alpha_old[j]) 
                    
            elif theta[j] > 0 and alpha[j] == np.inf: # j is not a parent but should be
                alpha[j] = (s[j]**2)/((q[j]**2)-s[j])
                alpha_diff[j] = abs(alpha[j]-alpha_old[j]) 
                A[j] = 1
                list_.append(j)
                basis = np.append(basis,e_j,axis=0)
                
            elif theta[j] <= 0 and alpha[j] < np.inf and sum(A) > 1: # j already a parent
                # cant have empty parent set
                A[j] = 0
                alpha[j] = np.inf
                alpha_diff[j] = np.inf 
                for k in range(len(basis)-1):
                    if basis[k,:].all() == e_j.all():
                        basis = np.delete(basis, (k), axis=0) # deletes e_j from basis
                        break
                      
                for l in range(len(list_)):
                    if list_[l] == j:
                        del list_[l] # deletes j from list of genes in basis
                        break
            else:
                alpha_diff[j] = 0 # j not a parent and shouldn't be
   
            it += 1
            
            if it == 2000:
                print('Gene',i,'did not converge in 2000 iterations') 
                A = np.ones(p-1)*np.inf # infinity values indicate no convergence
                break
            
        W[other,i] = A # updates weight matrix for gene i's parents
           
    return W