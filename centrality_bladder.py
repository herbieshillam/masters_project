# centrality_bladder.py

import numpy as np

# Note that genes 22 and 24 have been removed from consideration due to the 
# corresponding genes having zero affect on the main component


# Closeness Centrality:

# Distance matrix - calculated using Dijkstra's algorithm
B = np.array([[0,2,4,2,5,5,2,4,2,3,3,2,1,3,2,3,2,3,2,4,3,1,3],
              [2,0,6,4,7,7,4,4,4,5,3,2,1,5,4,3,2,1,4,6,5,3,5],
              [4,6,0,2,7,7,4,8,4,5,7,6,5,5,2,7,6,7,4,6,1,3,5],
              [2,4,2,0,5,5,2,6,2,3,5,4,3,3,2,5,4,5,2,4,1,1,3],
              [5,7,7,5,0,2,5,9,5,2,8,7,6,5,5,8,7,8,3,1,6,4,4],
              [5,7,7,5,2,0,5,9,5,2,8,7,6,5,5,8,7,8,3,1,6,4,4],
              [2,4,4,2,5,5,0,6,2,3,5,4,3,1,2,5,4,5,2,4,3,1,2],
              [4,4,8,6,9,9,6,0,6,7,1,2,3,7,6,3,4,5,6,8,7,5,7],
              [2,4,4,2,5,5,2,6,0,3,5,4,3,3,2,5,4,5,2,4,3,1,2],
              [3,5,5,3,2,2,3,7,3,0,6,5,4,3,3,6,5,6,1,1,4,2,2],
              [3,3,7,5,8,8,5,1,5,6,0,1,2,6,5,2,3,4,5,7,6,4,6],
              [2,2,6,4,7,7,4,2,4,5,1,0,1,5,4,1,2,3,4,6,5,3,5],
              [1,1,5,3,6,6,3,3,3,4,2,1,0,4,3,2,1,2,3,5,4,2,4],
              [3,5,5,3,5,5,1,7,3,3,6,5,4,0,3,6,5,6,2,4,4,2,1],
              [2,4,2,2,5,5,2,6,2,3,5,4,3,3,0,5,4,5,2,4,1,1,3],
              [3,3,7,5,8,8,5,3,5,6,2,1,2,6,5,0,3,4,5,7,6,4,6],
              [2,2,6,4,7,7,4,4,4,5,3,2,1,5,4,3,0,3,4,6,5,3,5],
              [3,1,7,5,8,8,5,5,5,6,4,3,2,6,5,4,3,0,5,7,6,4,6],
              [2,4,4,2,3,3,2,6,2,1,5,4,3,2,2,5,4,5,0,2,3,1,1],
              [4,6,6,4,1,1,4,8,4,1,7,6,5,4,4,7,6,7,2,0,5,3,3],
              [3,5,1,1,6,6,3,7,3,4,6,5,4,4,1,6,5,6,3,5,0,2,4],
              [1,3,3,1,4,4,1,5,1,2,4,3,2,2,1,4,3,4,1,3,2,0,2],
              [3,5,5,3,4,4,2,7,2,2,6,5,4,1,3,6,5,6,1,3,4,2,0]])


N = len(B) # number of genes
cen = np.zeros([N,1]) # initialise centrality vector

for i in range(len(B)):
    cen[i] = ((1/(N-1))*sum(B[i]))**(-1) # gives closeness centrality measure


# Spectral Centrality:

# Adjacency matrix
A = np.array([[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0],
              [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0],
              [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
              [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0],
              [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
              [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
              [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0],
              [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
              [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
              [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0],
              [0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
              [0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,0,0],
              [1,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0],
              [0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
              [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0],
              [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
              [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
              [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
              [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1],
              [0,0,0,0,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
              [0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
              [1,0,0,1,0,0,1,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0],
              [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0]])

v = np.ones([N,1])  # initial vector - mix of all eigenvectors
for i in range(100):
    v = A.dot(v)
    
v = v/max(v) # proportional to eigenvector that corresponds to leading eigenvalue,
             # gives weights of importance of all genes in network