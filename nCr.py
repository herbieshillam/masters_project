# nCr.py

import math

def nCr(n,k): 
    # Function for combinatorics 'n choose k'
    
    # n = value to choose from
    # k = value chosen
    
    choose = math.factorial(n) / (math.factorial(k)*math.factorial(n-k))   
    return choose