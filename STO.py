# STO.py

def STO(x,t):
    # Soft-threshold operator, used in lasso regression, results in convergence of betahat
    # STO(x,t) = sign(x)(|x|-t)_+
    
    if x > 0 and t < abs(x):
        return x-t
    elif x < 0 and t < abs(x):
        return x+t
    else:
        return 0