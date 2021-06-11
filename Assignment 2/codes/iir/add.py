import numpy as np

# This function adds two polynomials defined by vectors x and y
def add(x,y):
    m = len(x)
    n = len(y)
    if m == n:
        z = np.add(x,y)
    elif m > n:
        z = np.add(x,np.concatenate([np.zeros(m-n),y]))
    else:
        z = np.add(np.concatenate([np.zeros(n-m),x]),y)
    return z

