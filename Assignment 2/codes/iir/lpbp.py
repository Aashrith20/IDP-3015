import numpy as np

def lpbp(p,Omega0,B,Omega_p2):
    N = len(p)
    const = np.array([1,0,Omega0**2])
    v = const.copy()
    if N > 2:
        for i in range(1,N):
            M = len(v)
            v[M-i-1] = v[M-i-1]+(p[i]*B**i)
            if i < N-1:
                v = np.convolve(const,v)
        den = v
    elif N == 2:
        M = len(v)
        v[M-i-1] = v[M-i-1]+(p[N-1]*B)
        den = v
    else:
        den = p
    num = np.concatenate([np.array([1]),np.zeros(N-1)])
    G_bp = abs(np.polyval(den,1j*Omega_p2)/np.polyval(num,1j*Omega_p2))
    return num,den,G_bp

