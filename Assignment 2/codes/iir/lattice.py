import numpy as np


def lattice(c, v):
    # c = np.array([0,0.44,0.36,0.02])
    # v = np.array([1,0.4,0.18,-0.2])
    u = np.fliplr(v)
    m = len(v)
    K[m - 2] = v[m - 1]
    # C = np.zeros(m)
    C[m] = c[m]

    while m > 1 and K[m - 1] == 1:
        c = c - C[m] * u
        v = (v - K[m] * u) / (1 - K[m] ** 2)
        m = m - 1
        v = v[:m + 1]
        c = c[:m + 1]
        u = np.fliplr(v)
        if (m > 1):
            K[m - 1] = v(m)
        C[m] = c[m]
    return C, K