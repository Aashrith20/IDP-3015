import struct
import numpy as np


with open('dignum.dat', 'rb') as dat_file:
    dignum = struct.unpack('f'*9, dat_file.read())

with open('digden.dat', 'rb') as dat_file:
    digden = struct.unpack('f'*9, dat_file.read())

c = dignum
v = digden

u = np.flip(v)

m = len(v)
K = np.zeros(m)
K[m-2] = v[m-1]
C = np.zeros(m)
C[m-1]=c[m-1]

while m > 1 and K[m-2] == 1:
    c = c - C[m-1]*u
    v = (v - K[m-2]*u)/(1 - K[m-1]**2)
    m = m - 1
    v = v[0:m]
    c = c[0:m]
    u = np.flip(v)
    if m > 1:
        K[m-2] = v[m]

    C[m-1] = c[m-1]