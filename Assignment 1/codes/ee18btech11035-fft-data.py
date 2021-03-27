import soundfile as sf
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

def H(z,num,den):
	Num = np.polyval(num,z**(-1))
	Den = np.polyval(den,z**(-1))
	return Num/Den
	
x,fs = sf.read('Sound_Noise.wav')
order = 4
fc = 4000.0
Wn = 2*fc/fs

N_original = len(x)
temp = int(np.log2(N_original)) +1
N_fft = 2**temp
x_fft = np.pad(x, (0,N_fft-N_original), 'constant', constant_values=(0))
x_fft = np.array(x_fft).reshape(len(x_fft),1)

num,den = signal.butter(order,Wn,'low')

k = np.arange(len(x_fft))
w = 2*np.pi*k/len(x_fft)
z = np.exp(1j * w)
H_z = H(z,num,den)

np.savetxt("x_fft.dat",x_fft)

f = open("H_z.dat", "w")
for i in range(len(H_z)):
    f.write(str(H_z[i].real) + " " + str(H_z[i].imag) + "\n")
f.close()