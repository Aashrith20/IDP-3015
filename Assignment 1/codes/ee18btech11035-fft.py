#Code by P.Aditya
#21st March 2021

import soundfile as sf
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

#If using termux
import subprocess
import shlex

#Defining H(z) from numerator and denominator coefficents
def H(z,num,den):
	Num = np.polyval(num,z**(-1))
	Den = np.polyval(den,z**(-1))
	return Num/Den

#Implementing fft algorithm (Radix 2-FFT)
def fft(x):
	N = len(x)
	if(N==2):
		return np.vstack((x[0]+x[1],x[0]-x[1]))
	
	else:
		Xeven = fft(x[0::2])
		Xodd  = fft(x[1::2])
		
		dummy = np.ogrid[:int(N/2)]
		Wnu = np.exp(-1j*2*np.pi*dummy/N)
		Wnu = np.array(Wnu).reshape(len(Xodd),1)
		
		X = np.vstack((Xeven + np.multiply(Wnu,Xodd),Xeven - np.multiply(Wnu,Xodd)))
		
	return X
	
#Implementing ifft algorithm 
def ifft(Y):
	N = len(Y)
	if(N==2):
		return np.vstack(((Y[0]+Y[1]),(Y[0]-Y[1])))
	
	else:
		Yeven = ifft(Y[0::2])
		Yodd  = ifft(Y[1::2])
		
		dummy = np.ogrid[:int(N/2)]
		Wnu = np.exp(1j*2*np.pi*dummy/N)
		Wnu = np.array(Wnu).reshape(len(Yodd),1)
		Y1 = (Yeven + np.multiply(Wnu,Yodd))
		Y2 = (Yeven - np.multiply(Wnu,Yodd))
		y = np.vstack((Y1,Y2))
		
		
	return y

#Reading the soundfile 	
x,fs = sf.read('Sound_Noise.wav')
order = 4
fc = 4000.0
Wn = 2*fc/fs

#padding zeros to original signal to make is length as 2**n
N_original = len(x)
temp = int(np.log2(N_original)) +1
N_fft = 2**temp
x_fft = np.pad(x, (0,N_fft-N_original), 'constant', constant_values=(0))
x_fft = np.array(x_fft).reshape(len(x_fft),1)


#Passing butterworth filter
num,den = signal.butter(order,Wn,'low')

#Preparing own routine filter design
#Constructing H(z)
k = np.arange(len(x_fft))
w = 2*np.pi*k/len(x_fft)
z = np.exp(1j * w)
H_z = H(z,num,den)

#Computing DFT(FFT) of x(n)
X_fft = fft(x_fft)

#Multiplying X and H to give Y and computing IDFT(IFFT) to give y(n)
Y_fft = np.multiply(np.array(H_z).reshape(len(x_fft),1),X_fft)
y_fft = (1/len(Y_fft))*np.real(ifft(Y_fft))
y_fft = y_fft.flatten()

sf.write('7.1_Sound_With_ReducedNoise_using_python.wav',y_fft[0:len(x)],fs)

plt.figure(1)
plt.plot(y_fft[0:len(x)],'g')
plt.grid()
plt.title("Output from fft algorithm using python")
plt.savefig('../figs/ee18btech11034_1.eps')
#plt.show()
