import soundfile as sf
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np

# Function to compute H(e^jw) from b,a
def H(a,b,z):
    num = np.polyval(b[::-1],z**(-1))
    den = np.polyval(a[::-1],z**(-1))
    return num/den



def My_Filter(a,b,x,Wn):
    N = len(x)
    # fft of input signal
    X_k = np.fft.fft(x)
    # Calculating H(e^jw) from coefficients
    # Low pass filter can be observed from the plot
    omega = np.linspace(0,2*np.pi,N)
    H_k = H(a,b,np.exp(1j*omega))
    # Y(e^jw) = H(e^jw)*X(e^jw)
    Y_k = np.zeros(N)+1j*np.zeros(N)
    for k in range(N):
        Y_k[k] = X_k[k]*H_k[k]
    # Using ifft to compute y(n)
    y = np.fft.ifft(Y_k)
    # subplots
    plt.figure(3)
    plt.plot(np.abs(np.fft.fftshift(Y_k)))
    plt.grid()
    plt.ylabel("Y(K) with defined filter")

    plt.figure(4)
    plt.plot(y.real)
    plt.grid()
    plt.ylabel("Y(K) with defined filter")
    
    plt.show()
    return y

#read .wav file
input_signal,fs = sf.read('Sound_Noise.wav')

#sampling frequency of Input signal
sampl_freq=fs

#order of the filter
order=4

#cutoff frquency 4kHz
cutoff_freq=4000.0

#digital frequency
Wn=2*cutoff_freq/sampl_freq

# b and a are numerator and denominator polynomials respectively
b, a = signal.butter(order,Wn, 'low')


output_signal = signal.filtfilt(b, a, input_signal)
#output_signal = signal.lfilter(b, a, input_signal)

#write the output signal into .wav file
sf.write('Sound_With_ReducedNoise.wav', output_signal, fs)

plt.figure(1)
plt.plot(output_signal)
plt.grid()
plt.ylabel("y[n] with in-built filter")


plt.figure(2)
plt.plot(np.abs(np.fft.fftshift(np.fft.fft(output_signal))))
plt.grid()
plt.ylabel("Y(K) with in-built filter")

y = My_Filter(a,b,input_signal,Wn)
sf.write('7_1Sound_With_ReducedNoise.wav',np.abs(y), fs)

N = len(input_signal)