#Code by P.Aditya
#21st March 2021

#Importing the output of c code and constructing the audio file
import numpy as np
import matplotlib.pyplot as plt
import soundfile as sf

x,fs = sf.read('Sound_Noise.wav')

output_val = np.loadtxt("y_fft.dat")
output_val = output_val[:,0]

sf.write('7.1_Sound_With_ReducedNoise_using_c.wav',output_val[0:len(x)],fs)
plt.figure(1)
plt.plot(output_val[0:len(x)],'r')
plt.grid()
plt.title("Output from fft algorithm using c")
plt.savefig("../figs/ee18btech11035_c.eps")
plt.show()

