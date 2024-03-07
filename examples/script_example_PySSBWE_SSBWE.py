################################################################################
#HEADER

#This example script applies the State-Space Bandwidth Extrapolation package
#(PySSBWE) functions PySSBWE.SSBWE and PySSBWE.statespace_properties to a
#synthetic radar signal spectrum.

#The synthetic radar signal example (inspired by the WISDOM GPR of the ExoMars
#rover mission, Ciarletti et al. (2017)):
#   -A SFCW (Stepped Frequency Continuous Wave) radar working between 0.5 and 3
#    GHz measures a 1001 frequencies spectrum when sounding.
#   -Only the In-phase component (real part of the spectrum) is measured, the
#    Quadrature component (imaginary part of the spectrum) is reconstructed by
#    Hilbert transform.
#   -Two targets in free-space are seperated by 5 cm, slightly below the radar's
#    free-space resolution. With these targets' echoes (or complex sine waves in
#    frequency domain) are associated different complex amplitudes and different
#    decays in frequency-domain.
#   -The measured spectrum is corrupted by a white-noise of standard deviation
#    10X smaller than the complex sine-waves' amplitudes.

#The State-Space Bandwidth Extrapolation (SSBWE) is applied to this radar's
#signal using the PyBWE function_SSBWE:
#   -Most of the estimation errors when reconstructing a complex spectrum with
#    the Hilbert transform are on the far sides of the spectrum. For this reason,
#    we cut 5% of frequencies on each side of the spectrum before SSBWE. We
#    would do the same for a radar working in time-domain, as most errors in FFT
#    estimation are also on each side of the spectrum. This process is useless
#    for a radar working in the frequency-domain and measuring both In-Phase and
#    Quadrature components of the spectrum.
#   -We then fit a state-space model to the spectrum, first in the forward and
#    then the backward directions.
#   -We use these 2 models to extrapolate the spectrum on each side, to obtain a
#    bandwidth 3X larger (maximum extrapolation factor recommended by Cuomo
#    (1992)). A bandwidth X3 yields a resolution X3 better in time-domain.
#   -The extrapolated spectrum is eventually converted to a time-domain signal
#    by IFFT, with zero-padding to interpolate the signal X10. This
#    interpolation is purely aesthetic.

#From the state-space model of the spectrum, the properties of each echo
#(amplitude, time-delay, frequency-domain decay) can be estimated. This is done
#in the last part of this example.

#NB:
#   -In the following example the order of the model (= number of echoes) is
#    estimated using AIC. The order can also be forced by the user in the
#    PySSBWE functions statespace_model and SSBWE with an optional parameter
#    "order" if the number of echoes is known.
#   -2 methods can be used to fit a state-space model to a spectrum, using
#    either the observability or the controllability matrix. In this example we
#    display the results obtained with the 2 methods.

#References: Piou (1999)

################################################################################

#Libraries importation:---------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
from math import pi
from scipy.signal import hilbert

import PySSBWE

#Synthetic data generation:-----------------------------------------------------

#Generate a vector of 1001 frequencies between 0.5 and 3 GHz:
freq_vect = np.linspace(0.5e9,3e9,1001)

#Distance (m) between two radar targets in free-space (returning two echoes) and the radar:
dist_target1 = 1
dist_target2 = 1.07

#Amplitude associated with each target echo:
ampli_target1 = 1
ampli_target2 = 0.75

#Frequency decay associated with each target echo (1/Hz):
decay_target1 = 0.25e-9
decay_target2 = 0.5e-9

#Generate a sum of two complex sine-waves corresponding to the targets' echoes
#(each with an amplitude of 1):
spec_vect = (ampli_target1*np.exp((-decay_target1+(1j*4*pi*dist_target1/3e8))*freq_vect))+(ampli_target2*np.exp((-decay_target2+(1j*4*pi*dist_target2/3e8))*freq_vect))

#Only keep the real part of the spectrum (In-phase component):
spec_vect = np.real(spec_vect)

#Add a white noise (standard deviation 0.1) to this spectrum signal:
wn_vect = np.random.normal(0,0.1,spec_vect.shape)
spec_vect += wn_vect

#Reconstruct a complex signal with the Hilbert transform:
spec_vect = hilbert(spec_vect)[::2]
freq_vect = freq_vect[::2]

#SSBWE function application:----------------------------------------------------

#Extrapolation factor (=< 3 recommended by Cuomo (1992)):
extra_factor = 3
#Zero-padding factor (= time domain interpolation factor):
zp_factor = 10
#Cut 5% of samples on each side of the spectrum:
side_cut = True

#Calculate the spectrum's frequency step:
df = freq_vect[1]-freq_vect[0]

#Application of the BWE to the spectrum:
output_ssbwe_1, output_ssbwe_2, time_ssbwe_vect = PySSBWE.SSBWE(spec_vect,df,extra_factor,zp_factor,side_cut)

#Figure display:----------------------------------------------------------------

#Generate a time vector corresponding to the time-domain transform:
time_vect = np.linspace(0,1/df,zp_factor*len(spec_vect))

#Display the original radar sounding and the SSBWE version using method 1:
plt.plot(time_vect*1e9,abs(1.85*np.fft.fft(spec_vect*np.hamming(len(spec_vect)),zp_factor*len(spec_vect)))/len(spec_vect),'k-')
plt.plot(time_ssbwe_vect*1e9,abs(output_ssbwe_1),'r-')
plt.xlim([5,9])
plt.xlabel('Time delays (ns)')
plt.ylabel('Normalized amplitude')
plt.legend(['Original radar sounding','Radar sounding after SSBWE'], loc ='best')
plt.title('SSBWE application - method 1 (using the observability matrix)')
plt.grid()
plt.show()

#Display the original radar sounding and the SSBWE version using method 2:
plt.plot(time_vect*1e9,abs(1.85*np.fft.fft(spec_vect*np.hamming(len(spec_vect)),zp_factor*len(spec_vect)))/len(spec_vect),'k-')
plt.plot(time_ssbwe_vect*1e9,abs(output_ssbwe_2),'r-')
plt.xlim([5,9])
plt.xlabel('Time delays (ns)')
plt.ylabel('Normalized amplitude')
plt.legend(['Original radar sounding','Radar sounding after SSBWE'], loc ='best')
plt.title('SSBWE application - method 2 (using the controllability matrix)')
plt.grid()
plt.show()

#Echoes properties estimation:--------------------------------------------------

#Manually cut 5% of frequencies of each side of the spectrum:
spec_vect_cut = spec_vect[round(len(spec_vect)*0.05):len(spec_vect)-round(len(spec_vect)*0.05)]
freq_vect_cut = freq_vect[round(len(freq_vect)*0.05):len(freq_vect)-round(len(freq_vect)*0.05)]

#Fit a state-space model to the spectrum:
[A1,B1,C1,A2,B2,C2] = PySSBWE.statespace_model(spec_vect_cut)

#Retrieve the spectrum's frequency step and 1st frequency:
df = freq_vect_cut[1]-freq_vect_cut[0]
f1 = freq_vect_cut[0]

#Extract the echoes properties from the state-space model obtained with method 1:
[amp_1,td_1,dec_1] = PySSBWE.statespace_properties(A1,B1,C1,df,f1)

#Extract the echoes properties from the state-space model obtained with method 2:
[amp_2,td_2,dec_2] = PySSBWE.statespace_properties(A2,B2,C2,df,f1)

#Print the echoes amplitudes (module), range (m) and frequency decays (1/GHz)
#obtained using method 1:
print("---Echoes properties estimation (method 1: using observability matrix)---")
print("Module amplitude of echoes: "+str(abs(amp_1)))
print("Distance of targets from the radar (m): "+str(0.5*td_1*3e8))
print("Frequency-domain decay of echoes (1/GHz): "+str(dec_1*1e9))
print("-------------------------------------------------------------------------")

#Print the echoes amplitudes (module), range (m) and frequency decays (1/GHz)
#obtained using method 2:
print("---Echoes properties estimation (method 2: using controllability matrix)---")
print("Module amplitude of echoes: "+str(abs(amp_2)))
print("Distance of targets from the radar (m): "+str(0.5*td_2*3e8))
print("Frequency-domain decay of echoes (1/GHz): "+str(dec_2*1e9))
print("---------------------------------------------------------------------------")