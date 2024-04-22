################################################################################
#HEADER

#This example script applies the Bandwidth Extrapolation package (PyBWE)
#function PyBWE.BWE to a synthetic radar signal spectrum.

#The synthetic radar signal example (inspired by the WISDOM GPR of the ExoMars
#rover mission, Ciarletti et al. (2017)):
#   -A SFCW (Stepped Frequency Continuous Wave) radar working between 0.5 and 3
#    GHz measures a 1001 frequencies spectrum when sounding.
#   -Only the In-phase component (real part of the spectrum) is measured, the
#    Quadrature component (imaginary part of the spectrum) is reconstructed by
#    Hilbert transform.
#   -Two targets in free-space are seperated by 5 cm, slightly below the radar's
#    free-space resolution. These targets generate echoes of given complex
#    amplitudes in the radar's signal, or complex sine-waves in the measured
#    spectrum.
#   -The measured spectrum is corrupted by a white-noise of standard deviation
#    10X smaller than the complex sine-waves' amplitudes.

#The Bandwidth Extrapolation (BWE) is applied to this radar's signal using the
#PyBWE function_BWE:
#   -Most of the estimation errors when reconstructing a complex spectrum with
#    the Hilbert transform are on the far sides of the spectrum. For this reason,
#    we cut 5% of frequencies on each side of the spectrum before BWE. We would
#    do the same for a radar working in time-domain, as most errors in FFT
#    estimation are also on each side of the spectrum. This process is useless
#    for a radar working in the frequency-domain and measuring both In-Phase and
#    Quadrature components of the spectrum.
#   -We then fit an AR model to the spectrum, with an order equal to 1/3 of the
#    spectrum samples, as recommended by Cuomo (1992).
#   -We use this model to extrapolate the spectrum on each side, to obtain a
#    bandwidth 3X larger (maximum extrapolation factor recommended by Cuomo
#    (1992)). A bandwidth X3 yields a resolution X3 better in time-domain.
#   -The extrapolated spectrum is eventually converted to a time-domain signal
#    by IFFT, with zero-padding to interpolate the signal X10. This
#    interpolation is purely aesthetic.

#References: Bowling (1977), Cuomo (1992)

################################################################################

#Libraries importation:---------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
from math import pi
from scipy.signal import hilbert

import PyBWE

#Synthetic data generation:-----------------------------------------------------

#Generate a vector of 1001 frequencies between 0.5 and 3 GHz:
freq_vect = np.linspace(0.5e9,3e9,1001)

#Distance (m) between two radar targets in free-space (returning two echoes) and the radar:
dist_target1 = 1
dist_target2 = 1.07

#Amplitude of the 2 echoes corresponding to the 2 targets:
ampli_target1 = 1
ampli_target2 = 1

#Generate a sum of two complex sine-waves corresponding to the targets' echoes:
spec_vect = (ampli_target1*np.exp(1j*4*pi*dist_target1*freq_vect/3e8))+(ampli_target2*np.exp(1j*4*pi*dist_target2*freq_vect/3e8))

#Only keep the real part of the spectrum (In-phase component):
spec_vect = np.real(spec_vect)

#Add a white noise (standard deviation 0.1) to this spectrum signal:
wn_vect = np.random.normal(0,0.1,spec_vect.shape)
spec_vect += wn_vect

#Reconstruct a complex signal with the Hilbert transform:
spec_vect = hilbert(spec_vect)[::2]
freq_vect = freq_vect[::2]

#BWE function application:------------------------------------------------------

#Extrapolation factor (=< 3 recommended by Cuomo (1992)):
extra_factor = 3
#Order of the model as a ratio of the total number of samples in the spectrum
#(= 0.33 recommended by Cuomo (1992)):
model_order = 0.33
#Zero-padding factor (= time domain interpolation factor):
zp_factor = 10
#Cut 5% of samples on each side of the spectrum:
side_cut = True

#Calculate the original spectrum's frequency step:
df = freq_vect[1]-freq_vect[0]

#Application of the BWE to the spectrum:
output_bwe, time_bwe_vect = PyBWE.BWE(spec_vect,df,extra_factor,model_order,zp_factor,side_cut)

#Figure display:----------------------------------------------------------------

#Generate a time vector corresponding to the time-domain transform:
time_vect = np.linspace(0,1/df,zp_factor*len(spec_vect))

#Display the original radar sounding and the BWE version:
plt.figure(figsize=(10,5))
plt.plot(time_vect*1e9,abs(1.85*np.fft.fft(spec_vect*np.hamming(len(spec_vect)),zp_factor*len(spec_vect)))/len(spec_vect),'k-')
plt.plot(time_bwe_vect*1e9,abs(output_bwe),'r-')
plt.xlim([5,9])
plt.ylim([-0.05,1.3])
plt.xlabel('Time delays (ns)')
plt.ylabel('Normalized amplitude')
plt.legend(['Original radar sounding','Radar sounding after BWE'], loc ='best')
plt.title('BWE application')
plt.grid()
plt.tight_layout()
plt.show()