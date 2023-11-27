################################################################################
#HEADER

#This example script applies the Bandwidth Extrapolation package (PyBWE)
#functions PyBWE.burg and PyBWE.ar_extrapolation to a synthetic radar signal
#spectrum.

#The synthetic radar signal example (inspired by the WISDOM GPR of the ExoMars
#rover mission, Ciarletti et al. (2017)):
#   -A FMCW radar working between 0.5 and 3 GHz measures a 1001 frequencies
#    spectrum when sounding.
#   -Only the In-phase component (real part of the spectrum) is measured, the
#    Quadrature component (imaginary part of the spectrum) is reconstructed by
#    Hilbert transform.
#   -Two targets in free-space are seperated by 5 cm, slightly below the radar's
#    free-space resolution. These targets generate echoes of equal amplitudes in
#    the radar's signal, or complex sine-waves in the measured spectrum.
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

#Generate a sum of two complex sine-waves corresponding to the targets' echoes
#(each with an amplitude of 1):
spec_vect = np.exp(-1j*4*pi*dist_target1*freq_vect/3e8)+np.exp(-1j*4*pi*dist_target2*freq_vect/3e8)

#Only keep the real part of the spectrum (In-phase component):
spec_vect = np.real(spec_vect)

#Add a white noise (standard deviation 0.1) to this spectrum signal:
wn_vect = np.random.normal(0,0.1,spec_vect.shape)
spec_vect += wn_vect

#Reconstruct a complex signal with the Hilbert transform:
spec_vect = np.conjugate(hilbert(spec_vect))[::2]
freq_vect = freq_vect[::2]

#AR model fit:------------------------------------------------------------------

#Cutting 5% of frequencies on each side of the spectrum:
spec_vect = spec_vect[round(len(spec_vect)*0.05):len(spec_vect)-round(len(spec_vect)*0.05)]
freq_vect = freq_vect[round(len(freq_vect)*0.05):len(freq_vect)-round(len(freq_vect)*0.05)]

#Retrieve the number of samples in the spectrum:
N = len(spec_vect)

#Set the AR model's order as 1/3 of the number of samples:
model_order_nb = round(0.33*N)

#Fit an AR model to the spectrum using the Burg algorithm:
ar_coeff,ar_var,ar_rc = PyBWE.burg(spec_vect,model_order_nb)

#Spectrum extrapolation:--------------------------------------------------------

#Number of samples to be extrapolated on each side of the spectrum to get an
#extrapolation factor of 3:
Nextra = round(((3*N)-N)//2)+1

#Extrapolation on each side of the spectrum:
spec_vect_extra,spec_vect_forward,spec_vect_backward = PyBWE.ar_extrapolation(spec_vect,ar_coeff,Nextra,"both")

#Figure display:----------------------------------------------------------------

#Generate the expected spectrum after extrapolation for comparison:
df = freq_vect[1]-freq_vect[0]
freq_vect_expected = np.linspace(freq_vect[0]-(df*(Nextra-2)),freq_vect[-1]+(df*Nextra),N+(Nextra*2))
spec_vect_expected = np.exp(-1j*4*pi*dist_target1*freq_vect_expected/3e8)+np.exp(-1j*4*pi*dist_target2*freq_vect_expected/3e8)

#Display the expected spectrum after extrapolation:
plt.plot(freq_vect_expected*1e-9,np.real(spec_vect_expected),'r-')
plt.plot(freq_vect_expected*1e-9,abs(spec_vect_expected),'k-')
plt.axvspan((freq_vect[0]-(df*(Nextra-2)))*1e-9, (freq_vect[0]-df)*1e-9, facecolor='k', alpha=0.1)
plt.axvspan((freq_vect[-1]+df)*1e-9, (freq_vect[-1]+(df*Nextra))*1e-9, facecolor='k', alpha=0.1)
plt.xlim([(freq_vect[0]-(df*(Nextra-2)))*1e-9,(freq_vect[-1]+(df*Nextra))*1e-9])
plt.xlabel('Frequency (GHz)')
plt.ylabel('Amplitude')
plt.legend(['Expected spectrum after extrapolation (real part)','Expected spectrum after extrapolation (modulus)'], loc ='best')
plt.title('Expected spectrum after extrapolation')
plt.grid()
plt.show()

#Generate a vector of frequencies for the forward and backward extrapolations:
freq_vect_forward = np.linspace(freq_vect[-1]+df,freq_vect[-1]+(df*Nextra),Nextra)
freq_vect_backward = np.linspace(freq_vect[0]-(df*(Nextra-2)),freq_vect[0]-df,Nextra)

#Display the original radar spectrum and the extrapolation:
plt.plot(freq_vect*1e-9,np.real(spec_vect),'r-')
plt.plot(freq_vect*1e-9,abs(spec_vect),'k-')
plt.plot(freq_vect_backward*1e-9,np.real(spec_vect_backward),'y-')
plt.plot(freq_vect_backward*1e-9,abs(spec_vect_backward),'g-')
plt.plot(freq_vect_forward*1e-9,np.real(spec_vect_forward),'c-')
plt.plot(freq_vect_forward*1e-9,abs(spec_vect_forward),'b-')
plt.axvspan((freq_vect[0]-(df*(Nextra-2)))*1e-9, (freq_vect[0]-df)*1e-9, facecolor='k', alpha=0.1)
plt.axvspan((freq_vect[-1]+df)*1e-9, (freq_vect[-1]+(df*Nextra))*1e-9, facecolor='k', alpha=0.1)
plt.xlim([(freq_vect[0]-(df*(Nextra-2)))*1e-9,(freq_vect[-1]+(df*Nextra))*1e-9])
plt.xlabel('Frequency (GHz)')
plt.ylabel('Amplitude')
plt.legend(['Original radar spectrum (real part)','Original radar spectrum (modulus)','Backward extrapolation (real part)', 'Backward extrapolation (modulus)', 'Forward extrapolation (real part)', 'Forward extrapolation (modulus)'], loc ='best')
plt.title('Spectrum AR-extrapolation')
plt.grid()
plt.show()