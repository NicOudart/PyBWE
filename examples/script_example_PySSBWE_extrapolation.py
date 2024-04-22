################################################################################
#HEADER

#This example script applies the State-Space Bandwidth Extrapolation package
#(PySSBWE) functions PySSBWE.statespace_model and
#PySSBWE.statespace_extrapolation to a synthetic radar signal spectrum.

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

#The state-space extrapolation used in the SSBWE is applied manually to this
#spectrum, using the 2 possible methods:
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
#    (1992)).

#NB:
#   -In the following example the order of the model (= number of echoes) is
#    estimated using AIC. The order can also be forced by the user in the
#    PySSBWE functions statespace_model and SSBWE with an optional parameter
#    "order" if the number of echoes is known.
#   -2 methods can be used to fit a state-space model to a spectrum, using
#    either the observability or the controllability matrix. In this example we
#    display only the results obtained with the 1st method.

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

#State-space model fit:---------------------------------------------------------

#Cutting 5% of frequencies on each side of the spectrum:
spec_vect = spec_vect[round(len(spec_vect)*0.05):len(spec_vect)-round(len(spec_vect)*0.05)]
freq_vect = freq_vect[round(len(freq_vect)*0.05):len(freq_vect)-round(len(freq_vect)*0.05)]

#Flip for a backward version of the spectrum:
spec_vect_b = np.flip(spec_vect)

#Retrieve the numbers of samples in the spectrum:
N = len(spec_vect)

#Number of samples to be extrapolated on each side of the spectrum to get an
#extrapolation factor of 3:
Nextra = round(((3*N)-N)//2)+1

#Fit a state-space model to the spectrum for forward extrapolation:
[A1_f,B1_f,C1_f,A2_f,B2_f,C2_f] = PySSBWE.statespace_model(spec_vect,noise_type="white")

#Fit a state-space model to the flipped spectrum for backward extrapolation:
[A1_b,B1_b,C1_b,A2_b,B2_b,C2_b] = PySSBWE.statespace_model(spec_vect_b,noise_type="white")

#State-space extrapolation:-----------------------------------------------------

#Forward extrapolation of the spectrum (obtained with method 1 and 2):
spec_vect_extra_1_f = PySSBWE.statespace_extrapolation(spec_vect,A1_f,B1_f,C1_f,Nextra)
spec_vect_extra_2_f = PySSBWE.statespace_extrapolation(spec_vect,A2_f,B2_f,C2_f,Nextra)

#Backward extrapolation of the spectrum (obtained with method 1 and 2):
spec_vect_extra_1_b = PySSBWE.statespace_extrapolation(spec_vect_b,A1_b,B1_b,C1_b,Nextra)
spec_vect_extra_2_b = PySSBWE.statespace_extrapolation(spec_vect_b,A2_b,B2_b,C2_b,Nextra)

#Flip the backward extrapolation (obtained with method 1 and 2):
spec_vect_extra_1_b = np.flip(spec_vect_extra_1_b)
spec_vect_extra_2_b = np.flip(spec_vect_extra_2_b)

#Figure display:----------------------------------------------------------------

#Generate the expected spectrum after extrapolation for comparison:
df = freq_vect[1]-freq_vect[0]
freq_vect_expected = np.linspace(freq_vect[0]-(df*(Nextra-2)),freq_vect[-1]+(df*Nextra),N+(Nextra*2))
spec_vect_expected = (ampli_target1*np.exp((-decay_target1+(1j*4*pi*dist_target1/3e8))*freq_vect_expected))+(ampli_target2*np.exp((-decay_target2+(1j*4*pi*dist_target2/3e8))*freq_vect_expected))

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
plt.tight_layout()
plt.show()

#Generate a vector of frequencies for the forward and backward extrapolations:
freq_vect_forward = np.linspace(freq_vect[-1]+df,freq_vect[-1]+(df*Nextra),Nextra)
freq_vect_backward = np.linspace(freq_vect[0]-(df*(Nextra-2)),freq_vect[0]-df,Nextra)

#Display the original radar spectrum and the extrapolation using method 1:
plt.plot(freq_vect*1e-9,np.real(spec_vect),'r-')
plt.plot(freq_vect*1e-9,abs(spec_vect),'k-')
plt.plot(freq_vect_backward*1e-9,np.real(spec_vect_extra_1_b),'y-')
plt.plot(freq_vect_backward*1e-9,abs(spec_vect_extra_1_b),'g-')
plt.plot(freq_vect_forward*1e-9,np.real(spec_vect_extra_1_f),'c-')
plt.plot(freq_vect_forward*1e-9,abs(spec_vect_extra_1_f),'b-')
plt.axvspan((freq_vect[0]-(df*(Nextra-2)))*1e-9, (freq_vect[0]-df)*1e-9, facecolor='k', alpha=0.1)
plt.axvspan((freq_vect[-1]+df)*1e-9, (freq_vect[-1]+(df*Nextra))*1e-9, facecolor='k', alpha=0.1)
plt.xlim([(freq_vect[0]-(df*(Nextra-2)))*1e-9,(freq_vect[-1]+(df*Nextra))*1e-9])
plt.xlabel('Frequency (GHz)')
plt.ylabel('Amplitude')
plt.legend(['Original radar spectrum (real part)','Original radar spectrum (modulus)','Backward extrapolation (real part)', 'Backward extrapolation (modulus)', 'Forward extrapolation (real part)', 'Forward extrapolation (modulus)'], loc ='best')
plt.title('Spectrum state-space extrapolation - method 1 (using the observability matrix)')
plt.grid()
plt.tight_layout()
plt.show()

#Display the original radar spectrum and the extrapolation using method 2:
plt.plot(freq_vect*1e-9,np.real(spec_vect),'r-')
plt.plot(freq_vect*1e-9,abs(spec_vect),'k-')
plt.plot(freq_vect_backward*1e-9,np.real(spec_vect_extra_2_b),'y-')
plt.plot(freq_vect_backward*1e-9,abs(spec_vect_extra_2_b),'g-')
plt.plot(freq_vect_forward*1e-9,np.real(spec_vect_extra_2_f),'c-')
plt.plot(freq_vect_forward*1e-9,abs(spec_vect_extra_2_f),'b-')
plt.axvspan((freq_vect[0]-(df*(Nextra-2)))*1e-9, (freq_vect[0]-df)*1e-9, facecolor='k', alpha=0.1)
plt.axvspan((freq_vect[-1]+df)*1e-9, (freq_vect[-1]+(df*Nextra))*1e-9, facecolor='k', alpha=0.1)
plt.xlim([(freq_vect[0]-(df*(Nextra-2)))*1e-9,(freq_vect[-1]+(df*Nextra))*1e-9])
plt.xlabel('Frequency (GHz)')
plt.ylabel('Amplitude')
plt.legend(['Original radar spectrum (real part)','Original radar spectrum (modulus)','Backward extrapolation (real part)', 'Backward extrapolation (modulus)', 'Forward extrapolation (real part)', 'Forward extrapolation (modulus)'], loc ='best')
plt.title('Spectrum state-space extrapolation - method 2 (using the controllability matrix)')
plt.grid()
plt.tight_layout()
plt.show()