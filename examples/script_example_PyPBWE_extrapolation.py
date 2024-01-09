################################################################################
#HEADER

#This example script applies the Polarimetric Bandwidth Extrapolation package
#(PyPBWE) functions PyPBWE.polar_burg and PyBWE.polar_extrapolation to a
#synthetic radar signal spectrum.

#The synthetic radar signal example (inspired by the WISDOM GPR of the ExoMars
#rover mission, Ciarletti et al. (2017)):
#   -A SFCW (Stepped Frequency Continuous Wave) radar working between 0.5 and 3
#    GHz measures a 1001 frequencies spectrum when sounding.
#   -This radar can transmit and receive signals in two linear polarizations
#    named 0 and 1. This leads to 4 polirization channels: 2 co-polar named 00
#    and 11 (emission and reception in the same polarization), and 2 cross-polar
#    named 01 and 10 (emission and reception in different polarizations). For
#    simplicity, we will only consider the co-polar channels 00 and 11 in this
#    example.
#   -Only the In-phase component (real part of the spectrum) is measured, the
#    Quadrature component (imaginary part of the spectrum) is reconstructed by
#    Hilbert transform.
#   -Two targets in free-space are seperated by 5 cm, slightly below the radar's
#    free-space resolution. These targets generate echoes of equal amplitudes in
#    the radar's signal, or complex sine-waves in the measured spectrum. The two
#    targets are perfect reflectors, the 1st a perfect plate, the 2nd a perfect
#    dihedral corner. Echoes received by the 00 and 11 channels therefore have
#    the same scattering coefficients (1,1) for the 1st reflector, and different
#    coefficients for the 2nd (1,-1).
#   -The measured spectrum is corrupted by a white-noise of standard deviation
#    10X smaller than the complex sine-waves' amplitudes.

#The polarimetric extrapolation used in the PBWE is applied manually to this
#multi-channel spectrum:
#   -Most of the estimation errors when reconstructing a complex spectrum with
#    the Hilbert transform are on the far sides of the spectrum. For this reason,
#    we cut 5% of frequencies on each side of the spectrum before PBWE. We would
#    do the same for a radar working in time-domain, as most errors in FFT
#    estimation are also on each side of the multi-channel spectrum. This
#    process is useless for a radar working in the frequency-domain and
#    measuring both In-Phase and Quadrature components of the spectrum.
#   -We then fit an AR model to the multi-channel spectrum, with an order equal
#    to 1/3 of the spectrum samples, as recommended by Cuomo (1992).
#   -We use this model to extrapolate the spectrum on each side, to obtain a
#    bandwidth 3X larger (maximum extrapolation factor recommended by Cuomo
#    (1992)).

#NB: The polarimetric extrapolation is expected to perform better than the
#classic AR extrapolation only in cases where the targets have different radar
#responses in the different polarimetric channels of the radar.

#References: Suwa and Iwamoto (2003,2007)

################################################################################

#Libraries importation:---------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
from math import pi
from scipy.signal import hilbert

import PyPBWE

#Synthetic data generation:-----------------------------------------------------

#Generate a vector of 1001 frequencies between 0.5 and 3 GHz:
freq_vect = np.linspace(0.5e9,3e9,1001)

#Distance (m) between two radar targets in free-space (returning two echoes) and the radar:
dist_target1 = 1
dist_target2 = 1.07

#Amplitudes of the echoes for each polarimetric channel (00 and 11):
amp_target1_00 = 1
amp_target2_00 = 1
amp_target1_11 = 1
amp_target2_11 = -1

#Generate a sum of two complex sine-waves corresponding to the targets' echoes,
#for the 00 channel:
spec_vect_00 = (amp_target1_00*np.exp(-1j*4*pi*dist_target1*freq_vect/3e8))+(amp_target2_00*np.exp(-1j*4*pi*dist_target2*freq_vect/3e8))
#Generate a sum of two complex sine-waves corresponding to the targets' echoes
#for the 11 channel:
spec_vect_11 = (amp_target1_11*np.exp(-1j*4*pi*dist_target1*freq_vect/3e8))+(amp_target2_11*np.exp(-1j*4*pi*dist_target2*freq_vect/3e8))

#Only keep the real part of the spectrum (In-phase component) for the 2 channels:
spec_vect_00 = np.real(spec_vect_00)
spec_vect_11 = np.real(spec_vect_11)

#Add a white noise (standard deviation 0.1) to each spectrum signal:
wn_vect_00 = np.random.normal(0,0.1,spec_vect_00.shape)
wn_vect_11 = np.random.normal(0,0.1,spec_vect_11.shape)
spec_vect_00 += wn_vect_00
spec_vect_11 += wn_vect_11

#Reconstruct a complex signal with the Hilbert transform for the 2 channels:
spec_vect_00 = np.conjugate(hilbert(spec_vect_00))[::2]
spec_vect_11 = np.conjugate(hilbert(spec_vect_11))[::2]
freq_vect = freq_vect[::2]

#Assemble the 2 spectrum channels into a single matrix:
spec_mat = np.vstack((spec_vect_00,spec_vect_11))

#Multi-channel AR model fit:----------------------------------------------------

#Cutting 5% of frequencies on each side of the spectrum:
spec_mat = spec_mat[:,round(np.shape(spec_mat)[1]*0.05):np.shape(spec_mat)[1]-round(np.shape(spec_mat)[1]*0.05)]
freq_vect = freq_vect[round(len(freq_vect)*0.05):len(freq_vect)-round(len(freq_vect)*0.05)]

#Retrieve the number of polarimetric channels:
Npol = np.shape(spec_mat)[0]
#Retrieve the new number of samples in the spectrum:
M = np.shape(spec_mat)[1]

#Set the multichannel AR model's order as 1/3 of the number of samples:
model_order_nb = round(0.33*M)

#Fit a multichannel AR model to spec_mat:
[Thetaf,Thetab,err] = PyPBWE.polar_burg(spec_mat,model_order_nb)

#Multi-channel spectrum extrapolation:------------------------------------------

#Number of samples to be extrapolated on each side of the multichannel spectrum
#to get an extrapolation factor of 3:
Mextra = round(((3*M)-M)//2)+1

#Extrapolation of spec_mat on each side of the multichannel spectrum:
spec_mat_extra,spec_mat_extra_forward,spec_mat_extra_backward = PyPBWE.polar_extrapolation(spec_mat,Thetaf,Thetab,Mextra,"both")

#Figure display:----------------------------------------------------------------

#Generate the expected spectrum after extrapolation for comparison:
df = freq_vect[1]-freq_vect[0]
freq_vect_expected = np.linspace(freq_vect[0]-(df*(Mextra-2)),freq_vect[-1]+(df*Mextra),M+(Mextra*2))
spec_vect_expected_00 = np.exp(-1j*4*pi*dist_target1*freq_vect_expected/3e8)+np.exp(-1j*4*pi*dist_target2*freq_vect_expected/3e8)
spec_vect_expected_11 = np.exp(-1j*4*pi*dist_target1*freq_vect_expected/3e8)-np.exp(-1j*4*pi*dist_target2*freq_vect_expected/3e8)

#Display the expected spectrum after extrapolation for the 00 and 11
#channels:
plt.subplot(2,1,1)
plt.plot(freq_vect_expected*1e-9,np.real(spec_vect_expected_00),'r-')
plt.plot(freq_vect_expected*1e-9,abs(spec_vect_expected_00),'k-')
plt.axvspan((freq_vect[0]-(df*(Mextra-2)))*1e-9, (freq_vect[0]-df)*1e-9, facecolor='k', alpha=0.1)
plt.axvspan((freq_vect[-1]+df)*1e-9, (freq_vect[-1]+(df*Mextra))*1e-9, facecolor='k', alpha=0.1)
plt.xlim([(freq_vect[0]-(df*(Mextra-2)))*1e-9,(freq_vect[-1]+(df*Mextra))*1e-9])
plt.xlabel('Frequency (GHz)')
plt.ylabel('Amplitude')
plt.legend(['Expected spectrum after extrapolation (real part)','Expected spectrum after extrapolation (modulus)'], loc ='best')
plt.title('Expected spectrum after extrapolation - Polar channel 00')
plt.grid()

plt.subplot(2,1,2)
plt.plot(freq_vect_expected*1e-9,np.real(spec_vect_expected_11),'r-')
plt.plot(freq_vect_expected*1e-9,abs(spec_vect_expected_11),'k-')
plt.axvspan((freq_vect[0]-(df*(Mextra-2)))*1e-9, (freq_vect[0]-df)*1e-9, facecolor='k', alpha=0.1)
plt.axvspan((freq_vect[-1]+df)*1e-9, (freq_vect[-1]+(df*Mextra))*1e-9, facecolor='k', alpha=0.1)
plt.xlim([(freq_vect[0]-(df*(Mextra-2)))*1e-9,(freq_vect[-1]+(df*Mextra))*1e-9])
plt.xlabel('Frequency (GHz)')
plt.ylabel('Amplitude')
plt.legend(['Expected spectrum after extrapolation (real part)','Expected spectrum after extrapolation (modulus)'], loc ='best')
plt.title('Expected spectrum after extrapolation - Polar channel 11')
plt.grid()

plt.show()

#Generate a vector of frequencies for the forward and backward extrapolations:
freq_vect_forward = np.linspace(freq_vect[-1]+df,freq_vect[-1]+(df*Mextra),Mextra)
freq_vect_backward = np.linspace(freq_vect[0]-(df*(Mextra-2)),freq_vect[0]-df,Mextra)

#Display the original radar spectrum and the extrapolation for the 00 and 11
#channels:

plt.subplot(2,1,1)
plt.plot(freq_vect*1e-9,np.real(spec_mat[0,:]),'r-')
plt.plot(freq_vect*1e-9,abs(spec_mat[0,:]),'k-')
plt.plot(freq_vect_backward*1e-9,np.real(spec_mat_extra_backward[0,:]),'y-')
plt.plot(freq_vect_backward*1e-9,abs(spec_mat_extra_backward[0,:]),'g-')
plt.plot(freq_vect_forward*1e-9,np.real(spec_mat_extra_forward[0,:]),'c-')
plt.plot(freq_vect_forward*1e-9,abs(spec_mat_extra_forward[0,:]),'b-')
plt.xlabel('Frequency (GHz)')
plt.ylabel('Amplitude')
plt.legend(['Original radar spectrum (real part)','Original radar spectrum (modulus)','Backward extrapolation (real part)', 'Backward extrapolation (modulus)', 'Forward extrapolation (real part)', 'Forward extrapolation (modulus)'], loc ='best')
plt.title('Spectrum polarimetric extrapolation - Polar channel 00')
plt.grid()

plt.subplot(2,1,2)
plt.plot(freq_vect*1e-9,np.real(spec_mat[1,:]),'r-')
plt.plot(freq_vect*1e-9,abs(spec_mat[1,:]),'k-')
plt.plot(freq_vect_backward*1e-9,np.real(spec_mat_extra_backward[1,:]),'y-')
plt.plot(freq_vect_backward*1e-9,abs(spec_mat_extra_backward[1,:]),'g-')
plt.plot(freq_vect_forward*1e-9,np.real(spec_mat_extra_forward[1,:]),'c-')
plt.plot(freq_vect_forward*1e-9,abs(spec_mat_extra_forward[1,:]),'b-')
plt.xlabel('Frequency (GHz)')
plt.ylabel('Amplitude')
plt.legend(['Original radar spectrum (real part)','Original radar spectrum (modulus)','Backward extrapolation (real part)', 'Backward extrapolation (modulus)', 'Forward extrapolation (real part)', 'Forward extrapolation (modulus)'], loc ='best')
plt.title('Spectrum polarimetric extrapolation - Polar channel 11')
plt.grid()

plt.show()