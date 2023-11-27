################################################################################
#HEADER

#This example script applies the Polarimetric Bandwidth Extrapolation package
#(PyPBWE) function PyPBWE.PBWE to a synthetic radar signal spectrum.

#The synthetic radar signal example (inspired by the WISDOM GPR of the ExoMars
#rover mission, Ciarletti et al. (2017)):
#   -A FMCW radar working between 0.5 and 3 GHz measures a 1001 frequencies
#    spectrum when sounding.
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

#The Polarimetric Bandwidth Extrapolation (PBWE) is applied to this radar's
#signal using the PyBWE function_PBWE:
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
#    (1992)). A bandwidth X3 yields a resolution X3 better in time-domain.
#   -The extrapolated multi-channel spectrum is eventually converted to a
#    one time-domain signal per channel with IFFT, and zero-padding to
#    interpolate the signal X10. This interpolation is purely aesthetic.

#NB: The PBWE is expected to yield better performances than the classic BWE only
#in cases where targets have different scattering coefficients in the different
#polarimetric channels of the radar.

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

#Generate a sum of two complex sine-waves corresponding to the targets' echoes
#(each with an amplitude of 1), for the 00 channel:
spec_vect_00 = np.exp(-1j*4*pi*dist_target1*freq_vect/3e8)+np.exp(-1j*4*pi*dist_target2*freq_vect/3e8)
#Generate a sum of two complex sine-waves corresponding to the targets' echoes
#(with amplitudes of 1 and -1), for the 11 channel:
spec_vect_11 = np.exp(-1j*4*pi*dist_target1*freq_vect/3e8)-np.exp(-1j*4*pi*dist_target2*freq_vect/3e8)

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

#BWE function application:------------------------------------------------------

#Extrapolation factor (=< 3 recommended by Cuomo (1992)):
extra_factor = 3
#Order of the model as a ratio of the total number of samples in each spectrum
#(= 0.33 recommended by Cuomo (1992)):
model_order = 0.33
#Zero-padding factor (= time domain interpolation factor):
zp_factor = 10
#Cut 5% of samples on each side of the 2 spectra:
side_cut = True

#Application of the PBWE to the 2 polarimetric channels' spectra:
output_pbwe, time_pbwe_vect = PyPBWE.PBWE(spec_mat,freq_vect,extra_factor,model_order,zp_factor,side_cut)

#Figure display:----------------------------------------------------------------

#Calculate the spectrum's frequency step:
df = freq_vect[1]-freq_vect[0]
#Generate a time vector corresponding to the time-domain transform:
time_vect = np.linspace(0,1/df,zp_factor*np.shape(spec_mat)[1])

#Display the original radar sounding and the PBWE version for channel 00:
plt.subplot(1,2,1)
plt.plot(time_vect*1e9,abs(np.fft.fft(np.conjugate(spec_mat[0,:])*np.hamming(len(spec_mat[0,:])),zp_factor*len(spec_mat[0,:])))/len(spec_mat[0,:]),'k-')
plt.plot(time_pbwe_vect*1e9,abs(output_pbwe[0,:]),'r-')
plt.xlim([5,9])
plt.xlabel('Time delays (ns)')
plt.ylabel('Normalized amplitude')
plt.legend(['Original radar sounding','Radar sounding after PBWE'], loc ='best')
plt.title('PBWE application - Polar channel 00')
plt.grid()

#Display the original radar sounding and the PBWE version for channel 00:
plt.subplot(1,2,2)
plt.plot(time_vect*1e9,abs(np.fft.fft(np.conjugate(spec_mat[1,:])*np.hamming(len(spec_mat[1,:])),zp_factor*len(spec_mat[1,:])))/len(spec_mat[1,:]),'k-')
plt.plot(time_pbwe_vect*1e9,abs(output_pbwe[1,:]),'r-')
plt.xlim([5,9])
plt.xlabel('Time delays (ns)')
plt.ylabel('Normalized amplitude')
plt.legend(['Original radar sounding','Radar sounding after PBWE'], loc ='best')
plt.title('PBWE application - Polar channel 11')
plt.grid()

plt.show()