# PyPBWE package

## Description

The PyPBWE package is an "advanced" package of this library.
It allows you to apply the Polarimetric Bandwidth Extrapolation (PBWE) technique as defined by Suwa and Iwamoto (2003,2007) to any set of multi-polarimetric-channel radar spectra.
The PBWE is a super-resolution technique, as it yields a better range resolution in time-domain than classic spectral estimation techniques (FFT).
It is a multi-channel version of the classic Bandwidth Extrapolation (BWE), taking advantage of the polarimetric capabilities of a radar.

The PBWE works this way:

* The input is a set of frequency-domain spectra of the radar's reponse signals from a series of targets, received by different polarimetric channels of the instrument. In this situation, each channel signal should in theory be a sum of complex sine-waves.

* A set of autoregressive (AR) models is fitted to this set of spectra, using a multi-channel version of the Burg algorithm. The order of the model is a user-defined parameter.

* These models are used to extrapolate each spectrum forward and backward. The extrapolating factor is a user-defined parameter.

* The super-resolved time-domain soundings for each channel are obtained by IFFT. Zero-padding can be used for time-domain interpolation (this process is purely aesthetic).

Each spectrum's extrapolation factor is equal to the resolution enhancement factor after IFFT.

The multi-channel version of the Burg algorithm used by the PBWE accounts for the targets' echoes information present in each polarimetric channel: the principle of this technique consists in extrapolating each polarimetric channel with a weighted sum of previous/next samples not only using this channel but also the others.
For this reason, in the case where a target's scattering coefficients are different and non-zero for at least two polarimetric channels of the radar, we expect the PBWE to perform better than the BWE.

The PyPBWE package contains the **PyPBWE.PBWE** function, allowing you to apply the PBWE directly to a given set of multi-polarimetric-channel radar spectra and get a super-resolved sounding per channel.
This function calls several other functions from the package, that you can call independently if needed:

* **PyPBWE.polar_burg:** fits an AR model to each polarimetric channel of a spectrum, using a multi-channel version of the Burg algorithm.

* **PyPBWE.polar_extrapolation:** extrapolates forward, backward or both a multi-polarimetric-channel spectrum, given a set of AR models.

## Importation

Once the library is installed, the PyPBWE package can simply be imported:

~~~bash
import PyPBWE
~~~

## Recommendations

### The order of the AR model

In theory, the ideal order of an AR model depends on the number of complex sine-waves composing the spectrum.
However, in practical cases, a higher order is required to obtain good results in presence of noise.
For the BWE, Cuomo (1992) recommends an order equal to 1/3 of the number of samples in the spectrum. This is also what we recommend for PBWE.

### The extrapolation factor

In theory, a sum of sine-wave can be extrapolated infinitely.
However, in practical cases, errors accumulate with the extrapolation.
For the BWE, Cuomo (1992) recommends an extrapolation factor of maximum 3. 
Higher values can be tried for the PBWE, as the model is expected to be more robust, but the results must be interpreted cautiously.

### The effect of side-lobes

In the case of a radar signal spectrum obtained by FFT, 2 effects must be mitigated before PBWE:

* If a window was applied to any signal channel, it must be inverted.

* The application of the FFT generates side-lobes effects on the borders of a spectrum. For this reason, Cuomo (1992) recommends to remove 5% of frequencies on each side of a spectrum before PBWE. The same goes for PBWE.

In the case of a radar measuring only the real-part of each frequency spectrum, the imaginary part can be reconstructed by Hilbert transform.
The FFT being part of Hilbert transform implementations, errors can also be expected on the border of reconstructed spectra.
We thus also recommend to cut 5% of samples on each side of each polarimetric channel spectrum in this case.

The PyPBWE.PBWE function allows you to remove 5% of samples on each side of a spectrum before PBWE.

### The effect of distortions

The PBWE technique is based on the assumption that a target response signal will be a complex sine-wave of the same frequency (but potentially different amplitudes) in all polarimetric channels of a radar.
In many practical applications, a radar response signal can differ from this model, due to various distortion effects.
Such distortions can make the application of the direct application of the PBWE impossible.
We thus recommend to correct distortions prior to a PBWE application.

Different solutions are possible for different sources of distortion:

* **Effect of the antennas gain:**
No antenna can transmit all frequencies with the same gain. This will cause distortions of sine-waves in a radar spectrum.
You can try to compensate for this effect with a calibration radar measurement on a "perfect reflector" (a metallic sphere or a large metallic plate).
From the spectrum of the echo from such target, the antennas gain can be determined, and inverted in all other radar measurements.
This process is sometimes called "whitening" and is described in Oudart et al. (2021).
As different temperatures can lead to significantly different amplitudes over a radar's bandwidth, for a successful whitening process we recommend to apply first a temperature correction (as the one described in Hervé et al. (2020)).
If necessary, we recommend to apply this process for all polarimetric channels of the radar.

* **Effect of the ionosphere:**
The ionosphere of any planet has a frequency-dependent attenuation effect on electromagnetic waves. This will cause distortions of sine-waves in a radar spectrum. 
This effect can be compensated using a ionospheric model.
For instance Gambacorta et al. (2022) used a ionospheric Gamma model before applying the PBWE to MARSIS (Mars Express) soundings.
If necessary, we recommend to apply this process for all polarimetric channels of the radar.

* **Effect of subsurface attenuations:**
When a radar signal is transmitted through a sub-surface, it is affected by different frequency-dependent attenuations (absorption, scattering...).
It can thus be expected that deeper echoes in a high-loss subsurface will not be correctly reconstructed by PBWE.
We recommend to use the PySSBWE package, as the SSBWE (State-Space Bandwidth Extrapolation) technique accounts for such attenuations.

### The effect of noise

The presence of noise in a signal impact the quality of the AR model determined by the Burg algorithm.
For this reason, a high level of noise will impact the BWE results (echoes amplitudes and time-delays), as demonstrated in Oudart et al. (2021).
We expect this effect to be lower for the PBWE than for the BWE in the case of targets having different non-zero scattering coefficients in at least two polarimetric channels of the radar.

### Electromagnetic interferences

EM interferences can cause peaks in the measured radar spectrum. Such peaks will impair the AR modelling of the signal.
The effect can be mitigated by detecting outliers in the spectrum, and then replacing them by interpolation, as proposed by Raguso et al. (2018).
If necessary, we recommend to apply this process for all polarimetric channels of the radar using the multi-channel Burg algorithm and extrapolation.

### Cross-talk

If the emission and reception antennas / channels are not the same, a cross-talk phenomenon can generate unwanted signals in the radar response.
It is recommended to remove such effects prior to the application of the BWE.
This can for instance be done by subtracting a "free-space" (without any targets) calibration signal to any radar response.
We thus recommend to acquire a "free-space" measurement for each polarimetric channel of the radar, and to subtract it to any radar response coming from this channel.

## Functions

### PyPBWE.PBWE(spec_mat,df,extra_factor,model_order,zp_factor,side_cut=True)

The PyPBWE.PBWE function applies the Polarimetric Bandwidth Extrapolation to a radar polarimetric channels spectra.

**Inputs:**

* spec_mat: _complex 2D array_, each row containing the spectrum of a polarimetric channel to be extrapolated.

* df: _float_ frequency step corresponding to spec_mat [Hz].

* extra_factor: _float_ factor between the extrapolated and the original spectra's bandwidths.

* model_order: _float_ order of the AR models for extrapolation, expressed as a ratio of the original spectra's bandwidth.

* zp_factor: _float_ factor between the zero-padded and original spectra's bandwidth.

* side_cut: (optional) _boolean_, 5% of samples are cut on each side of each spectra if True.

**Outputs:**

* output_pbwe: _complex 2D array_, each row containing the time-domain radar signal of a polarimetric channel after PBWE.

* time_pbwe_vect: _float 1D array_ containing the time axis corresponding to output_pbwe [s].

**Notes:**
In this implementation, the time-domain output signal is scaled by the number of samples in each polarimetric channel's spectrum.
A Hamming window is applied to each spectrum before IFFT.

### PyPBWE.polar_burg(X,p)

The PyPBWE.polar_burg function fits an AR model to each polarimetric channel's spectrum, using a multi-channel version of the Burg algorithm.

**Inputs:**

* X: _complex 2D array_, each row containing the spectrum of a polarimetric channel to be modelled.

* p: _integer_ order of the AR models.

**Outputs:**

* Thetaf: _complex 2D array_ containing the coefficients of the multi-channel AR forward model.

* Thetab: _complex 2D array_ containing the coefficients of the multi-channel AR backward model.

* err: _float 1D array_ of prediction errors (forward/backward). 

### PyPBWE.polar_extrapolation(X,Thetaf,Thetab,Mextra,extra_mode='both')

The PyPBWE.polar_extrapolation function extrapolates a set of polarimetric channel spectra given their multi-channel AR model, forward, backward or both.

**Inputs:**

* X: _complex 2D array_, each row containing the spectrum of a polarimetric channel to be extrapolated.

* Thetaf: _complex 2D array_ containing the coefficients of the multi-channel AR forward model.

* Thetab: _complex 2D array_ containing the coefficients of the multi-channel AR backward model.

* Mextra: _integer_ number of samples to be extrapolated forward, backward or both for each channel.

* extra_mode: (optional) _string_ containing the extrapolation mode, "forward"/"backward"/"both".

**Outputs:**

* X_extra: _complex 2D array_, each row containing the extrapolated spectrum of a polarimetric channel.

* X_forward: _complex 2D array_, each row containing the forward extrapolation of a polarimetric channel.

* X_backward: _complex 2D array_, each row containing the backward extrapolation of a polarimetric channel.

## Example

2 example scripts are proposed with the same synthetic radar scenario for the PyBWE package:

* An application of PyPBWE.PBWE on a synthetic radar signal is proposed in _examples/script_example_PyPBWE_PBWE.py_.

* A manual use of PyPBWE.polar_burg and PyPBWE.polar_extrapolation on a synthetic radar signal is proposed in _examples/script_example_PyPBWE_extrapolation.py_.

In the following, these 2 scripts are explained as a single tutorial:

### Presentation of the scenario

The synthetic radar signal example (inspired by the WISDOM GPR of the ExoMars rover mission, Ciarletti et al. (2017)):

* A SFCW (Stepped Frequency Continuous Wave) radar working between 0.5 and 3 GHz measures a 1001 frequencies spectrum when sounding.

* This radar can transmit and receive signals in two linear polarizations named 0 and 1. This leads to 4 polirization channels: 2 co-polar named 00 and 11 (emission and reception in the same polarization), and 2 cross-polar named 01 and 10 (emission and reception in different polarizations). For simplicity, we will only consider the co-polar channels 00 and 11 in this example.

* Only the In-phase component (real part of the spectrum) is measured, the Quadrature component (imaginary part of the spectrum) is reconstructed by Hilbert transform.

* Two targets in free-space are seperated by 5 cm, slightly below the radar's free-space resolution. These targets generate echoes of given complex amplitudes in the radar's signal, or complex sine-waves in the measured spectrum.

* The measured spectrum is corrupted by a white-noise of standard deviation 10X smaller than the complex sine-waves' amplitudes.

The Polarimetric Bandwidth Extrapolation (PBWE) is applied to this radar's signal using the PyBWE function_PBWE:

* Most of the estimation errors when reconstructing a complex spectrum with the Hilbert transform are on the far sides of the spectrum. For this reason, we cut 5% of frequencies on each side of the spectrum before PBWE. We would do the same for a radar working in time-domain, as most errors in FFT estimation are also on each side of the multi-channel spectrum. This process is useless for a radar working in the frequency-domain and measuring both In-Phase and Quadrature components of the spectrum.

* We then fit an AR model to the multi-channel spectrum, with an order equal to 1/3 of the spectrum samples, as recommended by Cuomo (1992).

* We use this model to extrapolate the spectrum on each side, to obtain a bandwidth 3X larger (maximum extrapolation factor recommended by Cuomo (1992)). A bandwidth X3 yields a resolution X3 better in time-domain.

* The extrapolated multi-channel spectrum is eventually converted to a one time-domain signal per channel with IFFT, and zero-padding to interpolate the signal X10. This interpolation is purely aesthetic.

This scenario requires the following libraries:

~~~bash
import matplotlib.pyplot as plt
import numpy as np
from math import pi
from scipy.signal import hilbert

import PyPBWE
~~~

### The synthetic signal generation

In this section we describe how to generate analytically a synthetic radar sounding to test the PyPBWE package.

Generate a vector of 1001 frequencies between 0.5 and 3 GHz:
~~~bash
freq_vect = np.linspace(0.5e9,3e9,1001)
~~~

Distance (m) between two radar targets in free-space (returning two echoes) and the radar:
~~~bash
dist_target1 = 1
dist_target2 = 1.07
~~~

Amplitude associated with each target echo for each polarimetric channel:
~~~bash
amp_target1_00 = 1
amp_target2_00 = 1
amp_target1_11 = 1
amp_target2_11 = -1
~~~

Generate a sum of two complex sine-waves corresponding to the targets' echoes (each with an amplitude of 1), for the 00 channel:
~~~bash
spec_vect_00 = (amp_target1_00*np.exp(-1j*4*pi*dist_target1*freq_vect/3e8))+(amp_target2_00*np.exp(-1j*4*pi*dist_target2*freq_vect/3e8))
~~~

Generate a sum of two complex sine-waves corresponding to the targets' echoes (with amplitudes of 1 and -1), for the 11 channel:
~~~bash
spec_vect_11 = (amp_target1_11*np.exp(-1j*4*pi*dist_target1*freq_vect/3e8))+(amp_target2_11*np.exp(-1j*4*pi*dist_target2*freq_vect/3e8))
~~~

Only keep the real part of the spectrum (In-phase component) for the 2 channels:
~~~bash
spec_vect_00 = np.real(spec_vect_00)
spec_vect_11 = np.real(spec_vect_11)
~~~

Add a white noise (standard deviation 0.1) to each spectrum signal:
~~~bash
wn_vect_00 = np.random.normal(0,0.1,spec_vect_00.shape)
wn_vect_11 = np.random.normal(0,0.1,spec_vect_11.shape)
spec_vect_00 += wn_vect_00
spec_vect_11 += wn_vect_11
~~~

Reconstruct a complex signal with the Hilbert transform for the 2 channels:
~~~bash
spec_vect_00 = np.conjugate(hilbert(spec_vect_00))[::2]
spec_vect_11 = np.conjugate(hilbert(spec_vect_11))[::2]
freq_vect = freq_vect[::2]
~~~

Assemble the 2 spectrum channels into a single matrix:
~~~bash
spec_mat = np.vstack((spec_vect_00,spec_vect_11))
~~~

This is the matrix of signals to which we will apply the PBWE.

### The application of the PBWE:

The PBWE can be applied directly to the multi-channel spectrum using the PyPBWE.PBWE function.

Extrapolation factor (=< 3 recommended by Cuomo (1992)):
~~~bash
extra_factor = 3
~~~

Order of the model as a ratio of the total number of samples in each spectrum (= 0.33 recommended by Cuomo (1992)):
~~~bash
model_order = 0.33
~~~

Zero-padding factor (= time domain interpolation factor):
~~~bash
zp_factor = 10
~~~

Cut 5% of samples on each side of the 2 spectra:
~~~bash
side_cut = True
~~~

Calculate the spectrum's frequency step:
~~~bash
df = freq_vect[1]-freq_vect[0]
~~~

Application of the PBWE to the 2 polarimetric channels' spectra:
~~~bash
output_pbwe, time_pbwe_vect = PyPBWE.PBWE(spec_mat,df,extra_factor,model_order,zp_factor,side_cut)
~~~

output_bwe matrix contains the time-domain radar soundings after PBWE, ready to be displayed with time axis time_pbwe_vect.

### Display of the PBWE results:

In this section we display the PBWE results compared to the original radar soundings for polarimetric channels 00 and 11.

Generate a time vector corresponding to the time-domain transform:
~~~bash
time_vect = np.linspace(0,1/df,zp_factor*np.shape(spec_mat)[1])
~~~

Display the original radar sounding and the PBWE version for channel 00:
~~~bash
plt.subplot(1,2,1)
plt.plot(time_vect*1e9,abs(1.85*np.fft.fft(np.conjugate(spec_mat[0,:])*np.hamming(len(spec_mat[0,:])),zp_factor*len(spec_mat[0,:])))/len(spec_mat[0,:]),'k-')
plt.plot(time_pbwe_vect*1e9,abs(output_pbwe[0,:]),'r-')
plt.xlim([5,9])
plt.xlabel('Time delays (ns)')
plt.ylabel('Normalized amplitude')
plt.legend(['Original radar sounding','Radar sounding after PBWE'], loc ='best')
plt.title('PBWE application - Polar channel 00')
plt.grid()
~~~

Display the original radar sounding and the PBWE version for channel 00:
~~~bash
plt.subplot(1,2,2)
plt.plot(time_vect*1e9,abs(1.85*np.fft.fft(np.conjugate(spec_mat[1,:])*np.hamming(len(spec_mat[1,:])),zp_factor*len(spec_mat[1,:])))/len(spec_mat[1,:]),'k-')
plt.plot(time_pbwe_vect*1e9,abs(output_pbwe[1,:]),'r-')
plt.xlim([5,9])
plt.xlabel('Time delays (ns)')
plt.ylabel('Normalized amplitude')
plt.legend(['Original radar sounding','Radar sounding after PBWE'], loc ='best')
plt.title('PBWE application - Polar channel 11')
plt.grid()
plt.show()
~~~

(The X 1.85 factor in amplitude corresponds to the compensation for the Hamming window).

Here is an example of figure displayed with this command:

[<img src="https://github.com/NicOudart/PyBWE/raw/main/docs/figures/Figure_example_PyBWE_PBWE.png"/>](Figure_example_PyBWE_PBWE.png)

### Manual extrapolation:

If you want to model and extrapolate the multi-channel spectrum manually, it can also be done using the PyPBWE.polar_burg and PyPBWE.polar_extrapolation functions.

Cutting 5% of frequencies on each side of the spectrum:
~~~bash
spec_mat = spec_mat[:,round(np.shape(spec_mat)[1]*0.05):np.shape(spec_mat)[1]-round(np.shape(spec_mat)[1]*0.05)]
freq_vect = freq_vect[round(len(freq_vect)*0.05):len(freq_vect)-round(len(freq_vect)*0.05)]
~~~

Retrieve the number of polarimetric channels:
~~~bash
Npol = np.shape(spec_mat)[0]
~~~

Retrieve the new number of samples in the spectrum:
~~~bash
M = np.shape(spec_mat)[1]
~~~

Set the multichannel AR model's order as 1/3 of the number of samples:
~~~bash
model_order_nb = round(0.33*M)
~~~

Fit a multichannel AR model to spec_mat:
~~~bash
[Thetaf,Thetab,err] = PyPBWE.polar_burg(spec_mat,model_order_nb)
~~~

Number of samples to be extrapolated on each side of the multichannel spectrum to get an extrapolation factor of 3:
~~~bash
Mextra = round(((3*M)-M)//2)+1
~~~

Extrapolation of spec_mat on each side of the multichannel spectrum:
~~~bash
spec_mat_extra,spec_mat_extra_forward,spec_mat_extra_backward = PyPBWE.polar_extrapolation(spec_mat,Thetaf,Thetab,Mextra,"both")
~~~

spec_mat_extra is the extrapolated multichannel spectrum on both sides, spec_mat_extra_forward and spec_mat_extra_backward are the forward and backward extrapolations of this multichannel spectrum.

### Display of the extrapolation result:

In this section, we display the polarimetric extrapolated spectra (channels 00 and 11), along with the expected spectra (noise-free) after extrapolation for comparison.

Generate the expected spectrum after extrapolation for comparison:
~~~bash
df = freq_vect[1]-freq_vect[0]
freq_vect_expected = np.linspace(freq_vect[0]-(df*(Mextra-2)),freq_vect[-1]+(df*Mextra),M+(Mextra*2))
spec_vect_expected_00 = np.exp(-1j*4*pi*dist_target1*freq_vect_expected/3e8)+np.exp(-1j*4*pi*dist_target2*freq_vect_expected/3e8)
spec_vect_expected_11 = np.exp(-1j*4*pi*dist_target1*freq_vect_expected/3e8)-np.exp(-1j*4*pi*dist_target2*freq_vect_expected/3e8)
~~~

Display the expected spectrum after extrapolation for the 00 and 11 channels:
~~~bash
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
~~~

Here is an example of figure displayed with this command:

[<img src="https://github.com/NicOudart/PyBWE/raw/main/docs/figures/Figure_example_PyPBWE_extrapolation_1.png"/>](Figure_example_PyPBWE_extrapolation_1.png)

Generate a vector of frequencies for the forward and backward extrapolations:
~~~bash
freq_vect_forward = np.linspace(freq_vect[-1]+df,freq_vect[-1]+(df*Mextra),Mextra)
freq_vect_backward = np.linspace(freq_vect[0]-(df*(Mextra-2)),freq_vect[0]-df,Mextra)
~~~

Display the original radar spectrum and the extrapolation for the 00 and 11 channels:
~~~bash
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
~~~

Here is an example of figure displayed with this command:

[<img src="https://github.com/NicOudart/PyBWE/raw/main/docs/figures/Figure_example_PyPBWE_extrapolation_2.png"/>](Figure_example_PyPBWE_extrapolation_2.png)

## References

* [Cuomo (1992)](https://apps.dtic.mil/sti/tr/pdf/ADA258462.pdf)

* [Suwa and Iwamoto (2003)](https://doi.org/10.1109/IGARSS.2003.1295505)

* [Suwa and Iwamoto (2007)](https://doi.org/10.1109/TGRS.2006.885406)

* [Ciarletti et al. (2017)](https://doi.org/10.1089/ast.2016.1532)

* [Raguso et al. (2018)](https://doi.org/10.1109/MetroAeroSpace.2018.8453529)

* [Hervé et al. (2020)](https://doi.org/10.1016/j.pss.2020.104939)

* [Oudart et al. (2021)](https://doi.org/10.1016/j.pss.2021.105173)

* [Gambacorta et al. (2022)](https://doi.org/10.1109/TGRS.2022.3216893)
