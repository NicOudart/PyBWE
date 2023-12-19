# PySSBWE package

## Description

The PySSBWE package is an "advanced" package of this library.
It allows you to apply the State-Space Bandwidth Extrapolation (SSBWE) technique as defined by Piou (1999) to any to any radar spectrum.
The SSBWE is a super-resolution technique, as it yields a better range resolution in time-domain than classic spectral estimation techniques (FFT).
It is a state-space modelling version of the classic Bandwidth Extrapolation (BWE), this improved model making the technique more robust to noise and sub-surface attenuations.

The SSBWE works this way:

* The input is the frequency-domain spectrum of the radar's reponse signal from a series of targets, with sub-surface attenuations. In this situation, the signal should in theory be a sum of complex sine-waves, with distortions following decreasing exponentials.

* A state-space model is fitted to this spectrum, using a state-space identification technique. The order of the model is estimated by Akaike's Information Criterion (AIC).

* One model is used to extrapolate the spectrum forward, and another to extrapolate the spectrum backward. The extrapolating factor is a user-defined parameter.

* The super-resolved time-domain sounding is obtained by IFFT. Zero-padding can be used for time-domain interpolation (this process is purely aesthetic).

The spectrum's extrapolation factor is equal to the resolution enhancement factor after IFFT.

The PySSBWE package contains the **PySSBWE.SSBWE** function, allowing you to apply the SSBWE directly to a given radar spectrum and get a super-resolved sounding.
This function calls several other functions from the package, that you can call independently if needed:

* **PySSBWE.statespace_model:** fits a state-space model to a spectrum, using 2 state-space identification methods.

* **PySSBWE.AIC:** estimates the ideal order for a spectrum's state-space model using Akaike's Information Criterion.

* **PySSBWE.statespace_extrapolation:** extrapolates forward, backward or both a spectrum, given an AR model.

In addition, this package also contains a function that can be used to retrieve information on the different targets' echoes from a spectrum's state-space model:

* **PySSBWE.statespace_properties:** retrieve the amplitudes, time-delays, and frequency-decays estimations of each echoes from a spectrum's state-space model.

## Importation

Once the library is installed, the PySSBWE package can simply be imported:

~~~bash
import PySSBWE
~~~

## Recommendations

### The order of the AR model

The idea behind the state-space identification of the SSBWE technique is to seperate the singular values of signal from the singular values of noise in the input spectrum's Hankel matrix.
The ideal order of a state-space model is thus directly the number of complex sine-waves composing the spectrum, also equal to the number of targets sounded by the radar.
If you have an a priori knowledge of the number of targets, you directly set the order of the model.
If not, you can estimate it with Akaike's Information Criterion (AIC, see Akaike (1974)), as recommended by Piou (1999).

The stability of the model not being guaranteed in the case of the SSBWE, a wrong estimation of the ideal model order can lead to diverging extrapolations, and thus to very high errors in time-domain.

### The extrapolation factor

In theory, a sum of sine-wave can be extrapolated infinitely.
However, in practical cases, errors accumulate with the extrapolation.
For the BWE, Cuomo (1992) recommends an extrapolation factor of maximum 3. This is also the extrapolation factor used by Piou for the SSBWE (1999).

### The effect of side-lobes

In the case of a radar signal spectrum obtained by FFT, 2 effects must be mitigated before SSBWE:

* If a window was applied to the signal, it must be inverted.

* The application of the FFT generates side-lobes effects on the borders of the spectrum. For this reason, Cuomo (1992) recommends to remove 5% of frequencies on each side of a spectrum before SSBWE. The same goes for SSBWE.

In the case of a radar measuring only the real-part of a frequency spectrum, the imaginary part can be reconstructed by Hilbert transform.
The FFT being part of Hilbert transform implementations, errors can also be expected on the border of the spectrum.
We thus also recommend to cut 5% of samples on each side of a spectrum in this case.

The PySSBWE.SSBWE function allows you to remove 5% of samples on each side of a spectrum before SSBWE.

### The effect of distortions

The SSBWE technique is based on the assumption that the spectrum of a radar response signal from targets is a sum of complex sine-waves, with only decreasing exponential distortions.
In many practical applications, a radar response signal can differ from this model, due to various distortion effects.
Such distortions can make the application of the direct application of the BWE impossible.
We thus recommend to correct distortions prior to a SSBWE application.

Different solutions are possible for different sources of distortion:

* **Effect of the antennas gain:**
No antenna can transmit all frequencies with the same gain. This will cause distortions of sine-waves in a radar spectrum.
You can try to compensate for this effect with a calibration radar measurement on a "perfect reflector" (a metallic sphere or a large metallic plate).
From the spectrum of the echo from such target, the antennas gain can be determined, and inverted in all other radar measurements.
This process is sometimes called "whitening" and is described in Oudart et al. (2021).
As different temperatures can lead to significantly different amplitudes over a radar's bandwidth, for a successful whitening process we recommend to apply first a temperature correction (as the one described in Hervé et al. (2020)).

* **Effect of the ionosphere:**
The ionosphere of any planet has a frequency-dependent attenuation effect on electromagnetic waves. This will cause distortions of sine-waves in a radar spectrum. 
This effect can be compensated using a ionospheric model.
For instance Gambacorta et al. (2022) used a ionospheric Gamma model before applying the SSBWE to MARSIS (Mars Express) soundings.

* **Effect of subsurface attenuations:**
When a radar signal is transmitted through a sub-surface, it is affected by different frequency-dependent attenuations (absorption, scattering...).
It can thus be expected that deeper echoes in a high-loss subsurface will not be correctly reconstructed by the classic BWE.
The SSBWE accounts for this distortion effect, and the estimation of the level of attenuations at each target (the sum of absorption, scattering...) can even be extracted from the model.

### The effect of noise

The presence of noise in a signal impact the quality of the AR model determined by the Burg algorithm.
For this reason, a high level of noise will impact the classic BWE results (echoes amplitudes and time-delays), as demonstrated in Oudart et al. (2021).
The state-space model used by the SSBWE has a more robust noise model, and if the order of the model is chosen well, noise will impact much less the SSBWE than the BWE.

### Electromagnetic interferences

EM interferences can cause peaks in the measured radar spectrum. Such peaks will impair the AR modelling of the signal.
The effect can be mitigated by detecting outliers in the spectrum, and then replacing them by interpolation, as proposed by Raguso et al. (2018).

### Cross-talk

If the emission and reception antennas / channels are not the same, a cross-talk phenomenon can generate unwanted signals in the radar response.
It is recommended to remove such effects prior to the application of the SSBWE.
This can for instance be done by subtracting a "free-space" (without any targets) calibration signal to any radar response.

## Functions

### PySSBWE.AIC(sv,N,noise_type="white")

The PySSBWE.AIC function estimates the ideal state-space model order from the spectrum's Hankel matrix singular values, using AIC.
See Akaike (1974).

**Inputs:**

* sv: _float 1D array_ of the spectrum's Hankel matrix singular values, in decreasing order.

* N: _integer_ number of samples in the spectrum.

* noise_type: (optional) _string_, "white" if the spectrum is affected by a white-noise, "colored" if the spectrum is affected by a colored-noise.

**Outputs:**

* output_order: _integer_ ideal state-space model order according to AIC.

**Notes:**
The type of noise ("white" or "colored") is an optional input. The noise will be considered "white" by default.

### PySSBWE.SSBWE(spec,df,extra_factor,zp_factor,side_cut,order=0,noise_type="white")

The PySSBWE.SSBWE function applies the State-Space Bandwidth Extrapolation (SSBWE) to a radar signal's spectrum.
The State-Space model is estimated using two methods:

* Method 1: estimation of the state-space model from the observability matrix.

* Method 2: estimation of the state-space model from the controllability matrix.

Both models' results are returned by the function.

**Inputs:**

* spec: _complex 1D array_ containing the spectrum to which the SSBWE will be applied.

* df: _float_ frequency step corresponding to spec [Hz].

* extra_factor: _float_ factor between the extrapolated and the original spectrum's bandwidths.

* zp_factor: _float_ factor between the zero-padded and original spectrum's bandwidth.

* side_cut: _boolean_, 5% of samples are cut on each side of the spectrum if True.

* order: (optional) _integer_ order of the state-space model, if =0 the order is estimated using AIC.

* noise_type: (optional) _string_, "white" if the spectrum is affected by a white-noise, "colored" if the spectrum is affected by a colored-noise.

**Outputs:**

* output_ssbwe_1: _complex 1D array_ extrapolated spectrum using SSBWE method 1.

* output_ssbwe_2: _complex 1D array_ extrapolated spectrum using SSBWE method 2.

* time_ssbwe_vect: _complex 1D array_ time axis corresponding to output_ssbwe_1 and output_ssbwe_2 [s].

**Notes:**
Unlike the classic BWE (or the PBWE), where the Burg algorithm guaranties the stability of the model, the stability of the state-space model is not guaranteed.
The order of the model is by default estimated using AIC. It can also be set by the user.

### PySSBWE.statespace_extrapolation(y,A,B,C,Nextra)

The PySSBWE.statespace_extrapolation function extrapolates forward a spectrum given its state-space model.

**Inputs:**

* y: _complex 1D array_ containing the spectrum to be extrapolated.

* A: _complex 2D array_ containing the state-space model's "state matrix".

* B: _complex 1D array_ containing the state-space model's "input matrix".

* C: _complex 1D array_ containing the state-space model's "output matrix".

* Nextra: _integer_ number of samples to extrapolate from y.

**Outputs:**

* y_extra: _complex 1D array_ containing the forward extrapolation of y.

**Notes:**
This function only extrapolates a spectrum forward. 
To extrapolate a spectrum backward, you must flip the spectrum, fit a state-space model to it, extrapolate it forward, and eventually flip it again.

### PySSBWE.statespace_model(y,order=0,noise_type="white")

The PySSBWE.statespace_model function fits a state-space model to a spectrum.
The State-Space model is estimated using two methods:

* Method 1: estimation of the state-space model from the observability matrix.

* Method 2: estimation of the state-space model from the controllability matrix.

Both models are returned by the function.

**Inputs:**

* y: _complex 1D array_ containing the spectrum to be modelled.

* order: (optional) _integer_ order of the state-space model, if =0 the order is estimated using AIC.

* noise_type: (optional) _string_, "white" if the spectrum is affected by a white-noise, "colored" if the spectrum is affected by a colored-noise.

**Outputs:**

* A1: _complex 2D array_ containing the state-space model's "state matrix" obtained with method 1.

* B1: _complex 2D array_ containing the state-space model's "input matrix" obtained with method 1.

* C1: _complex 2D array_ containing the state-space model's "output matrix" obtained with method 1.

* A2: _complex 2D array_ containing the state-space model's "state matrix" obtained with method 2.

* B2: _complex 2D array_ containing the state-space model's "input matrix" obtained with method 2.

* C2: _complex 2D array_ containing the state-space model's "output matrix" obtained with method 2.

**Notes:**
Unlike the Burg algorithm which guaranties the stability of the model, the stability of the state-space model is not guaranteed.
The order of the model is by default estimated using AIC. It can also be set by the user.

### PySSBWE.statespace_properties(A,B,C,df,f1)

The PySSBWE.statespace_properties function retrieves the different echoes properties (complex amplitudes, time-delays, frequency-decays) from a spectrum's state-space model.

**Inputs:**

* A: _complex 2D array_ containing the spectrum state-space model's "state matrix".

* B: _complex 1D array_ containing the spectrum state-space model's "input matrix".

* C: _complex 1D array_ containing the spectrum state-space model's "output matrix".

* df: _float_ frequency step of the modelled spectrum [Hz].

* f1: _float_ 1st frequency of the modelled spectrum [Hz].

**Outputs:**

* amp: _complex 1D array_ containing the estimated echoes complex amplitudes in the signal.

* td: _float 1D array_ containing the estimated echoes time-delays [s].

* dec: _float_1D_array_ containing the estimated echoes frequency-decay [1/Hz].

## Example

2 example scripts are proposed with the same synthetic radar scenario for the PySSBWE package:

* An application of PySSBWE.SSBWE and PySSBWE.statespace_properties on a synthetic radar signal is proposed in _examples/script_example_PySSBWE_SSBWE.py_.

* A manual use of PySSBWE.statespace_model and PySSBWE.statespace_extrapolation on a synthetic radar signal is proposed in _examples/script_example_PySSBWE_extrapolation.py_.

In the following, these 2 scripts are explained as a single tutorial:

### Presentation of the scenario

The synthetic radar signal example (inspired by the WISDOM GPR of the ExoMars rover mission, Ciarletti et al. (2017)):

* A SFCW (Stepped Frequency Continuous Wave) radar working between 0.5 and 3 GHz measures a 1001 frequencies spectrum when sounding.

* Only the In-phase component (real part of the spectrum) is measured, the Quadrature component (imaginary part of the spectrum) is reconstructed by Hilbert transform.

* Two targets in free-space are seperated by 5 cm, slightly below the radar's free-space resolution. With these targets' echoes are associated different amplitudes and different decays in frequency-domain.

* The measured spectrum is corrupted by a white-noise of standard deviation 10X smaller than the complex sine-waves' amplitudes.

The State-Space Bandwidth Extrapolation (SSBWE) is applied to this radar's signal using the PyBWE function_SSBWE:

* Most of the estimation errors when reconstructing a complex spectrum with the Hilbert transform are on the far sides of the spectrum. For this reason, we cut 5% of frequencies on each side of the spectrum before SSBWE. We would do the same for a radar working in time-domain, as most errors in FFT estimation are also on each side of the spectrum. This process is useless for a radar working in the frequency-domain and measuring both In-Phase and Quadrature components of the spectrum.

* We then fit a state-space model to the spectrum, first in the forward and then the backward directions.

* We use these 2 models to extrapolate the spectrum on each side, to obtain a bandwidth 3X larger (maximum extrapolation factor recommended by Cuomo (1992)). A bandwidth X3 yields a resolution X3 better in time-domain.

* The extrapolated spectrum is eventually converted to a time-domain signal by IFFT, with zero-padding to interpolate the signal X10. This interpolation is purely aesthetic.

From the state-space model of the spectrum, the properties of each echo (amplitude, time-delay, frequency-domain decay) can be estimated. This is done in the last part of this example.

In the following example the order of the model (= number of echoes) is estimated using AIC. The order can also be forced by the user in the PySSBWE functions statespace_model and SSBWE with an optional parameter "order" if the number of echoes is known.

2 methods can be used to fit a state-space model to a spectrum:

*Method 1: using the observability matrix.

*Method 2: using the controllability matrix.

This scenario requires the following libraries:

~~~bash
import matplotlib.pyplot as plt
import numpy as np
from math import pi
from scipy.signal import hilbert

import PySSBWE
~~~

### The synthetic signal generation

In this section we describe how to generate analytically a synthetic radar sounding to test the PySSBWE package.

Generate a vector of 1001 frequencies between 0.5 and 3 GHz:
~~~bash
freq_vect = np.linspace(0.5e9,3e9,1001)
~~~

Distance (m) between two radar targets in free-space (returning two echoes) and the radar:
~~~bash
dist_target1 = 1
dist_target2 = 1.07
~~~

Amplitude associated with each target echo:
~~~bash
ampli_target1 = 1
ampli_target2 = 0.75
~~~

Frequency decay associated with each target echo (1/Hz):
~~~bash
decay_target1 = 0.25e-9
decay_target2 = 0.5e-9
~~~

Generate a sum of two complex sine-waves corresponding to the targets' echoes (each with an amplitude of 1):
~~~bash
spec_vect = (ampli_target1*np.exp(-(decay_target1+(1j*4*pi*dist_target1/3e8))*freq_vect))+(ampli_target2*np.exp(-(decay_target2+(1j*4*pi*dist_target2/3e8))*freq_vect))
~~~

Only keep the real part of the spectrum (In-phase component):
~~~bash
spec_vect = np.real(spec_vect)
~~~

Add a white noise (standard deviation 0.1) to this spectrum signal:
~~~bash
wn_vect = np.random.normal(0,0.1,spec_vect.shape)
spec_vect += wn_vect
~~~

Reconstruct a complex signal with the Hilbert transform:
~~~bash
spec_vect = np.conjugate(hilbert(spec_vect))[::2]
freq_vect = freq_vect[::2]
~~~

This is the signal to which we will apply the SSBWE.

### The application of the SSBWE:

The SSBWE can be applied directly to the spectrum using the PySSBWE.SSBWE function.

Extrapolation factor (=< 3 recommended by Cuomo (1992)):
~~~bash
extra_factor = 3
~~~

Zero-padding factor (= time domain interpolation factor):
~~~bash
zp_factor = 10
~~~

Cut 5% of samples on each side of the spectrum:
~~~bash
side_cut = True
~~~

Calculate the spectrum's frequency step:
~~~bash
df = freq_vect[1]-freq_vect[0]
~~~

Application of the BWE to the spectrum:
~~~bash
output_ssbwe_1, output_ssbwe_2, time_ssbwe_vect = PySSBWE.SSBWE(spec_vect,freq_vect,extra_factor,zp_factor,side_cut)
~~~

output_ssbwe_1 and output_ssbwe_2 contains the time-domain radar soundings after SSBWE (method 1 and 2), ready to be displayed with time axis time_ssbwe_vect.

### Display of the SSBWE results:

In this section we display the SSBWE results compared to the original radar sounding, for methods 1 and 2.

Generate a time vector corresponding to the time-domain transform:
~~~bash
time_vect = np.linspace(0,1/df,zp_factor*len(spec_vect))
~~~

Display the original radar sounding and the SSBWE version using method 1:
~~~bash
plt.plot(time_vect*1e9,abs(1.85*np.fft.fft(np.conjugate(spec_vect)*np.hamming(len(spec_vect)),zp_factor*len(spec_vect)))/len(spec_vect),'k-')
plt.plot(time_ssbwe_vect*1e9,abs(output_ssbwe_1),'r-')
plt.xlim([5,9])
plt.xlabel('Time delays (ns)')
plt.ylabel('Normalized amplitude')
plt.legend(['Original radar sounding','Radar sounding after SSBWE'], loc ='best')
plt.title('SSBWE application - method 1 (using the observability matrix)')
plt.grid()
plt.show()
~~~

Display the original radar sounding and the SSBWE version using method 2:
~~~bash
plt.plot(time_vect*1e9,abs(1.85*np.fft.fft(np.conjugate(spec_vect)*np.hamming(len(spec_vect)),zp_factor*len(spec_vect)))/len(spec_vect),'k-')
plt.plot(time_ssbwe_vect*1e9,abs(output_ssbwe_2),'r-')
plt.xlim([5,9])
plt.xlabel('Time delays (ns)')
plt.ylabel('Normalized amplitude')
plt.legend(['Original radar sounding','Radar sounding after SSBWE'], loc ='best')
plt.title('SSBWE application - method 2 (using the controllability matrix)')
plt.grid()
plt.show()
~~~

(The X 1.85 factor in amplitude corresponds to the compensation for the Hamming window).

### Extraction of the echoes properties:

In this section we show how to extract the echoes properties (amplitudes, time-delays, frequency-domain decays) from a spectrum's state-space model.

Manually cut 5% of frequencies of each side of the spectrum:
~~~bash
spec_vect_cut = spec_vect[round(len(spec_vect)*0.05):len(spec_vect)-round(len(spec_vect)*0.05)]
freq_vect_cut = freq_vect[round(len(freq_vect)*0.05):len(freq_vect)-round(len(freq_vect)*0.05)]
~~~

Fit a state-space model to the spectrum:
~~~bash
[A1,B1,C1,A2,B2,C2] = PySSBWE.statespace_model(spec_vect_cut)
~~~

Retrieve the spectrum's frequency step and 1st frequency:
~~~bash
df = freq_vect_cut[1]-freq_vect_cut[0]
f1 = freq_vect_cut[0]
~~~

Extract the echoes properties from the state-space model obtained with method 1:
~~~bash
[amp_1,td_1,dec_1] = PySSBWE.statespace_properties(A1,B1,C1,df,f1)
~~~

Extract the echoes properties from the state-space model obtained with method 2:
~~~bash
[amp_2,td_2,dec_2] = PySSBWE.statespace_properties(A2,B2,C2,df,f1)
~~~

Print the echoes amplitudes (module), range (m) and frequency decays (1/GHz) obtained using method 1:
~~~bash
print("---Echoes properties estimation (method 1: using observability matrix)---")
print("Module amplitude of echoes: "+str(abs(amp_1)))
print("Distance of targets from the radar (m): "+str(0.5*td_1*3e8))
print("Frequency-domain decay of echoes (1/GHz): "+str(dec_1*1e9))
print("-------------------------------------------------------------------------")
~~~

Print the echoes amplitudes (module), range (m) and frequency decays (1/GHz) obtained using method 2:
~~~bash
print("---Echoes properties estimation (method 2: using controllability matrix)---")
print("Module amplitude of echoes: "+str(abs(amp_2)))
print("Distance of targets from the radar (m): "+str(0.5*td_2*3e8))
print("Frequency-domain decay of echoes (1/GHz): "+str(dec_2*1e9))
print("---------------------------------------------------------------------------")
~~~

### Manual extrapolation:

If you want to model and extrapolate the spectrum manually, it can also be done using the PySSBWE.statespace_model and PySSBWE.statespace_extrapolation functions.

Cutting 5% of frequencies on each side of the spectrum:
~~~bash
spec_vect = spec_vect[round(len(spec_vect)*0.05):len(spec_vect)-round(len(spec_vect)*0.05)]
freq_vect = freq_vect[round(len(freq_vect)*0.05):len(freq_vect)-round(len(freq_vect)*0.05)]
~~~

Flip for a backward version of the spectrum:
~~~bash
spec_vect_b = np.flip(spec_vect)
~~~

Retrieve the numbers of samples in the spectrum:
~~~bash
N = len(spec_vect)
~~~

Number of samples to be extrapolated on each side of the spectrum to get an extrapolation factor of 3:
~~~bash
Nextra = round(((3*N)-N)//2)+1
~~~

Fit a state-space model to the spectrum for forward extrapolation:
~~~bash
[A1_f,B1_f,C1_f,A2_f,B2_f,C2_f] = PySSBWE.statespace_model(spec_vect,noise_type="white")
~~~

Fit a state-space model to the flipped spectrum for backward extrapolation:
~~~bash
[A1_b,B1_b,C1_b,A2_b,B2_b,C2_b] = PySSBWE.statespace_model(spec_vect_b,noise_type="white")
~~~

Forward extrapolation of the spectrum (obtained with method 1 and 2):
~~~bash
spec_vect_extra_1_f = PySSBWE.statespace_extrapolation(spec_vect,A1_f,B1_f,C1_f,Nextra)
spec_vect_extra_2_f = PySSBWE.statespace_extrapolation(spec_vect,A2_f,B2_f,C2_f,Nextra)
~~~

Backward extrapolation of the spectrum (obtained with method 1 and 2):
~~~bash
spec_vect_extra_1_b = PySSBWE.statespace_extrapolation(spec_vect_b,A1_b,B1_b,C1_b,Nextra)
spec_vect_extra_2_b = PySSBWE.statespace_extrapolation(spec_vect_b,A2_b,B2_b,C2_b,Nextra)
~~~

Flip the backward extrapolation (obtained with method 1 and 2):
~~~bash
spec_vect_extra_1_b = np.flip(spec_vect_extra_1_b)
spec_vect_extra_2_b = np.flip(spec_vect_extra_2_b)
~~~

spec_vect_extra_1_f and spec_vect_extra_1_b are the forward and backward spectrum extrapolations using method 1.
spec_vect_extra_2_f and spec_vect_extra_2_b are the forward and backward spectrum extrapolations using method 2.

### Display of the extrapolation result:

In this section, we display the state-space extrapolated spectra (methods 1 and 2), along with the expected spectra (noise-free) after extrapolation for comparison.

Generate the expected spectrum after extrapolation for comparison:
~~~bash
df = freq_vect[1]-freq_vect[0]
freq_vect_expected = np.linspace(freq_vect[0]-(df*(Nextra-2)),freq_vect[-1]+(df*Nextra),N+(Nextra*2))
spec_vect_expected = (ampli_target1*np.exp(-(decay_target1+(1j*4*pi*dist_target1/3e8))*freq_vect_expected))+(ampli_target2*np.exp(-(decay_target2+(1j*4*pi*dist_target2/3e8))*freq_vect_expected))
~~~

Display the expected spectrum after extrapolation:
~~~bash
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
~~~

Generate a vector of frequencies for the forward and backward extrapolations:
~~~bash
freq_vect_forward = np.linspace(freq_vect[-1]+df,freq_vect[-1]+(df*Nextra),Nextra)
freq_vect_backward = np.linspace(freq_vect[0]-(df*(Nextra-2)),freq_vect[0]-df,Nextra)
~~~

Display the original radar spectrum and the extrapolation using method 1:
~~~bash
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
plt.show()
~~~

Display the original radar spectrum and the extrapolation using method 2:
~~~bash
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
plt.show()
~~~

## References

* [Akaike (1974)](https://doi.org/10.1109/TAC.1974.1100705)

* [Cuomo (1992)](https://apps.dtic.mil/sti/tr/pdf/ADA258462.pdf)

* [Piou (1999):](https://apps.dtic.mil/sti/pdfs/ADA366105.pdf)

* [Ciarletti et al. (2017)](https://doi.org/10.1089/ast.2016.1532)

* [Raguso et al. (2018)](https://doi.org/10.1109/MetroAeroSpace.2018.8453529)

* [Hervé et al. (2020)](https://doi.org/10.1016/j.pss.2020.104939)

* [Oudart et al. (2021)](https://doi.org/10.1016/j.pss.2021.105173)

* [Gambacorta et al. (2022)](https://doi.org/10.1109/TGRS.2022.3216893)
