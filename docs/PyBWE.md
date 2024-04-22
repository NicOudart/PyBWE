# PyBWE package

## Description

The PyBWE package is the "plain vanilla" package of this library.
It allows you to apply the classic Bandwidth Extrapolation (BWE) technique as defined by Cuomo (1992) to any radar spectrum.
The BWE is a super-resolution technique, as it yields a better range resolution in time-domain than classic spectral estimation techniques (FFT).

The BWE works this way:

* The input is the frequency-domain spectrum of the radar's reponse signal from a series of targets. In this situation, the signal should in theory be a sum of complex sine-waves.

* An autoregressive (AR) model is fitted to this spectrum, using the Burg algorithm. The order of the model is a user-defined parameter.

* This model is used to extrapolate the spectrum forward and backward. The extrapolating factor is a user-defined parameter.

* The super-resolved time-domain sounding is obtained by IFFT. Zero-padding can be used for time-domain interpolation (this process is purely aesthetic).

The spectrum's extrapolation factor is equal to the resolution enhancement factor after IFFT.

The PyBWE package contains the **PyBWE.BWE** function, allowing you to apply the BWE directly to a given radar spectrum and get a super-resolved sounding.
This function calls several other functions from the package, that you can call independently if needed:

* **PyBWE.burg:** fits an AR model to a spectrum, using the Burg algorithm.

* **PyBWE.ar_extrapolation:** extrapolates forward, backward or both a spectrum, given an AR model.

## Importation

Once the library is installed, the PyBWE package can simply be imported:

~~~bash
import PyBWE
~~~

## Recommendations

### The order of the AR model

In theory, the ideal order of an AR model depends on the number of complex sine-waves composing the spectrum.
However, in practical cases, a higher order is required to obtain good results in presence of noise.
For the BWE, Cuomo (1992) recommends an order equal to 1/3 of the number of samples in the spectrum.

### The extrapolation factor

In theory, a sum of sine-wave can be extrapolated infinitely.
However, in practical cases, errors accumulate with the extrapolation.
For the BWE, Cuomo (1992) recommends an extrapolation factor of maximum 3.

### The effect of side-lobes

In the case of a radar signal spectrum obtained by FFT, 2 effects must be mitigated before BWE:

* If a window was applied to the signal, it must be inverted.

* The application of the FFT generates side-lobes effects on the borders of the spectrum. For this reason, Cuomo (1992) recommends to remove 5% of frequencies on each side of a spectrum before BWE.

In the case of a radar measuring only the real-part of a frequency spectrum, the imaginary part can be reconstructed by Hilbert transform.
The FFT being part of Hilbert transform implementations, errors can also be expected on the border of the spectrum.
We thus also recommend to cut 5% of samples on each side of a spectrum in this case.

The PyBWE.BWE function allows you to remove 5% of samples on each side of a spectrum before BWE.

### The effect of distortions

The BWE technique is based on the assumption that the spectrum of a radar response signal from targets is a sum of complex sine-waves.
In many practical applications, a radar response signal can differ from this model, due to various distortion effects.
Such distortions can make the application of the direct application of the BWE impossible.
We thus recommend to correct distortions prior to a BWE application.

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
For instance Gambacorta et al. (2022) used a ionospheric Gamma model before applying the BWE to MARSIS (Mars Express) soundings.

* **Effect of subsurface attenuations:**
When a radar signal is transmitted through a sub-surface, it is affected by different frequency-dependent attenuations (absorption, scattering...).
It can thus be expected that deeper echoes in a high-loss subsurface will not be correctly reconstructed by BWE.
We recommend to use the PySSBWE package, as the SSBWE (State-Space Bandwidth Extrapolation) technique accounts for such attenuations.

### The effect of noise

The presence of noise in a signal impact the quality of the AR model determined by the Burg algorithm.
For this reason, a high level of noise will impact the BWE results (echoes amplitudes and time-delays), as demonstrated in Oudart et al. (2021).
To mitigate the effect of noise, we recommend to use a more complex model, like the ones provided in PyPBWE and PySSBWE packages.

### Electromagnetic interferences

EM interferences can cause peaks in the measured radar spectrum. Such peaks will impair the AR modelling of the signal.
The effect can be mitigated by detecting outliers in the spectrum, and then replacing them by interpolation, as proposed by Raguso et al. (2018).

### Cross-talk

If the emission and reception antennas / channels are not the same, a cross-talk phenomenon can generate unwanted signals in the radar response.
It is recommended to remove such effects prior to the application of the BWE.
This can for instance be done by subtracting a "free-space" (without any targets) calibration signal to any radar response.

## Functions

### PyBWE.ar_extrapolation(x,ar_coeff,Nextra,extra_mode='both')

The PyBWE.ar_extrapolation function extrapolates a spectrum given its AR model, forward, backward or both.

**Inputs:**

* x: _complex 1D array_ containing the spectrum to be extrapolated.

* ar_coeff: _complex 1D array_ containing the coefficients of the spectrum's AR model.

* Nextra: _integer_ number of samples to be extrapolated forward, backward or both.

* extra_mode: (optional) _string_ containing the extrapolation mode, "forward"/"backward"/"both".

**Outputs:**

* x_extra: _complex 1D array_ containing the extrapolated spectrum.

* x_forward: _complex 1D array_ containing the forward extrapolation, empty if in "backward" mode.

* x_backward: _complex 1D array_ containing the backward extrapolation, empty if in "forward" mode.

### PyBWE.burg(x,p)

The PyBWE.burg function fits an AR model to a spectrum, using the Burg algorithm.

**Inputs:**

* x: _complex 1D array_ containing the spectrum to be modelled.

* p: _integer_ order of the AR model.

**Outputs:**

* a: _complex 1D array_ containing the AR model coefficients.

* e: _float_ estimated white-noise variance.

* rc: _complex 1D array_ containing the Burg algorithm's reflection coefficients.

### PyBWE.BWE(spec,df,extra_factor,model_order,zp_factor,side_cut=True)

The PyBWE.BWE function applies the Bandwidth Extrapolation to a radar signal's spectrum.

**Inputs:**

* spec: _complex 1D array_ containing the spectrum to which the BWE will be applied.

* df: _float_ frequency step corresponding to spec [Hz].

* extra_factor: _float_ factor between the extrapolated and the original spectrum's bandwidths.

* model_order: _float_ order of the AR model for extrapolation, expressed as a ratio of the original spectrum's bandwidth.

* zp_factor: _float_ factor between the zero-padded and original spectrum's bandwidth.

* side_cut: (optional) _boolean_, 5% of samples are cut on each side of the spectrum if True.

**Outputs:**

* output_bwe: _complex 1D array_ containing the time-domain radar signal after BWE.

* time_bwe_vect: _float 1D array_ containing the time axis corresponding to output_bwe [s].

**Notes:**
In this implementation, the time-domain output signal is scaled using the method described by Cuomo (1992).
A Hamming window is applied to the spectrum before IFFT.

## Example

2 example scripts are proposed with the same synthetic radar scenario for the PyBWE package:

* An application of PyBWE.BWE on a synthetic radar signal is proposed in _examples/script_example_PyBWE_BWE.py_.

* A manual use of PyBWE.burg and PyBWE.ar_extrapolation on a synthetic radar signal is proposed in _examples/script_example_PyBWE_extrapolation.py_.

In the following, these 2 scripts are explained as a single tutorial:

### Presentation of the scenario

The synthetic radar signal example (inspired by the WISDOM GPR of the ExoMars rover mission, Ciarletti et al. (2017)):

* A SFCW (Stepped Frequency Continuous Wave) radar working between 0.5 and 3 GHz measures a 1001 frequencies spectrum when sounding.

* Only the In-phase component (real part of the spectrum) is measured, the Quadrature component (imaginary part of the spectrum) is reconstructed by Hilbert transform.

* Two targets in free-space are seperated by 5 cm, slightly below the radar's free-space resolution. These targets generate echoes of given complex amplitudes in the radar's signal, or complex sine-waves in the measured spectrum.

* The measured spectrum is corrupted by a white-noise of standard deviation 10X smaller than the complex sine-waves' amplitudes.

The Bandwidth Extrapolation (BWE) is applied to this radar's signal using the PyBWE function_BWE:

* Most of the estimation errors when reconstructing a complex spectrum with the Hilbert transform are on the far sides of the spectrum. For this reason, we cut 5% of frequencies on each side of the spectrum before BWE. We would do the same for a radar working in time-domain, as most errors in FFT estimation are also on each side of the spectrum. This process is useless for a radar working in the frequency-domain and measuring both In-Phase and Quadrature components of the spectrum.

* We then fit an AR model to the spectrum, with an order equal to 1/3 of the spectrum samples, as recommended by Cuomo (1992).

* We use this model to extrapolate the spectrum on each side, to obtain a bandwidth 3X larger (maximum extrapolation factor recommended by Cuomo (1992)). A bandwidth X3 yields a resolution X3 better in time-domain.

* The extrapolated spectrum is eventually converted to a time-domain signal by IFFT, with zero-padding to interpolate the signal X10. This interpolation is purely aesthetic.

This scenario requires the following libraries:

~~~bash
import matplotlib.pyplot as plt
import numpy as np
from math import pi
from scipy.signal import hilbert

import PyBWE
~~~

### The synthetic signal generation

In this section we describe how to generate analytically a synthetic radar sounding to test the PyBWE package.

Generate a vector of 1001 frequencies between 0.5 and 3 GHz:
~~~bash
freq_vect = np.linspace(0.5e9,3e9,1001)
~~~

Distance (m) between two radar targets in free-space (returning two echoes) and the radar:
~~~bash
dist_target1 = 1
dist_target2 = 1.07
~~~

Amplitude of the 2 echoes corresponding to the 2 targets:
~~~bash
ampli_target1 = 1
ampli_target2 = 1
~~~

Generate a sum of two complex sine-waves corresponding to the targets' echoes (each with an amplitude of 1):
~~~bash
spec_vect = (ampli_target1*np.exp(1j*4*pi*dist_target1*freq_vect/3e8))+(ampli_target2*np.exp(1j*4*pi*dist_target2*freq_vect/3e8))
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
spec_vect = hilbert(spec_vect)[::2]
freq_vect = freq_vect[::2]
~~~

This is the signal to which we will apply the BWE.

### The application of the BWE:

The BWE can be applied directly to the spectrum using the PyBWE.BWE function.

Extrapolation factor (=< 3 recommended by Cuomo (1992)):
~~~bash
extra_factor = 3
~~~

Order of the model as a ratio of the total number of samples in the spectrum (= 0.33 recommended by Cuomo (1992)):
~~~bash
model_order = 0.33
~~~

Zero-padding factor (= time domain interpolation factor):
~~~bash
zp_factor = 10
~~~

Cut 5% of samples on each side of the spectrum:
~~~bash
side_cut = True
~~~

Calculate the original spectrum's frequency step:
~~~bash
df = freq_vect[1]-freq_vect[0]
~~~

Application of the BWE to the spectrum:
~~~bash
output_bwe, time_bwe_vect = PyBWE.BWE(spec_vect,df,extra_factor,model_order,zp_factor,side_cut)
~~~

output_bwe contains the time-domain radar sounding after BWE, ready to be displayed with time axis time_bwe_vect.

### Display of the BWE results:

In this section we display the BWE results compared to the original radar sounding.

Generate a time vector corresponding to the time-domain transform:
~~~bash
time_vect = np.linspace(0,1/df,zp_factor*len(spec_vect))
~~~

Display the original radar sounding and the BWE version:
~~~bash
plt.plot(time_vect*1e9,abs(1.85*np.fft.fft(spec_vect*np.hamming(len(spec_vect)),zp_factor*len(spec_vect)))/len(spec_vect),'k-')
plt.plot(time_bwe_vect*1e9,abs(output_bwe),'r-')
plt.xlim([5,9])
plt.xlabel('Time delays (ns)')
plt.ylabel('Normalized amplitude')
plt.legend(['Original radar sounding','Radar sounding after BWE'], loc ='best')
plt.title('BWE application')
plt.grid()
plt.tight_layout()
plt.show()
~~~

(The X 1.85 factor in amplitude corresponds to the compensation for the Hamming window).

Here is an example of figure displayed with this command:

[<img src="https://github.com/NicOudart/PyBWE/raw/main/docs/figures/Figure_example_PyBWE_BWE.png"/>](Figure_example_PyBWE_BWE.png)

### Manual extrapolation:

If you want to model and extrapolate the spectrum manually, it can also be done using the PyBWE.burg and PyBWE.ar_extrapolation functions.

Cutting 5% of frequencies on each side of the spectrum:
~~~bash
spec_vect = spec_vect[round(len(spec_vect)*0.05):len(spec_vect)-round(len(spec_vect)*0.05)]
freq_vect = freq_vect[round(len(freq_vect)*0.05):len(freq_vect)-round(len(freq_vect)*0.05)]
~~~

Retrieve the number of samples in the spectrum:
~~~bash
N = len(spec_vect)
~~~

Set the AR model's order as 1/3 of the number of samples:
~~~bash
model_order_nb = round(0.33*N)
~~~

Fit an AR model to the spectrum using the Burg algorithm:
~~~bash
ar_coeff,ar_var,ar_rc = PyBWE.burg(spec_vect,model_order_nb)
~~~

Number of samples to be extrapolated on each side of the spectrum to get an extrapolation factor of 3:
~~~bash
Nextra = round(((3*N)-N)//2)+1
~~~

Extrapolation on each side of the spectrum:
~~~bash
spec_vect_extra,spec_vect_forward,spec_vect_backward = PyBWE.ar_extrapolation(spec_vect,ar_coeff,Nextra,"both")
~~~

spec_vect_extra is the extrapolated spectrum on both sides, spec_vect_forward and spec_vect_backward are the forward and backward extrapolations of the spectrum.

### Display of the extrapolation result:

In this section, we display the AR-extrapolated spectrum, along with the expected spectrum (noise-free) after extrapolation for comparison.

Generate the expected spectrum after extrapolation for comparison:
~~~bash
df = freq_vect[1]-freq_vect[0]
freq_vect_expected = np.linspace(freq_vect[0]-(df*(Nextra-2)),freq_vect[-1]+(df*Nextra),N+(Nextra*2))
spec_vect_expected = np.exp(1j*4*pi*dist_target1*freq_vect_expected/3e8)+np.exp(1j*4*pi*dist_target2*freq_vect_expected/3e8)
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
plt.tight_layout()
plt.show()
~~~

Here is an example of figure displayed with this command:

[<img src="https://github.com/NicOudart/PyBWE/raw/main/docs/figures/Figure_example_PyBWE_extrapolation_1.png"/>](Figure_example_PyBWE_extrapolation_1.png)

Generate a vector of frequencies for the forward and backward extrapolations:
~~~bash
freq_vect_forward = np.linspace(freq_vect[-1]+df,freq_vect[-1]+(df*Nextra),Nextra)
freq_vect_backward = np.linspace(freq_vect[0]-(df*(Nextra-2)),freq_vect[0]-df,Nextra)
~~~

Display the original radar spectrum and the extrapolation:
~~~bash
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
plt.tight_layout()
plt.show()
~~~

Here is an example of figure displayed with this command:

[<img src="https://github.com/NicOudart/PyBWE/raw/main/docs/figures/Figure_example_PyBWE_extrapolation_2.png"/>](Figure_example_PyBWE_extrapolation_2.png)

## References

* [Cuomo (1992)](https://apps.dtic.mil/sti/tr/pdf/ADA258462.pdf)

* [Ciarletti et al. (2017)](https://doi.org/10.1089/ast.2016.1532)

* [Raguso et al. (2018)](https://doi.org/10.1109/MetroAeroSpace.2018.8453529)

* [Hervé et al. (2020)](https://doi.org/10.1016/j.pss.2020.104939)

* [Oudart et al. (2021)](https://doi.org/10.1016/j.pss.2021.105173)

* [Gambacorta et al. (2022)](https://doi.org/10.1109/TGRS.2022.3216893)
