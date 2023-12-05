[![DOI](https://zenodo.org/badge/724105763.svg)](https://zenodo.org/doi/10.5281/zenodo.10259072)

[<img src="https://github.com/NicOudart/PyBWE/raw/main/Logo_PyBWE.png" width="500"/>](Logo_PyBWE.png)

# PyBWE: Python tools for Bandwidth Extrapolation of radar signals

Resolution enhancement of radar signals with the super-resolution technique named "Bandwidth Extrapolation" (BWE).

## Introduction to PyBWE

### The Bandwidth Extrapolation technique

Range resolution enhancement is one of the main challenges in radar signal processing. It is limited by the frequency 
bandwidth of the instrument: the larger is the bandwidth, the better is the resolution. 

Classic Fourier transform techniques are efficient and robust for spectral estimation, but in practical applications 
the output resolution is limited. For this reason parametric spectral estimation techniques have been introduced (see 
Kay and Marple, 1981). Based on a signal model, assuming deterministic properties of the signal, the output of such 
techniques yields a better resolution. However, parametric techniques are less robust than classic Fourier transform 
ones, in particular in the presence of noise or distortions.

The **Bandwidth Extrapolation technique** (**BWE**) is a compromise between a classic Fourier transform and a 
parametric spectral estimation technique (Cuomo, 1992). A parametric model is fitted to the signal's spectrum, this 
model is then used to extrapolate this spectrum forward and backward, and the spectrum is eventually Fourier 
transformed using IFFT.
The extrapolation factor is equal to the resolution enhancement, and can be up to 3 in practical cases.

In the regular BWE, the signal is modelled by an autoregressive (AR) model, using the Burg algorithm (Burg, 1967). 
Several improvements to the BWE have been proposed: the **Polarimetric BWE** (**PBWE**) using the correlation between several 
polarimetric radar channels for an improved extrapolation (Suwa and Iwamoto, 2003,2007), or the **State-Space BWE** 
(**SSBWE**) using a State-Space model to account for noise and exponential distortions in the signal (Piou, 1999).

Also, the BWE can be used to fill gap between two spectra of multiband radars (Moore, 1997). This process is known as
**Bandwidth Interpolation** (**BWI**).

### Recent application to planetary radar sounders

If the BWE can be applied to any radar signals, it has been extensively applied to planetary radar sounders. 
Radar sounders unlocked a 3rd dimension in planetary studies, allowing scientists to discover the widely unknown 
subsurface of the Moon, Mars, Titan and soon Venus and Jupiter's moons. The development of space instruments being 
very constrained, the BWE is a mean to get as much information on the subsurface structure and composition from a 
given radar instrument.
 
Here is a non-exhaustive list of successful BWE applications in planetary science:

* The 1st bathymetry of a Titan sea using Cassini radar data (Mastrogiuseppe et al., 2014)

* The improvement of the stratigraphic analysis of Martian polar ice sheets using SHARAD (MRO) radar sounder data (Raguso et al., 2018)

* The improvement of the WISDOM (ExoMars) Ground Penetrating Radar soundings in preparation of the Rosalind Franklin rover mission (Oudart et al., 2021)

* The improvement of the MARSIS (Mars Express) radar sounder resolution by a factor of 6 using both BWE and BWI (Gambacorta et al., 2022)

* The BWE helped the estimation of attenuations in the Martian subsurface with the RIMFAX (Mars 2020) Ground Penetrating Radar data (Eide et al., 2022)

With the arrival of WISDOM (ExoMars) on Mars in 2028, as well as new radar sounders selected for the exploration of 
Venus (EnVision mission) and Jupiter's icy moons (Juice and Europa Clipper missions) in the 2030s, we expect the BWE 
techniques to be useful for future planetary science studies. Hence the publication of this Python library, intented 
to help the planetary science community use the BWE on any radar sounder data.
