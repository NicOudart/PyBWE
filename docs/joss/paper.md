---
title: "PyBWE: Python tools for Bandwidth Extrapolation of planetary radar signals"
tags:
 - radar
 - planetary science
 - bandwidth extrapolation
 - super-resolution
 - parametric spectral estimation
authors:
 - name: Nicolas Oudart
   orcid: 0000-0002-6691-808X
   affiliation: 1
affiliations:
 - name: LATMOS/IPSL, UVSQ Université Paris-Saclay, Sorbonne Université, CNRS, France
   index: 1
date: 5 January 2024
bibliography: paper.bib
---

# Summary

**PyBWE** is a Python library containing radar super-resolution methods known as "Bandwidth Extrapolation" (BWE).

Range resolution enhancement is one of the main challenges in radar signal processing. 
It is driven by the time resolution of radar soundings, and the speed of electromagnetic waves in the sounded material. 
Time resolution being limited by the frequency bandwidth of the instrument, the same applies to range resolution: the larger the bandwidth, the better the resolution. 

Fast Fourier Transform techniques are efficient and robust for spectral estimation, but their performances in time resolution are limited by the inverse of the bandwidth, and the necessary application of windowing to reduce the impact of side-lobes. 
For this reason parametric spectral estimation techniques have been introduced [@Kay:1981].
Based on a signal model, assuming deterministic properties of the signal, the output of such techniques yields a better resolution. 
However, parametric techniques are less robust than classic Fourier transform ones, in particular in the presence of noise or distortions.

The **Bandwidth Extrapolation technique** (**BWE**) is a compromise between a classic Fourier transform and a parametric spectral estimation technique [@Cuomo:1992]. 
A parametric model is fitted to the measured signal's spectrum, this model is then used to extrapolate this spectrum forward and backward, and the spectrum is eventually Fourier transformed using an IFFT.
The extrapolation factor is equal to the resolution enhancement, and can be up to 3 in practical cases. 

The extrapolation of a radar spectrum is indeed possible, as it can be modelled by a sum of complex sine-waves with different amplitudes, frequencies and phases. Each of these complex sine-waves corresponds to a target echo in time-domain. 
In practical cases, this deterministic signal is corrupted by various sources of noise, distorsions, and can be dampened by losses in the sounded material. Hence the limitation of a radar spectrum extrapolation.

An example of application on a synthetic planetary radar spectrum, inspired by the WISDOM (ExoMars, ESA) planetary radar [@Ciarletti:2017], is shown in \autoref{fig:example}.

In the regular BWE, the signal is modelled by an autoregressive (AR) model, using the Burg algorithm. 
Several improvements to the BWE have been proposed: the **Polarimetric BWE** (**PBWE**) using the correlation between several polarimetric radar channels for an improved extrapolation [@Suwa:2003] [@Suwa:2007], or the **State-Space BWE** (**SSBWE**) using a State-Space model accounting for noise and exponential distortions in the signal [@Piou:1999].

Also, the BWE can be used to fill a gap between two spectra of multiband radars [@Moore:1997]. This process is known as **Bandwidth Interpolation** (**BWI**).

This library contains 3 packages, each containing a different BWE technique, based on a different signal model:

* PyBWE: implementing the "classic" BWE technique, with an AR model determined by the Burg algorithm.

* PyPBWE: implementing the PBWE technique, with a multi-channel AR model determined by a multi-channel Burg algorithm.

* PySSBWE: implementing the SSBWE technique, with an ARMA model determined by a State-Space modelling approach.

Each package contains an integrated solution to directly apply the complete BWE process to a given radar spectrum, as well as all the individual functions for modelling and extrapolation.

This library relies on Numpy [@Harris:2020] for the linear algebra, on Scipy [@Virtanen:2020] for matrix operations and data processing functions used in test and examples, on Matplolib [@Hunter:2007] for figure display in examples, and on Pandas [@Pandas:2020] for performance tests data handling.

# Statement of need

If the BWE can be applied to any radar signal, it has been extensively applied to planetary radar sounders in the last decade. 

Radar sounders unlocked a 3rd dimension in planetary studies, unveiling the subsurface of the Moon, Mars, Titan and soon Venus and Jupiter's moons.
The design of planetary exploration instruments being highly constrained, the resolution performances of radar sounders are not only driven by the scientific objectives of a mission.
In this context, the BWE is a powerful tool to get further insights on the subsurface structure and composition from a given radar sounder.

Here is a non-exhaustive list of successful BWE applications in planetary science:

* The 1st bathymetry of a Titan sea using Cassini (NASA) radar data [@Mastrogiuseppe:2014]

* The improvement of the stratigraphic analysis of Martian polar ice sheets using SHARAD (MRO, NASA) radar sounder data [@Raguso:2018]

* The improvement of the WISDOM (ExoMars, ESA) Ground Penetrating Radar soundings in preparation of the Rosalind Franklin rover mission [@Oudart:2021]

* The improvement of the MARSIS (Mars Express, ESA) radar sounder resolution by a factor of 6 using both BWE and BWI [@Gambacorta:2022]

* The BWE helped the estimation of attenuations in the Martian subsurface with the RIMFAX (Mars 2020, NASA) Ground Penetrating Radar data [@Eide:2022]

However, the planetary science community has few radar experts, and to our knowledge no open-source integrated BWE solutions existed before the release of this library, limiting the planetary radar application of this technique.
For this reason, we propose in this library integrated BWE solutions for all planetary scientists, as well as the individual functions for planetary radar experts.

In addition, this library could also benefit other radar / sonar applications developed under harsh contraints, impacting the reachable range resolution.

![An example of BWE application to a synthetic planetary radar spectrum, with 2 targets in free-space separated by a distance of 7 cm, below the instrument's resolution limit with IFFT and Hamming windowing. This spectrum is generated with 2 complex sine-waves of amplitude 1, corrupted by a white-noise of standard deviation 0.1 (SNR = 20 dB), for frequencies between 0.5 and 3 GHz. (a) Extrapolation of this spectrum with a parametric model. (b) Corresponding time-domain signal obtained with classic FFT or BWE. \label{fig:example}](Figure_example.png)

# Acknowledgments

The author would like to thank Valérie Ciarletti, Alice Le Gall and Emile Brighi for their helpful comments and suggestions. 

# References