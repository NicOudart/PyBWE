[![DOI](https://zenodo.org/badge/724105763.svg)](https://zenodo.org/doi/10.5281/zenodo.10259072)

[![status](https://joss.theoj.org/papers/16898c7ed1ff84843966a673e06dd513/status.svg)](https://joss.theoj.org/papers/16898c7ed1ff84843966a673e06dd513)

[<img src="https://github.com/NicOudart/PyBWE/raw/main/Logo_PyBWE.png" width="500"/>](Logo_PyBWE.png)

# PyBWE: Python tools for Bandwidth Extrapolation of planetary radar signals

Resolution enhancement of radar signals with the super-resolution technique named "Bandwidth Extrapolation" (BWE).
Website: [https://nicoudart.github.io/PyBWE](https://nicoudart.github.io/PyBWE)

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
(**SSBWE**) using a State-Space model accounting for noise and exponential distortions in the signal (Piou, 1999).

Also, the BWE can be used to fill a gap between two spectra of multiband radars (Moore, 1997). This process is known as
**Bandwidth Interpolation** (**BWI**).

### Recent application to planetary radar sounders

If the BWE can be applied to any radar signals, it has been extensively applied to planetary radar sounders. 
Radar sounders unlocked a 3rd dimension in planetary studies, allowing scientists to discover the widely unknown 
subsurface of the Moon, Mars, Titan and soon Venus and Jupiter's moons. The development of space instruments being 
highly constrained, the BWE is a mean to get as much information on the subsurface structure and composition from a 
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

## Library installation and update

### Installation from PyPI

**Installation:**
~~~bash
pip install PyBWE
~~~

**Update:**
~~~bash
pip install PyBWE --upgrade
~~~

### Packages:

This library contains 3 packages:

* **PyBWE:** _Plain vanilla package_. It contains the BWE functions, as well as the sub-functions allowing you to model and extrapolate a radar spectrum.

* **PyPBWE:** _Advanced package_. It contains the PBWE function, as well as the sub-functions allowing you to model and extrapolate a polarimetric radar spectra.

* **PySSBWE:** _Advanced package_. It contains the SSBWE function, as well as the sub-functions allowing you to model and extrapolate a radar spectrum, as well as extracting echoes information from it.

PyBWE should be used by default. PyPBWE can be used only if the radar has more than one polarimetric channel. SSBWE 
can be used in high noise-level or in high attenuations situations, but the stability of the model is not guaranteed.
See the documentation for more information on the functions contained in each package.

Each package can be imported the following way:

~~~bash
import PyBWE
import PyPBWE
import PySSBWE
~~~

## Documentation

In **docs** can be found the full documentation of this library, with one file for each package, containing descriptions of each function in the package. 

In **examples** can be found one example scripts for each package, with applications on synthetic radar signals corrupted by a white-noise.

## Tests

In **test** can be found test scripts for each package, divided into unit, integration and performance tests.
Each test script can be launched on request with GitHub action **test-scripts.yml**. 

### Unit and integration tests

Each package (PyBWE, PyPBWE and PySSBWE) contains "unit" sub-functions that can be used independently, and "integrated" functions based on these sub-functions.  
We propose test scripts for the "unit functions", and then for the "integrated functions", to check if they behave as expected given a specific synthetic radar input.
Each script will print a message, with the test result ("OK"/"NOK").

### Performance tests

The performances of each package (PyBWE, PyPBWE and PySSBWE) are also tested, using a set of metrics (errors on the estimated distance between targets, errors on echoes amplitudes, etc.), on synthetic radar signals of different complexity.
For now, the performances of each package are evaluated for different distances between 2 radar targets (below the resolution limit with FFT), and for different SNR levels (considering a signal corrupted by a white-noise). Several noise cases are tested for a given SNR.
The performance test results are exported as Markdown reports.
Other performance test scenarios can be implemented in the future.

## Community guidelines

### Contribute

We are open to contributions to improve PyBWE.
Please fork the repository to make your updates.

For each new update, please follow these guidelines:

* For any new functions, add unit, integration and if possible performance tests to the corresponding test scripts in **test**. If possible, new example scripts should also be added to **examples**.

* For any change to PyBWE, update the documentation accordingly in **docs**.

* Keep the code style (the naming of functions/scripts, the headers, sections...) and comment your updates.

### Report issues

Feel free to report any issues through GitHub _Issues_.
To help us fix the issues you want to report, please:

* Include in the title and/or in the text, names of the packages and functions for which the issue is observed.

* For reproducibility, add to the text your OS, your Python version, and your PyBWE version.

* Add to the text a description of the expected and the observed behaviours.

* Include any output, error message or screenshot relevant to understand the issue.
 
### Support

If you need support using PyBWE, you can either:

* Look at the documentation: [PyBWE website](https://nicoudart.github.io/PyBWE)

* Contact us directly: nicolas.oudart@latmos.ipsl.fr

## Credits

© Nicolas OUDART (2024)

LATMOS/IPSL, UVSQ Université Paris-Saclay, Guyancourt, France

## Citation

OUDART Nicolas (2024), PyBWE: Python tools for Bandwidth Extrapolation of radar signals (2024.02.2).
DOI: 10.5281/zenodo.10612861

## License

This package is released under a MIT open source license. See **LICENSE**.

## References

Parametric spectral estimation:

* [Kay and Marple (1981)](https://doi.org/10.1109/PROC.1981.12184)

BWE:

* [Bowling (1977)](https://apps.dtic.mil/sti/pdfs/ADA042817.pdf)

* [Cuomo (1992)](https://apps.dtic.mil/sti/tr/pdf/ADA258462.pdf)

* [Moore (1997)](https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=5e09b4f4c00f6f660495a38c279961031a376e59)

PBWE:

* [Suwa and Iwamoto (2003)](https://doi.org/10.1109/IGARSS.2003.1295505)

* [Suwa and Iwamoto (2007)](https://doi.org/10.1109/TGRS.2006.885406)

SSBWE:

* [Piou (1999)](https://apps.dtic.mil/sti/pdfs/ADA366105.pdf)

Planetary science applications:

* [Mastrogiuseppe et al. (2014)](https://doi.org/10.1002/2013GL058618)

* [Raguso et al. (2018)](https://doi.org/10.1109/MetroAeroSpace.2018.8453529)

* [Oudart et al. (2021)](https://doi.org/10.1016/j.pss.2021.105173)

* [Gambacorta et al. (2022)](https://doi.org/10.1109/TGRS.2022.3216893)

* [Eide et al. (2022)](https://doi.org/10.1029/2022GL101429)
