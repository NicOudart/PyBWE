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