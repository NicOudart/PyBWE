# Unit tests

This library proposes unit tests for all unit functions of packages PyBWE, PyPBWE and PySSBWE, to ensure they run independently as expected.
These tests can be launched as GitHub actions (**test-script.yml**, jobs **unitTestBWE**, **unitTestPBWE** and **unitTestSSBWE**). The output is printed in the console.

The unit test scripts for the PyBWE, PyPBWE and PySSBWE packages in this library can be found at:

* test/script_unit_test_BWE.py

* test/script_unit_test_PBWE.py

* test/script_unit_test_SSBWE.py

## Test scenario

The unit tests for packages PyBWE, PyPBWE and PySSBWE are performed on synthetic radar data. All 3 are based on a similar scenario:

* 2 targets in sight of the radar, at distances 1 and 1.1 m, generate echoes in the radar signal of amplitudes 1 and 1.

* The radar output is a spectrum of 151 frequencies between 1.5 and 3 GHz.

* The spectrum is corrupted by a white-noise of a low standard deviation 1e-6 (to ensure a certain reproducibility of results).

In the case of PyPBWE, 2 signals are generated for 2 co-polar channels 00 and 11. The amplitudes of echoes are 1 and 1 for channel 00, 1 and -1 for channel 11.

In the case of PySSBWE, the echoes are affected by frequency-domain decays of 0.25 and 0.5 1/Hz. SSBWE methods 1 and 2 (using either the observability or the controlability matrix) are tested.

A reference for this spectrum is generated, for a larger frequency bandwidth and noise-free, for comparison.

## Test strategy

The unit test scripts for packages PyBWE, PyPBWE and PySSBWE follow different strategies depending on the function.

### Modelling functions

For the modelling functions PyBWE.burg, PyBWE.mcov, PyPBWE.polar_burg and PySSBWE.statespace_model, we follow this strategy:

* Check the error messages if the input order for the model is too high or too low.

* Compare the model's coefficients returned by the function to the expected values.

### Extrapolation functions

For the extrapolation functions PyBWE.ar_extrapolation, PyPBWE.polar_extrapolation and PySSBWE.statespace_extrapolation, we follow this strategy:

* Check the error messages if the model's coefficients arrays do not have dimensions consistent with the spectrum to extrapolate.

* Check the error messages if the number of extrapolated samples is not strictly positive.

* For all extrapolation modes (forward, backward, both), compare the extrapolation returned by the function to the reference (extrapolation factor 3).

### Echoes properties extraction functions

For the echoes properties extraction function PySSBWE.statespace_properties, we follow this strategy:

* Check the error messages if the model's coefficients arrays do not have dimensions consistent with the input spectrum.

* Compare the amplitudes of echoes extracted from the State-Space model to the expected values.

* Compare the targets distances estimated from the echoes time-delays, extracted from the State-Space model, to the expected values.

* Compare the frequency-domain decays of echoes extracted from the State-Space model to the expected values.

# Integration tests

This library proposes integration tests for the integrated methods PyBWE.BWE, PyPBWE.PBWE and PySSBWE.SSBWE, to ensure they run as expected.
These functions are indeed based on the sub-functions tested in the unit tests.
These tests can be launched as GitHub actions (**test-script.yml**, jobs **integrationTestBWE**, **integrationTestPBWE** and **integrationTestSSBWE**). The output is printed in the console.

## Test scenario

The unit tests for packages PyBWE, PyPBWE and PySSBWE are performed on synthetic radar data. All 3 are based on a similar scenario:

* 2 targets in sight of the radar, at distances 1 and 1.1 m, generate echoes in the radar signal of amplitudes 1 and 1.

* The radar output is a spectrum of 151 frequencies between 1.5 and 3 GHz.

* The spectrum is corrupted by a white-noise of a low standard deviation 1e-6 (to ensure a certain reproducibility of results).

In the case of PyPBWE, 2 signals are generated for 2 co-polar channels 00 and 11. The amplitudes of echoes are 1 and 1 for channel 00, 1 and -1 for channel 11.

SSBWE methods 1 and 2 (using either the observability or the controlability matrix) are tested.

## Test strategy

The integration test scripts for functions PyBWE.BWE, PyPBWE.PBWE and PySSBWE.SSBWE follow this strategy:

* Check if the function runs without failure with different parameters.

* Check the dimensions, the 1st and the last value of the time vector returned by the function.

* Check if the time vector and the radar sounding returned by the function have the same dimension.

* Check if 2 echoes are detected in the radar sounding returned by the function.

* Check if the amplitudes (module) of the 2 echoes in the radar sounding returned by the function are 1 and 1.

* Check if the distances of the 2 targets, estimated from the time-delays of echoes in the radar sounding returned by the function, are 1 and 1.1 m.
