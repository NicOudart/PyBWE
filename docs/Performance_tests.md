# Performance tests

This library proposes performance tests for the integrated methods PyBWE.BWE, PyPBWE.PBWE and PySSBWE.SSBWE, to ensure modifications will not cause regressions.

## Performance facing white-noise

This library proposes one white-noise performance test for each package's integrated method PyBWE.BWE, PyPBWE.PBWE and PySSBWE.SSBWE.
These tests are defined in the following scripts:

* test/script_test_performances_white_noise_BWE.py

* test/script_test_performances_white_noise_PBWE.py

* test/script_test_performances_white_noise_SSBWE.py

These tests follow a scenario similar to the example scripts of each package.

### The test scenario

Synthetic radar signals are generated (inspired by the WISDOM GPR of the ExoMars rover mission, Ciarletti et al. (2017)):

* A SFCW (Stepped Frequency Continuous Wave) radar working between 0.5 and 3 GHz measures a 1001 frequencies spectrum when sounding.

* Only the In-phase component (real part of the spectrum) is measured, the Quadrature component (imaginary part of the spectrum) is reconstructed by Hilbert transform.

* Two targets in free-space are seperated by a given distance, slightly below the radar's free-space resolution. These targets generate echoes of given amplitudes in the radar's signal, or complex sine-waves in the measured spectrum.

* The measured spectrum is corrupted by a white-noise, with a standard deviation leading to a given SNR.

Each script will apply the BWE, PBWE or SSBWE to synthetic radar signals generated for given distances between the targets, and given SNRs.
For each SNR, a given number of white-noise cases will be tested, with the same random seed to ensure the test's reproducibility.

### The test parameters

These test parameters can be modified by the user, but to follow the performances of the methods we recommend to keep them as defined in the script.

The list of distances between targets (in m) to be tested is given by:

~~~bash
list_dist_targets = [0.04,0.06,0.08,0.1,0.12]
~~~

The list of SNR values (in dB) to be tested is given by:

~~~bash
list_snr_levels = [6,14,20,40,60]
~~~

In a Monte-Carlo approach, the number of white-noise cases to be tested is given by:

~~~bash
nb_noise_case = 1000
~~~

1000 white-noise cases are proposed for the BWE and SSBWE. The multi-channel version of the Burg algorithm used in the PBWE being significantly longer to run, we propose 100 noise cases instead.

In addition to these test parameters, scenario parameters can also be modified if necessary.

### The scenario parameters

The amplitude of the 2 echoes corresponding to the 2 targets can be set:

~~~bash
amp_target1 = 1
amp_target2 = 1
~~~

In the case of the PBWE, the amplitudes can be set for each polarimetric channel 00 and 11:

~~~bash
amp_target1_00 = 1
amp_target2_00 = 1
amp_target1_11 = 1
amp_target2_11 = -1
~~~

The distance between the 1st target and the radar's antennas:

~~~bash
dist_target1 = 1
~~~

The peak detection threshold on amplitude, for an echo to be "detected":

~~~bash
detection_level = 0.5
~~~

The chosen value corresponds to a -6 dB criterion compared to the expected echoes amplitudes.

An additional parameter can be found for the SSBWE method, the order of the model:

~~~bash
param_order = 2
~~~

We set this parameter to 2, as we know 2 targets are in front of the antennas. This is of course an ideal case, and Akaike's Information Criterion (AIC) can be used instead by setting this parameter to zero.

For the BWE, PBWE and SSBWE functions, the default configuration will be tested.

### The test metrics

To evaluate the performances of the BWE, PBWE and SSBWE methods, we selected several metrics, corresponding to the objectives of radar range resolution enhancement:
Being able to detect the presence of the 2 echoes, to accurately measure the distance between them, and to accurately measure the amplitude of each echo.

One metric is returned for each couple distance between targets / SNR.

Here are the selected metrics:

* The percentage of times the 2 echoes are detected in the tested white-noise cases, according to the user defined threshold.

* The mean and std of the error on the distance between targets (m), estimated from the 2 echoes' time-delays measured in each white-noise cases.

* The mean and std of the relative error on the 1st echo's amplitude (%), measured in each white-noise case.

* The mean and std of the relative error on the 2nd echo's amplitude (%), measured in each white-noise case.

### The test report

For each method (BWE, PBWE and SSBWE), a Markdown test report will be exported.
Each section of this Markdown files will correspond to a specific test metric.

After running the tests, the reports can be found in:

* test/PyBWE_Report_test_performances_white_noise.md

* test/PyPBWE_Report_test_performances_white_noise.md

* test/PySSBWE_Report_test_performances_white_noise.md

In the case of the PBWE, as 2 co-polarimetric channels 00 and 11 are generated for the test, we return the metrics for each channel.

In the case of the SSBWE, as 2 methods can be used to determine the signal model (method 1 with the observability matrix, method 2 with the controllability matrix), we return the metrics for each method.

## References

The white-noise tests of the BWE proposed in this library are similar to the ones proposed in:

* [Oudart et al. (2021)](https://doi.org/10.1016/j.pss.2021.105173)

The white-noise tests of the BWE, PBWE and SSBWE proposed in this library are similar to the ones proposed in:

* [Oudart, thesis (2021)](https://theses.hal.science/tel-03364662)
