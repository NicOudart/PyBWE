################################################################################
#HEADER

#This script allows you to perform a test of the integrated functions in the
#PySSBWE package.

#The test is based on synthetic radar data, generated with the following scenario:
# - 2 targets in sight of the radar, at distances 1 and 1.1 m, generate
#   echoes in the radar signal of amplitudes 1 and 1.
# - The 2 echoes are affected by frequency-domain decays of 0.25 and 0.5 1/Hz.
# - The radar output is a spectrum of 151 frequencies between 1.5 and 3 GHz.
# - The spectrum is corrupted by a white-noise of standard deviation 1e-6.

################################################################################

#Libraries importation:---------------------------------------------------------

import numpy as np
from math import pi
from scipy.signal import find_peaks

import PySSBWE

#Test dataset generation:-------------------------------------------------------

#Generate a vector of 151 frequencies between 1.5 and 3 GHz:
freq_vect_test = np.linspace(1.5e9,3e9,151)

#Generate a sum of two complex sine-waves corresponding to 2 targets' echoes at
#the resolution limit, corrupted by a very low amplitude white-noise:
spec_vect_test = np.exp(1j*4*pi*freq_vect_test/3e8)+np.exp(1j*4*pi*1.1*freq_vect_test/3e8)+np.random.normal(0,1e-6,151)

#Test function SSBWE:-----------------------------------------------------------

#Check PySSBWE.SSBWE on the test dataset with different parameters:
try:
    output_ssbwe_1, output_ssbwe_2, time_ssbwe_vect = PySSBWE.SSBWE(spec_vect_test,freq_vect_test[1]-freq_vect_test[0],3,10,True,2) #Standard parameters
    output_ssbwe_1, output_ssbwe_2, time_ssbwe_vect = PySSBWE.SSBWE(spec_vect_test,freq_vect_test[1]-freq_vect_test[0],2,10,True,2)
    output_ssbwe_1, output_ssbwe_2, time_ssbwe_vect = PySSBWE.SSBWE(spec_vect_test,freq_vect_test[1]-freq_vect_test[0],3,5,True,2)
    output_ssbwe_1, output_ssbwe_2, time_ssbwe_vect = PySSBWE.SSBWE(spec_vect_test,freq_vect_test[1]-freq_vect_test[0],3,10,True,0) #AIC to estimate the order of the model
    output_ssbwe_1, output_ssbwe_2, time_ssbwe_vect = PySSBWE.SSBWE(spec_vect_test,freq_vect_test[1]-freq_vect_test[0],3,10,True,2)
    output_ssbwe_1, output_ssbwe_2, time_ssbwe_vect = PySSBWE.SSBWE(spec_vect_test,freq_vect_test[1]-freq_vect_test[0],3,10,False,2) #To be used for the rest of the test
except:
    assert False, 'PySSBWE.SSBWE NOK!'

#Application of a peak detection to the radar sounding (SSBWE method 1 and 2):
peaks_ssbwe_1 = find_peaks(abs(output_ssbwe_1),height=0.1)
peaks_ssbwe_2 = find_peaks(abs(output_ssbwe_2),height=0.1)
dist_peaks_1 = abs(time_ssbwe_vect[peaks_ssbwe_1[0]])*0.5*3e8
dist_peaks_2 = abs(time_ssbwe_vect[peaks_ssbwe_2[0]])*0.5*3e8

#Check the output time vector:
assert time_ssbwe_vect[-1]==1e-7, 'PySSBWE.SSBWE NOK!' #The last element
assert time_ssbwe_vect[0]==0, 'PySSBWE.SSBWE NOK!' #The 1st element
assert len(time_ssbwe_vect)==4550, 'PySSBWE.SSBWE NOK!' #The total length

#Check the output radar sounding:
assert (len(output_ssbwe_1)==len(time_ssbwe_vect))and(len(output_ssbwe_2)==len(time_ssbwe_vect)), 'PySSBWE.SSBWE NOK!' #The radar sounding has the right length
assert (len(peaks_ssbwe_1)==2)and(len(peaks_ssbwe_2)==2), 'PySSBWE.SSBWE NOK!' #Only 2 main echoes are detected in the radar sounding
assert (abs(abs(peaks_ssbwe_1[1]['peak_heights'][0])-1)<1e-2)and(abs(abs(peaks_ssbwe_1[1]['peak_heights'][1])-1)<1e-2)and(abs(abs(peaks_ssbwe_2[1]['peak_heights'][0])-1)<1e-2)and(abs(abs(peaks_ssbwe_2[1]['peak_heights'][1])-1)<1e-2), 'PySSBWE.SSBWE NOK!' #The echoes amplitudes are 1
assert (abs(abs(dist_peaks_1[0])-1)<1e-2)and(abs(abs(dist_peaks_1[1])-1.1)<1e-2)and(abs(abs(dist_peaks_2[0])-1)<1e-2)and(abs(abs(dist_peaks_2[1])-1.1)<1e-2), 'PySSBWE.SSBWE NOK!' #The distances of targets from the radar are 1 and 1.1 m

#Output:------------------------------------------------------------------------

print('PySSBWE integration tests OK!')