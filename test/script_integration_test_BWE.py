################################################################################
#HEADER

#This script allows you to perform a test of the integrated functions in the
#PyBWE package.

#The test is based on synthetic radar data, generated with the following scenario:
# - 2 targets in sight of the radar, at distances 1 and 1.1 m, generate
#   echoes in the radar signal of amplitudes 1 and 1.
# - The radar output is a spectrum of 151 frequencies between 1.5 and 3 GHz.
# - The spectrum is corrupted by a white-noise of standard deviation 1e-6.

################################################################################

#Libraries importation:---------------------------------------------------------

import numpy as np
from math import pi
from scipy.signal import find_peaks

import PyBWE

#Test dataset generation:-------------------------------------------------------

#Generate a vector of 151 frequencies between 1.5 and 3 GHz:
freq_vect_test = np.linspace(1.5e9,3e9,151)

#Generate a sum of two complex sine-waves corresponding to 2 targets' echoes at
#the resolution limit, corrupted by a very low amplitude white-noise:
spec_vect_test = np.exp(-1j*4*pi*freq_vect_test/3e8)+np.exp(-1j*4*pi*1.1*freq_vect_test/3e8)+np.random.normal(0,1e-6,151)

#Test function BWE:-------------------------------------------------------------

#Check PyBWE.BWE on the test dataset with different parameters:
try:
    output_bwe, time_bwe_vect = PyBWE.BWE(spec_vect_test,freq_vect_test[1]-freq_vect_test[0],3,0.33,10,True) #Standard parameters
    output_bwe, time_bwe_vect = PyBWE.BWE(spec_vect_test,freq_vect_test[1]-freq_vect_test[0],2,0.33,10,True)
    output_bwe, time_bwe_vect = PyBWE.BWE(spec_vect_test,freq_vect_test[1]-freq_vect_test[0],3,0.1,10,True)
    output_bwe, time_bwe_vect = PyBWE.BWE(spec_vect_test,freq_vect_test[1]-freq_vect_test[0],3,0.33,5,True)
    output_bwe, time_bwe_vect = PyBWE.BWE(spec_vect_test,freq_vect_test[1]-freq_vect_test[0],3,0.33,10,False) #To be used for the rest of the test
except:
    assert False, 'PyBWE.BWE NOK!'

#Application of a peak detection to the radar sounding:
peaks_bwe = find_peaks(abs(output_bwe),height=0.1)
dist_peaks = abs(time_bwe_vect[peaks_bwe[0]])*0.5*3e8

#Check the output time vector:
assert time_bwe_vect[-1]==1e-7, 'PyBWE.BWE NOK!' #The last element
assert time_bwe_vect[0]==0, 'PyBWE.BWE NOK!' #The 1st element
assert len(time_bwe_vect)==4550, 'PyBWE.BWE NOK!' #The total length

#Check the output radar sounding:
assert len(output_bwe)==len(time_bwe_vect), 'PyBWE.BWE NOK!' #The time vector and radar sounding have the same length
assert len(peaks_bwe)==2, 'PyBWE.BWE NOK!' #Only 2 main echoes are detected in the radar sounding
assert (abs(abs(peaks_bwe[1]['peak_heights'][0])-1)<1e-2)and(abs(abs(peaks_bwe[1]['peak_heights'][1])-1)<1e-2), 'PyBWE.BWE NOK!' #The echoes amplitudes are 1
assert (abs(abs(dist_peaks[0])-1)<1e-2)and(abs(abs(dist_peaks[1])-1.1)<1e-2), 'PyBWE.BWE NOK!' #The distances of targets from the radar are 1 and 1.1 m

#Output:------------------------------------------------------------------------

print('PyBWE integration tests OK!')