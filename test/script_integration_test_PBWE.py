################################################################################
#HEADER

#This script allows you to perform a test of the integrated functions in the
#PyPBWE package.

#The test is based on synthetic radar data, generated with the following scenario:
# - The radar acquires signals in 2 co-polarimetric channels 00 and 11.
# - 2 targets in sight of the radar, at distances 1 and 1.1 m, generate
#   echoes in the radar signal of amplitudes 1 and 1 in channel 00, 1 and -1 in
#   channel 11.
# - The radar output is 2 co-polarimetric spectra of 151 frequencies between
#   1.5 and 3 GHz.
# - The spectra are corrupted by a white-noise of standard deviation 1e-6.

################################################################################

#Libraries importation:---------------------------------------------------------

import numpy as np
from math import pi
from scipy.signal import find_peaks

import PyPBWE

#Test dataset generation:-------------------------------------------------------

#Generate a vector of 151 frequencies between 1.5 and 3 GHz:
freq_vect_test = np.linspace(1.5e9,3e9,151)

#Generate a sum of two complex sine-waves corresponding to 2 targets' echoes at
#the resolution limit, corrupted by a very low amplitude white-noise, for 2
#polarimetric channels (00 and 11):
spec_vect_00_test = np.exp(-1j*4*pi*freq_vect_test/3e8)+np.exp(-1j*4*pi*1.1*freq_vect_test/3e8)+np.random.normal(0,1e-6,151)
spec_vect_11_test = np.exp(-1j*4*pi*freq_vect_test/3e8)-np.exp(-1j*4*pi*1.1*freq_vect_test/3e8)+np.random.normal(0,1e-6,151)

#Concatenate the 2 polarimetric channels in a matrix:
spec_mat_test = np.vstack((spec_vect_00_test,spec_vect_11_test))

#Test function PBWE:------------------------------------------------------------

#Check PyPBWE.PBWE on the test dataset with different parameters:
try:
    output_pbwe, time_pbwe_vect = PyPBWE.PBWE(spec_mat_test,freq_vect_test[1]-freq_vect_test[0],3,0.33,10,True) #Standard parameters
    output_pbwe, time_pbwe_vect = PyPBWE.PBWE(spec_mat_test,freq_vect_test[1]-freq_vect_test[0],2,0.33,10,True)
    output_pbwe, time_pbwe_vect = PyPBWE.PBWE(spec_mat_test,freq_vect_test[1]-freq_vect_test[0],3,0.1,10,True)
    output_pbwe, time_pbwe_vect = PyPBWE.PBWE(spec_mat_test,freq_vect_test[1]-freq_vect_test[0],3,0.33,5,True)
    output_pbwe, time_pbwe_vect = PyPBWE.PBWE(spec_mat_test,freq_vect_test[1]-freq_vect_test[0],3,0.33,10,False) #To be used for the rest of the test
except:
    assert False, 'PyPBWE.PBWE NOK!'

#Application of a peak detection to the radar sounding:
peaks_pbwe_00 = find_peaks(abs(output_pbwe[0,:]),height=0.1)
peaks_pbwe_11 = find_peaks(abs(output_pbwe[1,:]),height=0.1)
dist_peaks_00 = abs(time_pbwe_vect[peaks_pbwe_00[0]])*0.5*3e8
dist_peaks_11 = abs(time_pbwe_vect[peaks_pbwe_11[0]])*0.5*3e8

#Check the output time vector:
assert time_pbwe_vect[-1]==1e-7, 'PyPBWE.PBWE NOK!' #The last element
assert time_pbwe_vect[0]==0, 'PyPBWE.PBWE NOK!' #The 1st element
assert len(time_pbwe_vect)==4550, 'PyPBWE.PBWE NOK!' #The total length

#Check the output radar sounding:
assert (len(output_pbwe)==2)and(len(output_pbwe[0])==len(time_pbwe_vect))and(len(output_pbwe[1])==len(time_pbwe_vect)), 'PyPBWE.PBWE NOK!' #The radar sounding has the right length
assert (len(peaks_pbwe_00)==2)and(len(peaks_pbwe_11)==2), 'PyPBWE.PBWE NOK!' #Only 2 main echoes are detected in the radar sounding
assert (abs(abs(peaks_pbwe_00[1]['peak_heights'][0])-1)<1e-2)and(abs(abs(peaks_pbwe_00[1]['peak_heights'][1])-1)<1e-2)and(abs(abs(peaks_pbwe_11[1]['peak_heights'][0])-1)<1e-2)and(abs(abs(peaks_pbwe_11[1]['peak_heights'][1])-1)<1e-2), 'PyPBWE.PBWE NOK!' #The echoes amplitudes are 1
assert (abs(abs(dist_peaks_00[0])-1)<1e-2)and(abs(abs(dist_peaks_00[1])-1.1)<1e-2)and(abs(abs(dist_peaks_11[0])-1)<1e-2)and(abs(abs(dist_peaks_11[1])-1.1)<1e-2), 'PyPBWE.PBWE NOK!' #The distances of targets from the radar are 1 and 1.1 m

#Output:------------------------------------------------------------------------

print('PyPBWE integration tests OK!')