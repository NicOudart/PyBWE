################################################################################
#HEADER

#This script allows you to perform a test of unit functions in the PyBWE package.

################################################################################

#Libraries importation:---------------------------------------------------------

import numpy as np
from math import pi,sqrt

import PyBWE

#Test dataset generation:-------------------------------------------------------

#Generate a vector of 151 frequencies between 1.5 and 3 GHz:
freq_vect_test = np.linspace(1.5e9,3e9,151)

#Generate a sum of two complex sine-waves corresponding to 2 targets' echoes at
#the resolution limit, corrupted by a very low amplitude white-noise:
spec_vect_test = np.exp(-1j*4*pi*freq_vect_test/3e8)+np.exp(-1j*4*pi*1.1*freq_vect_test/3e8)+np.random.normal(0,1e-6,151)

#Reference dataset generation:--------------------------------------------------

#Generate a vector of 451 frequencies between 0 and 2 GHz:
freq_vect_ref = np.linspace(0,4.5e9,451)

#Generate a sum of two complex sine-waves corresponding to 2 targets' echoes,
#but with a resolution 3X better:
spec_vect_ref = np.exp(-1j*4*pi*freq_vect_ref/3e8)+np.exp(-1j*4*pi*1.1*freq_vect_ref/3e8)

#Test function burg:------------------------------------------------------------

#Application of the burg function to the test dataset:
a,e,rc = PyBWE.burg(spec_vect_test,2)

#Check the AR model's coefficients:
assert (abs(abs(a[0])-2)<1e-3)and(abs(abs(a[1])-1)<1e-3), 'PyBWE.burg NOK!'

#Test function ar_extrapolation:------------------------------------------------

#Application of the ar_extrapolation function to the test dataset:
x_extra,x_forward,x_backward = PyBWE.ar_extrapolation(spec_vect_test,a,150,'both')

#Check the forward extrapolation:
assert sqrt(sum(abs(spec_vect_ref[301:]-x_forward)**2)/150)<1e-3, 'PyBWE.extrapolated NOK!'

#Check the backward extrapolation:
assert sqrt(sum(abs(spec_vect_ref[:150]-x_backward)**2)/150)<1e-3, 'PyBWE.extrapolated NOK!'

#Output:------------------------------------------------------------------------

print('PyBWE unit tests OK!')