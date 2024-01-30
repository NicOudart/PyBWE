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

#Check the error message if the model order is too low:
try:
    a,e,rc = PyBWE.burg(np.zeros(2),0)
    assert False, 'PyBWE.burg NOK!'
except ValueError:
    assert True

#Check the error message if the model order is too high:
try:
    a,e,rc = PyBWE.burg(np.zeros(2),3)
    assert False, 'PyBWE.burg NOK!'
except ValueError:
    assert True

#Application of the burg function to the test dataset:
a,e,rc = PyBWE.burg(spec_vect_test,2)

#Check the AR model's coefficients:
assert (abs(abs(a[0])-2)<1e-3)and(abs(abs(a[1])-1)<1e-3), 'PyBWE.burg NOK!'

#Test function ar_extrapolation:------------------------------------------------

#Check the error message if the number of coefficients is too high:
try:
    x_extra,x_forward,x_backward = PyBWE.ar_extrapolation(np.zeros(2),np.zeros(3),1,'both')
    assert False, 'PyBWE.extrapolated NOK!'
except ValueError:
    assert True

#Check the error message if the number of extrapolated samples is not strictly positive:
try:
    x_extra,x_forward,x_backward = PyBWE.ar_extrapolation(np.zeros(2),np.zeros(2),0,'both')
    assert False, 'PyBWE.extrapolated NOK!'
except ValueError:
    assert True

#Testing the different extrapolation modes:
for extra_mode in ['both','forward','backward']:

    #Application of the ar_extrapolation function to the test dataset:
    x_extra,x_forward,x_backward = PyBWE.ar_extrapolation(spec_vect_test,a,150,'both')

    #Check the forward extrapolation:
    if (extra_mode=='both')or(extra_mode=='forward'):
        assert sqrt(sum(abs(spec_vect_ref[301:]-x_forward)**2)/150)<1e-3, 'PyBWE.extrapolated NOK!'
        assert np.array_equal(x_forward,x_extra[301:]), 'PyBWE.extrapolated NOK!'

    #Check the backward extrapolation:
    if (extra_mode=='both')or(extra_mode=='backward'):
        assert sqrt(sum(abs(spec_vect_ref[:150]-x_backward)**2)/150)<1e-3, 'PyBWE.extrapolated NOK!'
        assert np.array_equal(x_backward,x_extra[:150]), 'PyBWE.extrapolated NOK!'

#Output:------------------------------------------------------------------------

print('PyBWE unit tests OK!')