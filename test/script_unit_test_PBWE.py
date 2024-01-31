################################################################################
#HEADER

#This script allows you to perform a test of unit functions in the PyPBWE package.

################################################################################

#Libraries importation:---------------------------------------------------------

import numpy as np
from math import pi,sqrt

import PyPBWE

#Test dataset generation:-------------------------------------------------------

#Generate a vector of 151 frequencies between 1.5 and 3 GHz:
freq_vect_test = np.linspace(1.5e9,3e9,151)

#Generate a sum of two complex sine-waves corresponding to 2 targets' echoes at
#the resolution limit, corrupted by a very low amplitude white-noise:
spec_vect_00_test = np.exp(-1j*4*pi*freq_vect_test/3e8)+np.exp(-1j*4*pi*1.1*freq_vect_test/3e8)+np.random.normal(0,1e-6,151)
spec_vect_11_test = np.exp(-1j*4*pi*freq_vect_test/3e8)-np.exp(-1j*4*pi*1.1*freq_vect_test/3e8)+np.random.normal(0,1e-6,151)

#Concatenate the 2 polarimetric channels in a matrix:
spec_mat_test = np.vstack((spec_vect_00_test,spec_vect_11_test))

#Reference dataset generation:--------------------------------------------------

#Generate a vector of 451 frequencies between 0 and 2 GHz:
freq_vect_ref = np.linspace(0,4.5e9,451)

#Generate a sum of two complex sine-waves corresponding to 2 targets' echoes,
#but with a resolution 3X better:
spec_vect_00_ref = np.exp(-1j*4*pi*freq_vect_ref/3e8)+np.exp(-1j*4*pi*1.1*freq_vect_ref/3e8)
spec_vect_11_ref = np.exp(-1j*4*pi*freq_vect_ref/3e8)-np.exp(-1j*4*pi*1.1*freq_vect_ref/3e8)

#Concatenate the 2 polarimetric channels in a matrix:
spec_mat_ref = np.vstack((spec_vect_00_ref,spec_vect_11_ref))

#Test function polar_burg:------------------------------------------------------

#Check the error message if the model order is too low:
try:
    Thetaf,Thetab,err = PyPBWE.polar_burg(np.zeros((2,2)),0)
    assert False, 'PyPBWE.polar_burg NOK!'
except ValueError:
    assert True

#Check the error message if the model order is too high:
try:
    Thetaf,Thetab,err = PyPBWE.polar_burg(np.zeros((2,2)),3)
    assert False, 'PyPBWE.polar_burg NOK!'
except ValueError:
    assert True

#Application of the polar_burg function to the test dataset:
Thetaf,Thetab,err = PyPBWE.polar_burg(spec_mat_test,1)

#Check the multichannel model's coefficients:
assert (abs(abs(Thetaf[0,0])-1)<1e-3)and(abs(abs(Thetaf[1,1])-1)<1e-3)and(abs(abs(Thetab[0,0])-1)<1e-3)and(abs(abs(Thetab[1,1])-1)<1e-3), 'PyPBWE.polar_burg NOK!'

#Test function polar_extrapolation:---------------------------------------------

#Check the error message if the number of coefficients is too high:
try:
    X_extra,X_forward,X_backward = PyPBWE.polar_extrapolation(np.zeros((2,3)),np.zeros((8,2)),np.zeros((8,2)),1,'both')
    assert False, 'PyPBWE.extrapolated NOK!'
except ValueError:
    assert True

#Check the error message if the 1st dimension of the coefficients matrices is not a multiple of the number of polarimetric channels:
try:
    X_extra,X_forward,X_backward = PyPBWE.polar_extrapolation(np.zeros((2,3)),np.zeros((7,2)),np.zeros((7,2)),1,'both')
    assert False, 'PyPBWE.extrapolated NOK!'
except ValueError:
    assert True

#Check the error message if the 2nd dimension of the coefficients matrices is not equal to the number of polarimetric channels:
try:
    X_extra,X_forward,X_backward = PyPBWE.polar_extrapolation(np.zeros((2,3)),np.zeros((6,3)),np.zeros((6,3)),1,'both')
    assert False, 'PyPBWE.extrapolated NOK!'
except ValueError:
    assert True

#Check the error message if the forward and backward coefficient matrices do not have the same dimensions:
try:
    X_extra,X_forward,X_backward = PyPBWE.polar_extrapolation(np.zeros((2,3)),np.zeros((6,2)),np.zeros((4,2)),1,'both')
    assert False, 'PyPBWE.extrapolated NOK!'
except ValueError:
    assert True

#Check the error message if the number of extrapolated samples is not strictly positive:
try:
    X_extra,X_forward,X_backward = PyPBWE.polar_extrapolation(np.zeros((2,3)),np.zeros((6,2)),np.zeros((6,2)),0,'both')
    assert False, 'PyPBWE.extrapolated NOK!'
except ValueError:
    assert True

#Check the 'both' extrapolation modes:
X_extra,X_forward,X_backward = PyPBWE.polar_extrapolation(spec_mat_test,Thetaf,Thetab,150,'both')
assert sqrt(np.sum(abs(spec_mat_ref[:,301:]-X_forward)**2)/300)<1e-3, 'PyPBWE.extrapolated NOK!'
assert sqrt(np.sum(abs(spec_mat_ref[:,:150]-X_backward)**2)/300)<1e-3, 'PyPBWE.extrapolated NOK!'
assert np.array_equal(X_forward,X_extra[:,301:]), 'PyPBWE.extrapolated NOK!'
assert np.array_equal(X_backward,X_extra[:,:150]), 'PyPBWE.extrapolated NOK!'

#Check the 'forward' extrapolation mode:
X_extra,X_forward,X_backward = PyPBWE.polar_extrapolation(spec_mat_test,Thetaf,Thetab,150,'forward')
assert sqrt(np.sum(abs(spec_mat_ref[:,301:]-X_forward)**2)/300)<1e-3, 'PyPBWE.extrapolated NOK!'
assert np.array_equal(X_forward,X_extra[:,151:]), 'PyPBWE.extrapolated NOK!'
assert np.sum(abs(X_backward))==0, 'PyPBWE.extrapolated NOK!'

#Check the 'backward' extrapolation mode:
X_extra,X_forward,X_backward = PyPBWE.polar_extrapolation(spec_mat_test,Thetaf,Thetab,150,'backward')
assert sqrt(np.sum(abs(spec_mat_ref[:,:150]-X_backward)**2)/300)<1e-3, 'PyPBWE.extrapolated NOK!'
assert np.array_equal(X_backward,X_extra[:,:150]), 'PyPBWE.extrapolated NOK!'
assert np.sum(abs(X_forward))==0, 'PyPBWE.extrapolated NOK!'

#Output:------------------------------------------------------------------------

print('PyPBWE unit tests OK!')