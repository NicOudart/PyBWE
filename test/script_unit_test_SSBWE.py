################################################################################
#HEADER

#This script allows you to perform a test of unit functions in the PySSBWE package.

#The test is based on synthetic radar data, generated with the following scenario:
# - 2 targets in sight of the radar, at distances 1 and 1.1 m, generate
#   echoes in the radar signal of amplitudes 1 and 1.
# - The 2 echoes are affected by frequency-domain decays of 0.25 and 0.5 1/Hz.
# - The radar output is a spectrum of 151 frequencies between 1.5 and 3 GHz.
# - The spectrum is corrupted by a white-noise of standard deviation 1e-6.

################################################################################

#Libraries importation:---------------------------------------------------------

import numpy as np
from math import pi,sqrt

import PySSBWE

#Test dataset generation:-------------------------------------------------------

#Generate a vector of 151 frequencies between 1.5 and 3 GHz:
freq_vect_test = np.linspace(1.5e9,3e9,151)

#Generate a sum of two complex sine-waves corresponding to 2 targets' echoes at
#the resolution limit, corrupted by a very low amplitude white-noise:
spec_vect_test = np.exp(((1j*4*pi/3e8)-0.25e-9)*freq_vect_test)+np.exp(((1j*4*pi*1.1/3e8)-0.5e-9)*freq_vect_test)+np.random.normal(0,1e-6,151)

#Reference dataset generation:--------------------------------------------------

#Generate a vector of 451 frequencies between 1.5 and 4.5 GHz:
freq_vect_ref = np.linspace(1.5e9,4.5e9,301)

#Generate a sum of two complex sine-waves corresponding to 2 targets' echoes,
#but with a resolution 3X better:
spec_vect_ref = np.exp(((1j*4*pi/3e8)-0.25e-9)*freq_vect_ref)+np.exp(((1j*4*pi*1.1/3e8)-0.5e-9)*freq_vect_ref)

#Test function statespace_model:------------------------------------------------

#Check the error message if the model order is too low:
try:
    A1,B1,C1,A2,B2,C2 = PySSBWE.statespace_model(np.zeros(2),-1)
    assert False, 'PySSBWE.statespace_model NOK!'
except ValueError:
    assert True

#Check the error message if the model order is too high:
try:
    A1,B1,C1,A2,B2,C2 = PySSBWE.statespace_model(np.zeros(6),4)
    assert False, 'PySSBWE.statespace_model NOK!'
except ValueError:
    assert True

#Check on the test dataset the order estimation with AIC for different noise types:
try:
    A1,B1,C1,A2,B2,C2 = PySSBWE.statespace_model(spec_vect_test,0,"white")
    A1,B1,C1,A2,B2,C2 = PySSBWE.statespace_model(spec_vect_test,0,"colored")
except:
    assert False, 'PySSBWE.statespace_model NOK!'

#Application of the statespace_model function to the test dataset:
A1,B1,C1,A2,B2,C2 = PySSBWE.statespace_model(spec_vect_test,order=2)

#Check the model's coefficients obtained with methods 1 and 2:
assert (abs(abs(A1[0,0])-1)<1e-2)and(abs(abs(A1[1,1])-1)<1e-2)and(abs(abs(A2[0,0])-1)<1e-2)and(abs(abs(A2[1,1])-1)<1e-2), 'PySSBWE.statespace_model NOK!'
assert (abs(abs(B1[0])-1.0425)<1e-3)and(abs(abs(B1[1])-0.3398)<1e-3)and(abs(abs(B2[0])-1.0425)<1e-3)and(abs(abs(B2[1])-0.3398)<1e-3), 'PySSBWE.statespace_model NOK!'
assert (abs(abs(C1[0,0])-1.0207)<1e-3)and(abs(abs(C1[0,1])-0.8164)<1e-3)and(abs(abs(C2[0])-1.0207)<1e-3)and(abs(abs(C2[1])-0.8164)<1e-3), 'PySSBWE.statespace_model NOK!'

#Test function statespace_extrapolation:----------------------------------------

#Check the error messages if the dimensions of the State-Space model matrices are not consistent:
try:
    y_extra = PySSBWE.statespace_extrapolation(np.zeros(3),np.zeros((2,3)),np.zeros((2,1)),np.zeros((1,2)),1)
    assert False, 'PySSBWE.statespace_extrapolation NOK!'
except ValueError:
    assert True
try:
    y_extra = PySSBWE.statespace_extrapolation(np.zeros(3),np.zeros((2,2)),np.zeros((3,1)),np.zeros((1,2)),1)
    assert False, 'PySSBWE.statespace_extrapolation NOK!'
except ValueError:
    assert True
try:
    y_extra = PySSBWE.statespace_extrapolation(np.zeros(3),np.zeros((2,2)),np.zeros((2,1)),np.zeros((1,3)),1)
    assert False, 'PySSBWE.statespace_extrapolation NOK!'
except ValueError:
    assert True

#Check the error message if the number of extrapolated samples is not strictly positive:
try:
    y_extra = PySSBWE.statespace_extrapolation(np.zeros(3),np.zeros((2,2)),np.zeros((2,1)),np.zeros((1,2)),0)
    assert False, 'PySSBWE.statespace_extrapolation NOK!'
except ValueError:
    assert True

#Check the statespace_extrapolation function on the test dataset models obtained with methods 1 and 2:
y_extra_1 = PySSBWE.statespace_extrapolation(spec_vect_test,A1,B1,C1,150)
y_extra_2 = PySSBWE.statespace_extrapolation(spec_vect_test,A2,B2,C2,150)
assert sqrt(sum(abs(spec_vect_ref[151:]-y_extra_1)**2)/150)<1e-3, 'PySSBWE.statespace_extrapolation NOK!'
assert sqrt(sum(abs(spec_vect_ref[151:]-y_extra_2)**2)/150)<1e-3, 'PySSBWE.statespace_extrapolation NOK!'

#Test function statespace_properties:-------------------------------------------

#Check the error messages if the dimensions of the State-Space model matrices are not consistent:
try:
    y_extra = PySSBWE.statespace_properties(np.zeros((2,3)),np.zeros((2,1)),np.zeros((1,2)),1,1)
    assert False, 'PySSBWE.statespace_properties NOK!'
except ValueError:
    assert True
try:
    y_extra = PySSBWE.statespace_properties(np.zeros((2,2)),np.zeros((3,1)),np.zeros((1,2)),1,1)
    assert False, 'PySSBWE.statespace_properties NOK!'
except ValueError:
    assert True
try:
    y_extra = PySSBWE.statespace_properties(np.zeros((2,2)),np.zeros((2,1)),np.zeros((1,3)),1,1)
    assert False, 'PySSBWE.statespace_properties NOK!'
except ValueError:
    assert True

#Check the statespace_properties results on the test dataset models obtained with methods 1 and 2:
amp_1,td_1,dec_1 = PySSBWE.statespace_properties(A1,B1,C1,freq_vect_test[1]-freq_vect_test[0],freq_vect_test[0])
amp_2,td_2,dec_2 = PySSBWE.statespace_properties(A2,B2,C2,freq_vect_test[1]-freq_vect_test[0],freq_vect_test[0])
dist_1 = td_1*0.5*3e8
dist_2 = td_2*0.5*3e8
dec_1 *= 1e9
dec_2 *= 1e9
assert (abs(abs(amp_1[0])-1)<1e-3)and(abs(abs(amp_1[1])-1)<1e-3)and(abs(abs(amp_2[0])-1)<1e-3)and(abs(abs(amp_2[1])-1)<1e-3), 'PySSBWE.statespace_properties NOK!'
assert (abs(abs(dist_1[0])-1)<1e-3)and(abs(abs(dist_1[1])-1.1)<1e-3)and(abs(abs(dist_2[0])-1)<1e-3)and(abs(abs(dist_2[1])-1.1)<1e-3), 'PySSBWE.statespace_properties NOK!'
(abs(abs(dec_1[0])-0.25)<1e-3)and(abs(abs(dec_1[1])-0.5)<1e-3)and(abs(abs(dec_2[0])-0.25)<1e-3)and(abs(abs(dec_2[1])-0.5)<1e-3), 'PySSBWE.statespace_properties NOK!'

#Output:------------------------------------------------------------------------

print('PySSBWE unit tests OK!')