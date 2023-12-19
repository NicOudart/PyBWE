################################################################################
#HEADER

#This function allows you to apply the super-resolution technique we named
#"State-Space Bandwidth Extrapolation" (SSBWE) to a signal. The frequency-domain
#spectrum of this signal is extrapolated (forward/backward) using an ARMA model
#with state-space representation, to obtain a bandwidth larger by a given factor
#N. The resolution of the time-domain signal after IFFT will be N times better.
#This technique accounts for white-noise and non-linear effects corrupting the
#signal.
#References: Piou (1999)

#-Inputs:
#   -spec: frequency spectrum to be extrapolated
#   -df: frequency step of spec
#   -extra_factor: factor by which spec will be extrapolated
#   -zp_factor: zero-padding factor, a 10 zp_factor will lead to a X10
#               interpolation in the time-domain output signal output_bwe. This
#               process is only aesthetic.
#   -side_cut: (boolean) if True, 5% of frequencies will be cut on each side of
#              the spectrum. This greatly improves result in the case of radars
#              working in time-domain, or in case the imaginary part of the
#              spectrum has been reconstructed by Hilbert transform
#   -order (optional): order of the model, by default the order will be
#                      estimated by AIC.
#   -noise_type (optional): noise corrupting the spectrum, "white" for a
#                white-noise only, "colored" if a colored noise is present,
#                "white" by default.

#-Outputs:
#   -output_ssbwe_1: extrapolated spectrum using SSBWE method 1
#   -output_ssbwe_2: extrapolated spectrum using SSBWE method 2
#   -time_ssbwe_vect: time axis corresponding to output_ssbwe_1 and
#                     output_ssbwe_2

#NB:
#   -Unlike the BWE, the SSBWE does not guaranty the stability of the model.
#   -The estimation of the state-space model from the observability or the
#    controllability matrices should yield very similar results. We name the 1st
#    estimation technique "method 1", and the 2nd "method 2".
#   -AIC = "Akaike's Information Criterion", see Akaike (1974).

################################################################################

#Libraries importation:---------------------------------------------------------

import numpy as np

from .function_statespace_model import statespace_model
from .function_statespace_extrapolation import statespace_extrapolation

#Function definition:-----------------------------------------------------------

def SSBWE(spec,df,extra_factor,zp_factor,side_cut,order=0,noise_type="white"):

    #Cutting 5% of frequencies on each side of the spectrum:
    spec = spec[round(len(spec)*0.05):len(spec)-round(len(spec)*0.05)]

    #Convert y to Numpy array:
    y = np.array(spec)
    #Backward version of spectrum y:
    y_b = np.flip(y)

    #Retrieve the numbers of samples:
    N = len(y) #In the spectrum
    Nextra = round(((extra_factor*N)-N)//2)+1 #To be extrapolated

    #Initialization of the extrapolated spectrum vector:
    y_extra_1 = np.zeros((N+(Nextra*2)),dtype=complex) #Method 1
    y_extra_2 = np.zeros((N+(Nextra*2)),dtype=complex) #Method 2

    #Add the original spectrum to the output:
    y_extra_1[Nextra:Nextra+N] = y #Method 1
    y_extra_2[Nextra:Nextra+N] = y #Method 2

    #Fit a state-space model to spectrum y for forward extrapolation:
    [A1_f,B1_f,C1_f,A2_f,B2_f,C2_f] = statespace_model(y,order,noise_type)

    #Fit a state-space model to spectrum y for backward extrapolation:
    [A1_b,B1_b,C1_b,A2_b,B2_b,C2_b] = statespace_model(y_b,order,noise_type)

    #Forward extrapolation of spectrum y:
    y_extra_1_f = statespace_extrapolation(y,A1_f,B1_f,C1_f,Nextra) #Method 1
    y_extra_2_f = statespace_extrapolation(y,A2_f,B2_f,C2_f,Nextra) #Method 2

    #Backward extrapolation of spectrum y:
    y_extra_1_b = statespace_extrapolation(y_b,A1_b,B1_b,C1_b,Nextra) #Method 1
    y_extra_2_b = statespace_extrapolation(y_b,A2_b,B2_b,C2_b,Nextra) #Method 2

    #Assemble the extrapolated spectrum:
    y_extra_1[:Nextra] = np.flip(y_extra_1_b) #Method 1
    y_extra_1[Nextra+N:] = y_extra_1_f #Method 1
    y_extra_2[:Nextra] = np.flip(y_extra_2_b) #Method 2
    y_extra_2[Nextra+N:] = y_extra_2_f #Method 2

    #Hamming windowing (a good compromise for resolution) and IFFT:
    y_extra_1 = np.hamming(len(y_extra_1))*y_extra_1 #Method 1
    y_extra_2 = np.hamming(len(y_extra_2))*y_extra_2 #Method 2
    output_ssbwe_1 = 1.85*np.fft.fft(np.conjugate(y_extra_1),round(zp_factor*len(y_extra_1)))/len(y_extra_1) #Method 1
    output_ssbwe_2 = 1.85*np.fft.fft(np.conjugate(y_extra_2),round(zp_factor*len(y_extra_2)))/len(y_extra_2) #Method 2

    #Create a new time vector for output_bwe:
    time_ssbwe_vect = np.linspace(0,1/df,zp_factor*len(y_extra_1))

    return output_ssbwe_1, output_ssbwe_2, time_ssbwe_vect