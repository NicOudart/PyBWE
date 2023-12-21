################################################################################
#HEADER

#This function allows you to apply the super-resolution technique named
#"Bandwidth Extrapolation" (BWE) to a signal. The frequency-domain spectrum of
#this signal is extrapolated (forward/backward) using an AR model to obtain a
#bandwidth larger by a given factor N. The resolution of the time-domain signal
#after IFFT will be N times better.
#References: Bowling (1977), Cuomo (1992)

#-Inputs:
#   -spec: frequency spectrum to be extrapolated
#   -df: frequency step of spec
#   -extra_factor: factor by which spec will be extrapolated
#   -model_order: order of the AR model expressed as a ratio of the total number
#                 of samples in spec
#   -zp_factor: zero-padding factor, a 10 zp_factor will lead to a X10
#               interpolation in the time-domain output signal output_bwe. This
#               process is only aesthetic.
#   -side_cut: (boolean) if True, 5% of frequencies will be cut on each side of
#              the spectrum. This greatly improves result in the case of radars
#              working in time-domain, or in case the imaginary part of the
#              spectrum has been reconstructed by Hilbert transform

#-Outputs:
#   -output_bwe: time-domain signal output after BWE application
#   -time_bwe_vect: time vector corresponding to output_bwe

#NB: The BWE assumes the input frequency-domain spectrum is a sum of complex
#sine-waves, and thus follows an autoregressive (AR) model. Any noise or
#non-linear effects will impact the quality of the extrapolation.

################################################################################

#Libraries importation:---------------------------------------------------------

import numpy as np
from statistics import variance as var

from .function_burg import burg
from .function_ar_extrapolation import ar_extrapolation

#Function definition:-----------------------------------------------------------

def BWE(spec,df,extra_factor,model_order,zp_factor,side_cut):

    #Cutting 5% of frequencies on each side of the spectrum:
    spec = spec[round(len(spec)*0.05):len(spec)-round(len(spec)*0.05)]

    #Retrieve the number of samples in the spectrum:
    N = len(spec)

    #Order of the model expressed as a number of samples:
    model_order_nb = round(model_order*N)

    #Fit an AR model to the spectrum using the Burg algorithm:
    ar_coeff,ar_var,ar_rc = burg(spec,model_order_nb)

    #Number of samples to be extrapolated on each side of the spectrum:
    Nextra = round(((extra_factor*N)-N)//2)+1

    #Extrapolation on each side of the spectrum:
    spec_extra,spec_forward,spec_backward = ar_extrapolation(spec,ar_coeff,Nextra,"both")

    #Hamming windowing (a good compromise for resolution) and IFFT:
    spec_extra = np.hamming(len(spec_extra))*spec_extra
    output_bwe = 1.85*np.fft.fft(np.conjugate(spec_extra),round(zp_factor*len(spec_extra)))/len(spec_extra)

    #Create a new time vector for output_bwe:
    time_bwe_vect = np.linspace(0,1/df,zp_factor*len(spec_extra))

    return output_bwe, time_bwe_vect