################################################################################
#HEADER

#This function allows you to apply the super-resolution technique named
#"Polarimetric Bandwidth Extrapolation" (PBWE) to a set of signals from
#different polarimetric channels.
#The frequency-domain spectrum of this signal is extrapolated (forward/backward)
#using a multi-channels AR model, to obtain a bandwidth larger by a given factor
#N. The resolution of the time-domain signal after IFFT will be N times better.
#This technique model the signals accounting for the different polarimetric
#channels of the radar.
#References: Suwa and Iwamoto (2003,2007)

#-Inputs:
#   -spec_mat: matrix containing the multi-channel frequency spectrum to be
#              extrapolated
#   -df: frequency step of spec_mat
#   -extra_factor: factor by which spec_mat will be extrapolated
#   -model_order: order of the AR models expressed as a ratio of the total number
#                 of samples in each spec_mat channel.
#   -zp_factor: zero-padding factor, a 10 zp_factor will lead to a X10
#               interpolation in the time-domain output signal output_bwe. This
#               process is only aesthetic.
#   -side_cut: (boolean) if True, 5% of frequencies will be cut on each side of
#              the spectrum. This greatly improves result in the case of radars
#              working in time-domain, or in case the imaginary part of the
#              spectrum has been reconstructed by Hilbert transform

#-Outputs:
#   -output_pbwe: time-domain multi-channel signal output after PBWE application
#   -time_pbwe_vect: time vector corresponding to output_pbwe

################################################################################

#Libraries importation:---------------------------------------------------------

import numpy as np

from .function_polar_burg import polar_burg
from .function_polar_extrapolation import polar_extrapolation

#Function definition:-----------------------------------------------------------

def PBWE(spec_mat,df,extra_factor,model_order,zp_factor,side_cut):

    #Convert the input spec_mat into a Numpy array:
    spec_mat = np.array(spec_mat)

    #Cutting 5% of frequencies on each side of the spectrum:
    spec_mat = spec_mat[:,round(np.shape(spec_mat)[1]*0.05):np.shape(spec_mat)[1]-round(np.shape(spec_mat)[1]*0.05)]

    #Retrieve the number of polarimetric channels:
    Npol = np.shape(spec_mat)[0]
    #Retrieve the new number of samples in the spectrum:
    M = np.shape(spec_mat)[1]

    #Order of the model expressed as a number of samples:
    model_order_nb = round(model_order*M)

    #Fit a multichannel AR model to spec_mat:
    [Thetaf,Thetab,err] = polar_burg(spec_mat,model_order_nb)

    #Calculate the number of samples to extrapolate on each side of spec_mat:
    Mextra = round(((extra_factor*M)-M)//2)+1

    #Extrapolation of spec_mat in both directions:
    spec_mat_extra,spec_mat_extra_forward,spec_mat_extra_backward = polar_extrapolation(spec_mat,Thetaf,Thetab,Mextra,"both")

    #Hamming windowing (a good compromise for resolution) and IFFT for each
    #spectrum channel:
    output_pbwe = np.zeros((Npol,round(zp_factor*np.shape(spec_mat_extra)[1])),dtype=complex)
    for idx_channel in range(Npol):
        spec_mat_extra[idx_channel,:] = np.hamming(np.shape(spec_mat_extra)[1])*spec_mat_extra[idx_channel,:]
        output_pbwe[idx_channel,:] = 1.85*np.fft.fft(np.conjugate(spec_mat_extra[idx_channel,:]),round(zp_factor*np.shape(spec_mat_extra)[1]))/np.shape(spec_mat_extra)[1]

    #Create a new time vector for output_bwe:
    time_pbwe_vect = np.linspace(0,1/df,zp_factor*np.shape(spec_mat_extra)[1])

    return output_pbwe, time_pbwe_vect