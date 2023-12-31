################################################################################
#HEADER

#This function allows you to extrapolate a spectrum signal modelled by an
#autoregressive (AR) model.
#References: Bowling (1977), Cuomo (1992)

#-Inputs:
#   -x: spectrum data vector
#   -ar_coeff: coefficients of the AR model of x
#   -Nextra: number of samples to be extrapolated forward, backward or both
#   -extra_mode: extrapolation mode, "forward","backward" or "both"

#-Outputs:
#   -x_extra: extrapolated spectrum x
#   -x_forward: forward extrapolation of spectrum x
#   -x_backward: backward extrapolation of spectrum x

################################################################################

#Libraries importation:---------------------------------------------------------

import numpy as np

#Function definition:-----------------------------------------------------------

def ar_extrapolation(x,ar_coeff,Nextra,extra_mode):

    x_forward = np.zeros((Nextra),dtype=complex)
    x_backward = np.zeros((Nextra),dtype=complex)

    #Retrieve the numbers of samples in the spectrum:
    N = len(x)

    #Check if the number of coefficients is below the number of samples:
    if len(ar_coeff) > N:
        raise ValueError("The number of AR coefficients must be less than the number of samples")

    #Add the x vector data to x_extra:
    if (extra_mode=="both"):
        x_extra = np.zeros((N+(2*Nextra)),dtype=complex)
        x_extra[Nextra:Nextra+N] = x
    elif (extra_mode=="forward"):
        x_extra = np.zeros((N+Nextra),dtype=complex)
        x_extra[:N] = x
    else:
        x_extra = np.zeros((N+Nextra),dtype=complex)
        x_extra[Nextra:] = x

    #Extrapolation in both directions:
    if (extra_mode=="both"):

        for idx_extra in range(Nextra+N,(2*Nextra)+N):
            sumcoeff = 0
            for idx_arcoeff in range(len(ar_coeff)):
                sumcoeff = sumcoeff - (ar_coeff[idx_arcoeff]*x_extra[idx_extra-idx_arcoeff-1])
            x_extra[idx_extra] = sumcoeff

        for idx_extra in range(Nextra-1,-1,-1):
            sumcoeff = 0
            for idx_arcoeff in range(len(ar_coeff)):
                sumcoeff = sumcoeff - (np.conjugate(ar_coeff[idx_arcoeff])*x_extra[idx_extra+idx_arcoeff+1])
            x_extra[idx_extra] = sumcoeff

        x_forward = x_extra[Nextra+N:]
        x_backward = x_extra[:Nextra]

    #Extrapolation in the forward direction:
    if (extra_mode=="forward"):

        for idx_extra in range(N,Nextra+N):
            sumcoeff = 0
            for idx_arcoeff in range(len(ar_coeff)):
                sumcoeff = sumcoeff - (ar_coeff[idx_arcoeff]*x_extra[idx_extra-idx_arcoeff-1])
            x_extra[idx_extra] = sumcoeff

        x_forward = x_extra[N:]

    #Extrapolation in the backward direction:
    if (extra_mode=="backward"):

        for idx_extra in range(Nextra-1,-1,-1):
            sumcoeff = 0
            for idx_arcoeff in range(len(ar_coeff)):
                sumcoeff = sumcoeff - (np.conjugate(ar_coeff[idx_arcoeff])*x_extra[idx_extra+idx_arcoeff+1])
            x_extra[idx_extra] = sumcoeff

        x_backward = x_extra[:Nextra]

    return x_extra,x_forward,x_backward