################################################################################
#HEADER

#This function allows you to extrapolate forward a spectrum signal modelled
#by an ARMA model using a state-space representation.
#References: Piou (1999)

#-Inputs:
#   -y: spectrum data vector
#   -A: y state-space model's "state matrix"
#   -B: y state-space model's "input matrix"
#   -C: y state-space model's "output matrix"
#   -Nextra: number of samples to extrapolate from y

#-Outputs:
#   -y_extra: forward extrapolation from spectrum y

################################################################################

#Libraries importation:---------------------------------------------------------

import numpy as np
from numpy.linalg import matrix_power as mp

#Function definition:-----------------------------------------------------------

def statespace_extrapolation(y,A,B,C,Nextra):

    #Check the consistency of the State-Space model matrices dimensions:
    if (np.shape(A)[0]!=np.shape(A)[1]):
        raise ValueError("The A matrix must be square")
    if (np.shape(A)[0]!=max(np.shape(B))):
        raise ValueError("The dimension of B is not consistent with A")
    if (np.shape(A)[0]!=max(np.shape(C))):
        raise ValueError("The dimension of C is not consistent with A")

    #Retrieve the number of samples in the spectrum and initialize the
    #extrapolation:
    N = len(y)
    y_extra = np.zeros(Nextra,dtype=complex)

    #Extrapolate the spectrum in the forward direction using its state-space
    #model:
    for idx in range(Nextra):
        y_extra[idx] = np.matmul(np.matmul(C,mp(A,N+idx)),B)

    return y_extra