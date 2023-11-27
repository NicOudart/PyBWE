################################################################################
#HEADER

#This function allows you to retrieve the properties of each echo (time-delay,
#amplitude and frequency-domain decay) in a radar signal based on the ARMA model
#of its frequency spectrum, using state-space representation.
#References: Piou (1999)

#-Inputs:
#   -A: state-space model's "state matrix"
#   -B: state-space model's "input matrix"
#   -C: state-space model's "output matrix"
#   -df: frequency step in the modelled spectrum
#   -f1: 1st frequency in the modelled spectrum

#-Outputs:
#   -amp: vector containing the estimated complex amplitudes of echoes
#   -td: vector containing the estimated time-delays of echoes (s)
#   -dec: vector containing the estimated frequency-domain decays of echoes (1/Hz)

################################################################################

#Libraries importation:---------------------------------------------------------

import numpy as np
from numpy.linalg import eig
from numpy.linalg import inv
from math import pi

#Function definition:-----------------------------------------------------------

def statespace_properties(A,B,C,df,f1):

    #Retrieve the eigenvalues and eigenvector of A:
    eigval, eigvect = eig(A)
    #Calculate the inverse of A's eigenvector:
    inv_eigvect = inv(eigvect)

    #Initialize the output vectors:
    amp = np.zeros(len(eigval),dtype=complex)
    td = np.zeros(len(eigval),dtype=float)
    dec = np.zeros(len(eigval),dtype=float)

    #For each eigenvalue in A (= for each echo in the signal):
    for idx in range(len(eigval)):

        #Fill to the output vectors of amplitudes, time-delays and
        #frequency-domain decays:
        amp[idx] = np.matmul(C,eigvect[:,idx])*np.matmul(inv_eigvect[idx,:],B)/(eigval[idx]**(f1/df))
        td[idx] = -np.angle(eigval[idx])/(2*pi*df)
        dec[idx] = -np.log(abs(eigval[idx]))/df

    return amp,td,dec