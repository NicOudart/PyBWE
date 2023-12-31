################################################################################
#HEADER

#This function allows you to fit a ARMA model to a spectrum signal, using a
#state-space representation.
#References: Piou (1999)

#-Inputs:
#   -y: spectrum data vector
#   -order (optional): order of the model, by default the order will be
#                      estimated by AIC.
#   -noise_type (optional): noise corrupting the spectrum, "white" for a
#                white-noise only, "colored" if a colored noise is present,
#                "white" by default.

#-Outputs:
    #-A1: model's "state matrix" estimated from the observability matrix
    #-B1: model's "input matrix" estimated from the observability matrix
    #-C1: model's "output matrix" estimated from the observability matrix
    #-A2: model's "state matrix" estimated from the controllability matrix
    #-B2: model's "input matrix" estimated from the controllability matrix
    #-C2: model's "output matrix" estimated from the controllability matrix

#NB:
#   -Unlike the Burg algorithm, the stability of the model is not guaranteed.
#   -The estimation of the state-space model from the observability or the
#    controllability matrices should yield very similar results. We name the 1st
#    estimation technique "method 1", and the 2nd "method 2".
#   -AIC = "Akaike's Information Criterion", see Akaike (1974).

################################################################################

#Libraries importation:---------------------------------------------------------

from math import floor
import numpy as np
from numpy.linalg import svd
from numpy.linalg import pinv
from numpy.linalg import matrix_power as mp
from scipy.linalg import fractional_matrix_power as fmp

from .function_AIC import AIC

#Function definition:-----------------------------------------------------------

def statespace_model(y,order=0,noise_type="white"):

    #Retrieve the number of samples in the spectrum and convert it to a
    #numpy array:
    N = len(y)
    y = np.array(y)

    #Define the correlation window length:
    L = floor(2*N/3)

    #Check if the order of the model is below 1/3 the number of samples:
    if order > N-L+1:
        raise ValueError("The order must be less than 1/3 the number of samples")

    #Create the Hankel matrix from the spectrum:
    H = np.zeros((N-L+1,L),dtype=complex)
    for idx in range(N-L+1):
        H[idx,:] = y[idx:idx+L]

    #Singular value decomposition of the Hankel matrix:
    [U,S,V] = svd(H)

    #Estimate the order of the model using Akaike's Information Criterion:
    if order==0:
        order = AIC(S,N,noise_type)
    #Else use the order defined by the user.

    #Create a diagonal matrix from the singular values, transpose V:
    S = np.diag(S)

    #Seperate "signal" from "noise" with the estimated order:
    U = U[:,:order]
    V = V[:order,:]
    S = S[:order,:order]

    #Calculate the observability and controllability matrices
    #(balanced coordinate method):
    Omega = np.matmul(U,fmp(S,0.5))
    Theta = np.matmul(fmp(S,0.5),V)

    #Method 1: calculate A from the observability matrix:
    A1 = np.matmul(pinv(Omega[:-1,:],rcond=1e-20),Omega[1:,:])
    #Method 2: calculate A from the controllability matrix:
    A2 = np.matmul(Theta[:,1:],pinv(Theta[:,:-1],rcond=1e-20))

    #Method 1: calculate C from the observability matrix:
    C1 = Omega[0,:].reshape(1,-1)
    #Method 2: calculate B from the controllability matrix:
    B2 = Theta[:,0].reshape(-1,1)

    #For the state-space impulse response to match y:
    #-Method 1: Calculate B by a least-squares method with a new observability
    #           matrix OmegaN.
    #-Method 2: Calculate C by a least-squares method with a new controllability
    #           matrix ThetaN.
    OmegaN = C1
    ThetaN = B2
    for idx in range(1,N):
        OmegaN = np.vstack((OmegaN,np.matmul(C1,mp(A1,idx))))
        ThetaN = np.hstack((ThetaN,np.matmul(mp(A2,idx),B2)))
    B1 = np.matmul(pinv(OmegaN,rcond=1e-20),y.transpose())
    C2 = np.matmul(y,pinv(ThetaN,rcond=1e-20))

    return A1,B1,C1,A2,B2,C2