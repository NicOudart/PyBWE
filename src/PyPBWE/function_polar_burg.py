################################################################################
#HEADER

#This function allows you to fit an autoregressive (AR) linear model to a set of
#spectrum signals, corresponding to different polarimetric channels of a same
#radar. This algorithm is a polarimetric version of the Burg algorithm. The
#order of the model is an input of the algorithm, and depends on the number of
#complex sine-waves composing each signal.
#References: Suwa and Iwamoto (2003,2007)

#-Inputs:
#   -X: matrix containing the spectrum for each polarimetric channel (one row
#       per channel)
#   -p: order of the AR model

#-Outputs:
#   -Thetaf: matrix of the model's coefficient in the forward direction
#   -Thetab: matrix of the model's coefficient in the backward direction
#   -err: vector of prediction errors (forward/backward)

################################################################################

#Libraries importation:---------------------------------------------------------

import numpy as np
from numpy.linalg import inv
from numpy.linalg import norm

#Function definition:-----------------------------------------------------------

def polar_burg(X,p):

    #Convert the input into a Numpy array:
    X = np.array(X)

    #Retrieve the number of polarimetric channels:
    Npol = np.shape(X)[0]
    #Retrieve the number of samples in the spectrum:
    M = np.shape(X)[1]

    #Check if the order of the model is below the number of samples:
    if p > M:
        raise ValueError("The order must be less than the number of samples")

    #Initialize the error matrices (forward/backward):
    Fmat = np.zeros((Npol,M),dtype=complex)
    Bmat = np.zeros((Npol,M),dtype=complex)
    Fmat[:,1:] = X[:,1:]
    Bmat[:,:-1] = X[:,:-1]

    #Initialize the models' errors output:
    err = np.zeros(p,dtype=float)

    #Iterate on the order of the model:
    for idx_order in range(p):

        #Initialize error matrices for this order:
        B = np.zeros((Npol,Npol),dtype=complex)
        F = np.zeros((Npol,Npol),dtype=complex)
        D = np.zeros((Npol,Npol),dtype=complex)

        #Fill the error matrices for this order:
        for idx_sample in range(idx_order+1,M):

            B += np.matmul(Bmat[:,idx_sample-idx_order-1].reshape(-1,1),Bmat[:,idx_sample-idx_order-1].conj().reshape(1,-1))
            F += np.matmul(Fmat[:,idx_sample].reshape(-1,1),Fmat[:,idx_sample].conj().reshape(1,-1))
            D += np.matmul(Fmat[:,idx_sample].reshape(-1,1),Bmat[:,idx_sample-idx_order-1].conj().reshape(1,-1))

        #Calculate the models' coefficients for this order:
        Cb = np.matmul(D.conj().T,inv(F))
        Cf = np.matmul(D,inv(B))

        #If the order is 1, add the 1st set of coefficients to the matrix:
        if idx_order==0:

            Thetaf = Cf.T
            Thetab = Cb.T

        #Else, add the new coefficients to the matrix and update the others:
        else:

            Thetaf0 = Thetaf - np.matmul(Thetab,Cf.T)
            Thetab0 = Thetab - np.matmul(Thetaf,Cb.T)
            Thetaf = np.vstack((Thetaf0,Cf.T))
            Thetab = np.vstack((Cb.T,Thetab0))

        #Calculate the prediction errors corresponding to this order:
        for idx_sample in range(idx_order+1,M):

            Ym = X[:,idx_sample].reshape(-1,1)
            Ym_1 = X[:,idx_sample-1].reshape(-1,1)

            for idx in range(idx_order):
                Ym = np.vstack((Ym,X[:,idx_sample-idx-1].reshape(-1,1)))
                Ym_1 = np.vstack((Ym_1,X[:,idx_sample-idx-2].reshape(-1,1)))

            Bmat[:,idx_sample-idx_order-1] = X[:,idx_sample-idx_order-1] - np.matmul(Thetab.T,Ym)[:,0]
            Fmat[:,idx_sample] = X[:,idx_sample] - np.matmul(Thetaf.T,Ym_1)[:,0]

        #Add the new prediction error to the output errors vector:
        e = 0
        for idx_sample in range(idx_order+1,M):
            e += (norm(Fmat[:,idx_sample])**2)+(norm(Bmat[:,idx_sample-idx_order-1])**2)
        err[idx_order] = e/(2*(M-idx_order))

    return Thetaf,Thetab,err