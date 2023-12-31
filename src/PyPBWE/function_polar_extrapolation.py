################################################################################
#HEADER

#This function allows you to extrapolate a set of spectrum signals
#corresponding to different radar polarimetric channels, modelled by an
#autoregressive (AR) model.
#References: Suwa and Iwamoto (2003,2007)

#-Inputs:
#   -X: matrix containing the spectrum for each polarimetric channel (one row
#       per channel)
#   -Thetaf: coefficients matrix of the AR forward model of X
#   -Thetab: coefficients matrix of the AR backward model of X
#   -Mextra: number of samples to be extrapolated forward, backward or both
#   -extra_mode: extrapolation mode, "forward","backward" or "both"

#-Outputs:
#   -X_extra: matrix containing the extrapolated X in the desired direction(s)
#   -X_forward: matrix containing the forward extrapolation of X
#   -X_backward: matrix containing the backward extrapolation of X

################################################################################

#Libraries importation:---------------------------------------------------------

import numpy as np

#Function definition:-----------------------------------------------------------

def polar_extrapolation(X,Thetaf,Thetab,Mextra,extra_mode):

    #Retrieve the number of samples in the spectrum:
    M = np.shape(X)[1]

    #Retrieve the number of polarimetric channels:
    Npol = np.shape(X)[0]
    #Retrieve the order of the model:
    order = round((np.shape(Thetaf)[0])/Npol)

    #Check if the order of the model is below the number of samples:
    if order > M:
        raise ValueError("The order must be less than the number of samples")

    #Initialize the extrapolated spectrum matri and add X to it:
    if (extra_mode=="both"):
        X_extra = np.zeros((Npol,M+(2*Mextra)),dtype=complex)
        X_extra[:,Mextra:Mextra+M] = X
    elif (extra_mode=="forward"):
        X_extra = np.zeros((Npol,M+Mextra),dtype=complex)
        X_extra[:,:M] = X
    else:
        X_extra = np.zeros((Npol,M+Mextra),dtype=complex)
        X_extra[:,Mextra:] = X

    #Extrapolation in both directions:
    if (extra_mode=="both"):

        #Forward direction:
        for idx_extra in range(Mextra+M,(2*Mextra)+M):

            #Reshape X into a new matrix named Y for extrapolation:
            Y = X_extra[:,idx_extra-1].reshape(-1,1)
            for idx_sample in range(order-1):
                Y = np.vstack((Y,X_extra[:,idx_extra-idx_sample-2].reshape(-1,1)))

            #Extrapolation
            X_extra[:,idx_extra] = np.matmul(Thetaf.T,Y)[:,0]

        #Backward direction:
        for idx_extra in range(Mextra,-1,-1):

            #Reshape X into a new matrix named Y for extrapolation:
            Y = X_extra[:,idx_extra+order].reshape(-1,1)
            for idx_sample in range(order-1):
                Y = np.vstack((Y,X_extra[:,idx_extra+order-idx_sample-1].reshape(-1,1)))

            #Extrapolation:
            X_extra[:,idx_extra] = np.matmul(Thetab.T,Y)[:,0]

        X_forward = X_extra[:,Mextra+M:(2*Mextra)+M]
        X_backward = X_extra[:,:Mextra]

    #Extrapolation in the forward direction:
    if (extra_mode=="forward"):

        for idx_extra in range(M,Mextra+M):

            #Reshape X into a new matrix named Y for extrapolation:
            Y = X_extra[:,idx_extra-1].reshape(-1,1)
            for idx_sample in range(order-1):
                Y = np.vstack((Y,X_extra[:,idx_extra-idx_sample-2].reshape(-1,1)))

            #Extrapolation
            X_extra[:,idx_extra] = np.matmul(Thetaf.T,Y)[:,0]

        X_forward = X_extra[:,M:Mextra+M]
        X_backward = np.array([])

    #Extrapolation in the backward direction:
    if (extra_mode=="backward"):

        for idx_extra in range(Mextra,-1,-1):

            #Reshape X into a new matrix named Y for extrapolation:
            Y = X_extra[:,idx_extra+order].reshape(-1,1)
            for idx_sample in range(order-1):
                Y = np.vstack((Y,X_extra[:,idx_extra+order-idx_sample-1].reshape(-1,1)))

            #Extrapolation:
            X_extra[:,idx_extra] = np.matmul(Thetab.T,Y)[:,0]

        X_forward = np.array([])
        X_backward = X_extra[:,:Mextra]

    return X_extra,X_forward,X_backward