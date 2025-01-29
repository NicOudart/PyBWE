################################################################################
#HEADER

#This function allows you to fit an autoregressive (AR) linear model to a
#spectrum, using the Modified Covariance algorithm. The order of the model is 
#an input of the algorithm, and depends on the number of complex sine-waves 
#composing the signal.
#References: Marple (1987), Cokelaer et al. (2017)

#-Inputs:
#   -x: spectrum data vector
#   -p: order of the AR model

#-Outputs:
#   -a: AR model coefficient
#   -e: model fit errors

################################################################################

#Libraries importation:---------------------------------------------------------

import numpy as np
from scipy.linalg import toeplitz,lstsq

#Function definition:-----------------------------------------------------------

def mcov(x,p):

    #Input data retrieval:
    N = len(x) #number of elements
    x = np.array(x) #numpy array conversion of the input spectrum

    #Check if the order of the model is below the number of samples:
    if p > N:
        raise ValueError("The order must be less than the number of samples")

    #Check if the order of the model is above 0:
    if p < 1:
        raise ValueError("The order cannot be 0 or negative")

    mat_toep = toeplitz(x[p:],x[p::-1]) #Create a Toeplitz matrix from the input data
    
    #From this matrix create the 'covariance' matrix linking forward-backward prediction errors to the model coefficients:
    mat_fb = np.zeros((2*(N-p),p+1),dtype='complex')
    mat_fb[:N-p,:] = mat_toep #Upper part of the matrix (forward) 
    mat_fb[N-p:2*(N-p),:] = np.conjugate(np.fliplr(mat_toep)) #Lower part of the matrix (backward)
    mat_fb/=(2*(N-p))**0.5 #Scaling
            
    #Calculate the model coefficients and prediction error by solving this set of equations with the least-square method:
    a,_,_,_ = lstsq(-mat_fb[:,1:],mat_fb[:,0])
    
    #Estimate the white-noise variance
    e = np.real(np.matmul(np.conjugate(mat_fb[:,0].T),mat_fb[:,0]) + np.matmul(np.matmul(np.conjugate(mat_fb[:,0].T),mat_fb[:,1:]),a))

    return a,e