################################################################################
#HEADER

#This function allows you to estimate the order of a State Space BWE (SSBWE)
#model from the SVD of its Hankel matrix, based on the "two-line fit" approach.
#References: Brindise & Vlachos (2017).

#-Inputs:
#   -sv: vector of the Hankel matrix singular values, in decreasing order

#-Outputs:
#   -output_order: order of the model (number of sources) estimated by the
#                  "two-line fit" approach.

################################################################################

#Libraries importation:---------------------------------------------------------

import numpy as np
from math import sqrt
from scipy.stats import linregress

#Function definition:-----------------------------------------------------------

def two_line_fit(sv):

    #Ensure we have singular values expressed as modules:
    sv = abs(sv)

    #Retrieve the number of singular values:
    Q = len(sv)

    #Initialize the two-line fit criteria's list:
    TLF = np.zeros(Q-2)

    #Calculate the total sum of squares of singular values:
    SSTot = np.sum((sv-np.mean(sv))**2)

    #Retrieve the total number of singular values:
    N = np.arange(Q)

    #Initialize the singular values fit vector:
    sv_fit = np.zeros(Q)

    #For increasing model order fit 2 lines to the singular values:
    for idx_model in range(1,Q-1):

        #Separate the sigular values in "signal" and "noise" groups:
        sv_sig = sv[:idx_model+1]
        sv_noise = sv[idx_model:]

        #Select the indices corresponding corresponding to these groups:
        N_sig = N[:idx_model+1]
        N_noise = N[idx_model:]

        #Linear regression of each group:
        reg_sig = linregress(N_sig,sv_sig)
        reg_noise = linregress(N_noise,sv_noise)

        #Calculate linear fit values for each group:
        sv_sig_fit = reg_sig.intercept + reg_sig.slope*N_sig
        sv_noise_fit = reg_noise.intercept + reg_noise.slope*N_noise

        #Calculate the error between the linear fits and the singular values:
        error_fit = np.sum(np.abs(sv_sig-sv_sig_fit))+np.sum(np.abs(sv_noise-sv_noise_fit))
        #Calculate the R2 coefficient of the linear fits:
        R2_fit = 1 - (np.sum((sv_sig-sv_sig_fit)**2)+np.sum((sv_noise-sv_noise_fit)**2))/SSTot

        #If the error is not zero add R2/error to the :
        if error_fit!=0:
            TLF[idx_model-1] = R2_fit/error_fit

    #Choose the order of the model by maximizing the criterion:
    output_order = np.argmax(TLF)+1

    return output_order
