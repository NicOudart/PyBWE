################################################################################
#HEADER

#This function allows you to estimate the order of a State Space BWE (SSBWE)
#model from the SVD of its Hankel matrix, based on Akaike's Information
#Criterion (AIC).
#References: Akaike (1974)

#-Inputs:
#   -sv: vector of the Hankel matrix singular values, in decreasing order
#   -N: number of samples in the spectrum corresponding to the Hankel matrix

#-Outputs:
#   -output_order: order of the model (number of sources) estimated by AIC.

################################################################################

#Libraries importation:---------------------------------------------------------

import numpy as np
from math import sqrt

#Function definition:-----------------------------------------------------------

def AIC(sv,N):

    #Initialize the information criteria's list:
    IC = []

    #Retrieve the number of singular values:
    Q = len(sv)

    #Retrieve the module of the singular values:
    lam = abs(sv)

    #For increasing model order calculate AIC:
    for ord in range(1,Q):

        #Calculate the maximum likelihood for this order:
        num = 1
        for idx in range(ord,Q):
            num *= lam[idx]**(1/(Q-ord))
        den = (1/(Q-ord))*sum(lam[ord:Q])
        lik = num/den

        #Calculate AIC for this order:
        crit = (-2*N*(Q-ord)*np.log(lik))+(2*ord*(2*Q-ord))

        #Add this AIC to the criteria's list:
        IC += [crit]

    #Choose the order of the model by minimizing AIC:
    output_order = np.argmin(IC)+1

    return output_order