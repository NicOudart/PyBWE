################################################################################
#HEADER

#This function allows you to estimate the order of a State Space BWE (SSBWE)
#model from the SVD of its Hankel matrix, based on Minimum Description Length
#(MDL).
#References: Wax & Kailath (1985).

#-Inputs:
#   -sv: vector of the Hankel matrix singular values, in decreasing order
#   -N: number of samples in the spectrum corresponding to the Hankel matrix

#-Outputs:
#   -output_order: order of the model (number of sources) estimated by MDL.

################################################################################

#Libraries importation:---------------------------------------------------------

import numpy as np
from math import sqrt

#Function definition:-----------------------------------------------------------

def MDL(sv,N):

    #Retrieve the number of singular values:
    Q = len(sv)

    #Initialize the information criteria's list:
    IC = np.zeros(Q-1)

    #Retrieve the module of the singular values:
    lam=abs(sv)

    #For increasing model order calculate MDL:
    for ord in range(1,Q):

        #Calculate the maximum likelihood for this order:
        num = 1
        for idx in range(ord,Q):
            num *= lam[idx]**(1/(Q - ord))
        den = (1/(Q - ord))*sum(lam[ord:Q])
        lik = num/den

        #Calculate MDL for this order:
        crit = -N*(Q - ord)*np.log(lik) + 0.5*ord*(2*Q - ord)*np.log(N)

        #Add this MDL to the criteria's list:
        IC[ord-1] = crit

    #Choose the order of the model by minimizing MDL:
    output_order = np.argmin(IC)+1

    return output_order
