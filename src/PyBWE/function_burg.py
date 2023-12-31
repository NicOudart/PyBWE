################################################################################
#HEADER

#This function allows you to fit an autoregressive (AR) linear model to a
#spectrum, using the Burg algorithm. The order of the model is an input of
#the algorithm, and depends on the number of complex sine-waves composing the
#signal.
#References: Burg (1967), Cokelaer et al. (2017)

#-Inputs:
#   -x: spectrum data vector
#   -p: order of the AR model

#-Outputs:
#   -a: AR model coefficient
#   -e: model fit errors
#   -rc: Burg algorithm's reflection coefficients

################################################################################

#Libraries importation:---------------------------------------------------------

import numpy as np

#Function definition:-----------------------------------------------------------

def burg(x,p):

    #Input data retrieval:
    N = len(x) #number of elements
    x = np.array(x) #numpy array conversion of the input spectrum

    #Check if the order of the model is below the number of samples:
    if p > N:
        raise ValueError("The order must be less than the number of samples")

    #Variables initialization:
    e = sum(abs(x)**2.) / float(N) #white-noise variance
    a = np.zeros(0, dtype=complex) #AR model coefficients vector
    rc = np.zeros(0, dtype=complex) #reflection coefficients vector
    f = x.astype(complex) #forward prediction errors vector
    b = x.astype(complex) #backward prediction errors vector
    rc_idx = 0 #current reflection coefficient

    #Loop of increasing model order:
    for idx_order in range(p):

        #Reflection coefficient calculation:
        num = sum([f[idx_i]*b[idx_i-1].conjugate() for idx_i in range(idx_order+1, N)])
        if idx_order==0:
            den = e*2.*N
        den = (1.-abs(rc_idx)**2)*den - abs(f[idx_order])**2 - abs(b[N-1])**2
        rc_idx = -2.* num/den
        #Add the new reflection coefficient to the list:
        rc.resize(rc.size+1)
        rc[idx_order] = rc_idx

        #Update of the estimated white-noise variance:
        e = e*(1.-abs(rc_idx)**2)

        #Check if the estimated white-noise variance is strictly positive:
        if e<0:
            print('[PyBWE.burg] Warning: Order too high - iteration stopped at '+str(idx_order))
            break

        #Add the new reflection coefficient to the AR model coefficients vector:
        a.resize(a.size+1)
        a[idx_order] = rc_idx
        #Update the AR model coefficients vector:
        if idx_order>0:
            for idx_i in range(0, (idx_order+1)//2):
                ap = a[idx_i]
                a[idx_i] = ap + rc_idx*a[idx_order-idx_i-1].conjugate()
                if idx_i != idx_order-idx_i-1:
                    a[idx_order-idx_i-1] = a[idx_order-idx_i-1] + rc_idx*ap.conjugate()

        #Update of the forward and backward prediction errors:
        for idx_i in range(N-1, idx_order, -1):
            f_old = f[idx_i]
            f[idx_i] = f_old + rc_idx * b[idx_i-1]
            b[idx_i] = b[idx_i-1] + rc_idx.conjugate()*f_old

    return a,e,rc