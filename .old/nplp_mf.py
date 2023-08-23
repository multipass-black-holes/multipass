import numpy as np
import math
mpiv = 30
normal_norm = (2*math.pi)**-0.5

def mf(mbh, parameters={'mmax':50, 'mum':30, 'sm':5, 'alpha':2, 'lambdap':.2, 'mmin':5, 'deltam':5}):
    """
    1g mass function for mbh
    """
    mmax, mum, sm, alpha, lp, mmin, dm = parameters['mmax'], parameters['mum'], parameters['sm'], parameters['alpha'], parameters['lambdap'], parameters['mmin'], parameters['deltam']
    lp = 10**lp if lp<0 else lp
    if (mbh < mmin or mbh > mmax):
        return 0
    else:
        dNdM = (1-lp)*mbh**(-alpha)*(1-alpha)/(mmax**(1-alpha)-mmin**(1-alpha)) + lp*math.exp(-0.5*((mbh - mum)/sm)**2)*normal_norm/sm
        if (mmin<mbh<mmin+dm):
            ssexparg = np.fmin(dm*(1/(mbh-mmin)+1/(mbh-mmin-dm)), 1e2)
            ss = 1/(math.exp(ssexparg)+1)
            dNdM *= ss
        return dNdM

def lmf(mbh, parameters={'mmax':50, 'mum':30, 'sm':5, 'alpha':2, 'lambdap':.2, 'mmin':5, 'deltam':5}):
    """
    1g mass function for all mbh supplied
    """
    mmax, mum, sm, alpha, lp, mmin, dm = parameters['mmax'], parameters['mum'], parameters['sm'], parameters['alpha'], parameters['lambdap'], parameters['mmin'], parameters['deltam']
    lp = 10**lp if lp<0 else lp
    out = np.zeros_like(mbh)
    mid = np.where((mbh > mmin) & (mbh < mmax))[0]
    out[mid] = (1-lp)*mbh[mid]**(-alpha)*(1-alpha)/(mmax**(1-alpha)-mmin**(1-alpha)) + lp*np.exp(-0.5*((mbh[mid] - mum)/sm)**2)*normal_norm/sm
    smid = np.where((mbh > mmin) & (mbh < mmin+dm))[0]
#     ssexparg = np.fmin(dm*(1/(mbh[smid]-mmin)+1/(mbh[smid]-mmin-dm)), 1e2)
#     ss = 1/(np.exp(ssexparg)+1)
#     out[smid] *= ss
    stargx2 = -1. + 2.*(mbh[smid]-mmin)/dm
    out[smid] *= 0.5*(1+np.tanh((stargx2+stargx2)/(1-stargx2*stargx2)))
    return out

def lpm2m1den(lm1, parameters={'mmax':50, 'mum':30, 'sm':5, 'alpha':2, 'lambdap':.2, 'mmin':5, 'deltam':5}):
    """
    returns the integral from mmin to m1 of mf_1g for all m1 in lm1; this is the denominator of the conditional probability for m2 given the assumption that m1 is in 1g
    """
    mmax, mum, sm, alpha, lp, mmin, dm = parameters['mmax'], parameters['mum'], parameters['sm'], parameters['alpha'], parameters['lambdap'], parameters['mmin'], parameters['deltam']
    lp = 10**lp if lp<0 else lp
    out=np.ones_like(lm1)
    low = np.where((lm1>mmin) & (lm1<mmin+dm))[0]
    out[low] = (lm1[low]-mmin)*lmf_1g(lm1[low], parameters=parameters)*0.5 # this is just an approximation
    mid = np.where((lm1>mmin+dm) & (lm1<mgap*mfudge))[0]
    lm1g, mdg = lm1[mid]/mgap, md/mgap
    t1 = (lm1[mid]**b1-md**b1)/b1
    t2 = sdmgtb(mdg, lm1g, b32, aa, 8)
    t2 *= 2*aa*aa*mgap**b1
    out[mid] = (t1+t2)/mpiv + mf_1g(mmin+dm, parameters=parameters)*dm*0.5
    return out

citations=['2010.14533']