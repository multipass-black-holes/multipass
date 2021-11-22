import numpy as np
import math
# from scipy.special import beta, hyp2f1, betainc, gamma
# import bincsdm as bs
mpiv = 30
mfudge = 0.99

def mf_1g(mbh, parameters={'mgap':50, 'a':.5, 'b':-2., 'd':-3., 'mmin':5, 'deltam':5}):
    """
    1g mass function for mbh
    """
    mgap, aa, bb, dd, mmin, dm = parameters['mgap'], parameters['a'], parameters['b'], parameters['d'], parameters['mmin'], parameters['deltam']
    if (mbh < mmin or mbh > mgap*mfudge):
        return 0
    else:
        dNdM = (mbh/mpiv)**bb*(1+2*aa*aa*math.sqrt(mbh/mgap)*(1-mbh/mgap)**(aa-1))
        if (mmin<mbh<mmin+dm):
            ssexparg = np.fmin(dm*(1/(mbh-mmin)+1/(mbh-mmin-dm)), 1e2)
            ss = 1/(math.exp(ssexparg)+1)
            dNdM *= ss
        return dNdM*mpiv**bb

def mf_2g(mbh, parameters={'mgap':50, 'a':.5, 'b':-2., 'd':-3., 'mmin':5, 'deltam':5}):
    """
    2g mass function for mbh
    """
    mgap, aa, bb, dd, mmin, dm = parameters['mgap'], parameters['a'], parameters['b'], parameters['d'], parameters['mmin'], parameters['deltam']
    if mbh < mgap*mfudge:
        return 0
    elif mbh < mgap + mmin + dm/2:
        return mpiv**bb
    else:
        return (mbh/(mgap + mmin + dm/2))**dd*mpiv**bb

def lmf_1g(mbh, parameters={'mgap':50, 'a':.5, 'b':-2., 'd':-3., 'mmin':5, 'deltam':5}):
    """
    1g mass function for all mbh supplied
    """
    mgap, aa, bb, dd, mmin, dm = parameters['mgap'], parameters['a'], parameters['b'], parameters['d'], parameters['mmin'], parameters['deltam']
    out = np.zeros_like(mbh)
    mid = np.where((mbh > mmin) & (mbh < mgap*mfudge))[0]
    out[mid] = (mbh[mid]/mpiv)**bb*(1+2*aa*aa*np.sqrt(mbh[mid]/mgap)*(1-mbh[mid]/mgap)**(aa-1))
    smid = np.where((mbh > mmin) & (mbh < mmin+dm))[0]
    ssexparg = np.fmin(dm*(1/(mbh[smid]-mmin)+1/(mbh[smid]-mmin-dm)), 1e2)
    ss = 1/(np.exp(ssexparg)+1)
    out[smid] *= ss
    return out*mpiv**bb

def lmf_2g(mbh, parameters={'mgap':50, 'a':.5, 'b':-2., 'd':-3., 'mmin':5, 'deltam':5}):
    """
    2g mass function for all mbh supplied
    """
    mgap, aa, bb, dd, mmin, dm = parameters['mgap'], parameters['a'], parameters['b'], parameters['d'], parameters['mmin'], parameters['deltam']
    out = np.ones_like(mbh)
    low = np.where(mbh < mmin)[0]
    out[low] = 0
    smid = np.where((mbh > mmin) & (mbh < mmin+dm))[0]
    ssexparg = np.fmin(dm*(1/(mbh[smid]-mmin)+1/(mbh[smid]-mmin-dm)), 1e2)
    ss = 1/(np.exp(ssexparg)+1)
    out[smid] *= ss
    high = np.where(mbh > mgap + mmin + dm/2)[0]
    out[high] = (mbh[high]/(mgap + mmin + dm/2))**dd
    return out*mpiv**bb

def sdmgtb(x1, x2, p, q, k):
    jarr = np.fromiter((1-q/(j+1) for j in range(k)), dtype=float)
    xindarr = np.append([1/p], np.fromiter((np.product(jarr[:kk])/(p+kk) for kk in range(1,k+1)), dtype=float)) #np.fromiter((gamma(kk+1-q)/math.factorial(kk) for kk in range(k+1)), dtype=float)/gamma(1-q)
#    x2p, x1p = x2**p, x1**p
    x2a = np.array([x2**kk for kk in range(k+1)])
    x1a = np.fromiter((x1**kk for kk in range(k+1)), dtype=float)
    xdeparr = x2**p*x2a - x1**p*x1a[:, None]
    return np.einsum('i,ij->j', xindarr, xdeparr)#np.dot(xdeparr.T, xindarr) np.sum(xindarr[:,None]*xdeparr, axis=0)

def lpm2m1den_m11g(lm1, parameters={'mgap':50, 'a':.3, 'b':-2.2, 'd':-3.2, 'mmin':5, 'deltam':5}):
    """
    returns the integral from mmin to m1 of mf_1g for all m1 in lm1; this is the denominator of the conditional probability for m2 given the assumption that m1 is in 1g
    """
    mgap, aa, bb, dd, mmin, dm = parameters['mgap'], parameters['a'], parameters['b'], parameters['d'], parameters['mmin'], parameters['deltam']
    if abs(bb)%1<1e-8:
        bb*=1.+1e-6
    if (aa%1)<1e-8:
        aa+=1e-8
#     x = np.linspace(mmin, np.fmin(lm1, mgap)).T
#     return np.array([np.trapz(lmf_1g(xx, parameters), xx) for xx in x])
    b1, b32, md = bb+1, bb+1.5, mmin+dm
    out=np.ones_like(lm1)
    low = np.where((lm1>mmin) & (lm1<mmin+dm))[0]
    out[low] = (lm1[low]-mmin)*lmf_1g(lm1[low], parameters=parameters)*0.5 # this is just an approximation
    mid = np.where((lm1>mmin+dm) & (lm1<mgap*mfudge))[0]
    lm1g, mdg = lm1[mid]/mgap, md/mgap
    t1 = (lm1[mid]**b1-md**b1)/b1
    t2 = sdmgtb(mdg, lm1g, b32, aa, 8) #pcagbinc, bs.pgtb
#     t2 =  bs.pcabinc(lm1g, b32, aa) - bs.pbinc(mdg, b32, aa)
    t2 *= 2*aa*aa*mgap**b1#*beta(b32, aa)
#     t2 = lm1g[mid]**b32*hyp2f1(b32, 1-aa, b32+1, lm1g[mid]) - mdg**b32*hyp2f1(b32, 1-aa, b32+1, mdg)
#     t2 *= 2*aa*mgap**b1
    out[mid] = (t1+t2)/mpiv + mf_1g(mmin+dm, parameters=parameters)*dm*0.5
    return out

def lpm2m1den_m12g(lm1, parameters={'mgap':50, 'a':.5, 'b':-2., 'd':-3., 'mmin':5, 'deltam':5}):
    """
    returns the integral from mmin to m1 of mf_2g for all m1 in lm1; this is the denominator of the conditional probability for m2 given the assumption that m1 is in 2g
    """
    mgap, aa, bb, dd, mmin, dm = parameters['mgap'], parameters['a'], parameters['b'], parameters['d'], parameters['mmin'], parameters['deltam']
#     x = np.linspace(mmin, lm1).T
#     return np.array([np.trapz(lmf_2g(xx, parameters), xx) for xx in x])
    out=np.ones_like(lm1)
    vlow, mmax = dm*0.5, mgap + mmin + dm/2
    low = np.where((lm1>mmin) & (lm1<mmin+dm))[0]
    out[low] = (lm1[low]-mmin)*0.5 # this is just an approximation
    mid = np.where((lm1>mmin+dm) & (lm1<mmax))[0]
    out[mid] = vlow+lm1[mid]
    high = np.where(lm1>mmax)[0]
    out[high] = vlow + mmax + ((lm1[high]/mmax)**(1+dd) - 1)*mmax/(1+dd)
    return out*mpiv**bb

citations=['2104.02685']