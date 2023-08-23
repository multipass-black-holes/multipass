import numpy as np
import nppisn_mf as npp
import nplp_mf as nplp
import h5py
from astropy.cosmology import Planck18_arXiv_v2 as cosmo
from scipy import interpolate

mass_fns = {'ppisn':{'mf1g':npp.lmf_1g, 'mf2g':npp.lmf_2g, 'pm2m1den_m11g':npp.lpm2m1den_m11g, 'pm2m1den_m12g':npp.lpm2m1den_m12g, 'citations':npp.citations}, 'plp':{'mf1g':nplp.lmf, 'pm2m1den_m11g':nplp.lpm2m1den, 'citations':nplp.citations}}
spin_fns = {'trivial':{'spin1g1g':(lambda *args : np.ones_like(args[0])), 'spin1g2g':(lambda *args : np.ones_like(args[0])), 'spin2g1g':(lambda *args : np.ones_like(args[0])), 'spin2g2g':(lambda *args : np.ones_like(args[0]))}}
redshift_fns = {'trivial':{'z1g1g':(lambda *args : np.ones_like(args[0])), 'z2g1g':(lambda *args : np.ones_like(args[0])), 'z1g2g':(lambda *args : np.ones_like(args[0])), 'z2g2g':(lambda *args : np.ones_like(args[0]))}}

def qflatprob(m1, mmin):
    return 1./(m1-mmin)

"""
things to do:
    make general likelihood for any given mass function
"""

def av_likelihood(observables, parameters={'mgap':50, 'a':.5, 'b':-2., 'd':-3., 'mmin':5, 'deltam':5}, mass_fn = 'ppisn', secondary_prob='physical', spin_fn = 'trivial', redshift_fn = 'trivial', norms = None, change_priors = None):
    """
    parameters should be ordered [mgap, aa, bb, dd, mmin, dm] for nppisn_mf
    observables should be lists of all the same size
    """
    m1, m2, chieff, redshift = observables
    mmin, dm = parameters['mmin'], parameters['deltam']
    L_out = mass_fns[mass_fn]['mf1g'](m1, parameters)*spin_fns[spin_fn]['spin1g1g'](chieff, parameters)*redshift_fns[redshift_fn]['z1g1g'](redshift, parameters)
    if secondary_prob=='physical':
        L_out*=mass_fns[mass_fn]['mf1g'](m2, parameters)/mass_fns[mass_fn]['pm2m1den_m11g'](m1, parameters)
    else:
        if secondary_prob=='qflat':
            L_out*=qflatprob(m1, mmin)
        elif secondary_prob=='qpow':
            k = parameters['k']
            L_out*=(m2/m1)**k
        mid2 = np.where((m2 > mmin) & (m2 < mmin+dm))[0]
        stargx2 = -1. + 2.*(m2[mid2]-mmin)/dm
        L_out[mid2] *= 0.5*(1+np.tanh((stargx2+stargx2)/(1-stargx2*stargx2)))
    if ((norms is not None) & (mass_fn == 'ppisn')):
        l10lam21, l10lam12, l10lam22 = norms
        L_21 = 10**l10lam21*mass_fns[mass_fn]['mf2g'](m1, parameters)*spin_fns[spin_fn]['spin2g1g'](chieff, parameters)*redshift_fns[redshift_fn]['z2g1g'](redshift, parameters)
        L_12 = 10**l10lam12*mass_fns[mass_fn]['mf1g'](m1, parameters)*spin_fns[spin_fn]['spin1g2g'](chieff, parameters)*redshift_fns[redshift_fn]['z1g2g'](redshift, parameters)
        L_22 = 10**l10lam22*mass_fns[mass_fn]['mf2g'](m1, parameters)*spin_fns[spin_fn]['spin2g2g'](chieff, parameters)*redshift_fns[redshift_fn]['z2g2g'](redshift, parameters)
        if secondary_prob=='physical':
            L_21*=mass_fns[mass_fn]['mf1g'](m2, parameters)/mass_fns[mass_fn]['pm2m1den_m12g'](m1, parameters)
            L_12*=mass_fns[mass_fn]['mf2g'](m2, parameters)/mass_fns[mass_fn]['pm2m1den_m11g'](m1, parameters)
            L_22*=mass_fns[mass_fn]['mf2g'](m2, parameters)/mass_fns[mass_fn]['pm2m1den_m12g'](m1, parameters)
        elif secondary_prob=='qflat':
            L_21*=qflatprob(m1, mmin)
            L_12*=qflatprob(m1, mmin)
            L_22*=qflatprob(m1, mmin)
        L_out += L_21 + L_12 + L_22
    if change_priors is not None:
        L_out /= change_priors(observables)
    return np.nanmean(L_out)

def injection_extraction(ifar_find=1, with_indices=False):
    with h5py.File("o3a_bbhpop_inj_info.hdf",'r') as f:
        found_indices = np.where((np.array(f['injections/ifar_gstlal'])>ifar_find) & (np.array(f['injections/ifar_pycbc_full'])>ifar_find) & (np.array(f['injections/ifar_pycbc_bbh'])>ifar_find))[0]
        m1s = np.array(f['injections/mass1_source'])
        m2s = np.array(f['injections/mass2_source'])
        s1z = np.array(f['injections/spin1z'])
        s2z = np.array(f['injections/spin2z'])
        redshift = np.array(f['injections/redshift'])
    if with_indices:
        return np.array([m1s, m2s, (m1s*s1z+m2s*s2z)/(m1s+m2s), redshift]).T[found_indices].T, found_indices
    else:
        return np.array([m1s, m2s, (m1s*s1z+m2s*s2z)/(m1s+m2s), redshift]).T[found_indices].T

def lvc_log_xi(ifar_find=None, found_injection_data=None, parameters={'mgap':50, 'a':.5, 'b':-2., 'd':-3., 'mmin':5, 'deltam':5}, mass_fn = 'ppisn', secondary_prob='physical', spin_fn = 'trivial', redshift_fn = 'trivial', norms=None):
    """
    parameters should be ordered [mgap, aa, bb, dd, mmin, dm] for nppisn_mf
    observables should be lists of all the same size
    """
    if found_injection_data is not None:
        fid = found_injection_data
    elif ifar_find is not None:
        fid = injection_extraction(ifar_find=ifar_find, with_indices=False)
    else:
        print('gotta find some injections, bro')
    return np.log(av_likelihood(fid, parameters=parameters, mass_fn = mass_fn, secondary_prob=secondary_prob, spin_fn = spin_fn, redshift_fn = redshift_fn, norms=norms, change_priors = lambda x: x[0]**-4.35*x[1]**2))

z_table = np.linspace(0, 15, 3000)
d_L_table = cosmo.luminosity_distance(z_table)
z_func = interpolate.interp1d(d_L_table, z_table)

def lvc_data(gw_tel_number, cm=True):
    try:
        with h5py.File("all_posterior_samples/GW"+gw_tel_number+('_comoving' if cm else '')+".h5",'r') as f:
            return np.array(f['PublicationSamples/posterior_samples'])
    except Exception:
        with h5py.File("all_posterior_samples/GW"+gw_tel_number+"_GWTC-1.hdf5",'r') as ff:
            f = ff['Overall_posterior']
            dists = f['luminosity_distance_Mpc']
            m1, m2, s1, s2, c1, c2, redshifts = f['m1_detector_frame_Msun'], f['m2_detector_frame_Msun'], f['spin1'], f['spin2'], f['costilt1'], f['costilt2'], z_func(dists)
            if cm:
                m1/=(1+redshifts)
                m2/=(1+redshifts)
            return {'mass_1_source':m1, 'mass_2_source':m2, 'chi_eff':(s1*m1*c1+s2*m2*c2)/(m1+m2), 'redshift':redshifts}

def gw_event_log_likelihood(ts0 = None, gw_tel_number=None, parameters={'mgap':50, 'a':.5, 'b':-2., 'd':-3., 'mmin':5, 'deltam':5}, mass_fn = 'ppisn', secondary_prob='physical', spin_fn = 'trivial', redshift_fn = 'trivial', norms=None):
    """
    parameters should be ordered [mgap, aa, bb, dd, mmin, dm] for nppisn_mf
    observables should be lists of all the same size
    """
    if ts0 is None:
        ts0 = lvc_data(gw_tel_number=gw_tel_number)
    res = av_likelihood([ts0['mass_1_source'], ts0['mass_2_source'], ts0['chi_eff'], ts0['redshift']], parameters=parameters, mass_fn = mass_fn, secondary_prob=secondary_prob, spin_fn = spin_fn, redshift_fn = redshift_fn, norms=norms)
    if res>0:
        return np.log(res)
    else:
        return -1e6

def weighted_log_likelihood(gw_dat_list, parameters={'mgap':50, 'a':.5, 'b':-2., 'd':-3., 'mmin':5, 'deltam':5}, mass_fn = 'ppisn', secondary_prob='physical', spin_fn = 'trivial', redshift_fn = 'trivial', norms=None, ifar_find=None, found_injection_data=None):
    llist = np.array([gw_event_log_likelihood(ts0=gwd, gw_tel_number=None, parameters=parameters, mass_fn = mass_fn, secondary_prob=secondary_prob, spin_fn = spin_fn, redshift_fn = redshift_fn, norms=norms) for gwd in gw_dat_list])
    Nlxi = len(gw_dat_list)*lvc_log_xi(ifar_find=ifar_find, found_injection_data=found_injection_data, parameters=parameters, mass_fn = mass_fn, secondary_prob=secondary_prob, spin_fn = spin_fn, redshift_fn = redshift_fn, norms=norms)
    return -Nlxi+np.sum(llist)