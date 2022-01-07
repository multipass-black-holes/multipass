import h5py
import numpy as np
import struct
import os
from astropy.cosmology import Planck18_arXiv_v2 as cosmo
import scipy


z_table = np.linspace(0, 15, 3000)
d_L_table = cosmo.luminosity_distance(z_table)
zfunc = scipy.interpolate.interp1d(d_L_table, z_table)


def write_record(fp, typ, content):
    if type(content) == np.ndarray:
        content = list(content.flatten('F'))

    body = b''.join(
        struct.pack('<'+typ, i) for i in content
    )
    fp.write(struct.pack('<I', len(body)))
    fp.write(body)
    fp.write(struct.pack('<I', len(body)))


def convert_injection(fi='../o3a_bbhpop_inj_info.hdf', fo='inj.rec'):
    with h5py.File(fi, 'r') as f:
        dat = np.column_stack((
            f['injections/ifar_gstlal'],
            f['injections/ifar_pycbc_full'],
            f['injections/ifar_pycbc_bbh'],
            f['injections/mass1_source'],
            f['injections/mass2_source'],
            f['injections/spin1z'],
            f['injections/spin2z'],
            f['injections/redshift']
        ))

    with open(fo, 'wb') as fp:
        write_record(fp, 'd', dat)


def load_gw(base='../all_posterior_samples/', v=2, cm=True):
    if v == 1:
        pat = re.compile('GW([\d_]*)_GWTC-1.hdf5')
    elif v == 2:
        pat = re.compile('GW([\d_]*)_comoving.h5')

    files = [pat.match(i) for i in os.listdir(base)]
    files = [base + i.group(0) for i in files if i]

    m1 = np.array([])
    m2 = np.array([])
    rs = np.array([])
    ce = np.array([])

    offsets = []

    if v == 1:
        files = [files[0]]

    for fn in files:
        with h5py.File(fn, 'r') as f:
            if v == 1:
                d = f['Overall_posterior']

                rsi = zfunc(d['luminosity_distance_Mpc'])
                m1i = d['m1_detector_frame_Msun']
                m2i = d['m2_detector_frame_Msun']

                if cm:
                    m1i /= (1+rsi)
                    m2i /= (1+rsi)

                m1 = np.concatenate((m1, m1i))
                m2 = np.concatenate((m2, m2i))
                rs = np.concatenate((rs, rsi))
                ce = np.concatenate((ce,
                    (
                        d['spin1'] * m1i * d['costilt1']
                       +d['spin2'] * m2i * d['costilt2']
                    ) / (m1 + m2)
                ))

            elif v == 2:
                d = f['PublicationSamples/posterior_samples']
                m1 = np.concatenate((m1, d['mass_1_source']))
                m2 = np.concatenate((m2, d['mass_2_source']))
                rs = np.concatenate((rs, d['redshift']))
                ce = np.concatenate((ce, d['chi_eff']))

            offsets.append(len(m1))

    return np.column_stack((m1, m2, rs, ce)), np.array(offsets)

