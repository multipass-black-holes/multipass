import h5py
import numpy as np
import struct
import os
from astropy.cosmology import Planck18_arXiv_v2 as cosmo
import scipy.interpolate
import re


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


def convert_injection(ifar_find=1, fi='../endo3_mixture-LIGO-T2100113-v12.hdf5', fo='inj.rec'):
    with h5py.File(fi, 'r') as f:
        mask = [
            np.array(f['injections/ifar_gstlal']) > ifar_find,
            np.array(f['injections/ifar_pycbc_bbh']) > ifar_find
        ]
        if 'injections/ifar_pycbc_full' in f:
            mask.append(
                np.array(f['injections/ifar_pycbc_full']) > ifar_find
            )
        mask = np.all(mask, axis=0)
        m1 = f['injections/mass1_source'][mask]
        m2 = f['injections/mass2_source'][mask]
        rs = f['injections/redshift'][mask]
        s1 = f['injections/spin1z'][mask]
        s2 = f['injections/spin2z'][mask]

        dat = np.column_stack((
            m1, m2, rs,
            (m1*s1+m2*s2)/(m1+m2)
        ))

    with open(fo, 'wb') as fp:
        write_record(fp, 'i', [len(dat)])
        write_record(fp, 'd', dat)


def load_gw(veto, base='../tmp/', v=2, cm=True):
    if v == 1:
        pat = re.compile(r'GW([\d_]*)_GWTC-1.hdf5')
    elif v == 2:
        pat = re.compile(r'GW([\d_]*)_comoving.h5')
    elif v == 3:
        pat = re.compile(r'IGWN-GWTC3p0-v\d*-GW([\d_]*)_PEDataRelease_mixed_cosmo.h5')

    files = [pat.match(i) for i in os.listdir(base)]
    files = [
        base + i.group(0)
        for i in files
        if i and i.group(1) not in veto
    ]

    m1 = np.array([])
    m2 = np.array([])
    rs = np.array([])
    ce = np.array([])

    offsets = []

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
                ce = np.concatenate((
                    ce,
                    (
                        d['spin1'] * m1i * d['costilt1']
                        + d['spin2'] * m2i * d['costilt2']
                    ) / (m1i + m2i)
                ))

            elif v == 2:
                d = f['PublicationSamples/posterior_samples']
                m1 = np.concatenate((m1, d['mass_1_source']))
                m2 = np.concatenate((m2, d['mass_2_source']))
                rs = np.concatenate((rs, d['redshift']))
                ce = np.concatenate((ce, d['chi_eff']))
            elif v == 3:
                d = f['C01:Mixed/posterior_samples']
                m1 = np.concatenate((m1, d['mass_1_source']))
                m2 = np.concatenate((m2, d['mass_2_source']))
                rs = np.concatenate((rs, d['redshift']))
                ce = np.concatenate((ce, d['chi_eff']))

            offsets.append(len(m1))

    return np.column_stack((m1, m2, rs, ce)), np.array(offsets)


def convert_gw(fo='data.rec', base='../tmp/', cm=True):
    veto = set()

    # From 2104.02685
    veto.update([
        "170817",
        "190521",
        "190425",
        "190814",
        "190909_114149",
        "190719_215514",
        "190426_152155"
    ])

    # 2108.01045
    veto.update([
        "190425_081805",
        "190707_093326",
        "190720_000836",
        "190725_174728",
        "190728_064510",
        "190814_211039",
        "190924_021846",
        "190930_133541"
    ])

    # 2111.03634
    veto.update([
        "200105_162426",
        "200115_042309",
        "190426_152155"
    ])

    # 2111.03634
    veto.update([
        "170817",
        "190425",
        "200105",
        "200115",
        "190426",
        "190917"
    ])

    d1, o1 = load_gw(veto, base, v=1, cm=cm)

    d2, o2 = load_gw(veto, base, v=2, cm=cm)
    o2 += len(d1)

    d3, o3 = load_gw(veto, base, v=3, cm=cm)
    o3 += len(d1) + len(d2)

    d = np.concatenate((d1, d2, d3))
    o = np.concatenate((o1, o2, o3))

    with open(fo, 'wb') as fp:
        write_record(fp, 'i', [len(d), len(o)])
        write_record(fp, 'i', o)
        write_record(fp, 'd', d)


if __name__ == '__main__':
    convert_injection()
    convert_gw()
