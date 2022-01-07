import h5py
import numpy as np
import struct
import os


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
    if v == 2:
        pat = re.compile('GW([\d_]*)_comoving.h5')

    files = [pat.match(i) for i in os.listdir(base)]
    files = [base + i.group(0) for i in files if i]

    m1 = np.array([])
    m2 = np.array([])
    rs = np.array([])
    ce = np.array([])

    offsets = []

    for fn in files:
        with h5py.File(fn, 'r') as f:
            if v == 2:
                d = f['PublicationSamples/posterior_samples']
                m1 = np.concatenate((m1, d['mass_1_source']))
                m2 = np.concatenate((m2, d['mass_2_source']))
                rs = np.concatenate((rs, d['redshift']))
                ce = np.concatenate((ce, d['chi_eff']))

            offsets.append(len(m1))

    return np.column_stack((m1, m2, rs, ce)), np.array(offsets)

