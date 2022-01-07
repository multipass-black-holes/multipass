import h5py
import numpy as np
import struct


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



