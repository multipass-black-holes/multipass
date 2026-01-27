import h5py
import numpy as np
import struct
import os
from astropy.cosmology import Planck15 as cosmo
import scipy.interpolate
import re
import numpy.typing
from typing import Any, IO, Optional

nda = numpy.typing.NDArray[np.float64]

z_table = np.linspace(0, 15, 3000)
d_L_table = cosmo.luminosity_distance(z_table)
zfunc = scipy.interpolate.make_interp_spline(d_L_table, z_table, k=1)
invzfunc = scipy.interpolate.make_interp_spline(z_table, d_L_table, k=1)


def write_record(
    fp: IO[bytes], typ: str, content: nda | list[int] | list[float]
) -> None:
    if type(content) == np.ndarray:
        content = list(content.flatten("F"))

    body = b"".join(struct.pack("<" + typ, i) for i in content)
    fp.write(struct.pack("<I", len(body)))
    fp.write(body)
    fp.write(struct.pack("<I", len(body)))


def read_record(fp: IO[bytes], typ: str) -> Any:
    l1 = struct.unpack("<I", fp.read(4))[0]
    body_raw = fp.read(l1)
    l2 = struct.unpack("<I", fp.read(4))[0]
    n = l1 // struct.calcsize(typ)
    body = struct.unpack("<" + typ * n, body_raw)
    return body


def extract_columns_default(d: dict[str, nda]) -> nda:
    return np.column_stack(
        (
            d["mass_1_source"],
            d["mass_2_source"],
            d["mass_1"],
            d["mass_2"],
            d["redshift"],
            d["luminosity_distance"],
            np.sqrt(d["spin_1x"] ** 2 + d["spin_1y"] ** 2 + d["spin_1z"] ** 2),
            np.sqrt(d["spin_2x"] ** 2 + d["spin_2y"] ** 2 + d["spin_2z"] ** 2),
        )
    )


def extract_columns_GWTC1(d: dict[str, nda]) -> nda:
    rsi = zfunc(d["luminosity_distance_Mpc"])
    m1i = d["m1_detector_frame_Msun"]
    m2i = d["m2_detector_frame_Msun"]

    return np.column_stack(
        (
            d["m1_detector_frame_Msun"] / (1 + rsi),
            d["m2_detector_frame_Msun"] / (1 + rsi),
            d["m1_detector_frame_Msun"],
            d["m2_detector_frame_Msun"],
            zfunc(d["luminosity_distance_Mpc"]),
            d["luminosity_distance_Mpc"],
            d["spin1"],
            d["spin2"],
        )
    )


def load_all_files(base: str = "../tmp") -> dict[str, dict[str, nda]]:
    regexprs = [
        (
            re.compile(r"IGWN-GWTC.p.-(v\d*)-GW([\d_]*)_PEDataRelease_mixed_cosmo.h5"),
            "C01:Mixed/posterior_samples",
            extract_columns_default,
            lambda m: m.groups(),
        ),
        (
            re.compile(r"IGWN-GWTC.p.-(v\d*)-GW([\d_]*)_PEDataRelease.h5"),
            "PrecessingSpinIMRHM/posterior_samples",
            extract_columns_default,
            lambda m: m.groups(),
        ),
        (
            re.compile(r"GW([\d_]*)_comoving.h5"),
            "PublicationSamples/posterior_samples",
            extract_columns_default,
            lambda m: ("v0", m.group(1)),
        ),
        (
            re.compile(r"GW([\d_]*)_GWTC-1.hdf5"),
            "Overall_posterior",
            extract_columns_GWTC1,
            lambda m: ("v0", m.group(1)),
        ),
    ]

    events: dict[str, dict[str, nda]] = {"v0": {}, "v1": {}, "v2": {}}
    for i in os.listdir(base):
        if i.endswith(".tar"):
            continue

        for pat, key, extract, label in regexprs:
            if m := pat.match(i):
                with h5py.File(base + "/" + i, "r") as f:
                    try:
                        d = f[key]
                    except KeyError:
                        print("no posterior", i, f.keys())
                        break
                    group, ev = label(m)
                    if ev in events[group]:
                        print("Collision", i)
                    events[group][ev] = extract(d)

                    break
        else:
            print("bad file", i)

    return events


def load_rec(path: str) -> tuple[nda, list[nda]]:
    fp = open(path, "rb")
    lend, leno = read_record(fp, "i")
    o = np.array([0] + list(read_record(fp, "i")))
    _d = read_record(fp, "d")
    assert leno == len(o) - 1
    assert 8 * lend == len(_d)
    fp.close()

    d = np.array(_d).reshape((lend, 8), order="F")

    assert np.all(d[:, 0] > d[:, 1])
    assert np.all(d[:, 2] > d[:, 3])

    return o, [d[i:j] for i, j in zip(o[:-1], o[1:])]


def match_old_new(
    oldpath: str, events: Optional[dict[str, dict[str, nda]]] = None
) -> list[tuple[Optional[str], Optional[str]]]:
    if events is None:
        events = load_all_files()

    o, eventsold = load_rec(oldpath)
    old_events: list[tuple[Optional[str], Optional[str]]] = [(None, None)] * len(
        eventsold
    )
    diffs = list(np.diff(o))
    for campaign, e in events.items():
        for i, j in e.items():
            try:
                n = diffs.index(j.shape[0])
                if old_events[n][0]:
                    print("Collision", n, old_events[n], i)
                old_events[n] = i, campaign
                assert np.all(np.isclose(j, eventsold[n]))
            except ValueError:
                pass
    return old_events


def get_veto() -> set[str]:
    veto = set()

    # From 2104.02685
    veto.update(
        [
            "170817",
            "190521",
            "190425",
            "190814",
            "190909_114149",
            "190719_215514",
            "190426_152155",
        ]
    )

    # 2108.01045
    veto.update(
        [
            "190425_081805",
            "190707_093326",
            "190720_000836",
            "190725_174728",
            "190728_064510",
            "190814_211039",
            "190924_021846",
            "190930_133541",
        ]
    )

    # 2111.03634
    veto.update(["200105_162426", "200115_042309", "190426_152155"])

    # 2111.03634
    veto.update(
        ["170817", "190425", "200105", "200115", "190426", "190426_190642", "190917"]
    )

    # others
    veto.update(
        [
            "190917_114630",
            "191219_163120",
            # "200210_092254"
        ]
    )

    return veto


DEFAULT_LIST = [
    "170729",
    "170814",
    "170608",
    "170809",
    "151012",
    "170823",
    "150914",
    "151226",
    "170104",
    "170818",
    "190707_093326",
    "190803_022701",
    "190512_180714",
    "190929_012149",
    "190513_205428",
    "190602_175927",
    "190701_203306",
    "190930_133541",
    "190706_222641",
    "190828_063405",
    "190521",
    "190708_232457",
    "190412",
    "190517_055101",
    "190620_030421",
    "190719_215514",
    "190413_052954",
    "190910_112807",
    "190720_000836",
    "190731_140936",
    "190413_134308",
    "190924_021846",
    "190728_064510",
    "190503_185404",
    "190521_074359",
    "190408_181802",
    "190828_065509",
    "190519_153544",
    "190915_235702",
    "190527_092055",
    "190421_213856",
    "190630_185205",
    "190727_060333",
    "190925_232845",
    "190805_211137",
    "190725_174728",
    "191127_050227",
    "191109_010717",
    "191105_143521",
    "191129_134029",
    "191204_171526",
    "200216_220804",
    "200311_115853",
    "200225_060421",
    "191222_033537",
    "200112_155838",
    "200224_222234",
    "200302_015811",
    "191215_223052",
    "191216_213338",
    "200129_065458",
    "191230_180458",
    "200208_130117",
    "200209_085452",
    "191103_012549",
    "200219_094415",
    "200202_154313",
    "200316_215756",
    "200128_022011",
]


def convert_events(
    fo: str,
    events: Optional[dict[str, nda]] = None,
    all_events: Optional[dict[str, dict[str, nda]]] = None,
    lst: list[str] = DEFAULT_LIST,
) -> tuple[list[int], nda]:
    if isinstance(all_events, dict):
        # favour newer events
        myevents = {}
        myevents.update(all_events["v0"])
        myevents.update(all_events["v1"])
        myevents.update(all_events["v2"])
    elif isinstance(events, dict):
        myevents = events
    else:
        raise KeyError

    o = []
    d = np.zeros((0, 8))
    for name in lst:
        d = np.concatenate((d, myevents[name]))
        o.append(len(d))

    with open(fo, "wb") as fp:
        write_record(fp, "i", [len(d), len(o)])
        write_record(fp, "i", o)
        write_record(fp, "d", d)

    return o, d


def convert_injection(
    ifar_find: float = 1,
    fi: str = "../o1+o2+o3_mixture_real+semianalytic-LIGO-T2100377-v2.hdf5",
    fo: str = "inj.rec",
    version: int = 4,
    auto_sampling_pdf: bool = True,
):
    with h5py.File(fi, "r") as f:
        if version == 2:
            mask = (
                (np.array(f["injections/ifar_gstlal"]) > ifar_find)
                & (np.array(f["injections/ifar_pycbc_full"]) > ifar_find)
                & (np.array(f["injections/ifar_pycbc_bbh"]) > ifar_find)
            )

            s1 = f["injections/spin1z"][mask]
            s2 = f["injections/spin2z"][mask]

            m1D = f["injections/mass1_source"][mask] * (
                1 + f["injections/redshift"][mask]
            )
            m2D = f["injections/mass2_source"][mask] * (
                1 + f["injections/redshift"][mask]
            )
            ld = f["injections/distance"][mask]

            if auto_sampling_pdf:
                raise ValueError("not supported")
            else:
                pdf = m1D**-4.35 * m2D**2
        elif version == 3:
            mask = [
                np.array(f["injections/ifar_gstlal"]) > ifar_find,
                np.array(f["injections/ifar_pycbc_bbh"]) > ifar_find,
            ]
            if "injections/ifar_pycbc_full" in f:
                mask.append(np.array(f["injections/ifar_pycbc_full"]) > ifar_find)
            mask = np.all(mask, axis=0)

            s1 = np.sqrt(
                f["injections/spin1x"][mask] ** 2
                + f["injections/spin1y"][mask] ** 2
                + f["injections/spin1z"][mask] ** 2
            )
            s2 = np.sqrt(
                f["injections/spin2x"][mask] ** 2
                + f["injections/spin2y"][mask] ** 2
                + f["injections/spin2z"][mask] ** 2
            )

            m1D = f["injections/mass1"][mask]
            m2D = f["injections/mass2"][mask]
            ld = f["injections/distance"][mask]

            if auto_sampling_pdf:
                pdf = f["injections/sampling_pdf"][mask]
            else:
                pdf = m1D**-4.35 * m2D**2

        elif version == 4:
            injO1 = np.array(f['injections/name']) == b'o1'
            injO2 = np.array(f['injections/name']) == b'o2'
            injO3 = np.array(f['injections/name']) == b'o3'
            snr_cut = np.array(f['injections/optimal_snr_net'])>6
            far_cut = np.any([
                np.array(f["injections/ifar_cwb"]) > ifar_find,
                np.array(f["injections/ifar_gstlal"]) > ifar_find,
                np.array(f["injections/ifar_mbta"]) > ifar_find,
                np.array(f["injections/ifar_pycbc_bbh"]) > ifar_find,
                np.array(f["injections/ifar_pycbc_hyperbank"]) > ifar_find,
            ], axis=0)
            mask = np.any([
                np.all([np.any([injO1, injO2], axis=0), snr_cut, far_cut], axis=0),
                np.all([injO3, far_cut], axis=0)
            ], axis=0)

            s1 = np.sqrt(
                f["injections/spin1x"][mask] ** 2
                + f["injections/spin1y"][mask] ** 2
                + f["injections/spin1z"][mask] ** 2
            )
            s2 = np.sqrt(
                f["injections/spin2x"][mask] ** 2
                + f["injections/spin2y"][mask] ** 2
                + f["injections/spin2z"][mask] ** 2
            )

            m1D = f["injections/mass1_source"][mask] * (
                1 + f["injections/redshift"][mask]
            )
            m2D = f["injections/mass2_source"][mask] * (
                1 + f["injections/redshift"][mask]
            )
            ld = invzfunc(f["injections/redshift"][mask])

            pdf = f["injections/sampling_pdf"][mask]

        m1 = f["injections/mass1_source"][mask]
        m2 = f["injections/mass2_source"][mask]
        rs = f["injections/redshift"][mask]
        pdf /= np.abs(invzfunc.derivative()(f["injections/redshift"]))
        # pdf = m1**-4.35 * m2**2

        dat = np.column_stack((m1, m2, m1D, m2D, rs, ld, s1, s2, pdf))

    with open(fo, "wb") as fp:
        write_record(fp, "i", [len(dat)])
        write_record(fp, "d", dat)


if __name__ == "__main__":
    convert_injection(fi="../o1+o2+o3_bbhpop_real+semianalytic-LIGO-T2100377-v2.hdf5")
    o, d = convert_events("data.rec", all_events=load_all_files())
