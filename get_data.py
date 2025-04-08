import requests
import re
import time
import subprocess
import tarfile
import os


def list_catalog(name):
    print(f"Fetching catalog {name}...", end="")
    alias = {
        'gwtc-1': 'GWTC-1-confident',
        'gwtc-2': 'GWTC-2',
        'gwtc-2.1': 'GWTC-2.1-confident',
        'gwtc-3': 'GWTC-3-confident',
        'all': 'GWTC'
    }
    url = f"https://www.gw-openscience.org/eventapi/jsonfull/{alias[name]}"
    catalog = requests.get(url).json()['events']
    print(f" found {len(catalog)} events")

    for k, v in catalog.items():
        p = re.match(r"(GW[\d_]*)-v\d*", k)
        if p:
            if v['far'] > 1:
                print(p.group(1), v['far'])
            yield (p.group(1), v['jsonurl'])


def parse_date(date):
    return time.mktime(time.strptime(date, "%Y-%m-%d"))


def get_download(url):
    m = requests.get(url).json()['events']
    assert len(m) == 1
    m = next(iter(m.values()))

    parameters = [
        (k, parse_date(v['date_added']))
        for k,v in m['parameters'].items()
        if 'pe' in k
    ]
    assert len(parameters) > 0
    parameters = sorted(parameters, key=lambda a: a[1], reverse=True)
    return [
        m['parameters'][k[0]]['data_url']
        for k in parameters
    ]


def download_file(url, cwd="tmp/"):
    if url.endswith("/content"):
        filename = url.split("/")[-2]
    elif url.endswith("?download=1"):
        filename = url.split("/")[-2][:-11]
    elif url.endswith(".h5") or url.endswith(".hdf5") or url.endswith(".hdf") or url.endswith('.tar'):
        filename = url.split("/")[-1]
    else:
        print(url)

    if not os.path.exists(f"{cwd}{filename}"):
        args = ['wget', url, '-O', filename]
        print(args)
        # subprocess.Popen(args, cwd=cwd).wait()


def untar(fn):
    t = tarfile.open(fn)
    pat = re.compile(r'.*/(GW[\d_]*_comoving.h5)')
    for i in t:
        p = pat.match(i.name)
        if p:
            with open("tmp/" + p.group(1), "wb") as fp:
                fp.write(t.extractfile(i).read())
            return
    raise KeyError


def list_catalogs(lst):
    for name in lst:
        yield from list_catalog(name)


if __name__ == "__main__":
    failed_downloads = []
    for k, url in list_catalogs(['all', 'gwtc-2']):
        try:
            print(f"Getting metadata for {k}...", end="")
            urls = get_download(url)
            print(f" done")
            for url in urls:
                download_file(url)
        except KeyboardInterrupt:
            break
        except:
            import traceback
            failed_downloads.append(k)
            print("fail in",k)
            traceback.print_exc()

    print(failed_downloads)

    # for i in os.listdir("tmp/"):
    #     if i.endswith(".tar"):
    #         print(f"Untarring {i}...", end="")
    #         untar("tmp/" + i)
    #         print(" done")

    download_file("https://zenodo.org/records/5546676/files/endo3_mixture-LIGO-T2100113-v12.hdf5?download=1", cwd=".")
    download_file("https://zenodo.org/records/7890398/files/o1+o2+o3_mixture_real+semianalytic-LIGO-T2100377-v2.hdf5?download=1", cwd=".")
    download_file("https://dcc-llo.ligo.org/public/0168/P2000217/002/o3a_bbhpop_inj_info.hdf", cwd=".")
