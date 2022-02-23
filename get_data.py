import requests
import re
import time
import subprocess


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
    k = sorted(parameters, key=lambda a: a[1], reverse=True)[0][0]
    return m['parameters'][k]['data_url']


def download_file(url):
    subprocess.Popen(['wget', url], cwd='tmp/').wait()


if __name__ == "__main__":
    urls = []
    for k, url in list_catalog('all'):
        try:
            print(f"Getting metadata for {k}...", end="")
            url = get_download(url)
            print(f" done")
            urls.append(url)
            download_file(url)
        except:
            print(" failed!")