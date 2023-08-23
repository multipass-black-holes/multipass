import sys
import re
try:
    import getdist
    import getdist.plots
except ImportError:
    sys.path.append("../bhmf/lib/python3.10/site-packages")
    sys.path.append("../bhmf/lib/python3.11/site-packages")
    import getdist
    import getdist.plots
import numpy as np
from matplotlib.pyplot import *
import interface

rc('text', usetex=True)


parameters = {
    "plp+flat+trivial+trivial": [
        ("m_{min}\\ [M_{\\odot}]", "mmin"),
        ("\\delta_{m}\\ [M_{\\odot}]", "dm"),
        ("m_{max}\\ [M_{\\odot}]", "mmax"),
        ("\\mu_{m}\\ [M_{\\odot}]", "mu"),
        ("\\sigma_{m}\\ [M_{\\odot}]", "sigma"),
        ("\\alpha", "alpha"),
        ("\\lambda_{p}", "lp")
    ],
    "plp+pow+trivial+trivial": [
        ("m_{min}\\ [M_{\\odot}]", "mmin"),
        ("\\delta_{m}\\ [M_{\\odot}]", "dm"),
        ("m_{max}\\ [M_{\\odot}]", "mmax"),
        ("\\mu_{m}\\ [M_{\\odot}]", "mu"),
        ("\\sigma_{m}\\ [M_{\\odot}]", "sigma"),
        ("\\alpha", "alpha"),
        ("\\lambda_{p}", "lp"),
        ("k", "k")
    ],
    "ppisn+flat+trivial+trivial": [
        ("m_{min}\\ [M_{\\odot}]", "mmin"),
        ("\\delta_{m}\\ [M_{\\odot}]", "dm"),
        ("m_{gap}\\ [M_{\\odot}]", "mgap"),
        ("a", "a"),
        ("b", "b"),
        ("d", "d"),
        ("\\log_{10}\\lambda_{21}", "lam21"),
        ("\\log_{10}\\lambda_{12}", "lam12"),
    ],
    "ppisn+trivial+trivial": [
        ("m_{min}\\ [M_{\\odot}]", "mmin"),
        ("\\delta_{m}\\ [M_{\\odot}]", "dm"),
        ("m_{gap}\\ [M_{\\odot}]", "mgap"),
        ("a", "a"),
        ("b", "b"),
        ("d", "d"),
        ("\\log_{10}\\lambda_{21}", "lam21"),
        ("\\log_{10}\\lambda_{12}", "lam12"),
    ]
}

def parseInfo(root):
    with open(root + ".timing") as fp:
        txt = fp.read()
    model = re.findall("model: (.*)", txt)[0]
    np = len(parameters[model])
    lower = [
        float(i)
        for i in re.findall(
            r"prior range, lower bounds" + r" *([\d\.-]+)"*np,
            txt
        )[0]
    ]
    upper = [
        float(i)
        for i in re.findall(
            r"Prior range, upper bounds" + r" *([\d\.-]+)"*np,
            txt
        )[0]
    ]
    return model, list(zip(lower, upper))



def loadMC(root, model="plp+pow+trivial+trivial"):
    dat = np.loadtxt(root + '.txt')
    samples = getdist.MCSamples(
        samples=dat[:, 2:],
        weights=dat[:, 0],
        names=[j for _, j in parameters[model]],
        labels=[j for j, _ in parameters[model]]
    )
    samples.dat = dat
    return samples


def triangle_plot(samples, lim, model="plp+pow+trivial+trivial"):
    g = getdist.plots.get_subplot_plotter()
    g.triangle_plot(
        samples,
        param_limits={n: l for ((_, n), l) in zip(parameters[model], lim)},
        filled=True
    )
    return gcf()


def best_fit(samples):
    return samples.dat[samples.dat[:, 1].argmin(), 2:]


def plot_bestfit_m1(samples, model="plp+pow+trivial+trivial", fig=None, prescale=1, col=None):
    h = 1e-6
    n = 1000

    bestfit = best_fit(samples)

    x = np.linspace(0, 100, n)
    y = prescale * interface.pyinterface(model, 'm1', bestfit, x)

    plots = parameters[model]

    J = np.zeros((n, len(plots)))

    mat = (1 + h * np.identity(len(plots))) * bestfit
    for i in range(len(plots)):
        yy = prescale * interface.pyinterface(model, 'm1', mat[i], x)
        J[:, i] = (yy - y) / h / bestfit[i]

    tot = J @ samples.cov() @ J.T
    err = np.sqrt(tot.diagonal())

    if fig == None:
        fig = figure()

    if col:
        plot(x, y, color=col)
        fill_between(x, y-err/2, y+err/2, color=col, alpha=0.2)
    else:
        plot(x, y)
        fill_between(x, y-err/2, y+err/2, alpha=0.2)
    yscale('log')

    ylabel(r"${\rm d}R/{\rm d}m_1$")
    xlabel("$m_1/M_\odot$")
    xlim(0,65)
    ylim(3e-4, 1)
    return fig


def central_values(root, model):
    mod = parameters[model]
    with open(root + 'stats.dat') as fp:
        for line in fp.readlines():
            if m := re.match(" *(\d) +([\d.E+-]+) +([\d.E+-]+) *", line):
                ind, y, e = m.groups()
                _, para = mod[int(ind)-1]
                print(f"{para:8s} {float(y):7.2f} +- {float(e):7.2f}")


def main(roots):
    figc = figure()
    l = []
    logz = {}

    for root in roots:
        print(root)
        model, lim = parseInfo(root)

        samples = loadMC(root, model)

        fig = triangle_plot(samples, lim, model)
        fig.savefig(root + model + "triangle.pdf")

        fig = plot_bestfit_m1(samples, model)
        fig.savefig(root + model + "fitm1.pdf")

        plot_bestfit_m1(samples, model, figc)
        l.append(model)
        l.append('$1\sigma$')

        central_values(root, model)
        with open(root+'stats.dat') as fp:
            s = fp.read()
        model = root
        logz[model] = [float(i) for i in re.findall(
            '([\d\.E\+-]+) *\+/- *([\d\.E\+-]+)',
            s.splitlines()[0]
        )[0]]

    legend(l)

    for i, ki in enumerate(logz.keys()):
        for j, kj in enumerate(logz.keys()):
            if i < j:
                bf = (
                    logz[ki][0] - logz[kj][0],
                    np.sqrt(logz[ki][1]**2 + logz[kj][1]**2)
                )
                if bf[0] > 0:
                    kf = kj
                else:
                    kf = ki
                print(f"BF {ki} and {kj} = {bf[0]} +- {bf[1]} favours {kf}")

    return figc


if __name__ == "__main__":
    main(sys.argv[1:])