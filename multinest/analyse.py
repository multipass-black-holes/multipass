import sys
try:
    import getdist
    import getdist.plots
except ImportError:
    sys.path.append("../bhmf/lib/python3.10/site-packages")
    import getdist
    import getdist.plots
import numpy as np
from matplotlib.pyplot import *
import interface

rc('text', usetex=True)


plots = [
    ((2, 10), "m_{min}\\ [M_{\\odot}]", "mmin"),
    ((0, 10), "\\delta_{m}\\ [M_{\\odot}]", "dm"),
    ((30, 100), "m_{max}\\ [M_{\\odot}]", "mmax"),
    ((20, 50), "\\mu_{m}\\ [M_{\\odot}]", "mu"),
    ((1, 10), "\\sigma_{m}\\ [M_{\\odot}]", "sigma"),
    ((-4, 12), "\\alpha", "alpha"),
    ((0, 1), "\\lambda_{p}", "lp"),
    ((0, 10), "k", "k")
]


def loadMC(root):
    dat = np.loadtxt(root + '.txt')
    samples = getdist.MCSamples(
        samples=dat[:, 2:],
        weights=dat[:, 0],
        names=[j for _, _, j in plots],
        labels=[j for _, j, _ in plots]
    )
    samples.dat = dat
    return samples


def triangle_plot(samples):
    g = getdist.plots.get_subplot_plotter()
    g.triangle_plot(
        samples,
        #param_limits={i: j for j, _, i in plots},
        filled=True
    )
    return gcf()


def best_fit(samples):
    return samples.dat[samples.dat[:, 1].argmin(), 2:]


def plot_bestfit_m1(samples, model="plp+pow+trivial+trivial"):
    h = 1e-6
    n = 1000

    bestfit = best_fit(samples)

    x = np.linspace(0, 100, n)
    y = interface.pyinterface(model, 'm1', bestfit, x)

    J = np.zeros((n, len(plots)))

    mat = (1 + h * np.identity(len(plots))) * bestfit
    for i in range(len(plots)):
        yy = interface.pyinterface(model, 'm1', mat[i], x)
        J[:, i] = (yy - y) / h / bestfit[i]

    tot = J @ samples.cov() @ J.T
    err = np.sqrt(tot.diagonal())

    fig = figure()
    plot(x, y)
    fill_between(x, y-err/2, y+err/2, alpha=0.2)
    yscale('log')

    ylabel(r"${\rm d}R/{\rm d}m_1$")
    xlabel("$m_1/M_\odot$")
    xlim(0,65)
    ylim(3e-4, 1)
    return fig

if __name__ == "__main__":
    root = sys.argv[1] + "/test-"
    samples = loadMC(root)

    fig = triangle_plot(samples)
    fig.savefig(root + "triangle.pdf")

    fig = plot_bestfit_m1(samples)
    fig.savefig(root + "fitm1.pdf")
