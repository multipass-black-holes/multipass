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
    return getdist.MCSamples(
        samples=dat[:, 2:],
        weights=dat[:, 0],
        names=[j for _, _, j in plots],
        labels=[j for _, j, _ in plots]
    )


def triangle_plot(samples):
    g = getdist.plots.get_subplot_plotter()
    g.triangle_plot(
        samples,
        param_limits={i: j for j,_,i in plots},
        filled=True
    )
    return gcf()

if __name__ == "__main__":
    root = sys.argv[1] + "/test-"
    samples = loadMC(root)

    fig = triangle_plot(samples)
    fig.savefig(root + "triangle.pdf")

