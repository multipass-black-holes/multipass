import sys
import numpy as np
import scipy.stats
import scipy.optimize
from matplotlib.pyplot import *

rc('text', usetex=True)

alpha = 1-scipy.stats.chi2.cdf([4,1], 1)


def get_pdf2(x, y,post, nbins=70, limits=None):
    if limits == None:
        limits = np.array([
            [x.min(), x.max()],
            [y.min(), y.max()]
        ])
    pdf, _, _ = np.histogram2d(x, y, nbins, limits, weights=post)
    pdf /= pdf.sum()
    return pdf, limits


def critical_density(pdf, alpha):
    r"""
     \int_{p > p_\text{critical}} p(x, y) dx dy = 1. - \alpha
    """
    def prob(d):
        return pdf[pdf > d].sum()
    return scipy.optimize.bisect(
        lambda d : prob(d) - (1-alpha), 0, pdf.max()
    )


#root = 'testrun/test-'
#root = 'chains2/test-'
root = sys.argv[1] + '/test-'
dat = np.loadtxt(root+'.txt')


def mkplot2(i, j, limits=None):
    post = dat[:,0]
    chis = dat[:,1]
    xdat = dat[:,i+2]
    ydat = dat[:,j+2]

    pdf, limits = get_pdf2(xdat, ydat, post, limits=limits)

    levels = [critical_density(pdf, i) for i in alpha] + [pdf.max()]

    gca().contourf(
        pdf.T,
        levels,
        extent=[i  for j in limits for i in j]
    )

    bx = xdat[chis.argmin()]
    by = ydat[chis.argmin()]

    gca().plot([bx],[by],'r*')


def get_pdf1(x,post, nbins=70, limits=None):
    if limits == None:
        limits = np.array([
            x.min(), x.max()
        ])
    pdf, bin_edges = np.histogram(x, nbins, limits, weights=post)
    pdf /= pdf.sum()

    bin_centers = (bin_edges[:-1] + bin_edges[1:]) * 0.5
    return bin_centers, pdf, limits


def mkplot1(i, limits=None):
    post = dat[:,0]
    xdat = dat[:,i+2]

    x, pdf, limits = get_pdf1(xdat, post, limits=limits)
    gca().plot(x, pdf)


plots = [
    ((2, 10), "m_{min}\\ [M_{\\odot}]"),
    ((0, 10), "\\delta_{m}\\ [M_{\\odot}]"),
    ((30, 100), "m_{max}\\ [M_{\\odot}]"),
    ((20, 50), "\\mu_{m}\\ [M_{\\odot}]"),
    ((1, 10), "\\sigma_{m}\\ [M_{\\odot}]"),
    ((-4, 12), "\\alpha"),
    ((0, 1), "\\lambda_{p}"),
    ((0, 10), "k")
]

ax = np.empty((8,8), dtype=object)
for i in range(len(plots)):
    ax[i,i] = subplot2grid(
        (8,8), (i,i)
    )
    mkplot1(i)
    for j in range(i+1, len(plots)):
        ax[j,i] = subplot2grid(
            (8,8), (j,i),
            sharex=ax[i,i],
            sharey=ax[j,0]
        )
        print(i,j, plots[i][1])
        mkplot2(
            i, j,
            [plots[i][0], plots[j][0]]
        )
        xlim(plots[i][0])
        ylim(plots[j][0])

        if i == 0:
            ylabel("$"+plots[j][1]+"$")
        if j == 7:
            xlabel("$"+plots[i][1]+"$")
    if i == 7:
        xlabel("$"+plots[i][1]+"$")

gcf().subplots_adjust(hspace=0, wspace=0)
gcf().savefig(root+'plot.pdf')
