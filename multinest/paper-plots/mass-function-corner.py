import sys

sys.path.append("..")
from analyse import *

samples_ppisn = loadMC("../ppisn-cut/long", "ppisn+trivial+trivial")

fig, axs = subplots(
    7,
    7,
    gridspec_kw={"hspace": 0, "wspace": 0},
    figsize=(1.5 * 7.84, 1.5 * 5.9),
    # sharex="col",
    # sharey="row",
)

transformations = [
    ("mgap", lambda a: a, "linear"),
    ("a", lambda a: a, "linear"),
    ("b", lambda a: a, "linear"),
    ("d", lambda a: a, "linear"),
    ("lam21", lambda a: 10**a, "log"),
    ("lam12", lambda a: 10**a, "log"),
    ("bq0", lambda a: a, "linear"),
]

sharedx = axs[0,0].get_shared_x_axes()
sharedy = axs[0,0].get_shared_y_axes()


for i2, (name2, t2, s2) in enumerate(transformations):
    for i1, (name1, t1, s1) in enumerate(transformations):
        if i1 == i2:
            density = samples_ppisn.get1DDensityGridData(name1).normalize()
            axs[i1, i2].plot(t1(density.x), density.P)
            axs[i1, i2].get_yaxis().set_ticks([])
            axs[i1, i2].set_xscale(s1)

        elif i1 > i2:
            # axs[i1, i2]._shared_axes["y"].remove(axs[i1, i1])
            density = samples_ppisn.get2DDensityGridData(name2, name1).normalize()
            axs[i1, i2].contourf(
                *np.meshgrid(t2(density.x), t1(density.y)),
                density.P,
                levels=density.getContourLevels((0.95, 0.68, 1e-10))
            )

            axs[i1, i2].set_xscale(s2)
            axs[i1, i2].set_yscale(s1)
            if i2 > 0:
                axs[i1, i2].get_yaxis().set_ticklabels([])
        else:
            axs[i1, i2].set_axis_off()


axs[1,0].set_ylabel("$a$")
axs[2,0].set_ylabel("$b$")
axs[3,0].set_ylabel("$d$")
axs[4,0].set_ylabel("$\\lambda_{21}$")
axs[5,0].set_ylabel("$\\lambda_{12}$")
axs[6,0].set_ylabel("$\\beta_q$")

axs[6,0].set_xlabel(r"$M_{\rm BHMG}\,/\,M_{\odot}$")
axs[6,1].set_xlabel("$a$")
axs[6,2].set_xlabel("$b$")
axs[6,3].set_xlabel("$d$")
axs[6,4].set_xlabel("$\\lambda_{21}$")
axs[6,5].set_xlabel("$\\lambda_{12}$")
axs[6,6].set_xlabel("$\\beta_q$")

lims = [
    (14.,143.),
    (-0.14,0.64),
    (-4.48, 0.34),
    (-13.6, 3.60),
    (2e-9, 25),
    (2e-9, 25),
    (-3.6, 14.7),
]

for i in range(7):
    for j, l in enumerate(lims):
        axs[i,j].set_xlim(l)

        if i < j:
            axs[j,i].set_ylim(l)

fig.savefig("mass-function-corner.pdf")
