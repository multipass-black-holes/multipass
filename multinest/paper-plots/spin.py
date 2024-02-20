import sys
sys.path.append("..")
from analyse import *

matplotlib.rcParams.update({'font.size': 14})

samples = loadMC('../ppisn-spin-mass/medium-narrow', 'ppisn+trivial+beta-mass')

fig, axs = subplots(
    4,
    4,
    gridspec_kw={"hspace": 0, "wspace": 0},
    # figsize=(1.5 * 7.84, 1.5 * 5.9),
    # sharex="col",
    # sharey="row",
)

vars = parameters["ppisn+trivial+beta-mass"]

for i2, (latex, name2) in enumerate(vars):
    for i1, (_, name1) in enumerate(vars):
        if i1 == i2:
            density = samples.get1DDensityGridData(name1).normalize()
            axs[i1, i2].plot(density.x, density.P)
            axs[i1, i2].get_yaxis().set_ticks([])

        elif i1 > i2:
            # axs[i1, i2]._shared_axes["y"].remove(axs[i1, i1])
            density = samples.get2DDensityGridData(name2, name1).normalize()
            axs[i1, i2].contourf(
                *np.meshgrid(density.x, density.y),
                density.P,
                levels=density.getContourLevels((0.95, 0.68, 1e-10))
            )

            if i2 > 0:
                axs[i1, i2].get_yaxis().set_ticklabels([])
        else:
            axs[i1, i2].set_axis_off()


for i, (latex, _) in enumerate(vars[1:]):
    axs[i+1,0].set_ylabel(f"${latex}$")

for i, (latex, _) in enumerate(vars):
    axs[-1,i].set_xlabel(f"${latex}$")

lims = [
    (-.7,11.1),
    (-.7,11.1),
    (-.7,11.1),
    (-.7,11.1),
]

for i in range(4):
    for j, l in enumerate(lims):
        axs[i,j].set_xlim(l)

        if i < j:
            axs[j,i].set_ylim(l)

fig.savefig("spin-corner.pdf")
