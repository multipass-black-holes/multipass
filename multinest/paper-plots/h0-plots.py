import sys
sys.path.append("..")
from analyse import *

samples_plp = loadMC('../plp-h0/long.card', 'plp+plp+planck+trivial')
samples_ppisn = loadMC('../ppisn-h0/long', 'ppisn+planck+trivial')

fig = figure(figsize=(3,2))

density = samples_plp.get1DDensityGridData('h0').normalize()
plot(density.x, density.P, label='PLP')

density = samples_ppisn.get1DDensityGridData('h0').normalize()
plot(density.x, density.P, label='PPISN')
legend()

gca().set_xlabel(r"$H_0\,/\,{\rm km}\,{\rm s}^{-1}\,{\rm Mpc}^{-1}$")
fig.tight_layout()

fig.savefig("h0.pdf")


# corner plot, H0-mgap
densityH0 = samples_ppisn.get1DDensityGridData('h0').normalize()
densityMG = samples_ppisn.get1DDensityGridData('mgap').normalize()
density=samples_ppisn.get2DDensityGridData('mgap','h0').normalize()

fig, axs = subplots(2,2, sharex='col', sharey='row', gridspec_kw={"hspace": 0, "wspace": 0})
axs[1,1].plot(densityH0.P, densityH0.x)
axs[0,0].plot(densityMG.x, densityMG.P)
axs[1,0].contourf(*np.meshgrid(density.x,density.y),density.P, levels=density.getContourLevels((0.95,0.68,1e-10)))
axs[1,0].set_ylabel(r"$H_0\,/\,{\rm km}\,{\rm s}^{-1}\,{\rm Mpc}^{-1}$")
axs[1,0].set_xlabel(r"$M_{\rm BHMG}\,/\,M_{\odot}$")
axs[0,1].set_axis_off()
axs[0,0].get_yaxis().set_ticks([])
axs[1,1].get_xaxis().set_ticks([])

fig.savefig("h0-mgap.pdf")
