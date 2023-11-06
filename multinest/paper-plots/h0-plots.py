import sys
sys.path.append("..")
from analyse import *

samples_plp = loadMC('../plp-h0/long.card', 'plp+plp+planck+trivial')
samples_ppisn = loadMC('../ppisn-h0/long', 'ppisn+planck+trivial')

fig = figure(figsize=(3,2))

density = samples_plp.get1DDensityGridData('h0')
density.normalize()
plot(density.x, density.P, label='PLP')

density = samples_ppisn.get1DDensityGridData('h0')
density.normalize()
plot(density.x, density.P, label='PPISN')
legend()

gca().set_xlabel(r"$H_0\,/\,{\rm km}\,{\rm s}^{-1}\,{\rm Mpc}^{-1}$")
fig.tight_layout()

fig.savefig("h0.pdf")
