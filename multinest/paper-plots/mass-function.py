import sys

sys.path.append("..")
from analyse import *

samples_plp = loadMC("../plp/long.card", "plp+plp+trivial+trivial")
samples_ppisn = loadMC("../ppisn/long", "ppisn+trivial+trivial")

fig = figure()
plot_bestfit_m1(samples_plp, "plp+plp+trivial+trivial", fig=fig, label="PLP")
plot_bestfit_m1(samples_ppisn, "ppisn+trivial+trivial", fig=fig, label="PPISN")

legend()

xlim(0, 90)
ylim(2e-5, 2e-1)

fig.tight_layout()

fig.savefig("mass-function-best-fit.pdf")
