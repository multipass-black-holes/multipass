import sys

sys.path.append("..")
from analyse import *

samples_plp = loadMC("../plp/long.card", "plp+plp+trivial+trivial")
samples_ppisn = loadMC("../ppisn/long", "ppisn+trivial+trivial")

fig = figure()
fill_between(*best_fit_band(samples_plp, "plp+plp+trivial+trivial", CL=0.68), alpha=0.2, label="PLP")
fill_between(*best_fit_band(samples_ppisn, "ppisn+trivial+trivial", prescale=15, CL=0.68), alpha=0.2, label="PPISN")

yscale('log')
ylabel(r"${\rm d}R/{\rm d}m_1$")
xlabel("$m_1/M_\odot$")
legend()

xlim(0, 90)
ylim(2e-5, 2e-1)

fig.tight_layout()

fig.savefig("mass-function-best-fit.pdf")
