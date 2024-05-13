import sys

sys.path.append("..")
from analyse import *

samples_plp = loadMC("../plp/long.card", "plp+plp+trivial+trivial")
samples_ppisn = loadMC("../ppisn-cut/long", "ppisn+trivial+trivial")

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



samples_ppisnL = loadMC("../ppisn-cut-bsplit/long-low", "ppisn+trivial+trivial")
samples_ppisnH = loadMC("../ppisn-cut-bsplit/long-high", "ppisn+trivial+trivial")

fig = figure()
fill_between(*best_fit_band(samples_plp, "plp+plp+trivial+trivial", CL=0.68), alpha=0.2, label="PLP")
fill_between(*best_fit_band(samples_ppisnL, "ppisn+trivial+trivial", prescale=15, CL=0.68), alpha=0.2, label="PPISN $-4<b<-2$")
fill_between(*best_fit_band(samples_ppisnH, "ppisn+trivial+trivial", prescale=2, CL=0.68), alpha=0.2, label="PPISN $-2<b<0$")
fill_between(*best_fit_band(samples_ppisn, "ppisn+trivial+trivial", prescale=4, CL=0.68), alpha=0.1, label="PPISN $-4<b<0$", color='grey')

yscale('log')
ylabel(r"${\rm d}R/{\rm d}m_1$")
xlabel("$m_1/M_\odot$")
legend()

xlim(0, 90)
ylim(2e-5, 2e-1)

fig.tight_layout()

fig.savefig("mass-function-best-fit-bsplit.pdf")
