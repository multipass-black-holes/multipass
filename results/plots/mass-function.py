import sys

sys.path.append("../")
from analyse import *

samples_plp = loadMC("../plp/long", "plp+plp+trivial+trivial")
samples_ppisn = loadMC("../ppisn/long", "ppisn+trivial+trivial")

m1, yaplp, ybplp, ycplp = best_fit_band_norms(samples_plp, "plp+plp+trivial+trivial", CL=0.68)
m1, yappisn, ybppisn, ycppisn = best_fit_band_norms(samples_ppisn, "ppisn+trivial+trivial", CL=0.68)

fig = figure()
fill_between(m1, yaplp, ybplp, alpha=0.2, label="PLP")
plot(m1, ycplp, 'C0')
fill_between(m1, 1.5e-2*yappisn, 1.5e-2*ybppisn, alpha=0.2, label="PPISN")
plot(m1, 1.5e-2*ycppisn, 'C1')

yscale('log')
ylim(2e-10,5e3)
ylabel(r"${\rm d}R/{\rm d}m_1$")
xlabel("$m_1/M_\odot$")
legend()

# mmin,dm,mgap,_,_,_,_,_,_=best_fit(samples_ppisn)
# axvline(mmin, color='black', linewidth=0.4, zorder=1)
# axvline(mmin+dm, color='black', linewidth=0.4, zorder=1)
# axvline(mgap, color='black', linewidth=0.4, zorder=1)

fig.savefig("mass-function-best-fit.pdf")
