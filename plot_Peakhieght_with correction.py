import numpy as np
import pylab as pl
import sys
sys.path.append("../../aum/install-3.9/lib/python3.9/site-packages/")
import cosmology as cc

cosmo = "pla"
ax = pl.subplot(111)

a = 0.8081
z = 1/a - 1

if cosmo=="pla":
   aa = cc.cosmology(0.315,0.0,-1.0,0.0,0.02222/0.6731**2,0.6731,2.726,0.829,0.9655,np.log10(8.0),1.0)
   Om0 = 0.315
   h = 0.6731
elif cosmo=="bol":
   aa = cc.cosmology(0.27,0.0,-1.0,0.0,0.0469,0.7,2.726,0.82,0.95,np.log10(8.0),1.0)
   Om0 = 0.27
   h = 0.7
else:
   print(f"Cosmology {cosmo} not supported yet.")
   exit(11)

Om = aa.Omega(z)

mass = np.logspace(10, 15, 10000)
var = mass * 0.0
for ii in range(mass.size):
   var[ii] = aa.varM_TH_num(mass[ii], z)
            
stdev = np.sqrt(var)
nu = (1.686*pow(Om,0.0055))/stdev

nus = np.array([0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25])

for i in nus:
   idx = np.argwhere(np.diff(np.sign(i - nu))).flatten()
   ax.plot(mass[idx], nu[idx], 'ro')

    
ax.plot(mass, nu)
#ax.plot(mass[idx], nu[idx], 'ro')
ax.axhline(y=1)
ax.set_xscale("log")
ax.set_xlabel(r"M (${\rm h^{-1}Msun}$)")
ax.set_ylabel(r"v")
ax.set_ylim(0.3,3.2)
ax.grid()
#ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
pl.savefig(f"../Plots/PeakhieghtVsMass_a={a}.png")

