import numpy as np
import pylab as pl
import sys
sys.path.append("../../aum/install-3.9/lib/python3.9/site-packages/")
import cosmology as cc

#cosmo = "pla"
#a = 1


def nutomass(cosmo, a, nus):
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
	
	z = 1/a - 1
	Om = aa.Omega(z)
	mass = np.logspace(10, 15, 10000)
	var = mass * 0.0
	for ii in range(mass.size):
	   var[ii] = aa.varM_TH_num(mass[ii], z)
		    
	stdev = np.sqrt(var)
	nu = (1.686*pow(Om,0.0055))/stdev

	#nus = np.array([0.75, 1, 1.25, 1.5, 1.75, 2, 2.25])
	#nus = np.array([0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25])

	mval = []
	for i in nus:
		idx = np.argwhere(np.diff(np.sign(i - nu))).flatten()
		mval.append(mass[idx][0])
	
	#print (mval)
	mval = np.array(mval)
	varcheck = mval * 0.0
	for jj in range(mval.size):
	   varcheck[jj] = aa.varM_TH_num(mval[jj], z)
	   
	stdevcheck = np.sqrt(varcheck)
	nucheck = (1.686*pow(Om,0.0055))/stdevcheck
	print(nucheck)
	
	return mval#, nus
	
#print(nutomass("pla",0.97))
   
