import numpy as np
import pylab as pl
import sys
sys.path.append("../../aum/install-3.9/lib/python3.9/site-packages/")
import cosmology as cc
import glob
import math
import mass_from_peakhieght as mfp


colordict={'0_1':'tab:pink','1_2':'tab:cyan', '2_3':'tab:olive', '3_4':'tab:green', '4_5':'tab:red', '5_6':'tab:blue', '6_7':'tab:purple'} 
ptsdict={500:'10', 250:'8', 125:'4'}
dirname = sys.argv[1]
flist = np.sort(np.array(glob.glob(f"{dirname}/moria*powspec*")))
cosmo = dirname.split("Data/c")[1].split("_")[0]
if cosmo=="pla":
    aa = cc.cosmology(0.315,0.0,-1.0,0.0,0.02222/0.6731**2,0.6731,2.726,0.829,0.9655,np.log10(8.0),1.0)
elif cosmo=="bol":
    aa = cc.cosmology(0.27,0.0,-1.0,0.0,0.0469,0.7,2.726,0.82,0.95,np.log10(8.0),1.0)
else:
    print(f"Cosmology {cosmo} not supported yet.")
    exit(11)
    
massbins, nubin = mfp.nutomass(cosmo)


Lbox = float(dirname.split("_l")[-1].split("_")[0])
ax = pl.subplot(111)
for i, fname in enumerate(flist):
    if Lbox==500:
    	i=i+1
    fin = open(fname.replace("_powspec", ""))
    Npart = int(fin.readline())
    fin.close()
    print(Lbox, Npart)
    
    kn, Pkrough, wN, dPk = np.loadtxt(fname, usecols=(0, 3, 4, 5), unpack=True)
    Pkmeas = Pkrough - wN/Npart  #removed shot noise
    dPkmeas = Pkmeas*dPk         #statistical error
    
    kk = kn*2*np.pi/Lbox         #normalised k 
    
    
    #delete values where k=0
    index=[]
    for ii in range(kk.size):
         if kk[ii] == 0:
            index.append(ii)
    kk = np.delete(kk, index)
    Pkmeas = np.delete(Pkmeas, index)
    dPkmeas = np.delete(dPkmeas, index)
    
    #remove values where kk>0.2
    for jj in range(kk.size):       
    	if kk[jj]>0.21:
            break
    
    kk = np.delete(kk, range(jj,kk.size))
    Pkmeas = np.delete(Pkmeas, range(jj,Pkmeas.size))
    dPkmeas = np.delete(dPkmeas, range(jj, dPkmeas.size))
    print(kk.size,Pkmeas.size,dPkmeas.size)
    
    
    #remove values where Pkmeas is -ve
    index=[]
    print(kk.size,Pkmeas.size,dPkmeas.size)
    for ii in range(kk.size):
         if Pkmeas[ii] < 0:
            index.append(ii)
    kk = np.delete(kk, index)
    Pkmeas = np.delete(Pkmeas, index)
    dPkmeas = np.delete(dPkmeas, index)
    
    #find mattermatter power spectrum
    Pkmatter = kk * 0.0          
    for ii in range(Pkmatter.size):
          Pkmatter[ii] = aa.Delta2_NL_num(kk[ii], 0.0)*(2*np.pi**2/kk[ii]**3)
        #ax.plot(kk[1:], Pktheory[1:], marker =".")
    
    #bias = root(Pkdark/Pkmatter)
    Pkdark = Pkmeas*Lbox**3
    bias = (Pkdark/Pkmatter)**0.5
    dbias = (dPkmeas*Lbox**3)/(2*(Pkdark*Pkmatter)**0.5)
    m=0
    a=int(ptsdict[Lbox])
    c, var = np.polyfit(kk[:a], bias[:a], 0, w=1/dbias[:a], cov=True)
    print(f"no. of points fitted = {a}")
    err = var**0.5
    
    col=colordict[fname.split('_powspec_')[-1].split('.dat')[0]]
    
    #ax.errorbar(kk[:10], bias[:10], dbias[:10], fmt=".", label=f"%.2e_%.2e _c=%.2f res=%.2f"%(massbins[i],massbins[i+1],c,res), color=col)
    #ax.plot(kk,c+kk*m,color=col)
    #ax.set_xscale("log")
    nu = (nubin[i]+nubin[i+1])/2
    x = np.linspace(nubin[i], nubin[i+1], 10)
    #ax.plot(x, c+x*0, color='r')
    ax.errorbar(nu, c, err[0], fmt='_', label=f"err=%.3f"%err, color=col)
    #ax.set_xlim(None, 0.2)
    ax.set_xlabel(r"ν")
    ax.set_ylabel("bias")
    ax.set_title("Lbox = %i" % Lbox)

ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
pl.xticks(nubin)
pl.savefig(f"{dirname}/Planck{Lbox}_bias_nu_cov2.png", bbox_inches='tight')
