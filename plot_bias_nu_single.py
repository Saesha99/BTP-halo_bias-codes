import numpy as np
import pylab as pl
import sys
sys.path.append("../../aum/install-3.9/lib/python3.9/site-packages/")
import cosmology as cc
import glob
import math
import mass_from_peakhieght as mfp

f=['.','_','x']
#colordict={'0_1':'cyan','1_2':'olive', '2_3':'green', '3_4':'red', '4_5':'blue', '5_6':'pink', '6_7':'brown'} 
colordict={'0.5_0.75':'tab:pink','0.75_1.0':'tab:cyan', '1.0_1.25':'tab:olive', '1.25_1.5':'tab:green', '1.5_1.75':'tab:red', '1.75_2.0':'tab:blue', '2.0_2.25':'tab:purple'}
ptsdict={500:'10', 250:'8', 125:'4'}

ax = pl.subplot(111)
for j in range(1,4):
	dirname = sys.argv[j]
	def_type=dirname.split("1024/")[-1].split("/a")[0]
	print(def_type)
	sf = float(dirname.split("/a_")[-1])
	flist = np.sort(np.array(glob.glob(f"{dirname}/moria*powspec_subhr*")))
	cosmo = dirname.split("Data/c")[1].split("_")[0]
	if cosmo=="pla":        
	    aa = cc.cosmology(0.315,0.0,-1.0,0.0,0.02222/0.6731**2,0.6731,2.726,0.829,0.9655,np.log10(8.0),1.0)
	elif cosmo=="bol":
	    aa = cc.cosmology(0.27,0.0,-1.0,0.0,0.0469,0.7,2.726,0.82,0.95,np.log10(8.0),1.0)
	else:
	    print(f"Cosmology {cosmo} not supported yet.")
	    exit(11)


	Lbox = float(dirname.split("_l")[-1].split("_")[0])
	
	for i, fname in enumerate(flist):
	    
	    fin = open(fname.replace("_powspec", ""))
	    Npart = int(fin.readline())
	    fin.close()
	    print(Lbox, Npart)
	    
	    vmin=float(fname.split("subhr_")[-1].split("_")[0])
	    vmax=float(fname.split("subhr_")[-1].split("_")[-1].split(".da")[0])
	    
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
	    	if kk[jj]>0.22:
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
	    
	    #col=colordict[fname.split('_powspec_')[-1].split('.dat')[0]]
	    col=colordict[fname.split('_powspec_subhr_')[-1].split('.dat')[0]]
	    
	    #ax.errorbar(kk[:10], bias[:10], dbias[:10], fmt=".", label=f"%.2e_%.2e _c=%.2f res=%.2f"%(massbins[i],massbins[i+1],c,res), color=col)
	    #ax.plot(kk,c+kk*m,color=col)
	    #ax.set_xscale("log")
	    #nu = (nubin[i]+nubin[i+1])/2
	    ax.errorbar((vmin+vmax)/2, c, err[0], fmt=f[j-1], label=f"{Lbox}_ b=%.2f err=%.3f"%(c,err), color=col)
	    #ax.set_xlim(None, 0.2)
	    ax.set_xlabel(r"ν")
	    ax.set_ylabel("bias")
	    ax.set_title(f"a={sf}")

ax.legend(loc='upper left', bbox_to_anchor=(1, 1))
pl.xticks([0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25])
pl.savefig(f"../Plots/Planck_{def_type}_{sf}_subhr_bias_nu_cov_single_19-4.png", bbox_inches='tight')
