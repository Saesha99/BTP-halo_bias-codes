import h5py
import numpy as np
import pandas
import sys
from scipy.spatial import KDTree
import mass_from_peakhieght as mfp

fname = sys.argv[1]

vmin = float(sys.argv[2])
vmax = float(sys.argv[3])
cosmo = fname.split("_l")[0].split("/c")[-1]
a = float(fname.split("ria_")[-1].split(".hdf")[0])
massval = mfp.nutomass(cosmo, a, np.array([vmin, vmax]))

#mmin = int(sys.argv[2])
#mmax = int(sys.argv[3])
#print(mmin,mmax)


def remove_subhalos(pos, mass, radius, box):
     sort_key = mass.argsort()
     
     m_rs = mass[sort_key[::-1]]
     r_rs = radius[sort_key[::-1]] /1000   #convert to kpc
     pos_rs = pos[sort_key[::-1]] * a      #convert to physical units
     
     length = len(m_rs)
     print(length)
     T = KDTree(pos_rs, boxsize=box)
     index = np.array([True]*length)
     for i in range(0,length):
	#print(i)
	#if m_rs[i] == True and r_rs[i] == True and pos_rs[i] == True:
        if index[i] == True:
           idx = T.query_ball_point(pos_rs[i],r=r_rs[i], return_sorted=True)
           index[idx[1:]] = False
           
     mass = m_rs[index]
     position = pos_rs[index]
     radius = r_rs[index]
     print(len(mass))
     	   
     return(position, mass)    
     
     
method = 1
     
with h5py.File(fname, "r") as fin:

     Lbox = float(fname.split("_l")[-1].split("_")[0])
     #print(fin.keys())
     print(Lbox,a)

     pos = np.array(fin["x"])
     print(pos)
     mkey = "Msp-apr-mn" #"M200m_bnd_cat" 
     rkey = "Rsp-apr-mn" #R200m_bnd_cat"
     print(mkey,rkey)
     
     mass = np.array(fin[mkey]) #"M200m_bnd_cat"]) 
     radius = np.array(fin[rkey]) #"R200m_bnd_cat"]) #[
     status = np.array(fin["status_sparta"])
     
     if method == 0:
        print("removing subhalos using kdtree")
        posr, massr = remove_subhalos(pos, mass, radius, Lbox)
        posr = posr/a

     elif method == 1:
        print("removing subhalos using sparta flag")
        length = len(mass)
        index = np.array([True]*length)
        for i in range(0,length):
            if status[i] in (20,21,22,23,24):
                index[i] = False
	    
        massr = mass[index]
        posr = pos[index]
        print(len(posr))
        print(len(massr))
        print(posr)
        
     else:
     	print("subhalos not removed")
     	exit(11)
     	 

     #idx = (mass>10**mmin) & (mass<10**mmax)
     #for i in range(0,1):
     #print(i)
     idx = (massr>massval[0]) & (massr<massval[1])
     x = posr[ : , 0][idx]/Lbox
     y = posr[ : , 1][idx]/Lbox
     z = posr[ : , 2][idx]/Lbox
     massr = massr[idx]

     fout = open(f"{fname}_subhr_{vmin}_{vmax}.dat", "w")
     fout.write(f"{np.sum(idx)}\n")
     for ii in range(x.size):
        fout.write("%.5f %.5f %.5f %.5e\n" % (x[ii], y[ii], z[ii], massr[ii]))
         
         
      

