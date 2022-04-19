import h5py
import numpy as np
import pandas
import sys
import mass_from_peakhieght as mfp
from scipy.spatial import KDTree




fname = sys.argv[1]
counttype = int(sys.argv[2])
cosmo = fname.split("Data/c")[1].split("_")[0]
if cosmo=="pla":
   Om0 = 0.315
   h = 0.6731
elif cosmo=="bol":
   Om0 = 0.27
   h = 0.7
else:
   print(f"Cosmology {cosmo} not supported yet.")
   exit(11)
print(cosmo)
Lbox = float(fname.split("_l")[-1].split("_")[0])
a = float(fname.split("ria_")[-1].split(".hdf")[0])
print(Lbox,a)
     

def remove_subhalos(pos, mass, radius):
     sort_key = mass.argsort()
     
     m_rs = mass[sort_key[::-1]]
     r_rs = radius[sort_key[::-1]] /1000
     pos_rs = pos[sort_key[::-1]] * a 
     
     length = len(m_rs)
     print(length)
     T = KDTree(pos_rs)
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
     
     
     

mp = Om0 * Lbox**3 * 2.77 * 10**11 / 1024**3
#mp = round(Om0 * Lbox**3 * 2.77 / 1024**3, 2) * 10**11 
print("mp = %.2e" % mp)
mres = mp*100
print("mass resolution = %.2e" %mres)
#print(counttype)

with h5py.File(fname, "r") as fin:

     #print(fin.keys())
     pos = np.array(fin["x"])
     
     mkey = "Msp-apr-mn" #"M200m_bnd_cat" 
     rkey = "Rsp-apr-mn" #R200m_bnd_cat"
     print(mkey,rkey)
     
     mass = np.array(fin[mkey]) #"M200m_bnd_cat"]) 
     radius = np.array(fin[rkey]) #"R200m_bnd_cat"]) #[
     status = np.array(fin["status_sparta"])
     
     print(len(mass))
     #pos, mass = remove_subhalos(pos, mass, radius)
     #pos = pos/a

     print("removing subhalos using sparta flag")
     length = len(mass)
     index = np.array([True]*length)
     for i in range(0,length):
          if status[i] in (20,21,22,23,24):
                index[i] = False
	    
     mass = mass[index]
     pos = pos[index]
     print(len(pos))
     print(len(mass))   
     
     
     if counttype == 1:
     	print("counting in mass resolution bins")
     	idx = (mass>mres) #& (mass<10**mmax)
     	print("halos above mres: %i" % np.sum(idx))
     
     	idx = (mass>mres) & (mass<mres*10)
     	print("0_1: %i" %np.sum(idx))
     
     	idx = (mass>mres) & (mass<mres*10**0.5)
     	print("0_0.5: %i" %np.sum(idx))
     
     	idx = (mass>mres*10**0.5) & (mass<mres*10**1)
     	print("0.5_1: %i" %np.sum(idx))
     
     	idx = (mass>mres*10**1) & (mass<mres*10**2)
     	print("1_2: %i" %np.sum(idx))
     
     	idx = (mass>mres*10**2) & (mass<mres*10**3)
     	print("2_3: %i" %np.sum(idx))
     
     elif counttype == 2:
     	print("counting in peakhieght bins")
     	mval= mfp.nutomass(cosmo, a, np.array([0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25]))
     	print(mval)
     	print(len(mass))
     	for i in range(0,mval.size-1):
     	   idx = (mass>mval[i]) & (mass<mval[i+1])
     	   print(np.sum(idx))
     	
'''idx = (mass>mval[0]) & (mass<mval[1])
	     	print(np.sum(idx))
	     	
	     	idx = (mass>mval[1]) & (mass<mval[2])
	     	print(np.sum(idx))
	     	
	     	idx = (mass>mval[2]) & (mass<mval[3])
	     	print(np.sum(idx))
	     	
	     	idx = (mass>mval[3]) & (mass<mval[4])
	     	print(np.sum(idx))
	     	
	     	idx = (mass>mval[4]) & (mass<mval[5])
	     	print(np.sum(idx))
	     	
	     	idx = (mass>mval[5]) & (mass<mval[6])
	     	print(np.sum(idx))'''
     	
     	
     	

    

     
