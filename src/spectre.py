#Test spectre en puissance

import matplotlib.pyplot as plt
import healpy as hp
import numpy as np
import healpy.pixelfunc as px
from sys import exit
from numpy.linalg import inv
from healpy.sphtfunc import anafast

##Recup_carte

cmb_maps=['256Norm1P']
map_masque = hp.read_map("data/HFI_2048_convol_mask_25_percent.fits")
ind_masque=np.array(np.where(map_masque==hp.UNSEEN)[0])

for i in range(len(cmb_maps)):

	#Recupere la carte, et la met dans un tableau (une colonne)
	maps = hp.read_map("data/HFI_"+str(cmb_maps[i])+".fits")
	if i==0:
		size = hp.get_map_size(maps)
		array_maps=np.zeros((len(cmb_maps),size))
		array_maps[0,:]=maps
	else:
		array_maps[i,:]=maps
	
	array_masque=np.copy(array_maps)
	array_masque[i,:][ind_masque]=hp.UNSEEN
	
	##hp.mollview(array_masque[i,:], title="Cartes avec masque", coord=['G'], unit=r"$K_{CMB}$", norm="hist", min=min(array_masque[i,:]), max=max(array_masque[i,:]))
	##plt.show()
	##exit()

##Power spectrum
	CL=anafast(array_masque[i,:])#,gal_cut=17)
	l=np.arange(len(CL))
	plt.plot(l,CL*(l*(l+1))/(2*np.pi),'*',label=cmb_maps[i])
#	plt.xscale('log')
	plt.legend(loc='best')

plt.show()

