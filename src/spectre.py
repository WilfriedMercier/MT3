#Test spectre en puissance

import matplotlib.pyplot as plt
import healpy as hp
import numpy as np
import healpy.pixelfunc as px
from sys import exit
from numpy.linalg import inv
from healpy.sphtfunc import anafast

##Recup_carte

cmb_maps=[]

for i in range(len(cmb_maps)):

	#Recupere la carte, et la met dans un tableau (une colonne)
	maps = hp.read_map("data/HFI_"+str(cmb_maps[i])+".fits")
	if i==0:
		size = hp.get_map_size(maps)
		array_maps=np.zeros((len(cmb_maps),size))
		array_maps[0,:]=maps
	else:
		array_maps[i,:]=maps

##Power spectrum
	CL=anafast(array_maps[i,:])
	l=np.arange(len(CL))
	plt.plot(l,CL*(l*(l+1))/(2*np.pi),'*',label=cmb_maps[i])
	plt.xscale('log')
	plt.legend(loc='best')

plt.show()
