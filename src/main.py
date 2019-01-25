#Programme principal
#Massonneau/Mercier

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import inv
import healpy as hp
import healpy.pixelfunc as px

##Recup_carte

freq=[100,143,217,353,545,857]

for i in range(len(freq)):

	#Recupere la carte, le nb de pixels et le NSIDE min acceptable
	maps = hp.read_map("/home/massonne/p1_massonne_mercier/data/HFI_256_"+str(freq[i])+".fits")
	if i==0:
		size = hp.get_map_size(maps)
		array_maps=np.zeros((size,len(freq)))
		array_maps[:,0]=maps
		print(maps)
	else:
		array_maps[:,i]=maps


	#
	#minval = px.get_min_valid_nside(size)
	#NSIDE  = px.get_nside(map)

	#Affichage
	#print("Taille : ", size, "\nNSDIDE minimum : ", minval, "\nCurrent : ", NSIDE)
	#hp.mollview(map, title="Reduced, freq = "+str(i) + "GHz", coord=['G', 'E'], unit="mK", norm="hist", min=-1, max=1)#, xsize=2000)

print(array_maps)



##Covariance

cov_vraie=np.zeros((len(freq),len(freq)))

for i in range(size):
	cov_vraie+=np.array([array_maps[i]]).T.dot(np.array([array_maps[i]]))
		
cov_vraie/=size



##Weights

sp_emis=[1]*len(freq) #le vecteur f_c dans MILCA qui vaut 1 pour CMB_unique

weights=inv(cov_vraie).dot(np.array([sp_emis]).T).dot(inv(np.array([sp_emis]).dot(inv(cov_vraie)).dot((np.array([sp_emis]).T))))
print(np.shape(weights))

Var=weights.T.dot(cov_vraie).dot(weights)
print(Var)

CMB_unique=(array_maps.dot(weights)).T


hp.mollview(CMB_unique[0], title="Premier CMB Elliptique", norm="hist")










plt.show()
