#Programme principal
#Massonneau/Mercier

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import inv
import healpy as hp
import healpy.pixelfunc as px
from sys import exit

##Recup_carte

freq=[100,143,217,353,545,857]

for i in range(len(freq)):

	#Recupere la carte, et la met dans un tableau (une colonne)
	maps = hp.read_map("data/HFI_256_"+str(freq[i])+"_convol.fits")
	if i==0:
		size = hp.get_map_size(maps)
		array_maps=np.zeros((size,len(freq)))
		array_maps[:,0]=maps
	else:
		array_maps[:,i]=maps


##Masque unique galaxie
#création du masque depuis la carte ayant la plus haute fréquence (cad celle ou la galaxie est le plus visible) et application de ce dernier sur toutes les autres cartes

T_gal=20.5 #a changer pour trouver la bonne valeur
masque_gal=np.copy(array_maps[:,5])
ind_masque=np.array(np.where(masque_gal>T_gal)[0])
#masque_gal[ind_masque]=min(array_maps[:,5])-1.
#hp.mollview(masque_gal, title="Masque galactique", coord=['G'], unit=r"$K_{CMB}$", norm="hist", min=min(masque_gal), max=max(masque_gal))

for i in range(len(freq)):
	array_maps[:,i][ind_masque]=min(array_maps[:,i])-0.05
	#hp.mollview(array_maps[:,i], title="Cartes avec le meme masque", coord=['G'], unit=r"$K_{CMB}$", norm="hist", min=min(array_maps[:,i]), max=max(array_maps[:,i]))
	#print(array_maps[:,i][ind_masque])



#for i in ind_masque:
#	np.delete(array_maps,(i),axis=0)

##Covariance (qui est la moyenne de chaque valeur de chaque carte pour chaque pixel)

cov_vraie=np.zeros((len(freq),len(freq)))

test=min(array_maps[:,0])
for i in range(size):
	if array_maps[i,0]!=test:
		cov_vraie+=np.array([array_maps[i]]).T.dot(np.array([array_maps[i]]))
		
cov_vraie/=size


##Weights

sp_emis=[1]*len(freq) #le vecteur f_c dans MILCA qui vaut 1 pour CMB_unique

weights=inv(cov_vraie).dot(np.array([sp_emis]).T).dot(inv(np.array([sp_emis]).dot(inv(cov_vraie)).dot((np.array([sp_emis]).T))))
print(np.shape(weights))

Var=weights.T.dot(cov_vraie).dot(weights)
print(Var)


CMB_unique=(array_maps.dot(weights)).T


hp.mollview(CMB_unique[0], title="CMB Galactique", coord=['G'], unit=r"$K_{CMB}$", norm="hist", min=min(CMB_unique[0]), max=max(CMB_unique[0]))
#hp.mollview(CMB_unique[0], title="CMB Ecliptique", coord=['G','E'], unit=r"$K_{CMB}$", norm="hist", min=min(CMB_unique[0]), max=max(CMB_unique[0]))

#Ecriture nouvelle carte
#hp.write_map("data/HFI_rien.fits", CMB_unique[0])








plt.show()

