#Programme principal
#Massonneau/Mercier

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import inv
import healpy as hp
import healpy.pixelfunc as px

freq=[100,143,217,353,545,857]
sp_emis=[1]*len(freq) #le vecteur f_c dans MILCA qui vaut 1 pour CMB_unique

def recup_carte_et_copie(freq_cartes):
	for i in range(len(freq_cartes)):

	maps = hp.read_map("data/HFI_256_"+str(freq_cartes[i])+"_convol.fits")

<<<<<<< HEAD
	#Recupere la carte, et la met dans un tableau (une colonne)
	maps = hp.read_map("data/HFI_256_"+str(freq[i])+"_convol.fits")

=======
>>>>>>> cf81467c5194c075f0357eeeec742ab4977fdac4
	if i==0:
		size = hp.get_map_size(maps)
		array_maps=np.zeros((size,len(freq_cartes)))
		array_maps[:,0]=maps
		print(maps)
	else:
		array_maps[:,i]=maps

	return array_maps

<<<<<<< HEAD
##Covariance (qui est la moyenne de chaque valeur de chaque carte pour chaque pixel)

cov_vraie=np.zeros((len(freq),len(freq)))

for i in range(size):
	cov_vraie+=np.array([array_maps[i]]).T.dot(np.array([array_maps[i]]))

cov_vraie/=size
=======
def covariance(donnees_cartes):
	cov=np.zeros((len(donnees_cartes[0,:]),len(donnees_cartes[0,:])))
	for i in range(size):
		cov+=np.array([array_maps[i]]).T.dot(np.array([array_maps[i]]))
		
	cov/=len(donnees_cartes[:,0])
>>>>>>> cf81467c5194c075f0357eeeec742ab4977fdac4

	return cov

def poids(f_c,cov):
	weights=inv(cov).dot(np.array([f_c]).T).dot(inv(np.array([f_c]).dot(inv(cov)).dot((np.array([f_c]).T))))

	return weights


###

array_maps = recup_carte_et_copie(freq)
cov = covariance(array_maps)
weights = poids(sp_emis,cov)

Var=weights.T.dot(cov).dot(weights)
print(Var)

CMB_unique=(array_maps.dot(weights)).T


hp.mollview(CMB_unique[0], title="CMB Galactique", coord=['G'], unit=r"$K_{CMB}$", norm="hist", min=min(CMB_unique[0]), max=max(CMB_unique[0]))
hp.mollview(CMB_unique[0], title="CMB Ecliptique", coord=['G','E'], unit=r"$K_{CMB}$", norm="hist", min=min(CMB_unique[0]), max=max(CMB_unique[0]))

#Ecriture nouvelle carte
hp.write_map("data/HFI_colbol.fits", CMB_unique[0])








plt.show()
