#Programme principal
#Massonneau/Mercier

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import inv
import healpy as hp
import healpy.pixelfunc as px
from sys import exit, argv

freq = [100,143,217,353,545,857]
sp_emis = [1]*len(freq) #le vecteur f_c dans MILCA qui vaut 1 pour CMB_unique
nside = 2048

'''
taille = len(argv)
if (taille==1):
	chemin = "full"
	print("coucou")
else:
	if (int(argv[1])==1):
		chemin = "halfmission-1"
	else:
		chemin = "halfmission-2"
'''
def recup_carte_et_copie(freq_cartes,chemin):
	for i in range(len(freq_cartes)):

		maps = hp.read_map("data/"+str(nside)+"/HFI_"+str(nside)+"_"+str(freq_cartes[i])+"_"+chemin+".fits")

		if i==0:
			size = hp.get_map_size(maps)
			array_maps=np.zeros((size,len(freq_cartes)))
			array_maps[:,0]=maps
			print(maps)
		else:
			array_maps[:,i]=maps

	return array_maps

def covariance(donnees_cartes):
	cov=np.zeros((len(donnees_cartes[0,:]),len(donnees_cartes[0,:])))
	for i in range(len(donnees_cartes[:,0])):
		if donnees_cartes[i,0]!=hp.UNSEEN:
			cov+=np.array([donnees_cartes[i]]).T.dot(np.array([donnees_cartes[i]]))
		
	cov/=len(donnees_cartes[:,0])

	return cov

def poids(f_c,cov):
	weights=inv(cov).dot(np.array([f_c]).T).dot(inv(np.array([f_c]).dot(inv(cov)).dot((np.array([f_c]).T))))

	return weights


def masque(donnees_cartes,chemin,voir_masque=False):
	#création du masque depuis la carte ayant la plus haute fréquence (cad celle ou la galaxie est le plus visible) et application de ce dernier sur toutes les autres cartes

	T_gal=0.0005*max(donnees_cartes[:,5]) #a changer pour trouver la bonne valeur
	masque_gal=np.copy(donnees_cartes[:,5])
	ind_masque=np.array(np.where(masque_gal>T_gal)[0])
	print(ind_masque)
	array_masque=np.copy(donnees_cartes)

	for i in range(len(freq)):
		array_masque[:,i][ind_masque]=hp.UNSEEN
		#hp.mollview(array_maps[:,i], title="Cartes avec le meme masque", coord=['G'], unit=r"$K_{CMB}$", norm="hist", min=min(array_maps[:,i]), max=max(array_maps[:,i]))
		#print(array_maps[:,i][ind_masque])

	hp.write_map("data/"+str(nside)+"/HFI_"+str(nside)+"_"+chemin+"_mask_"+str(np.around(len(ind_masque)/len(array_masque[:,0]),3))+"_percent.fits", array_masque[:,5])
	#hp.mollview(array_masque[:,5], title="Cartes avec le meme masque", coord=['G'], unit=r"$K_{CMB}$", norm="hist", min=min(array_masque[:,5]), max=max(array_masque[:,5]))
	#plt.show()
	#exit()

	return array_masque

###

chemin="halfmission-2"
array_maps = recup_carte_et_copie(freq,chemin)
array_masque = masque(array_maps,chemin)
cov = covariance(array_masque)
weights = poids(sp_emis,cov)

Var=weights.T.dot(cov).dot(weights)
print(Var)


CMB_unique=(array_maps.dot(weights)).T


hp.mollview(CMB_unique[0], title="CMB Galactique", coord=['G'], unit=r"$K_{CMB}$", norm="hist", min=min(CMB_unique[0]), max=max(CMB_unique[0]))
#hp.mollview(CMB_unique[0], title="CMB Ecliptique", coord=['G','E'], unit=r"$K_{CMB}$", norm="hist", min=min(CMB_unique[0]), max=max(CMB_unique[0]))

#Ecriture nouvelle carte
hp.write_map("data/"+str(nside)+"/HFI_"+str(nside)+"_"+chemin+".fits", CMB_unique[0])



plt.show()

