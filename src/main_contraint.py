#Programme principal
#Massonneau/Mercier

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import inv

#CONSTANTES
from scipy.constants import Planck as CONST_h
from scipy.constants import Boltzmann as CONST_k
CONST_TCMB = 2.72548 #K

import healpy as hp
import healpy.pixelfunc as px
from sys import exit, argv

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

def poids(a, b, cov):
	a = np.array([a]).T
	b = np.array([b]).T

	bT       = b.T
	aT       = a.T
	invCov   = inv(cov)
	aTInvCov = aT.dot(invCov)
	bTInvCov = bT.dot(invCov)

	#Numerateur
	fnum = bTInvCov.dot(b)
	snum = aTInvCov.dot(b)
	num  = fnum.dot(aTInvCov) - snum.dot(bTInvCov)

	#Denominateur
	tnum = aTInvCov.dot(a)
	den  = tnum.dot(fnum) - snum.dot(snum)

	#Ancienne formule
	#weights=inv(cov).dot(np.array([f_c]).T).dot(inv(np.array([f_c]).dot(inv(cov)).dot((np.array([f_c]).T))))

	#print("test", num/den, np.shape((num/den).T))
	return (num/den).T


def masque(donnees_cartes,chemin,voir_masque=False):
	#création du masque depuis la carte ayant la plus haute fréquence (cad celle ou la galaxie est le plus visible) et application de ce dernier sur toutes les autres cartes

	T_gal=0.005*max(donnees_cartes[:,5]) #a changer pour trouver la bonne valeur
	masque_gal=np.copy(donnees_cartes[:,5])
	ind_masque=np.array(np.where(masque_gal>T_gal)[0])
	print(ind_masque)
	array_masque=np.copy(donnees_cartes)

	for i in range(len(freq)):
		array_masque[:,i][ind_masque]=hp.UNSEEN
		#hp.mollview(array_maps[:,i], title="Cartes avec le meme masque", coord=['G'], unit=r"$K_{CMB}$", norm="hist", min=min(array_maps[:,i]), max=max(array_maps[:,i]))
		#print(array_maps[:,i][ind_masque])

	hp.write_map("data/"+str(nside)+"/HFI_"+str(nside)+"_"+chemin+"_mask_"+str(np.around(len(ind_masque)/len(array_masque[:,0]),3))+"_percent.fits", array_masque[:,5])
	print(np.around(len(ind_masque)/len(array_masque[:,0]),3))
	exit()
	#hp.mollview(array_masque[:,5], title="Cartes avec le meme masque", coord=['G'], unit=r"$K_{CMB}$", norm="hist", min=min(array_masque[:,5]), max=max(array_masque[:,5]))
	#plt.show()
	#exit()

	return array_masque

###

freq       = [100,143,217,353,545,857]
sZ         = 1.0/np.array([-0.24815, -0.35923, 5.152, 0.161098, 0.06918, 0.0380])

sp_emis = np.zeros(len(freq))+1 #le vecteur f_c dans MILCA qui vaut 1 pour CMB_unique
nside   = 2048

chemin="halfmission-2"
array_maps = recup_carte_et_copie(freq,chemin)
array_masque = masque(array_maps,chemin)
cov = covariance(array_masque)
weights = poids(sp_emis, sZ, cov)

Var=weights.T.dot(cov).dot(weights)
print(Var)


CMB_unique=(array_maps.dot(weights)).T


hp.mollview(CMB_unique[0], title="CMB Galactique", coord=['G'], unit=r"$K_{CMB}$", norm="hist", min=min(CMB_unique[0]), max=max(CMB_unique[0]))
#hp.mollview(CMB_unique[0], title="CMB Ecliptique", coord=['G','E'], unit=r"$K_{CMB}$", norm="hist", min=min(CMB_unique[0]), max=max(CMB_unique[0]))

#Ecriture nouvelle carte
hp.write_map("data/"+str(nside)+"/HFI_"+str(nside)+"_"+chemin+"_CM208.fits", CMB_unique[0])



plt.show()

