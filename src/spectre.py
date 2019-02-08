#Test spectre en puissance

import matplotlib.pyplot as plt
import healpy as hp
import numpy as np
import healpy.pixelfunc as px
from sys import exit, argv
from numpy.linalg import inv
from healpy.sphtfunc import anafast, pixwin, gauss_beam


<<<<<<< HEAD
cmb_maps=['2048_convol']
map_masque = hp.read_map("data/HFI_2048_convol_mask_25_percent.fits")
ind_masque=np.array(np.where(map_masque==hp.UNSEEN)[0])
=======
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
chemin0 = "full"
chemin1 = "halfmission-1"
chemin2 = "halfmission-2"
cmb_maps = [str(nside)+"_"+chemin0,str(nside)+"_"+chemin1,str(nside)+"_"+chemin2]
pourcentage = 0.2#a completer
map_masque = hp.read_map("data/"+str(nside)+"/HFI_"+str(nside)+"_"+chemin0+"_mask_"+str(pourcentage)+"_percent.fits")
ind_masque = np.array(np.where(map_masque==hp.UNSEEN)[0])

>>>>>>> cf81467c5194c075f0357eeeec742ab4977fdac4

for i in range(len(cmb_maps)):

	#Recupere la carte, et la met dans un tableau (une colonne)
	maps = hp.read_map("data/"+str(nside)+"/HFI_"+str(cmb_maps[i])+".fits")

	if i==0:
		size = hp.get_map_size(maps)
		array_maps = np.zeros((len(cmb_maps),size))
		array_masque = np.copy(array_maps)
		array_maps[0,:] = maps
	else:
		array_maps[i,:] = maps
	
	array_masque[i,:] = np.copy(array_maps[i,:])
	array_masque[i,:][ind_masque]=hp.UNSEEN
	
	##hp.mollview(array_masque[i,:], title="Cartes avec masque", coord=['G'], unit=r"$K_{CMB}$", norm="hist", min=min(array_masque[i,:]), max=max(array_masque[i,:]))
	##plt.show()
	##exit()


##Bruit blanc
array_white_noise = (array_masque[1,:] - array_masque[2,:])/2


##Power spectrum
<<<<<<< HEAD
	CL=anafast(array_masque[i,:], lmax=8192)#,gal_cut=17)

	CL /= pixwin(2048)**2 * (1-len(ind_masque)/len(array_maps[i,:])) * gauss_beam(4.33*np.pi/(180*60), lmax=8192)

	l=np.arange(len(CL))
	plt.plot(l,CL*(l*(l+1))/(2*np.pi),'*',label=cmb_maps[i])
=======
CL_CMB = anafast(array_masque[0,:],lmax=2500)#4*nside)#,gal_cut=17)
CL_WN = anafast(array_white_noise,lmax=2500)#4*nside)
CL_cross = anafast(array_masque[0,:],map2=array_white_noise,lmax=2500)#4*nside)
Pix = pixwin(nside)
Pix = Pix[0:2501]
Gauss_sq = gauss_beam(4.33*np.pi/(180*60),lmax=2500)**2
CL_vrai = CL_CMB / (Pix**2*Gauss_sq*(1-pourcentage)) - CL_WN / (Pix**2) - 2*CL_cross / (Pix)
l=np.arange(len(CL_vrai))
plt.plot(l,CL_vrai*(l*(l+1))/(2*np.pi),'*',label="CMB_Power_Spectrum")
>>>>>>> cf81467c5194c075f0357eeeec742ab4977fdac4
#	plt.xscale('log')
plt.legend(loc='best')
#	hp.write_cl("data/spectre_vrai2_2048.fits",CL_vrai,overwrite=True)
plt.show()

