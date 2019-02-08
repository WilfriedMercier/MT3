#Reduit les cartes et les ecrit dans data
#Massonneau/Mercier

import matplotlib.pyplot as plt
import healpy as hp
import numpy as np
import healpy.pixelfunc as px
from sys import exit, argv

taille = len(argv)
if (taille==1):
	chemin = "full"
	print("coucou")
else:
	if (int(argv[1])==1):
		chemin = "halfmission-1"
	else:
		chemin = "halfmission-2"

nside           = 2048
freq            = [100,143,217,353,545,857]
FWHM_gauss_conv = np.array([9.65,7.25,4.99,4.82,4.68,4.33])*(np.pi/(180*60))

#Facteurs correction en unite K_CMB
ToCMB = 4*[1.0] + [1.0/58.04, 1.0/2.27]

#Facteurs correction en unite MJy/sr
FrCMB = [244.1, 371.74, 483.69, 287.45, 1.0, 1.0]

for i, corr, conv in zip(freq, ToCMB, FWHM_gauss_conv):

	#Recupere la carte, le nb de pixels et le NSIDE min acceptable
	map    = hp.read_map("/home/vbonjean/PLANCK/HFI_SkyMap_" + str(i) + "_2048_R3.01_" + chemin + ".fits")
	print("/home/vbonjean/PLANCK/HFI_SkyMap_" + str(i) + "_2048_R3.01_" + chemin + ".fits")

	#Nettoyage professionnel
	map[np.array(np.where(abs(map)==abs(hp.UNSEEN))[0])]=0.

	hp.mollview(map, title="CMB Galactique " + str(i) + "GHz", coord=['G'], unit=r"$K_{CMB}$", norm="hist", min=min(map), max=max(map))
	size   = hp.get_map_size(map)
	minval = px.get_min_valid_nside(size)
	NSIDE  = px.get_nside(map)
	print("Taille : ", size, "\nNSDIDE minimum : ", minval, "\nCurrent : ", NSIDE)

	#Convertit en unite K_CMB
	map *= corr

	#Convolution pour degrader la resolution
	FWHM=np.sqrt(conv**2-FWHM_gauss_conv[-1]**2)
	if i!=freq[-1]:
		map=hp.sphtfunc.smoothing(map,FWHM)

	#Degradation de carte
	newmap = px.ud_grade(map, nside)

	#Ecriture nouvelle carte
	hp.write_map("data/HFI_"+str(nside)+"_" + str(i) + "_" + chemin + ".fits", newmap, overwrite=True)
	print("Wrote "+str(i))

<<<<<<< HEAD
	exit()



aze
q
dsd
s
dfs
df
sd
fsd
f
=======
plt.show()
>>>>>>> c9842b4fca2b314b1bb3b1f43ac120f592b90c7d
