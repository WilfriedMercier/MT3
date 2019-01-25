#Reduit les cartes et les ecrit dans data
#Massonneau/Mercier

import matplotlib.pyplot as plt
import healpy as hp
import numpy as np
import healpy.pixelfunc as px

nside = 256
freq  = [100,143,217,353,545,857]

#Facteurs correction en unite K_CMB
ToCMB = 4*[1.0] + [1.0/58.04, 1.0/2.27]

#Facteurs correction en unite MJy/sr
FrCMB = [244.1, 371.74, 483.69, 287.45, 1.0, 1.0]

for i, corr in zip(freq, ToCMB):

	#Recupere la carte, le nb de pixels et le NSIDE min acceptable
	map    = hp.read_map("/home/vbonjean/PLANCK/HFI_SkyMap_"+str(i)+"_2048_R3.01_full.fits")
	size   = hp.get_map_size(map)
	minval = px.get_min_valid_nside(size)
	NSIDE  = px.get_nside(map)
	print("Taille : ", size, "\nNSDIDE minimum : ", minval, "\nCurrent : ", NSIDE)

	#Convertit en unite K_CMB
	map *= corr

	#Degradation de carte
	newmap = px.ud_grade(map, nside)

	#Ecriture nouvelle carte
	hp.write_map("/home/massonne/p1_massonne_mercier/data/HFI_"+str(nside)+"_"+str(i)+".fits", newmap)
	print("Wrote "+str(i))

