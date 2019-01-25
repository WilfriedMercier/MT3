#Reduit les cartes et les ecrit dans data
#Massonneau/Mercier

import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import healpy.pixelfunc as px

nside = 256
freq  = [100,143,217,353,545,857]

for i in freq:

	#Recupere la carte, le nb de pixels et le NSIDE min acceptable
	map    = hp.read_map("/home/vbonjean/PLANCK/HFI_SkyMap_"+str(i)+"_2048_R3.01_full.fits")
	size   = hp.get_map_size(map)
	minval = px.get_min_valid_nside(size)
	NSIDE  = px.get_nside(map)
	print("Taille : ", size, "\nNSDIDE minimum : ", minval, "\nCurrent : ", NSIDE)

	#Degradation de carte
	newmap = px.ud_grade(map, nside)

	#Ecriture nouvelle carte
	hp.write_map("/home/massonne/p1_massonne_mercier/data/HFI_"+str(nside)+"_"+str(i)+".fits", newmap)
	print("Wrote "+str(i))

