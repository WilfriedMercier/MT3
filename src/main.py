#Programme principal
#Massonneau/Mercier

import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import healpy.pixelfunc as px

NSIDE = 20

freq=[100]#,143,217,353,545,857]

for i in freq:

	#Recupere la carte, le nb de pixels et le NSIDE min acceptable
	map    = hp.read_map("/home/vbonjean/PLANCK/HFI_SkyMap_"+str(i)+"_2048_R3.01_full.fits")
	size   = hp.get_map_size(map)
	minval = px.get_min_valid_nside(size)
	NSIDE  = px.get_nside(map)
	print("Taille : ", size, "\nSDIDE minimum : ", minval, "\nCurrent : ", NSIDE)

	#Degradation de carte
	newmap = px.ud_grade(map, 512)

	#Affichage
	hp.mollview(map, coord=['G','E'], title='Full, freq = ' + str(i) +' GHz', unit='mK', norm='hist', min=-1, max=1, xsize=2000)
	hp.mollview(newmap, coord=['G','E'], title='Version degradee, freq = ' + str(i) +' GHz', unit='mK', norm='hist', min=-1, max=1, xsize=2000)

	print("V. deg. -> NSIDE = ", px.get_nside(newmap))

	#Ecriture nouvelle carte
#	hp.write_map("testout.fits", rmap)

plt.show()
