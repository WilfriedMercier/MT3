#Programme principal
#Massonneau/Mercier

import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import healpy.pixelfunc as px

freq=[100,143,217,353,545,857]

for i in freq:

	#Recupere la carte, le nb de pixels et le NSIDE min acceptable
	map    = hp.read_map("data/HFI_256_"+str(i)+".fits")
	size   = hp.get_map_size(map)
	minval = px.get_min_valid_nside(size)
	NSIDE  = px.get_nside(map)

	#Affichage
	print("Taille : ", size, "\nNSDIDE minimum : ", minval, "\nCurrent : ", NSIDE)
	hp.mollview(map, title="Reduced, freq = "+str(i) + "GHz", coord=['G', 'E'], unit="mK", norm="hist", min=-1, max=1, xsize=2000)


plt.show()
