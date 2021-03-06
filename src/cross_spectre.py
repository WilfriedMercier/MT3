#Test spectre en puissance

import matplotlib.pyplot as plt
import healpy as hp
import numpy as np
import healpy.pixelfunc as px
from sys import exit, argv
from numpy.linalg import inv
from healpy.sphtfunc import anafast, pixwin, gauss_beam

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
cmb_maps = [str(nside)+"_"+chemin1, str(nside)+"_"+chemin2]
pourcentage = 0.211 #a completer
map_masque = hp.read_map("data/"+str(nside)+"/HFI_"+str(nside)+"_"+chemin0+"_mask_"+str(pourcentage)+"_percent.fits")
ind_masque = np.array(np.where(map_masque==hp.UNSEEN)[0])


for i in range(len(cmb_maps)):

	#Recupere la carte, et la met dans un tableau (une colonne)
	maps = hp.read_map("data/" + str(nside)+"/HFI_" + str(cmb_maps[i]) + "_M" + str(int(1000*pourcentage)) + ".fits")

	if i==0:
		size = hp.get_map_size(maps)
		array_maps = np.zeros((len(cmb_maps),size))
		array_masque = np.copy(array_maps)
		array_maps[0,:] = maps
	else:
		array_maps[i,:] = maps
	
	array_masque[i,:] = np.copy(array_maps[i,:])
	array_masque[i,:][ind_masque]=hp.UNSEEN
	

##Power spectrum
CL_CMB = anafast(array_masque[0,:], array_masque[1, :], lmax=2500) #4*nside)#,gal_cut=17)
Pix = pixwin(nside)
Pix = Pix[0:2501]
Gauss_sq = gauss_beam(9.65*np.pi/(180*60),lmax=2500)**2

CL_vrai = CL_CMB / (Pix**2*Gauss_sq*(1-pourcentage))
l=np.arange(len(CL_vrai))

#	plt.xscale('log')
f = plt.figure()
ax = f.add_subplot(111)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.tick_params(which='both', direction='in', labelsize=13)

plt.plot(l, CL_vrai*(l*(l+1))/(2*np.pi), '*', label="CMB_Power_Spectrum")

plt.title("Cross-correlation : NSIDE=2048 et masque de " + str(round(100*pourcentage, 1)) + "%")
plt.xlabel("l", size=18)
plt.ylabel(r'$C_l l(l+1)/(2\pi)$', size=18)
plt.grid()
plt.legend(loc='best', fontsize=13)
plt.show()

hp.write_cl("data/spectre/crossSpectreCMB_2048_M" + str(int(1000*pourcentage)) + ".fits", CL_vrai, overwrite=False)
