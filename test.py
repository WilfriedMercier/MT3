import matplotlib.pyplot as plt
import healpy as hp

NSIDE = 258



freq=[100,143]#,217,353,545,857]

for i in freq:

	map = hp.read_map("/home/vbonjean/PLANCK/HFI_SkyMap_"+str(i)+"_2048_R3.01_full.fits")
	hp.mollview(map, coord=['G','E'], title='Histogram equalized Ecliptic', unit='mK', norm='hist', min=-1,max=1, xsize=2000)
	hp.graticule()

plt.show()
