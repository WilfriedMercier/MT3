#Programme principal
#Massonneau/Mercier

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import inv
import healpy as hp
import healpy.pixelfunc as px

maps_rien   = hp.read_map("data/HFI_rien.fits")
maps_corr   = hp.read_map("data/HFI_corr.fits")
maps_convol = hp.read_map("data/HFI_convol.fits")

hp.mollview(maps_convol-maps_corr, title="CMB Convol - Corr", coord=['G'], unit=r"$K_{CMB}$", norm="hist", min=min(maps_convol-maps_corr), max=max(maps_convol-maps_corr))
hp.mollview(maps_convol-maps_rien, title="CMB Convol - Rien", coord=['G'], unit=r"$K_{CMB}$", norm="hist", min=min(maps_convol-maps_rien), max=max(maps_convol-maps_rien))
hp.mollview(maps_corr-maps_rien, title="CMB Corr - Rien", coord=['G'], unit=r"$K_{CMB}$", norm="hist", min=min(maps_corr-maps_rien), max=max(maps_corr-maps_rien))

plt.show()
