import matplotlib.pyplot as plt
import healpy as hp


freq  = [100,143,217,353,545,857]
for i in freq:
	maps    = hp.read_map("data/HFI_256_"+str(i)+".fits")
	hp.mollview(maps, title="galactique", coord=['G'], unit="K", norm="hist", xsize=400)
	hp.mollview(maps, title="ecliptique", coord=['G', 'E'], unit="K", norm="hist", xsize=400)

plt.show()
