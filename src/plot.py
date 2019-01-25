import matplotlib.pyplot as plt
import healpy as hp


freq  = [100] #,143,217,353,545,857]
for i in freq:
	map    = hp.read_map("data/HFI_256_"+str(i)+".fits")
	hp.mollview(map, title="galactique", coord=['G'], unit="mK", norm="hist", xsize=400)
	hp.mollview(map, title="galactique", coord=['G', 'E'], unit="mK", norm="hist", xsize=400)
plt.show()
