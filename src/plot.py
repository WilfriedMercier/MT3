import matplotlib.pyplot as plt
import healpy as hp


freq  = [100] #,143,217,353,545,857]
for i in freq:
<<<<<<< HEAD
	maps = hp.read_map("/home/massonne/p1_massonne_mercier/data/HFI_256_"+str(i)+".fits")
	hp.mollview(maps, title="galactique", coord=['G'], unit=r"$K_{CMB}$", norm="hist")#, min=-1, max=1, xsize=200)
=======
	map    = hp.read_map("data/HFI_256_"+str(i)+".fits")
	hp.mollview(map, title="galactique", coord=['G'], unit="mK", norm="hist", xsize=400)
	hp.mollview(map, title="galactique", coord=['G', 'E'], unit="mK", norm="hist", xsize=400)
>>>>>>> aab4853e1915d111ed447efaff808a691f0a0b61
plt.show()
