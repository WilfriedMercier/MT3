import numpy as np
from scipy.stats import chi2_contingency, chisquare
import matplotlib.pyplot as plt
import healpy as hp
import healpy.fitsfunc as ft
from classy import Class

f = plt.figure()
ax = f.add_subplot(111)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.tick_params(which='both', direction='in', labelsize=13)

#Recup spectre
sp = hp.read_cl("data/spectre/crossSpectreCMB_2048_M211.fits")
sp = sp[0:1500]

l = np.arange(len(sp))

#Fonction calcule A a Omega_b, Omega_CDM fixes
def calc_amplitude(y, f):
	return np.sum(f*y)/np.sum(f*f)

#Calcule le spectre theorique a Omega_b, Omega_CDM, A fixes
def calc_spectre_th(params):
	#Variables globales
	global cosmo

	cosmo.set(params)
	cosmo.compute()
	cl = cosmo.lensed_cl(len(sp)-1)
	return cl['tt']

##MAIN##
#Classe class
cosmo = Class()

#parametres (Omega*h^2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)
params = {'output': 'tCl lCl',
         'l_max_scalars': 2000,
         'lensing': 'yes',
         'A_s': 2.3e-9*(2.2758**2),
         'n_s': 0.9624, 
         'h': 0.6711,
         'omega_b': 0.02230,
         'omega_cdm': 0.1188}

#list1 = np.arange(0.1, 1, 0.1)
list2 = np.geomspace(0.01, 2, 100)*1e-8
plist = np.copy(list2)*0-1

for i, j in zip(list2, range(len(plist))):
	print(j)
	params['A_s'] = i

	cl       = calc_spectre_th(params)
	plist[j] = np.max(cl*l*(l+1)/(2*np.pi))


#plt.plot(list2, plist, "r.")
#ax.set_xscale('log')
#plt.show()

diff = np.abs(plist-np.max(sp*l*(l+1)/(2*np.pi)))
closest = np.where(diff == np.min(diff))[0]
print("clos", closest)

A    = list2[closest]
print("A", A)
params['A_s'] = A[0]
cl       = calc_spectre_th(params)


plt.plot(cl*l*(l+1)/(2*np.pi), "r", label="Best fit : A = " + np.format_float_scientific(A, exp_digits=2, precision=4))
plt.plot(sp*l*(l+1)/(2*np.pi), "b.", label="Donnees : A = " + np.format_float_scientific(np.max(sp*l*(l+1)/(2*np.pi)), exp_digits=2, precision=4))

plt.xlabel("Multipole " + r"$l$", size=18)
plt.ylabel(r'$C_l l (l+1) / (2 \pi)$', size=18)
plt.legend(loc='best', fontsize=10)
plt.grid()
plt.show()
