#Elimine la galaxie en la renormalisant et cherche la meilleur normalisation par minimisation de la variance
#Massonneau/Mercier

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import inv
from scipy.optimize import minimize, differential_evolution
import healpy as hp
import healpy.pixelfunc as px
from sys import exit

PI = np.pi

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------
# Fonction a minimiser -
#-----------------------
def Variance(minprop):
	#################################################################################
	# IN:										#
	#	limit	:	fraction de la carte dans laquelle on cherche la galaxy	#
	#	minprop	:	ecart relatif au maximum de la carte en temperature	#
	# OUT:										#
	#	Var	:	Variance des cartes					#
	#################################################################################

	#Tableaux et variables globales
	global normArray, weights


	#Pas encore utilise
	limit = 1.0

	#Normalisation de la zone galactique avec les parametres choisis
	normalizeGalaxy(limit, minprop)
        #Calcule de la matrice de covariance
	cov_vraie = computeCov(normArray)
        #Calcule des poids
	weights = computeWeights(cov_vraie)
        #Calcul de la variance
	Var = weights.T.dot(cov_vraie).dot(weights)

	print("VAR = ", Var, minprop)
	return Var


#--------------------------------------
# Manipulation des pixels galactiques -
#--------------------------------------
def normalizeGalaxy(limit, minprop):
	#################################################################################################################################
	# IN:															      	#
	# 	limit		:	Proportion de latitude dans laquelle on veut normaliser les pixels (0=rien, 1=toute la carte)	#
	# 	minprop		:	Proportion du maximum sur la bande centrale que l'on souhaite normaliser                      	#
	#################################################################################################################################

	#Tableaux et variables globales
	global array_maps, normArray, weights, rgSz, sz

	# Declaration tableaux utiles pour la suite
	maxVal = np.zeros(3)		#Valeur max de la bande centrale pour chaque carte
	minVal = np.copy(maxVal)        #Idem mais min
	meanIn = np.copy(maxVal)        #Valeur moyenne de la zone galactique pour chaque carte
	meanOu = np.copy(maxVal)        #Idem mais pour la zone extra-galactique
	rgMax  = np.arange(3, 6)
	rgLs   = np.arange(0, 3)

	#Calcul du max de la bande centrale pour chaque carte
	mask   = theta==0
	for i, j in zip(rgLs, rgMax):
		tmp       = array_maps[:, j]
		tmp       = tmp[mask]
		maxVal[i] = np.max(tmp)
		minVal[i] = np.min(array_maps[:, j])
	#Pondere par le poids voulu
	maxVal -= minprop*(maxVal-minVal)

	#Calcule les valeurs moyennes dans chaque zone (galactique et hors galactique)
	for i, j in zip(rgLs, rgMax):
		tmp       = array_maps[:,j]
		inside    = tmp[tmp>=maxVal[i]]
		outside   = tmp[tmp<maxVal[i]]
		meanIn[i] = np.sum(inside)/np.size(inside)
		meanOu[i] = np.sum(outside)/np.size(outside)

	normArray = array_maps
	#Normalisation de la bande centrale
	for i, j in zip(rgLs, rgMax):
		normArray[:,j] = np.where(array_maps[:,j]>=maxVal[i], array_maps[:,j]/(1000*meanIn[i])*meanOu[i], array_maps[:,j])


#----------------------------------------------------------------------------------------------------------------
# Covariance (moyenne  sur l'ensemble des pixels des matrices de covariance sur les frequences de chaque pixel) -
#----------------------------------------------------------------------------------------------------------------
def computeCov(array_in):
	########################################################################
	# IN:								       #
	# 	array_in	:	matrice contenant les cartes modifiees #
	# OUT:								       #
	#	cov_vraie	:	matrice de covariance		       #
	########################################################################

	#Variables globales
	global sz, size

	#Initialisation
	cov_vraie = np.zeros((sz, sz))

	#Calcule la matrice de covariance de chaque pixel et l'ajoute aux precedentes
	for i in range(size):
		tmp        = np.array([array_in[i]])
		cov_vraie += tmp.T.dot(tmp)
	cov_vraie /= size
	return cov_vraie


#--------------------
# Calcule des poids -
#--------------------
def computeWeights(cov_vraie):
	#######################################################
	# IN :						      #
	# 	cov_vraie	:	matrice de covariance #
	# OUT:						      #
	#	weights		:	vecteur des poids     #
	#######################################################

	#Variables globales
	global sz

	#le vecteur f_c dans MILCA qui vaut 1 pour CMB_unique
	sp_emis  = np.array([np.zeros(sz)+1])
	Tsp_emis = sp_emis.T

	#Inverse de la comatrice
	inv_cov = inv(cov_vraie)

	weights  = inv_cov.dot(Tsp_emis).dot(inv(sp_emis.dot(inv_cov).dot(Tsp_emis)))
	return weights


#------------------------------------------
#- Determination des meilleurs parametres -
#------------------------------------------
def bestFit():
        differential_evolution(Variance, bounds=[(0,0.3)], popsize=1, maxiter=1)
#	minimize(Variance, 0.1, method="Nelder-Mead", tol=0.001)

#------------
# Affichage -
#------------
def showCMB(CMB_unique):
	hp.mollview(CMB_unique, title="CMB unique Galactique", coord=['G'], unit=r"$K_{CMB}$", norm="hist", min=min(CMB_unique), max=max(CMB_unique))
	plt.show()
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#MAIN

freq = [100,143, 217,353,545,857]
sz   = len(freq)
rgSz = np.arange(sz)		#Indices des cartes utilisees

#--------------------------------------------
# Declaration tableaux utiles pour la suite -
#--------------------------------------------
maxVal = np.zeros(sz)		#Valeur max de la bande centrale pour chaque carte
minVal = np.copy(maxVal)	#Idem mais min
meanIn = np.copy(maxVal)	#Valeur moyenne de la zone galactique pour chaque carte
meanOu = np.copy(maxVal) 	#Idem mais pour la zone extra-galactique

for i, fr in zip(rgSz, freq):

	#Recupere la carte, et la met dans un tableau (une colonne)
	maps = hp.read_map("data/HFI_256_"+str(fr)+"_convol.fits")

	#Cree le tableau complet lorsque la premiere carte est ouverte
	if (i==0):
		size            = hp.get_map_size(maps)
		array_maps      = np.zeros((size, sz))
		normArray       = np.copy(array_maps)
		nside 		= hp.get_nside(maps)
		nbPix		= px.nside2npix(nside)


		#Inutile pour le moment
		#Positions des pixels dans la carte
		posPix		= np.arange(nbPix)
		#Angles theta, phi des pixels des cartes (en RING)
		theta, phi	= px.pix2ang(nside, posPix, lonlat=True)

	#Stocke chaque carte dans la bonne colonne
	array_maps[:,i] = maps

#Reconstruction CMB
bestFit()
#a = Variance(0.0)
CMB_unique = (normArray.dot(weights)).T
showCMB(CMB_unique[0])

#Ecriture nouvelle carte
hp.write_map("data/HFI_256Norm1P.fits", CMB_unique[0])

