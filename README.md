**Utiliser python 3 ou plus.**

## Dossier carte Planck
/home/vbonjean/PLANCK/

## Librairies pour lire les cartes Planck
Pour l'affichage:

import matplotlib as plt

Pour l'utilisation des cartes

import healpy as hp

import healpy.pixelfunc as px


# Fonctions healpy utiles

## Ouvrir une carte
carte = hp.read_map("nom_de_la_carte")

## Afficher une carte
Avec **ipython** utiliser **plt.ion()** pour le mode interactif.

Pour afficher la carte avec [hp.mollview](https://healpy.readthedocs.io/en/latest/generated/healpy.visufunc.mollview.html "page d'information"):

**hp.mollview(carte, title="titre", coord=['G', 'E'], unit="unite", norm="hist", min=-1, max=1, xsize=2000)**

Pour afficher une grille par-dessus avec [hp.graticule](https://healpy.readthedocs.io/en/latest/generated/healpy.visufunc.graticule.html#healpy.visufunc.graticule):

**hp.graticule()**

## Recuperer des infos sur les cartes
Pour avoir le nombre de pixels dans la carte avec [hp.get_map_size](https://healpy.readthedocs.io/en/latest/generated/healpy.pixelfunc.get_map_size.html#healpy.pixelfunc.get_map_size):

size = hp.get_map_size(map)

Pour avoir la valeur minimum du NSIDE recommandee avec [px.get_min_valid_nside](https://healpy.readthedocs.io/en/latest/generated/healpy.pixelfunc.get_min_valid_nside.html#healpy.pixelfunc.get_min_valid_nside):

min = px.get_min_valid_nside(size)

Pour avoir la valeur du NSIDE pour la carte avec [px.get_nside(map)](https://healpy.readthedocs.io/en/latest/generated/healpy.pixelfunc.get_nside.html#healpy.pixelfunc.get_nside):

nside = px.get_nside(map)

## Degrader ou augmenter la resolution d'une carte
Pour changer la resolution vers une valeur NSIDE avec [px.ud_grade](https://healpy.readthedocs.io/en/latest/generated/healpy.pixelfunc.ud_grade.html#healpy.pixelfunc.ud_grade):

**newmap = px.ud_grade(map, NSIDE)**

