**Utiliser python 3 ou plus.**

## Dossier carte Planck
/home/vbonjean/PLANCK/

## Librairies pour lire les cartes Planck
import healpy as hp

import matplotlib as plt

## Ouvrir une carte
carte = hp.read_map("nom_de_la_carte")

## Afficher une carte
Avec **ipython** utiliser **plt.ion()** pour le mode interactif.

Pour afficher la carte avec [mollview](https://healpy.readthedocs.io/en/latest/generated/healpy.visufunc.mollview.html "page d'information"):

**hp.mollview(carte, title="titre", coord=['G', 'E'], unit="unite", norm="hist", min=-1, max=1, xsize=2000)**

Pour afficher une grille par-dessus:

**hp.graticule()**
