**Utiliser python 3 ou plus.**

## Dossier carte Planck
/home/vbonjean/PLANCK/

**WARNING** : les cartes 100-353GHz sont en unite K_CMB, les cartes 545-857GHz sont en MJy/sr.

**WARNING** : la carte a 100GHz contient un pixel UNSEEN.

La conversion est donnee en Table 6 (p16) de l'article [Planck IX](https://arxiv.org/abs/1303.5070).

## Librairies pour lire les cartes Planck
Pour l'affichage:

```python
import matplotlib as plt
```

Pour l'utilisation des cartes

```python
import healpy as hp
import healpy.pixelfunc as px
```

# Ouvrir le spectre Planck

## Fichier spectre Planck (fichier fits)

/home/vbonjean/Planck

## Librairie pour ouvrir le fichier

```python
from astropy.io import fits
```

## Ouvrir les parties basses et hautes frequences

```python
loc = "/home/vbonjean/Planck"
lfi = fits.getdata(loc, 1) #basses frequences
hfi = fits.getdata(loc, 7) #hautes frequences
```

## Afficher le nom des colonnes

```python
lfi.columns
hfi.columns
```

## Afficher les spectre complet

```python
plt.plot(np.concatenate(lfi['Ell'], hfi['Ell']), np.concatenate(lfi['D_Ell'], hfi['D_Ell']))
```

# Fonctions healpy utiles

## Ouvrir une carte

```python
carte = hp.read_map("nom_de_la_carte")
```

## Afficher une carte
Avec **ipython** utiliser `plt.ion()` pour le mode interactif.

Pour afficher la carte avec [hp.mollview](https://healpy.readthedocs.io/en/latest/generated/healpy.visufunc.mollview.html "page d'information"):

```pyton
hp.mollview(carte, title="titre", coord=['G', 'E'], unit="unite", norm="hist", min=-1, max=1, xsize=2000)
```

Pour afficher une grille par-dessus avec [hp.graticule](https://healpy.readthedocs.io/en/latest/generated/healpy.visufunc.graticule.html#healpy.visufunc.graticule):

```python
hp.graticule()
```

## Recuperer des infos sur les cartes
Pour avoir le nombre de pixels dans la carte avec [hp.get_map_size](https://healpy.readthedocs.io/en/latest/generated/healpy.pixelfunc.get_map_size.html#healpy.pixelfunc.get_map_size):

```python
size = hp.get_map_size(map)
```

Pour avoir la valeur minimum du NSIDE recommandee avec [px.get_min_valid_nside](https://healpy.readthedocs.io/en/latest/generated/healpy.pixelfunc.get_min_valid_nside.html#healpy.pixelfunc.get_min_valid_nside):

```python
minval = px.get_min_valid_nside(size)
```

Pour avoir la valeur du NSIDE pour la carte avec [px.get_nside(map)](https://healpy.readthedocs.io/en/latest/generated/healpy.pixelfunc.get_nside.html#healpy.pixelfunc.get_nside):

```python
nside = px.get_nside(map)
```

## Degrader ou augmenter la resolution d'une carte
Pour changer la resolution vers une valeur NSIDE avec [px.ud_grade](https://healpy.readthedocs.io/en/latest/generated/healpy.pixelfunc.ud_grade.html#healpy.pixelfunc.ud_grade):

```python
newmap = px.ud_grade(map, NSIDE)
```

