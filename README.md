# Projet CGDI : Courbes caractéristiques

## Librairies python nécéssaires

Les librairies utilisées dans le code se situent dans `requirements.txt`. A installer avec :

```pip install requirements.txt```

## Structure du code

Une classe `Surface` est définie dans `src/make_objects.py` qui va servir à manipuler les objets, et les afficher dans polyscope.
L'essentiel des algorithmes se situent dans `src/path.py`.
Pour reproduire les figures du rapport, il suffit d'exécuter les fonction dans `src/fig_gallery.py`.
Pour calculer la plus courte base homotopique sur une surface :

1. Placer une surface en `.obj` dans `data/`
2. Exécuter `python src/main.py`

Pour avoir plus de détails sur les structures intermédiaires utilisées, plutôt exécuter `fig4()` dans `src/fig_gallery.py` en y mettant le fichier `.obj` de votre choix.

## Rapport

Le rapport en pdf se situe dans le dossier principal et le code $\LaTeX$ dans le dossier `rapport`.

## Maillages exemples

Divers mailages sont présents dans `data`, réalisés à la main ou à l'aide d'Autodesk Fusion.

