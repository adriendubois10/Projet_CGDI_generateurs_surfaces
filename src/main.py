import polyscope
import numpy as np
from pygel3d import hmesh
from make_objects import *
from path import *

# Placer le .obj dans data
objet = "2-tore.obj"
S = Surface(obj=objet)
source = 1 # Indice du point source

# Calcul de la plus courte base homotopique
loops = [l[0] for l in shortest_loops(S.mesh, 1)]

# Affichage
S.view(
    point=source, 
    edges_sets=[(np.array([S.mesh.positions()[v] for v in loops[i]]), 'line', "Boucle n°"+str(i), colors[i]) for i in range(len(loops))]
)