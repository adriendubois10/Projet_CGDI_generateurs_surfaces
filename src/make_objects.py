import numpy as np
import random as rd
from pygel3d import hmesh
import polyscope

colors = [(205,92,92), (255,69,0), (255,140,0), (189,183,107), (0,128,128), (138,43,226), (218,112,214), (127,255,212), (0,191,255), (128,128,128)]
colors = [(a/255, b/255, c/255) for (a, b, c) in colors] # liste des couleurs pour les cycles

# Structure
class Surface:

    def __init__(self, name="", vertices=[], faces=[], mesh=None, obj=None):
        self.edges = None
        self.mesh = None
        if obj != None:
            mesh = hmesh.load("data/"+obj)
        if mesh == None:
            self.name = name
            self.vertices = np.array(vertices) # liste de position des sommets
            self.faces = faces # liste de liste d'indices de sommets 
        else:
            self.mesh = mesh
            self.name = "my_mesh" if obj==None else obj
            self.vertices = mesh.positions()
            self.faces = [[v for v in mesh.circulate_face(f)] for f in mesh.faces()]
        self.nb_vertex = len(self.vertices)
        self.nb_edge = sum([len(face) for face in self.faces]) // 2 # pour une surface fermée, chaque arête est présente exactement dans 2 faces différentes
        self.nb_face = len(self.faces)
        self.euler_carac = self.nb_vertex - self.nb_edge + self.nb_face
        self.g = (2 - self.euler_carac) // 2
    
    def __str__(self):
        s = "Name: " + str(self.name) + "\n   Vertices: " + str(self.nb_vertex) +  "\n   Edges: "
        s += str(self.nb_edge) + "\n   Faces: " + str(self.nb_face) + "\n   Euler Caracteristic: "
        s += str(self.euler_carac) + "\n   Genus: " + str(self.g)
        return s
        
    def save_as_obj(self, path="data/"):
        with open(path+self.name+".obj", 'w') as f:
            for vertex in self.vertices:
                f.write(f"v {' '.join(map(str, vertex+1))}\n")
            for face in self.faces:
                f.write(f"f {' '.join(map(str, [v+1 for v in face]))}\n")

    def dist(self, x, y):
        return np.sqrt(sum([(self.vertices[x][k]-self.vertices[y][k])**2 for k in [0,1,2]]))

    def len_path(self, path):
        assert self.mesh != None
        if self.adj_dist == None:
            self.compute_dist()
        return sum([self.dist[path[i]][path[i+1]] for i in range(len(path)-1)])
    
    def view(self, point=None, paths=[], edges_sets=[], cloudpoint=False, faces=None):
        polyscope.init()
        if cloudpoint:
            polyscope.register_point_cloud("data", self.vertices)
        ps = polyscope.register_surface_mesh(self.name, self.mesh.positions(), self.faces, transparency=0.9, enabled=True)
        for (positions, edges, name, c) in edges_sets:
            highlighted_edges = polyscope.register_curve_network(
                name,
                positions,
                edges, 
                color = c
            )
        if point!=None:
            highlighted_point = polyscope.register_curve_network("Source", np.array([self.vertices[point]]),'line', color = (1,0,0))
        if faces!=None:
            ps = polyscope.register_surface_mesh(
                "fg", 
                self.mesh.positions(), 
                [[v for v in self.mesh.circulate_face(f)] for f in faces[0]], 
                transparency=0.9,
                enabled=True,
                color=(0.8, 0, 0)
            )
            ps = polyscope.register_surface_mesh(
                "fd", 
                self.mesh.positions(), 
                [[v for v in self.mesh.circulate_face(f)] for f in faces[1]], 
                transparency=0.9,
                enabled=True,
                color=(0, 0.8, 0)
            )
        polyscope.show()


# Exemples divers

def sphere(center=(0,0,0), radius=1, n=10):
    v , f = [[0,0,1]], []
    # Parcours de la sphere en coordonées polaire theta=k*pi/n et phi = 2k'*pi/n
    # Cas ni au nord ni au sud (k!=0, k!=n)
    for k in range(1, n):
        theta = k * np.pi / n
        for l in range(n):
            phi = 2 * l * np.pi / n
            x = np.sin(theta) * np.cos(phi)
            y = np.sin(theta) * np.sin(phi)
            z = np.cos(theta)
            v.append([radius * x + center[0], radius * y + center[1], radius * z + center[2]])
            if k<(n-1):
                f.append([(k-1)*n+2+(l % n), k*n+2+(l%n), k*n+2+((l+1)%n), (k-1)*n+2+((l+1)%n)])
    # Cas k=0 et k=n
    for l in range(n):
        f.append([0, 1+(l%n), 1+((l+1)%n)])
        f.append([n*(n-1)+1, n*(n-1)-((l+1)%n), n*(n-1)-(l%n)])
    v.append([0,0,-1])
    return Surface("sphere", reversed(v), f)

# Tore avec équation paramétrique utilisant les "B-Lines"
# Résultat : lignes du tore diagonales
def tore_blines(a=10, b=1, m=100, n=20):
    vertices, faces = [], []
    for k in range(m):
        u = 2 * k * np.pi / m
        cu, su = np.cos(u), np.sin(u)
        for l in range(n):
            v = 2 * l * np.pi / n
            cv, sv = np.cos(v), np.sin(v)
            x = np.sqrt(a**2-b**2)*sv*cu - (b+a*cv)*su
            y = np.sqrt(a**2-b**2)*sv*su + (b+a*cv)*cu
            z = b*sv
            vertices.append([x, y, z])
            faces.append([n*(k%m)+l+1, n*(k%m)+((l+1)%n)+1, n*((k+1)%m)+((l+1)%n)+1, n*((k+1)%m)+l+1])
    return Surface("tore", vertices, faces)

# Tore avec équation paramétrique déterminée à la main
# Résultat : lignes droitezs
def tore(a=2, b=1, m=100, n=40):
    vertices, faces = [], []
    for k in range(m):
        theta = 2 * k * np.pi / m
        ct, st = np.cos(theta), np.sin(theta)
        for l in range(n):
            phi = 2 * l * np.pi / n
            cphi, sphi = np.cos(phi), np.sin(phi)
            x = (a + b*cphi) * ct
            y = (a + b*cphi) * st
            z = b * sphi
            vertices.append([x, y, z])
            faces.append([n*(k%m)+l+1, n*(k%m)+((l+1)%n)+1, n*((k+1)%m)+((l+1)%n)+1, n*((k+1)%m)+l+1])
    return Surface("tore_smallsmall", vertices, faces)

# Surface de genre 2 : tore avec 2 trous (incomplet)
def tore2(a=2, b=1, m=100, n=40, eps=0):
    t1 = tore(a, b, m, n)
    t2 = tore(a, b, m, n)
    # Décalage du 2e tore en x (eps=0 collés)
    for i in range(t2.nb_vertex):
        t2.vertices[i][0] += 2*(a+b) + eps
    for j in range(t2.nb_face):
        t2.faces[j] = [v + t1.nb_vertex for v in t2.faces[j]]
    return Surface(name="2-tore_try", vertices=t1.vertices + t2.vertices, faces=t1.faces + t2.faces)

tore(a=2, b=1, m=10, n=4).save_as_obj()