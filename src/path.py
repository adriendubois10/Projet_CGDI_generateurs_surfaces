from make_objects import Surface
from pygel3d import hmesh
import random as rd
import numpy as np
from time import time
from tqdm import tqdm

def example_path(mesh, nb=200, origin=1):
    l = [origin]
    for _ in range(nb):
        lm = list(mesh.circulate_vertex(l[-1]))
        l.append(lm[rd.randint(0,len(lm)-1)])
    return l

# Arbre couvrant de poids (distance) min
def dichotomy_insertion(e, file, dist):
    index_a, index_b = 0, len(file)
    while (index_b - index_a) > 0:
        # print(index_a, index_b)
        index = (index_a + index_b) // 2
        if file[index][2] == dist:
            break
        elif file[index][2] > dist:
            index_b = index
        else:
            index_a = index + 1
    return file[:index_a] + [(*e, dist)] + file[index_a:]

def short_tree(source, mesh0, adjacent_distances=None):
    file = [(source, source, 0)]
    pos = mesh0.positions()
    visited = [False for _ in mesh0.vertices()]
    tree = [source for _ in mesh0.vertices()]
    while len(file) > 0 and False in visited:
        x, ancestor_x, dist = file.pop(0)
        if visited[x]:
            continue
        visited[x] = True
        tree[x] = (ancestor_x, dist)
        for y in mesh0.circulate_vertex(x):
            if not visited[y]:
                if adjacent_distances!=None:
                    file = dichotomy_insertion((y, x), file, dist + adjacent_distances[x][y])
                else:
                    file = dichotomy_insertion((y, x), file, dist + np.sqrt(sum([(pos[x][k]-pos[y][k])**2 for k in [0,1,2]])) )
    return tree, source

def edges_from_tree(tree):
    edges = []
    for y in range(len(tree)):
        if y < tree[y][0]:
            edges.append([y, tree[y][0]])
        if y > tree[y][0]:
            edges.append([tree[y][0], y])
    return edges

# Eppstein tree-cotree decompositon

## Graphe dual à partir d'un maillage
def common_edge(mesh, f1, f2):
    f1_vertices = set(mesh.circulate_face(f1))
    common_vertices = [v for v in mesh.circulate_face(f2) if v in f1_vertices]
    return common_vertices if common_vertices else None

def common_edge_knowing_u(mesh, f1, f2, u):
    lf1_set = set(mesh.circulate_face(f1))
    lf2_list = list(mesh.circulate_face(f2))
    for i in range(len(lf2_list)):
        if lf2_list[i] == u:
            v = lf2_list[(i+1) % len(lf2_list)]
            w = lf2_list[(i-1) % len(lf2_list)]
            if v in lf1_set:
                return v
            elif w in lf1_set:
                return w
            else:
                return None
    return None

def dual_G(mesh):
    edges = [] # liste des e*, les arêtes du graphe dual comme paires de faces
    crossed_edges = [] # liste des arêtes du graphe normal traversées par l'arête e* correspondante
    for f1 in mesh.faces():
        for u in mesh.circulate_face(f1, mode='v'):
            for f2 in mesh.circulate_vertex(u, mode='f'):
                if f2==f1:
                    continue
                cv = common_edge_knowing_u(mesh, f2, u)
                if cv != None:
                    edges.append((f1, f2))
                    crossed_edges.append(cv)
    return edges, crossed_edges              

## Dual de G privé de son arbre couvrant directement (Ceci est le reduce cut locus de G en x, la source de l'arbre couvrant T)
def GminusT_adjacence(Gmesh, T, adjacent_distance=None):
    t0 = time()
    t_ce = 0
    face_adj = [[] for _ in Gmesh.faces()]
    for f in Gmesh.faces():
        for fadj in Gmesh.circulate_face(f, mode='f'):
            tcheck = time()
            u, v = common_edge(Gmesh, f, fadj)
            t_ce += time()-tcheck
            if T[u][0]==v or T[v][0]==u:
                continue
            face_adj[f].append((fadj, (u,v)))
            face_adj[fadj].append((f, (u,v)))
    t_total = time()-t0
    # print("Common edges:", round(100 * t_ce/t_total, 4))
    return face_adj

def face_in_GminusT(mesh, f, GmT_adj):
    l_vertices = list(mesh.circulate_face(f))
    for i in range(len(l_vertices)):
        if l_vertices[(i+1) % len(l_vertices)] not in GmT_adj[l_vertices[i]]:
            return False
    return True

# def dual_GminusT_fail(Gmesh, T): # Version ou l'on calcule le dual de (G\T) littéralement et non 
#     # Gmesh est le maillage correspondant au graphe G et T matrice d'adjacence de l'arbre couvrant
#     edges = [] # liste des e*, les arêtes du graphe dual comme paires de faces
#     crossed_edges = [] # liste des arêtes du graphe normal traversées par l'arête e* correspondante
#     GminusT_adj = GminusT_adjacence(Gmesh, T) # attention, adjacence de faces
#     for f1 in Gmesh.faces():
#         if not face_in_GminusT(Gmesh, f1, GminusT_adj): # on cherche une face de G qui soit dans G\T, pour en faire un sommet dans (G\T)*
#             continue
#         for u in Gmesh.circulate_face(f1, mode='v'):
#             for f2 in Gmesh.circulate_vertex(u, mode='f'):
#                 if f2==f1 or not face_in_GminusT(Gmesh, f2, GminusT_adj):
#                     continue
#                 cv = common_edge(Gmesh, f1, f2)
#                 if cv != None:
#                     edges.append((f1, f2))
#                     crossed_edges.append(cv)
#     return edges, crossed_edges

def sigma_e(T, source, e):
    # Trouve la boucle la plus courte passant par x (source de T) et uen arête e dans G\T
    # Le chemin le plus court est connu car c'est celui qui relie les extrémités de e à x dans T
    u, v = e
    x = source
    path_u_x, path_v_x = [u], [v]
    while u != x:
        u = T[u][0]
        path_u_x.append(u)
    while v != x:
        v = T[v][0]
        path_v_x.append(v)
    return list(reversed(path_u_x)) + path_v_x

def len_path(Gmesh, path):
    pos = Gmesh.positions()
    return sum([np.sqrt(sum((pos[path[i]][k]-pos[path[i+1]][k])**2 for k in [0,1,2])) for i in range(len(path)-1)])  

def dual_GminusT(Gmesh, T, adjacent_distance=None):
    num_faces = len(Gmesh.faces())
    dual = [[] for _ in range(num_faces)]
    for f in Gmesh.faces():
        for fadj in Gmesh.circulate_face(f, mode='f'):
            if f < fadj:
                u, v = common_edge(Gmesh, f, fadj)
                if T[u][0]==v or T[v][0]==u:
                    continue
                dist_e0 = T[u][1]  # distance de la source de T à e[0]
                dist_e1 = T[v][1]  # idem
                len_e = len_path(Gmesh, [u,v]) if adjacent_distance==None else adjacent_distance[u][v]
                dual[f].append((fadj, (u,v), dist_e0 + len_e + dist_e1))
                dual[fadj].append((f, (u,v), dist_e0 + len_e + dist_e1))
    return dual


def faces_dual(Gmesh, T, dual=None):
    is_in_face = [False for _ in Gmesh.faces()]
    if dual==None:
        dual = dual_GminusT(Gmesh, T)
    for f in Gmesh.faces():
        for (fadj, _, _) in dual[f]:
            is_in_face[fadj] = True
    return [f for f in Gmesh.faces() if is_in_face[f]]

def barycenter(mesh, indexes_points):
    return (1 / len(indexes_points)) * sum([mesh.positions()[i] for i in indexes_points])

def mesh_dual(Gmesh, T, dual=None):
    if dual==None:
        dual = dual_GminusT(Gmesh, T)
    pos = [barycenter(Gmesh, list(Gmesh.circulate_face(f))) for f in Gmesh.faces()]
    edges = []
    for f in range(len(dual)):
        for (fadj, _, _) in dual[f]:
            if f < fadj and not [f, fadj] in edges:
                edges.append([f, fadj])
    return pos, edges

## Arbre couvrant de poids max sur (G\T)*

def compute_tstar(Gmesh, T, dual=None, adjacent_distance=None):
    # Calcul des structures (G\T)* avec les poids
    if dual==None:
        dual = dual_GminusT(Gmesh, T, adjacent_distance=adjacent_distance)
    # Arbre couvrant de poids max (poids min avec poids négatifs)
    fd = faces_dual(Gmesh, T, dual=dual)
    source = fd[0] # source au hasard
    file = [(source, source, 0)]
    visited = [True for _ in Gmesh.faces()] 
    for f in fd:
        visited[f] = False
    tree_star = [[] for _ in Gmesh.faces()]

    while len(file) > 0 and False in visited:
        f, ancestor_f, _ = file.pop(0)
        if visited[f]:
            continue
        visited[f] = True
        tree_star[f] = (ancestor_f, 0)
        for (fadj, _, dist) in dual[f]:
            if not visited[fadj]:
                file = dichotomy_insertion((fadj, f), file, -dist)
    return tree_star

# Boucles les plus courtes
def precompute_adjacent_distances(mesh):
    adjacent_distances = {}
    positions = mesh.positions()
    for v in mesh.vertices():
        adjacent_distances[v] = {}
        for neighbor in mesh.circulate_vertex(v):
            distance = np.linalg.norm(np.array(positions[v]) - np.array(positions[neighbor]))
            adjacent_distances[v][neighbor] = distance
    return adjacent_distances

def shortest_loops(Gmesh, source, T=None, adjacent_distance=None):
    if T==None:
        T, source = short_tree(source, Gmesh, adjacent_distances=adjacent_distance)
    t = time()
    dual = dual_GminusT(Gmesh, T, adjacent_distance=adjacent_distance)
    tdual = time()-t
    t = time()
    Tstar = compute_tstar(Gmesh, T, dual=dual, adjacent_distance=adjacent_distance)
    t_tstar = time()-t
    t = time()
    inter_edges, inter_faces = [], []
    for f in Gmesh.faces():
        for (fadj, e, dist) in dual[f]:
            if fadj <= f:
                continue
            if Tstar[f][0]!=fadj and Tstar[fadj][0]!=f:
                if (f, fadj) not in inter_faces:
                    inter_faces.append((f, fadj))
                    inter_edges.append((e, dist))
    t_inter = time()-t
    #print([round(t, 3) for t in [tdual, t_tstar, t_inter]])
    return [(sigma_e(T, source, e), dist) for (e, dist) in inter_edges]

# Homolgies

# Liste des plus courtes boucles obtenues pour depuis toutes les sources
def all_shortest_loops(Gmesh):
    candidates_loops = [] # liste (boucle, longueur)
    adj_dist = precompute_adjacent_distances(Gmesh)
    for x in tqdm(Gmesh.vertices()):
        candidates_loops += shortest_loops(Gmesh, x, adjacent_distance=adj_dist)
    print("Sorting {} elements...".format(len(candidates_loops)))
    candidates_loops.sort(key = lambda x: x[1])
    print("Sort done")
    return candidates_loops

def is_connected(vertices, adj, pr=False):
    start_vertex = vertices[0]
    visited = [False for _ in vertices]
    pile = [start_vertex]
    visited[start_vertex] = True
    while pile:
        current_vertex = pile.pop()
        for neighbor in adj[current_vertex]:
            if not visited[neighbor]:
                pile.append(neighbor)
                visited[neighbor] = True
    if pr:
        print(visited)
    return all(visited)

def adjacent_faces_edge(mesh, u, v):
    fadj_u = set(mesh.circulate_vertex(u, mode='f'))
    return tuple([f for f in mesh.circulate_vertex(v, mode='f') if f in fadj_u and f!=-1])

# Couper un graphe selon le chemin cyclique path, 
# liste de sommets avec même sommet au départ et à l'arrivée
# On considère le graph comme une adjacence de ses faces, 
# cela marche car on ne coupe qu'avec des cycles
def cut_graph(mesh, face_adj, path):
    for i in range(len(path)-1):
        e = path[i], path[i+1]
        tf = tuple([f for f in adjacent_faces_edge(mesh, *e) if f!=-1])
        if len(tf)>1:
            f, fadj = tf
            if fadj in face_adj[f]:
                face_adj[f].remove(fadj)
            if f in face_adj[fadj]:
                face_adj[fadj].remove(f)
    return face_adj

def genus(mesh):
    nb_vertex, nb_face = len(mesh.vertices()), len(mesh.faces())
    nb_edge = sum([len(list(mesh.circulate_face(f))) for f in mesh.faces()]) // 2
    euler_carac = nb_vertex - nb_edge + nb_face
    return (2 - euler_carac) // 2

def equality_cycles(c1, c2):
    if len(c1)!=len(c2):
        return False
    for v in c1:
        if not v in c2:
            return False
    return True

def is_cycle_in_basis(basis, cycle):
    for c in basis:
        if equality_cycles(c, cycle):
            return True
    return False

def shortest_cycle(mesh):
    candidates_cycles = all_shortest_loops(mesh)
    print(len(candidates_cycles))
    faces = mesh.faces()
    face_adj = [list(mesh.circulate_face(f)) for f in faces]
    g = genus(mesh)
    basis = []
    print(is_connected(faces, face_adj, pr=True))
    for (cycle, _) in candidates_cycles:
        if is_cycle_in_basis(basis, cycle):
            continue
        face_adj_copy = face_adj[:]
        cut_graph(mesh, face_adj_copy, cycle)
        if is_connected(faces, face_adj_copy):
            basis.append(cycle)
            face_adj = face_adj_copy
        if len(basis) == 2*g:
            break
    return basis

        















    