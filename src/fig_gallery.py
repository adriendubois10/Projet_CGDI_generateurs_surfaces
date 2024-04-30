from make_objects import *
from path import *
from time import time

# Fig 1.a: Arbre couvrant tore
def fig1():
    S = Surface(obj="tore_droit.obj")
    print(S)
    T, source = short_tree(1, S.mesh)
    tree_edges = np.array(edges_from_tree(T))
    S.view(
        point=source, 
        edges_sets=[
            (S.mesh.positions(), tree_edges, "Arbre couvrant", (0.8, 0.03, 0.1))
        ]
    )

# Fig XX: Arbre couvrant cube percé
def fig2():
    S = Surface(obj="cube_perce.obj")
    print(S)
    T, source = short_tree(1, S.mesh)
    tree_edges = np.array(edges_from_tree(T))
    S.view(
        point=source, 
        edges_sets=[
            (S.mesh.positions(), tree_edges, "Arbre couvrant", (0.8, 0.03, 0.1))
        ]
    )

# Fig 1.a 1.b: Arbre couvrant et/ou Reduced cut locus sur tore
def fig3():
    S = Surface(obj="tore_droit.obj")
    print(S)
    T, source = short_tree(1, S.mesh)
    tree_edges = np.array(edges_from_tree(T))
    pos_dual, edges_dual = mesh_dual(S.mesh, T)
    S.view(point=source, 
           edges_sets=[
            (np.array(pos_dual), np.array(edges_dual), "(G\T)*", (1,1,1)),
            (S.mesh.positions(), tree_edges, "Arbre couvrant", (0.8, 0.03, 0.1))
            ]
    )

# Fig 1.a 1.b 1.c: Arbre couvrant tore ET reduced cut locus ET Tstar
def fig4():
    S = Surface(obj="tore_droit.obj")
    print(S)
    T, source = short_tree(1, S.mesh)
    tree_edges = np.array(edges_from_tree(T))
    pos_dual, edges_dual = mesh_dual(S.mesh, T)
    tstar = compute_tstar(S.mesh, T)
    treestar_edges = np.array(edges_from_tree(tstar))
    S.view(point=source, 
           edges_sets=[
            (np.array(pos_dual), np.array(edges_dual), "(G\T)*", (1,1,1)),
            (S.mesh.positions(), tree_edges, "Arbre couvrant", (0.8, 0.03, 0.1)),
            (np.array(pos_dual), treestar_edges, "T*", (0,0.9,0)),
            ]
    )

# Fig 1.d: Shortest loops [tore]
def fig5():
    S = Surface(obj="tore_droit.obj")
    print(S)
    T, source = short_tree(1, S.mesh)
    loops = [l[0] for l in shortest_loops(S.mesh, 1)]
    S.view(point=source, 
           edges_sets=[(np.array([S.mesh.positions()[v] for v in loops[i]]), 'line', "Boucle n°"+str(i), colors[i]) for i in range(len(loops))]
    )

# Fig 2: Shortest loops sur autres sturctures
def fig6():
    S = Surface(obj="2-tore.obj")
    print(S)
    t = time()
    T, source = short_tree(1, S.mesh)
    tree_time = time()-t
    t = time()
    adj_dist = precompute_adjacent_distances(S.mesh)
    loops = [l[0] for l in shortest_loops(S.mesh, 1, T=T, adjacent_distance=adj_dist)]
    tf = time()-t 
    print("Computing time =", round(tf, 1), "sec")
    print("Estimated time for homology =", round(4000*tf, 1), "sec")
    S.view(point=source, 
           edges_sets=[(np.array([S.mesh.positions()[v] for v in loops[i]]), 'line', "Boucle n°"+str(i), colors[i % len(colors)]) for i in range(len(loops))]
    )

# Fig XX: Shortest cycles
# Est correct sur le principe mais ne fonctionne pas à cause d'un détail d'implémentation
# La coupure du graphe "cut_graph(...)" et le test de connexité "is_connect" ne fonctionnent pas 
# malgré les efforts investis pour le corriger

def fig7():
    S = Surface(obj="2-tore.obj")
    print(S)
    loops = shortest_cycle(S.mesh)
    print(len(loops))
    S.view(
        # edges_sets=[(np.array([S.mesh.positions()[v] for v in list(e)]), 'line', "Arete", colors[0])],
        edges_sets=[(np.array([S.mesh.positions()[v] for v in loops[i]]), 'line', "Cycle n°"+str(i), colors[i%len(colors)]) for i in range(len(loops))]
        # faces=[f]
    )

fig7()
