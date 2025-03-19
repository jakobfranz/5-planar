# issue
# K_6
# Edges as vertices, complete
# => K_12

# Boundary exists => vertex has deg >= 3

# maximize edges - crossings * 1/100000


import pulp

# from functional_chord_solver import basic_ilp, issue_chords
from chord_solver import Vertex, Edge, crossings, ConnectingSolution
from itertools import combinations
from enum import Enum


def make_edge(a: Vertex, b: Vertex, n: int = 12) -> Edge:
    return tuple(sorted((a % n, b % n)))


class EdgeType(Enum):
    PLANAR = 0
    CHORD = 1


n = 6
k = 5

n_cubed = (2 * n) ** 3

vertices: list[Vertex] = [0, 2, 4, 6, 8, 10]
edge_vertices: list[Edge] = [1, 3, 5, 7, 9, 11]

c3_chords = [(0, 6), (2, 8), (4, 10)]

real_chords = list(combinations(vertices, 2))

# the pairwise crossing edges exist in any case
for c3_chord in c3_chords:
    real_chords.remove(c3_chord)

outer_edges = [make_edge(v, v + 2) for v in vertices]
# remove outer edges from real_chords to avoid duplicates
for o_edge in outer_edges:
    real_chords.remove(o_edge)

short_sticks = sum(
    [[make_edge(v - 3, v), make_edge(v, v + 3)] for v in edge_vertices], []
)

middle_segments = [make_edge(v, v + 2) for v in edge_vertices[:-1]]

# each short_stick and middle segment exists up to 3 times
# outer edges may be real outer edges (planar) or chords of the graph
edges = (
    real_chords
    + sum(
        [[(a, b, i) for i in range(3)] for (a, b) in short_sticks + middle_segments], []
    )
    + sum([[(a, b, typ) for typ in EdgeType] for (a, b) in outer_edges], [])
)


prob = pulp.LpProblem("Replace_c3", pulp.LpMaximize)
ilp_vars = pulp.LpVariable.dicts(
    "edges", edges, lowBound=0, upBound=1, cat=pulp.LpInteger
)

# any edge may be crossed at most k times if it exists
for edge in edges:
    prob += (
        pulp.lpSum([ilp_vars[other_edge] for other_edge in crossings(edge, edges)])
        + len(crossings(edge, c3_chords))
        + n_cubed * ilp_vars[edge]
        <= k + n_cubed
    )

# the same goes for c3_chords, although they always exist
for edge in c3_chords:
    prob += (
        pulp.lpSum([ilp_vars[other_edge] for other_edge in crossings(edge, edges)])
        <= k - 2
    )

# every (pseudo)-outer edge is either planar and real or a chord and crossed at least 3 times
for a, b in outer_edges:
    chord = (a, b, EdgeType.CHORD)
    planar = (a, b, EdgeType.PLANAR)
    prob += (
        pulp.lpSum([ilp_vars[other_edge] for other_edge in crossings(chord, edges)])
        + n_cubed
        >= 3 + n_cubed * ilp_vars[chord]
    )
    prob += (
        pulp.lpSum([ilp_vars[other_edge] for other_edge in crossings(planar, edges)])
        + n_cubed * ilp_vars[planar]
        <= 0 + n_cubed
    )
    prob += pulp.lpSum([ilp_vars[planar], ilp_vars[chord]]) == 1

# optimal solution with all 18 chords is still possible at this point

# Every vertex needs chord degree >= 3. If the 2-hop containing v does not exist,
# v needs chord degree > 3. Including the outer edges (which must exist, otherwise immediately <= 17) this results in minimum degree 5 or 6
for vertex in vertices:
    prob += pulp.lpSum(
        [
            ilp_vars[incident_edge]
            for incident_edge in edges
            if vertex in incident_edge[:2]
        ]
    ) - ilp_vars[
        make_edge(vertex - 2, vertex + 2)
    ] + 2 * n_cubed >= 5 + n_cubed * pulp.lpSum(
        [
            ilp_vars[(*planar_edge, EdgeType.PLANAR)]
            for planar_edge in (
                make_edge(vertex - 2, vertex),
                make_edge(vertex, vertex + 2),
            )
        ]
    )


# target function, maximize number of chords
prob += (
    pulp.lpSum([ilp_vars[edge] for edge in edges])
    + 3
    + 1 / n_cubed * pulp.lpSum([ilp_vars[real_chord] for real_chord in real_chords])
)

prob.solve()
# analyse solution
edges_in_solution = [edge for edge in edges if ilp_vars[edge].value() == 1]
print(len(edges_in_solution) + 3)
print(edges_in_solution)
for i in range(len(edges_in_solution)):
    if len(edges_in_solution[i]) > 2:
        edges_in_solution[i] = edges_in_solution[i][:2]

solution = ConnectingSolution(
    k,
    12,
    edges_in_solution + c3_chords,
    edge_vertices,
    draw_outer_edges=False,
)

solution.visualize()
