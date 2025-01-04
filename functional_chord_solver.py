from pulp import LpProblem, LpMaximize, LpInteger, pulp
from chord_solver import ChordSolution, ConnectingSolution, ChordSolver, Vertex, Edge
from persistent_cache import persistent_cache
import os

FILE = os.path.basename(__file__)


def issue_chords(vertices: list[Vertex]) -> list[Edge]:
    n = len(vertices)
    chords = []
    for i in range(0, n - 1):
        for j in range(i + 1, n):
            chords.append((vertices[i], vertices[j]))
    return chords


def basic_ilp(k: int, chords: list[Edge]) -> pulp.LpProblem:
    prob = LpProblem("Chords in Polygon", LpMaximize)
    n = len(chords)
    n_squared = n * n
    ilp_variables = pulp.LpVariable.dicts(
        "chords", chords, lowBound=0, upBound=1, cat=LpInteger
    )
    for chord in chords:
        prob += (
            pulp.lpSum(
                [
                    ilp_variables[otherChord]
                    for otherChord in chords
                    if ChordSolver.cross(chord, otherChord)
                ]
            )
            + n_squared * ilp_variables[chord]
            <= k + n_squared
        )
    return prob, ilp_variables


@persistent_cache(0, FILE)
def solve_connectable(
    k: int,
    n: int,
    connecting_edges: tuple[int] | int = 0,
    force_connections: tuple[tuple[int, int]] = (),
) -> ChordSolution:
    if type(connecting_edges) is not tuple:
        connecting_edges = tuple((connecting_edges,))
    real_vertices = list(range(n))
    connector_vertices = [edge + 0.5 for edge in connecting_edges]
    vertices = sorted(real_vertices + connector_vertices)

    edges = issue_chords(vertices)
    (
        outer_edges,
        real_chords,
        connecting_chords,
        through_chords,
    ) = ConnectingSolution.split_edges(edges, connector_vertices)

    ilp_prob, ilp_variables = basic_ilp(k, edges)

    # target function
    ilp_prob += pulp.lpSum(
        [ilp_variables[real_chord] for real_chord in real_chords]
        + [
            0.5 * ilp_variables[shared_edge]
            for shared_edge in connecting_chords + outer_edges
        ]
        + [0.01 * ilp_variables[through_chord] for through_chord in through_chords]
    )

    # force number of connections
    for connecting_edge, number_of_connections in force_connections:
        ilp_prob += pulp.LpSum(
            [
                ilp_variables[chord]
                for chord in connecting_chords
                if connecting_edge + 0.5 in chord
            ]
            == number_of_connections
        )

    ilp_prob.solve()

    # analyse solution
    edges_in_solution = [edge for edge in edges if ilp_variables[edge].value() == 1]

    solution = ConnectingSolution(
        k,
        n,
        edges_in_solution,
        connector_vertices,
        draw_outer_edges=False,
    )
    return solution
