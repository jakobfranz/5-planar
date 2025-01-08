from pulp import LpProblem, LpMaximize, LpInteger, pulp
from chord_solver import ConnectingSolution, ChordSolution, Vertex, Edge, crossings
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


def basic_ilp(
    k: int, chords: list[Edge]
) -> tuple[pulp.LpProblem, dict[Edge, pulp.LpVariable]]:
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
                    ilp_variables[crossed_chord]
                    for crossed_chord in crossings(chord, chords)
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
    force_connections: tuple[tuple[int, int]] | int = (),
) -> ConnectingSolution:
    if type(connecting_edges) is not tuple:
        connecting_edges = tuple((connecting_edges,))

    n_squared = n * n

    real_vertices = list(range(n))
    connector_vertices = [edge + 0.5 for edge in connecting_edges]

    if type(force_connections) is not tuple:
        force_connections = tuple(
            ((con_edge, force_connections) for con_edge in connecting_edges)
        )

    force_connections_dict = {
        connecting_edge + 0.5: max_crossings
        for connecting_edge, max_crossings in force_connections
    }

    def max_number_of_crossings(connector_chord: Edge) -> int:
        if (
            connector_chord[0] in force_connections_dict.keys()
            and connector_chord[1] in force_connections_dict.keys()
        ):
            return min(
                force_connections_dict[connector_chord[0]],
                force_connections_dict[connector_chord[1]],
            )
        if connector_chord[0] in force_connections_dict.keys():
            return min(force_connections_dict[connector_chord[0]], k)
        if connector_chord[1] in force_connections_dict.keys():
            return min(force_connections_dict[connector_chord[1]], k)
        return k

    real_chords = issue_chords(real_vertices)

    connection_chords: list[Edge] = []
    # create connectable chords
    # issue chords
    connection_chords_template = [
        tuple(sorted((connector_vertex, real_vertex)))
        for real_vertex in real_vertices
        for connector_vertex in connector_vertices
        if abs(connector_vertex - real_vertex) >= 1
        and abs(connector_vertex - real_vertex) < n - 1
    ]

    connection_chords_grouped: list[list[Edge]] = [
        [(a, b, i) for i in range(max_number_of_crossings((a, b)) + 1)]
        for a, b in connection_chords_template
    ]
    connection_chords: list[Edge] = [
        chord for chord_group in connection_chords_grouped for chord in chord_group
    ]

    edges = real_chords + connection_chords
    (
        outer_edges,
        real_chords,
        connecting_chords,
        through_chords,
    ) = ConnectingSolution.split_edges(edges, connector_vertices)

    ilp_prob, ilp_variables = basic_ilp(k, edges)

    # limit number of crossings for connection chords according to i
    for connection_chord in connection_chords:
        if connection_chord[2] < k:
            ilp_prob += (
                pulp.lpSum(
                    [
                        ilp_variables[crossed_edge]
                        for crossed_edge in crossings(connection_chord, edges)
                    ]
                )
                + n_squared * ilp_variables[connection_chord]
                <= connection_chord[2] + n_squared
            )

    # force only one of identical chords to exist
    for connection_chord_group in connection_chords_grouped:
        ilp_prob += (
            pulp.lpSum(
                [
                    ilp_variables[parallel_chord]
                    for parallel_chord in connection_chord_group
                ]
            )
            <= 1
        )

    # target function
    ilp_prob += pulp.lpSum(
        [ilp_variables[real_chord] for real_chord in real_chords]
        + [0.5 * ilp_variables[outer_edge] for outer_edge in outer_edges]
        + [
            (0.5 + (k / 2 - connection_chord[2]) * 0.05)
            * ilp_variables[connection_chord]
            for connection_chord in connection_chords
        ]
        + [0.01 * ilp_variables[through_chord] for through_chord in through_chords]
    )

    ilp_prob.solve()

    # analyse solution
    edges_in_solution = [edge for edge in edges if ilp_variables[edge].value() == 1]

    print(ilp_prob.objective.value())
    print([con_chord for con_chord in edges_in_solution if len(con_chord) > 2])

    for i in range(len(edges_in_solution)):
        if len(edges_in_solution[i]) > 2:
            edges_in_solution[i] = edges_in_solution[i][:2]

    solution = ConnectingSolution(
        k,
        n,
        edges_in_solution,
        connector_vertices,
        draw_outer_edges=False,
    )
    return solution


def double(k, n) -> tuple[ConnectingSolution, ChordSolution]:
    half_gon = solve_connectable(k, n, n - 1, 3)
    chords = half_gon.chords

    chords2 = [tuple(sorted((a + n - 1, (b + n - 1) % (2 * n - 2)))) for a, b in chords]

    connecting_chords = half_gon.connecting_chords

    connected_chords = [
        tuple(sorted((a, (a + n - 1) % (2 * n - 2)))) for a, b in connecting_chords
    ]

    full_gon = ChordSolution(
        5, 2 * n - 2, chords + chords2 + connected_chords + [(0, n - 1)]
    )

    return half_gon, full_gon
