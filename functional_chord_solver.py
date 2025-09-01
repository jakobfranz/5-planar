from pulp import LpProblem, LpMaximize, LpInteger, pulp
from chord_solver import (
    ConnectingSolution,
    ChordSolution,
    Vertex,
    Edge,
    crossings,
)
from persistent_cache import persistent_cache
import os
from typing import Callable

FILE = os.path.basename(__file__)


def edge(a: Vertex, b: Vertex) -> Edge:
    """Make an edge connecting ```a``` and ```b```

    Edges are, by convention here, tuples of length >= 2
    with the vertices as first and second argument
    sorted ascending.

    Args:
        a (Vertex):
            A vertex
        b (Vertex):
            Another vertex (not necessarily different)
    Returns:
        Edge connecting ```a``` and ```b``` according to
        internal convention. ```len(edge) == 2```."""
    return tuple(sorted([a, b]))


def issue_chords(vertices: list[Vertex], issue_boundary: bool = False) -> list[Edge]:
    """Returns all chords of the complete graph K_n on ```vertices```.

    Order of ```vertices``` is the order in which the vertices lie on the outer face.
    Chords are all edges, that do not connect consecutive vertices.

    Optionally, the boundary edges can be included as well.

    Args:
        vertices (list[Vertex]):
            List of Vertices
        issue_boundary (bool):
            Boolean flag to determine, whether boundary edges should be
            included in the return. By default ```False```.
    Returns:
        List with all chords or all edges of the complete graph on ```vertices```.
        The number of returned edges is ```n choose 2 - n``` or ```n choose 2``` depending on ```issue_boundary```.
    """
    n = len(vertices)
    if issue_boundary:
        offset = 1
    else:
        offset = 2
    chords = []
    for i in range(0, n - offset):
        for j in range(i + offset, n):
            chords.append((vertices[i], vertices[j]))  #
    if not issue_boundary and n >= 3:
        chords.remove((vertices[0], vertices[n - 1]))
    return chords


def exists_var_name(edge: Edge) -> str:
    """Get name of the binary ILP variable that describes if an edge exists in the solution.

    If ```edge``` is possibly a multi-edge in the solution of an ILP,
    we need two variables for that edge. One to determine how often it
    exists, and one to determine whether it exists at least once.

    The function ```basic_ilp``` returnes a dictionary of all ilp
    variables it created. In it the first variable type (the
    multiplicity type) has ```edge``` as its key, while the second
    (existence) uses the name created by this function.

    Args:
        edge (Edge):
            An Edge
    Returns:
        name of the ilp-variable describing the existence of ```edge```.
    """
    return f"exists_{edge}"


def basic_ilp(
    k: int,
    chords: list[Edge],
    optimizer: int = LpMaximize,
    multiplicity: Callable[[Edge], int] = lambda chord: 1,
) -> tuple[pulp.LpProblem, dict[Edge, pulp.LpVariable]]:
    """Create a basic ILP encodig crossing constraints of outer-k-planar graphs.

    A pulp ILP is created. Each chord is represented by one or two
    variables. In a valid solution of that ILP each existing chord
    is crossed by at most ```k``` other chords.

    Both the pulp ILP ```prob``` and a dictionary with the chord
    variables ```ilp_variables``` is returned.

    Chords may be possible multi-edges, determined by ```multiplicity```. ```ilp_variables[chord]``` describes, how often
    '```chord``` is in the solution of the ILP. If chord may in fact be
    a multi-edge, we need an auxiliary variable describing whether
    '```chord``` exists at least once. This variable can be accessed
    through ```ilp_variables[exists_var_name(chord)]```.

    The chords crossing a given a chord are determined by ```chord_solver.crossings()```.

    The returned ILP can be expanded upon. Especially a target function
    must be set.

    Args:
        k (int):
            Maximal number of crossings per edge
        chords (list[Edge]):
            List of the chords that may exist in the solution
        optimizer (int):
            '```pulp.LpMaximize``` or ```pulp.LpMinimize``` depending on
            whether the target function should be maximized or minimized.
            By default ```LpMaximize```.
        multiplicity (Callable[[Edge], int]):
            A function that, given a chord, returns how often this chord
            may be in the solution at most. Must be an integer.
    Returns:
        A tuple containing

        - prob (pulp.LpProblem): A pulp ILP encoding the crossing constraint of outer-k-planar graphs for the given chords.
        - ilp_variables (dict[Edge, pulp.LpVariable]): A dictionary with the used variables in the ILP
    """
    prob = LpProblem("Chords in Polygon", optimizer)
    n = len(chords)

    # Big M
    n_squared = n * n

    # create multiplicity variables
    ilp_variables = {
        chord: pulp.LpVariable(
            f"chord_{chord}", lowBound=0, upBound=multiplicity(chord), cat=LpInteger
        )
        for chord in chords
    }

    for chord in chords:
        mult = multiplicity(chord)
        if mult == 1:
            chord_exist_var = ilp_variables[chord]
        else:
            # if a chord has a higher multiplicity we need an
            # auxiliary  binary variable (existence variable)
            chord_exist_var = pulp.LpVariable(
                f"exists_{chord}", lowBound=0, upBound=1, cat=LpInteger
            )
            ilp_variables[exists_var_name(chord)] = chord_exist_var
            # exist_var == 1 iff multiplicity >= 1
            prob += 1 / (mult + 1) * ilp_variables[chord] <= chord_exist_var

        # crossing constraint
        prob += (
            pulp.lpSum(
                [
                    ilp_variables[crossed_chord]
                    for crossed_chord in crossings(chord, chords)
                ]
            )
            + n_squared * chord_exist_var
            <= k + n_squared
        )
    return prob, ilp_variables


def is_incident(edge: Edge, vertices: Vertex | list[Vertex]) -> bool:
    """Decides if an edge is incident to one or more of the given vertices.

    Args:
        edge (Edge):
            An edge
        vertices (Vertex | list[Vertex]):
            A single vertex or a list of vertices
    Returns:
        Bool ```True``` if there is a vertex in ```vertices``` that is incident to ```edge```
    """
    if type(vertices) is Vertex:
        vertices = [vertices]
    return len(set(edge[:2]).intersection(set(vertices))) > 0


def incident_edges(vertices: Vertex | list[Vertex], edges: list[Edge]) -> list[Edge]:
    """Returns all edges in ```edges``` that are incident to a vertex of ```vertices```.

    Incidence is determined by ```is_incident```.

    Args:
        vertices (Vertex | list[Vertex]):
            A single or a list of vertices.
        edges (list[Edge]):
            A list of edges
    Returns:
        Sublist of ```edges``` with all edges that are incident to at least on vertex in ```vertices```.
    """
    return [edge for edge in edges if is_incident(edge, vertices)]


# TODO documentation and revision
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

    real_chords = issue_chords(real_vertices, issue_boundary=True)

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
