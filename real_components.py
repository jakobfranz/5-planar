# Solve some cases when partitioning an outer-k-planar graph
# into a H-Block and real components.

from pulp import LpMinimize, LpInteger, LpVariable, lpSum
from chord_solver import ChordSolution, Edge, crossings, ConnectingSolution
from functional_chord_solver import (
    issue_chords,
    basic_ilp,
    is_incident,
    incident_edges,
    exists_var_name,
    edge,
)

from persistent_cache import persistent_cache
import os
import math
import networkx as nx

FILE = os.path.basename(__file__)


class RCSolution(ChordSolution):
    def __init__(self, k, n, chords, rc_vertices, draw_outer_edges=True):
        super().__init__(k, n, chords, draw_outer_edges)
        self.rc_vertices = rc_vertices

    def _point_positions(self):
        return super()._point_positions() | {
            rc_vertex: (
                math.sin((2 * math.pi / self.n) * rc_vertex),
                math.cos((2 * math.pi / self.n) * rc_vertex),
            )
            for rc_vertex in self.rc_vertices
        }

    # def to_nx_graph(self):
    #     G = super().to_nx_graph()
    #     vertices = list(range(self.n)) + self.rc_vertices
    #     boundary_edges = [
    #         (
    #             *edge(vertices[i], vertices[(i + 1) % len(vertices)]),
    #             {"multiplicity": 1, "boundary": True},
    #         )
    #         for i in range(len(vertices) - 1)
    #         if not (vertices[i] % 1 == 1 / 4 or vertices[i + 1] % 1 == 3 / 4)
    #     ]
    #     G.add_edges_from(boundary_edges)
    #     return G

    def visualize(self):
        G = self.to_nx_graph()
        colors = {True: "grey", False: "black"}
        edge_colors = [colors[edge.boundary] for edge in G.edges]
        nx.draw(G, self._point_positions())


def non_isomorphic_configurations(
    q: int, rc_types: list[int] = [0, 1, 2, 3]
) -> list[list[int]]:
    """Generate all non-isomorph real component configurations.

    A real component configuration is represented by a C_2q graph
    (Cycle of length 2q), where the vertices represent V_D and
    the edges are assigned the rc-type which sits between the
    adjacent V_D-vertices. Type `0` means, that there is a boundary
    edge instead of a real component.

    Real component configurations are isomorph to each other,
    if one can be transformed into the other through
    rotation or mirroring.

    We label the vertices in the order they appear on C_2q. We order
    the edges in the same order. Of all the vertex labelings we
    define that labeling to be canonical where the sequence of
    rc-types of the edges is lexicographiccaly the highest.

    There may be actually multiple canonical labelings for an rc-
    configuration. But then the corresponding rc-type-sequences
    are identical so the difference does not matter.

    The elimination of isomorphic rc-configurations reduces the
    number of configurations for `q=3` from 1300 to 430 and for
    `q=4` from 65536 to 4435.

    Args:
        q (int):
            Number of edges in `D`, important for the number of
            sides in C_2q
        rc_types (list[int]):
            Possible types of real components. Type `0` should
            mean, that the rc-type is non-existent.

    Returns:
        A list with all non-isomorphic real component
        configurations
    """

    def _cases(remaining_sides: int, rc_types: list[int]) -> list[list[int]]:
        """Generate possible rc-sequences recursively.

        Some pruning is done to the search-space."""
        if remaining_sides <= 0:
            return []
        elif remaining_sides == 1:
            return [[rc_type] for rc_type in rc_types]
        else:
            rek_cases = _cases(remaining_sides - 1, rc_types)
            cases = [
                case + [rc_type]
                for rc_type in rc_types
                for case in rek_cases
                if rc_type <= case[0]
                # if rc_type > case[0] the canonical labeling would
                # rather put rc_type then case[0] as the first
                # element of the rc-sequence
            ]
            return cases

    cases = _cases(2 * q, rc_types)
    # Remove non-canonical labellings
    cases = [
        case
        for case in cases
        if max(
            # These are the rc-sequences for all vertex labelings
            # that are ordered by the appearence on C_2q
            [
                case[i::direction] + case[:i:direction]
                for i in range(2 * q)
                for direction in (1, -1)
            ]
        )
        # only if case if the lexicographically highest rc-sequence
        # the vertex labeling was canonical and the rc-sequence
        # needs to be considered.
        == case
    ]

    return cases


def analyse_real_component_configuration(
    k: int, q: int, rc_configuration: list[int], bound_to_show: tuple[float, float]
) -> tuple[bool, int, RCSolution]:
    number_of_real_components = len([rc for rc in rc_configuration if rc != 0])

    V_D = list(range(2 * q))
    D = [(V_D[i], V_D[i + q]) for i in range(q)]

    bound_linear, bound_constant = bound_to_show

    Delta_min = (
        number_of_real_components * (2 * bound_linear - bound_constant - 1)
        + q * (k - q - 2 * bound_linear + 4)
        + bound_constant
    )

    min_Delta = {
        0: 0,
        1: 3 * bound_linear - bound_constant - 3,
        2: 4 * bound_linear - bound_constant - 6,
        3: 0,
    }

    # if the sum of the min_Deltas is greater than Delta_min
    # we can stop here
    min_delta_sum = sum(
        [min_Delta[real_component_type] for real_component_type in rc_configuration]
    )
    if min_delta_sum >= Delta_min:
        print(f"Delta >= {min_delta_sum} >= Delta_min = {Delta_min}")
        print(f"therefore upper bound is shown for configuration {rc_configuration}")
        # TODO Alternative to RCSolution for these min_Delta cases
        return (min_delta_sum, None)

    # Create ILP
    rc_vertices = [
        [
            (V_D[i] + ii / (rc_configuration[i] + 1)) % (2 * q)
            for ii in range(rc_configuration[i] + 2)
        ]
        for i in range(2 * q)
    ]
    vertices = sum([rcv[:-1] for rcv in rc_vertices], [])

    chords = issue_chords(vertices, False)
    # chords within a real component
    # including for shared edge between rc and V_C
    inner_rc_chords: list[Edge] = [
        issue_chords(real_component_vertices, issue_boundary=False)
        + [edge(real_component_vertices[0], real_component_vertices[-1])]
        for real_component_vertices in rc_vertices
    ]

    # in the case |V_i| >= 5 for any i, the middle vertex represents an arbitrary number of vertices in G_j.
    # As such, this chords incident to the middle vertex may represent several chords
    multi_vertices = [vi[2] for vi in rc_vertices if len(vi) == 5]
    multi_chords = [chord for chord in chords if is_incident(chord, multi_vertices)]

    def multiplicity(chord):
        if chord in multi_chords:
            # C-edges have multiplicity <= k-q+1
            # if non C-multi-edges have multiplicity >= k then no C-edge runs over them.
            return k
        else:
            return 1

    prob, chord_vars = basic_ilp(k, chords, LpMinimize, multiplicity)

    BIG = len(vertices) ** 2

    for diagonal in D:
        prob += chord_vars[diagonal] == 1

    # chords crossing D, <= 12 for k=6
    C = {
        chord
        for diagonal in D
        for chord in crossings(diagonal, chords)
        if chord not in D
    }

    max_number_of_C_chords = q * (k - q + 1)

    Delta_C = lpSum([max_number_of_C_chords] + [-chord_vars[chord] for chord in C])

    Delta_i_vars = LpVariable.dicts("Delta_", range(2 * q), lowBound=0, cat=LpInteger)

    for i in range(2 * q):
        if rc_configuration[i] == 0:
            prob += Delta_i_vars[i] == 0
        elif rc_configuration[i] == 1:
            prob += (
                Delta_i_vars[i]
                >= 3 * bound_linear
                - bound_constant
                - 2
                - chord_vars[inner_rc_chords[i][0]]
            )
        elif rc_configuration[i] == 2:
            prob += Delta_i_vars[i] >= 4 * bound_linear - bound_constant - 3 - lpSum(
                [chord_vars[chord] for chord in inner_rc_chords[i]]
            )
        elif rc_configuration[i] == 3:
            v_fo1, v_so1, v_to, v_so2, v_fo2 = rc_vertices[i]
            not_rc_i_chords = [
                chord for chord in chords if chord not in inner_rc_chords
            ]
            # splitting C-edge
            for splitting_edge in incident_edges(v_to, C):
                prob += (
                    Delta_i_vars[i]
                    >= -bound_linear
                    + bound_constant
                    - k
                    + chord_vars[exists_var_name(splitting_edge)] * BIG  # l start
                    - BIG
                    + lpSum(
                        [
                            chord_vars[crossing_chord]
                            for crossing_chord in crossings(
                                splitting_edge,
                                not_rc_i_chords,
                            )
                        ]
                    )  # l end
                    - lpSum(
                        [
                            chord_vars[exists_var_name(subgraph_boundary)]
                            for subgraph_boundary in [
                                edge(v_fo1, v_to),
                                edge(v_fo2, v_to),
                            ]
                        ]
                    )
                    + 2
                )

            # enclosing C-edge
            second_order_vertices = [v_so1, v_so2]
            for iter in range(2):
                incident_sov = second_order_vertices[iter]
                if iter == 0:
                    not_enclosed_fo_vertex = v_fo2
                else:
                    not_enclosed_fo_vertex = v_fo1

                for enclosing_edge in incident_edges(incident_sov, C):
                    prob += (
                        Delta_i_vars[i]
                        >= bound_linear
                        - k
                        + chord_vars[enclosing_edge] * BIG  # l start
                        - BIG
                        + lpSum(
                            [
                                chord_vars[crossing_chord]
                                for crossing_chord in crossings(
                                    enclosing_edge,
                                    not_rc_i_chords,
                                )
                            ]
                        )  # l end
                        - chord_vars[edge(incident_sov, not_enclosed_fo_vertex)]
                    )

            # real coponent chord configuration
            prob += (
                Delta_i_vars[i]
                >= 2 * bound_linear
                - lpSum([chord_vars[inner_chord] for inner_chord in inner_rc_chords[i]])
                - 1
            )

    # target function
    prob += Delta_C + lpSum(
        [Delta_i_vars[i] for i in range(2 * q) if rc_configuration[i] > 0]
    )

    prob.solve()

    # analyse solution
    chords_in_solution = [
        chord
        # (*chord, {"multiplicity": chord_vars[chord].value(), "boundary": False})
        for chord in chords
        if chord_vars[chord].value() >= 1
    ]

    print(f"Case {rc_configuration} solved with value {prob.objective.value()}")
    print(f"Delta_min = {Delta_min}")
    print(f"Delta_C = {Delta_C.value()}")
    for i in range(number_of_real_components):
        print(f"Delta_i = {Delta_i_vars[i].value()}")

    return (
        Delta_min <= prob.objective.value(),
        prob.objective.value(),
        ConnectingSolution(
            k,
            2 * q,
            chords_in_solution,
            [v for v in vertices if v not in V_D],
            verbose=False,
        ),
        # RCSolution(
        #     k,
        #     2 * q,
        #     chords_in_solution,
        #     rc_vertices=[v for v in vertices if v not in V_D],
        # ),
    )


def q3_bound(
    k: int, bound: tuple[float, float], ignore_failures: bool = True
) -> tuple[bool, dict[tuple[int, ...], tuple[bool, int, RCSolution]]]:
    q = 3
    results = dict()
    success = True
    first_critical = None
    for real_component_positions in non_isomorphic_configurations(q, [0, 1, 2, 3]):
        analysis = analyse_real_component_configuration(
            k, q, real_component_positions, bound
        )
        results[tuple(real_component_positions)] = analysis
        success &= analysis[0]
        if not ignore_failures and not analysis[0]:
            print(f"""
================================================
       Result of automatic bound proving:
       ----------------------------------
Upper bound of {bound[0]}n - {bound[1]} could
NOT be proven for k={k} with q=3.
Critical configuration was {real_component_positions[0]}
with objective value {analysis[1]}.
================================================
""")
            return results
        if not analysis[0] and first_critical is None:
            first_critical = real_component_positions

    if success:
        print(
            f"""
              
================================================
       Result of automatic bound proving:
       ++++++++++++++++++++++++++++++++++
Upper bound of {bound[0]}n - {bound[1]} has been
successfully proven for k={k} with q=3.
================================================

"""
        )
    else:
        print(
            f"""
              
================================================
       Result of automatic bound proving:
       ----------------------------------
Upper bound of {bound[0]}n - {bound[1]} could
NOT be proven for k={k} with q=3.
Critical configuration was {first_critical[0]}
with objective value {results[first_critical][1]}.
================================================

"""
        )
    return results


# def case_1rc():
#     return analyse_real_components(k=6, position_of_real_components=(0,))


# def case_2rc_connected():
#     return analyse_real_components(k=6, position_of_real_components=(0, 1))


# def case_2rc_adjacent():
#     return analyse_real_components(k=6, position_of_real_components=(0, 2))


# case_1rc()
# case_2rc_connected()

# q3_bound(6, (4, 7))

uncrit = analyse_real_component_configuration(6, 3, [1, 0, 0, 0, 0, 0], (4, 7))
