# Solve some cases when partitioning an outer-k-planar graph
# into a H-Block and real components.

from pulp import LpMinimize, LpInteger, LpVariable, lpSum
from chord_solver import Edge, crossings
from functional_chord_solver import (
    issue_chords,
    basic_ilp,
    is_incident,
    incident_edges,
    exists_var_name,
    edge,
)
from rc_solution import RCSolution

from persistent_cache import persistent_cache
import os

FILE = os.path.basename(__file__)


def non_isomorphic_configurations(
    q: int, rc_types: list[int] = [0, 1, 2, 3]
) -> list[tuple[int]]:
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
        tuple(case)
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


@persistent_cache((0, 3), FILE)
def analyse_real_component_configuration(
    k: int, q: int, rc_configuration: tuple[int], bound_to_show: tuple[float, float]
) -> tuple[bool, int, RCSolution]:
    """Try proving the upper bound for the given real component configuration.

    We fix `q` pairwise crossing edges, call them `D` and add
    real components to the boundary of this subgraph. Call the
    chords crossing `D` `C`. Calculate a lower bound on the
    missing edges in C and the real components compared to their
    theoretical maximum (Delta_C and Delta_i). If the sum of those
    is larger than Delta_min, the upper bound is shown for this
    configuration.

    Args:
        k (int):
            Maximum number of crossings per edge
        q (int):
            Assumption of the existence of `q` pairwise crossing
            edges, which are called `D` and are central to this
            analysis approach
        rc_configuration (list[int]):
            The types of real components in the order they appear
            around the D-region. Allowed types are 0,1,2 or 3 and
            stands for the number of additional vertices between
            two V_D vertices.
        bound_to_show (tuple[float, float]):
            Linear Upper bound that is to be shown. Takes the form
            (a,b) for a bound of an - b

    Returns:
        A tuple
        - success (bool): wether or not the boundary was shown
        - delta (int): lower bound for Delta in this config
        - rc_graph (RCSolution): The critical case for this config
    """
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

            # missing shared edge
            prob += Delta_i_vars[i] >= 1 - chord_vars[edge(v_fo1, v_fo2)]

    # target function
    prob += Delta_C + lpSum(
        [Delta_i_vars[i] for i in range(2 * q) if rc_configuration[i] > 0]
    )

    prob.solve()

    # analyse solution
    chords_in_solution = [
        (*chord, int(chord_vars[chord].value()))
        for chord in chords
        if chord_vars[chord].value() >= 1
    ]

    print(f"Case {rc_configuration} solved with value {prob.objective.value()}")
    print(f"Delta_min = {Delta_min}")
    print(f"Delta_C = {Delta_C.value()}")
    for i in range(2 * q):
        if rc_configuration[i] > 0:
            print(f"Delta_{i} = {Delta_i_vars[i].value()}")
    deltas = [Delta_min, prob.objective.value(), Delta_C.value()] + [
        delta_i.value() for delta_i in Delta_i_vars.values()
    ]

    return (
        Delta_min <= prob.objective.value(),
        prob.objective.value(),
        # ConnectingSolution(
        #     k,
        #     2 * q,
        #     chords_in_solution,
        #     [v for v in vertices if v not in V_D],
        #     verbose=False,
        # ),
        RCSolution(
            k,
            vertices,
            chords_in_solution,
            deltas,
        ),
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
