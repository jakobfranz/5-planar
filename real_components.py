# Solve some cases when partitioning an outer-k-planar graph
# into a H-Block and real components.

from pulp import LpMinimize, LpInteger, LpVariable, lpSum
from chord_solver import ChordSolution, Edge, crossings
from functional_chord_solver import (
    issue_chords,
    basic_ilp,
    incident_edges,
    exists_var_name,
    edge,
)

from persistent_cache import persistent_cache
import os
import math
import itertools
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


def all_cases(options: list[int], length: int, filter: list[int]) -> list[list[int]]:
    if length == 1:
        return [[option] for option in options]
    else:
        remaining_options = all_cases(options, length - 1, filter)
        if length - 1 in filter and length - 1 != filter[0]:
            compare_to_index = filter[filter.index(length - 1) - 1]
            return [
                remaining_option + [option]
                for option in options
                for remaining_option in remaining_options
                if option >= remaining_option[compare_to_index]
            ]
        else:
            return [
                remaining_option + [option]
                for option in options
                for remaining_option in remaining_options
            ]


def congruent_filter(position_of_real_components: list[int], q: int) -> list[int]:
    number_of_real_components = len(position_of_real_components)

    # signature: distances to other real components
    rc_signature = [
        [
            (
                position_of_real_components[(i + ii) % number_of_real_components]
                - position_of_real_components[i]
            )
            % (2 * q)
            for ii in range(1, number_of_real_components)
        ]
        for i in range(number_of_real_components)
    ]

    congruent_rc = []

    def _choices():
        yield [position_of_real_components[0]]
        if number_of_real_components >= 2:
            for comb in itertools.combinations(position_of_real_components, 2):
                yield comb
        if number_of_real_components >= 3:
            for comb in itertools.combinations(position_of_real_components, 3):
                yield comb

    def _signature_eq(signature1, signature2) -> bool:
        return signature1 == signature2 or signature1 == list(
            reversed([2 * q - diff for diff in signature2])
        )

    for con_cand in _choices():
        if all(
            [
                _signature_eq(
                    rc_signature[position_of_real_components.index(con_cand[i])],
                    rc_signature[position_of_real_components.index(con_cand[i + 1])],
                )
                for i in range(len(con_cand) - 1)
            ]
        ):
            # normal signatures are equal
            # check candidate specific signature

            cand_len = len(con_cand)
            cand_signatures = [
                [
                    (con_cand[(i + ii) % cand_len] - con_cand[i]) % (2 * q)
                    for ii in range(1, cand_len)
                ]
                for i in range(cand_len)
            ]

            if all(
                [
                    _signature_eq(cand_signatures[i], cand_signatures[i + 1])
                    for i in range(len(con_cand) - 1)
                ]
            ):
                congruent_rc += [con_cand]

    congruent_rc.sort(key=lambda arr: len(arr), reverse=True)
    return [position_of_real_components.index(con_rc) for con_rc in congruent_rc[0]]


@persistent_cache((0,), FILE)
def analyse_real_components(
    k: int = 6,
    q: int = 3,
    position_of_real_components: tuple[int, ...] = (0,),
    bound_to_show: tuple[int, int] = (4, 7),
) -> dict[tuple[int, ...], tuple[int, RCSolution]]:
    """Analyses missing edges in an optimal setting with the specified real components."""
    # for each real component G_i there are 3 options:
    # |V_i| is 3, 4 or >= 5   <=>   |V_i - V_C| = 1, 2 or >= 3
    # the latter is used to describe the option
    number_of_real_components = len(position_of_real_components)

    bound_linear, bound_constant = bound_to_show

    # We only consider the case, if the rc types for the congruent
    # rc positions are order ascending
    rc_types_cases = [
        tuple(case)
        for case in all_cases(
            [1, 2, 3],
            number_of_real_components,
            congruent_filter(position_of_real_components, q),
        )
    ]

    V_D = list(range(2 * q))
    D = [(V_D[i], V_D[i + q]) for i in range(q)]
    results = dict()

    Delta_min = (
        number_of_real_components * (2 * bound_linear - bound_constant - 1)
        + q * (k - q - 2 * bound_linear + 4)
        + bound_constant
    )

    min_Delta = {
        1: 3 * bound_linear - bound_constant - 3,
        2: 4 * bound_linear - bound_constant - 6,
        3: 0,
    }

    for rc_types in rc_types_cases:
        # if the sum of the min_Deltas is greater than Delta_min
        # we can stop here
        min_delta_sum = sum(
            [min_Delta[real_component_type] for real_component_type in rc_types]
        )
        if min_delta_sum >= Delta_min:
            print(f"Delta >= {min_delta_sum} >= Delta_min")
            print(f"therefore upper bound is shown for configuration {rc_types}")
            results[rc_types] = (min_delta_sum, None)
            continue

        # vertices of real components (equally spaced between neighboring V_C vertices)
        V_i = [
            [
                position_of_real_components[i] + (h + 1) / (rc_types[i] + 1)
                for h in range(rc_types[i])
            ]
            for i in range(number_of_real_components)
        ]

        vertices = V_D + sum(V_i, [])
        vertices.sort()

        chords = issue_chords(vertices, False)
        # chords within a real component
        # including for shared edge between rc and V_C
        inner_rc_chords: list[Edge] = [
            issue_chords(
                sorted(vi + [math.floor(vi[0]), math.ceil(vi[0]) % (2 * q)]),
                issue_boundary=False,
            )
            + [edge(math.floor(vi[0]), math.ceil(vi[0]) % (2 * q))]
            for vi in V_i
        ]

        # in the case |V_i| >= 5 for any i, the middle vertex represents an arbitrary number of vertices in G_j.
        # As such, this chords incident to the middle vertex may represent several chords
        multi_vertices = [vi[1] for vi in V_i if len(vi) == 3]
        multi_chords = [
            chord
            for chord in chords
            if chord[0] in multi_vertices or chord[1] in multi_vertices
        ]

        def multiplicity(chord):
            if chord in multi_chords:
                # C-edges have multiplicity <= k-q+1
                # if non C-multi-edges have multiplicity >= k then no C-edge runs over them.
                return k
            else:
                return 1

        prob, chord_vars = basic_ilp(k, chords, LpMinimize, multiplicity)

        BIG = (2 * q * 4) ** 2

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

        Delta_i_vars = LpVariable.dicts(
            "Delta_", range(number_of_real_components), lowBound=0, cat=LpInteger
        )

        for i in range(number_of_real_components):
            if rc_types[i] == 1:
                prob += (
                    Delta_i_vars[i]
                    >= 3 * bound_linear
                    - bound_constant
                    - 2
                    - chord_vars[inner_rc_chords[i][0]]
                )
            elif rc_types[i] == 2:
                prob += Delta_i_vars[
                    i
                ] >= 4 * bound_linear - bound_constant - 3 - lpSum(
                    [chord_vars[chord] for chord in inner_rc_chords[i]]
                )
            elif rc_types[i] == 3:
                v_to = V_i[i][1]
                v_so1 = V_i[i][0]
                v_so2 = V_i[i][2]
                v_fo1 = math.floor(v_to)
                v_fo2 = math.ceil(v_to) % (2 * q)
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
                    - lpSum(
                        [chord_vars[inner_chord] for inner_chord in inner_rc_chords[i]]
                    )
                    - 1
                )

        # # edge between V_D and real components
        # shared_boundary_edges = [
        #     tuple(sorted([v, v + 1 % 6])) for v in position_of_real_components
        # ]

        # # bonus edges for each real component
        # # |E_i_inner|
        # # |V_i| >= 5  ->  6 - |E_i_inner|
        # # |V_i| == 4  ->  3 - |E_i_inner|
        # # |V_i| == 3  ->  0
        # base_rc_bonus = {
        #     1: 3 * bound_linear - bound_constant - 3,
        #     2: 4 * bound_linear - bound_constant - 4,
        #     3: 2 * bound_linear - 2,
        # }
        # neg_rc_bonus = [
        #     [-base_rc_bonus[case[i]]]
        #     + [chord_vars[inner_chord] for inner_chord in inner_rc_chords[i]]
        #     for i in range(number_of_real_components)
        # ]
        # neg_rc_bonus_flat = [item for rc_list in neg_rc_bonus for item in rc_list]

        # target function
        prob += Delta_C + lpSum(
            [Delta_i_vars[i] for i in range(number_of_real_components)]
        )

        prob.solve()

        # analyse solution
        chords_in_solution = [
            (*chord, {"multiplicity": chord_vars[chord].value(), "boundary": False})
            for chord in chords
            if chord_vars[chord].value() >= 1
        ]

        print(f"Case {rc_types} solved with value {prob.objective.value()}")
        print(f"Delta_min = {Delta_min}")
        print(f"Delta_C = {Delta_C.value()}")
        for i in range(number_of_real_components):
            print(f"Delta_i = {Delta_i_vars[i].value()}")

        results[rc_types] = (
            prob.objective.value(),
            RCSolution(k, 6, chords_in_solution, rc_vertices=sum(V_i, [])),
        )

    critical_case = sorted(results.values(), key=lambda item: item[0])[0]
    results["critical"] = critical_case

    results["delta_min"] = Delta_min

    print(f"solved all cases with critical case having value {critical_case[0]}")
    if critical_case[0] >= Delta_min:
        print(
            f"Therefore, upper bound is proven for real component configuration {position_of_real_components}."
        )
        results["success"] = True
    else:
        print(
            f"Upper bound could not be proven for real component configuration {position_of_real_components}."
        )
        results["success"] = False

    return results


def all_hexagon_position() -> list[tuple[int, ...]]:
    return [
        (0,),
        (0, 1),
        (0, 2),
        (0, 3),
        (0, 1, 2),
        (0, 1, 3),
        (0, 2, 4),
        (0, 1, 2, 3),
        (0, 1, 2, 4),
        (0, 1, 3, 4),
        (0, 1, 2, 3, 4),
        (0, 1, 2, 3, 4, 5),
    ]


def q3_bound(k: int, bound: tuple[float, float]):
    results = dict()
    for real_component_positions in all_hexagon_position():
        results[real_component_positions] = analyse_real_components(
            k,
            3,
            real_component_positions,
            bound,
        )
    success = all([results[rc_pos]["success"] for rc_pos in all_hexagon_position()])
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
        critical_configuration = sorted(
            list(results.items()), key=lambda kv_pair: kv_pair[1]["critical"][0]
        )[0]
        print(
            f"""
              
================================================
       Result of automatic bound proving:
       ----------------------------------
Upper bound of {bound[0]}n - {bound[1]} could
NOT be proven for k={k} with q=3.
Critical configuration was {critical_configuration[0]}
with objective value {critical_configuration[1]["critical"][0]}.
================================================

"""
        )
    return results


def case_1rc():
    return analyse_real_components(k=6, position_of_real_components=(0,))


def case_2rc_connected():
    return analyse_real_components(k=6, position_of_real_components=(0, 1))


def case_2rc_adjacent():
    return analyse_real_components(k=6, position_of_real_components=(0, 2))


# case_1rc()
# case_2rc_connected()

# q3_bound(6, (4, 7))

congruent_filter([0, 1], 3)
