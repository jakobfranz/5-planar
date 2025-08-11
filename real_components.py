# Solve some cases when partitioning an outer-k-planar graph
# into a H-Block and real components.

from pulp import LpMinimize, LpInteger, LpVariable, lpSum
from chord_solver import ConnectingSolution, Edge, crossings
from functional_chord_solver import issue_chords, basic_ilp

# from persistent_cache import persistent_cache
import os
import math

FILE = os.path.basename(__file__)


def all_cases(options: list[int], length: int) -> list[list[int]]:
    if length == 1:
        return [[option] for option in options]
    else:
        remaining_options = all_cases(options, length - 1)
        return [
            [option] + remaining_option
            for option in options
            for remaining_option in remaining_options
        ]


def analyse_real_components(
    k: int = 6, position_of_real_components: tuple[int, ...] = (0,)
) -> dict[tuple[int, ...], tuple[int, ConnectingSolution]]:
    """Analyses missing edges in an optimal setting with the specified real components."""
    # for each real component G_i there are 3 options:
    # |V_i| is 3, 4 or >= 5   <=>   |V_i - V_C| = 1, 2 or >= 3
    # the latter is used to describe the option
    number_of_real_components = len(position_of_real_components)

    cases = [tuple(case) for case in all_cases([1, 2, 3], number_of_real_components)]

    V_C = list(range(6))
    C = [(V_C[i], V_C[i + 3]) for i in range(3)]
    results = dict()

    for case in cases:
        # vertices of real components (equally spaced between neighboring V_C vertices)
        V_i = [
            [
                position_of_real_components[i] + (h + 1) / (case[i] + 1)
                for h in range(case[i])
            ]
            for i in range(number_of_real_components)
        ]

        vertices = V_C + sum(V_i, [])
        vertices.sort()

        chords = issue_chords(vertices, False)
        # chords within a real component
        # except for shared edge between rc and V_C
        # inner_rc_chords: list[Edge] = sum(
        #     [
        #         issue_chords(vi + [math.floor(vi[0]), math.ceil(vi[0])])
        #         - [(math.floor(vi[0]), math.ceil(vi[0]))]
        #         for vi in V_i
        #     ],
        #     [],
        # )
        inner_rc_chords: list[Edge] = [
            issue_chords(
                sorted(vi + [math.floor(vi[0]), math.ceil(vi[0])]), issue_boundary=False
            )
            for vi in V_i
        ]
        # for i in range(number_of_real_components):
        #     inner_rc_chords[i].remove(
        #         tuple(
        #             sorted(
        #                 [
        #                     position_of_real_components[i],
        #                     (position_of_real_components[i] + 1) % 6,
        #                 ]
        #             )
        #         )
        #     )
        # inner_rc_chords_flat = [
        #     chord for real_component in inner_rc_chords for chord in real_component
        # ]

        # chords = [chord for chord in chords if chord not in inner_rc_chords]

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
                return None
            else:
                return 1

        prob, chord_vars = basic_ilp(k, chords, LpMinimize, multiplicity)

        BIG = (6 * 4) ** 2

        for diagonal in C:
            prob += chord_vars[diagonal] == 1

        # chords crossing C, <= 12 for k=6
        C_cr = {
            chord
            for diagonal in C
            for chord in crossings(diagonal, chords)
            if chord not in C
        }

        # chords incident to a real component vertex neighboring a V_C-vertex
        # if they are crossed 5 times in H, the real component has one chord less
        # than the optimum
        # encapsulated_vertices = {v for v in position_of_real_components}.union(
        #     {v + 1 for v in position_of_real_components}
        # )
        # encapsulated_vertices_vars = LpVariable.dicts(
        #     "encapsulated_vertex", encapsulated_vertices, lowBound=0, cat=LpInteger
        # )

        # for enc_vertex in encapsulated_vertices:
        #     # find neighbors of enc_vertex that are part of a real component
        #     idx = vertices.index(enc_vertex)
        #     neighbors = [
        #         neighbor
        #         for neighbor in [vertices[idx - 1], vertices[idx + 1]]
        #         if abs(enc_vertex - neighbor)
        #         < 1  # if difference is 1, it's another V_C vertex
        #     ]

        #     encapsulating_chords = [
        #         chord
        #         for chord in chords
        #         if (chord[0] in neighbors or chord[1] in neighbors)
        #         and not math.floor(chord[0])
        #         == math.floor(chord[1])  # remove chords within the same real component
        #     ]

        #     for enc_chord in encapsulating_chords:
        #         prob += (
        #             BIG * chord_vars[enc_chord]
        #             + lpSum(
        #                 [
        #                     chord_vars[crossing_chord]
        #                     for crossing_chord in crossings(enc_chord, chords)
        #                 ]
        #             )
        #             - encapsulated_vertices_vars[enc_vertex]
        #             <= BIG + k - 2
        #         )

        # a C_cr chord cutting a real component into two regions with >= 3 vertices each (case |V_i| == 5)
        # and is crossed l times
        # gives l-2 bonus edges
        # multi_vertex_vars = LpVariable.dicts(
        #     "cut_component", multi_vertices, lowBound=0, cat=LpInteger
        # )
        # for multi_vertex in multi_vertices:
        #     for cutting_chord in [chord for chord in chords if multi_vertex in chord]:
        #         prob += (
        #             BIG * chord_vars[cutting_chord]
        #             + lpSum(
        #                 [
        #                     chord_vars[crossing_chord]
        #                     for crossing_chord in crossings(cutting_chord, chords)
        #                 ]
        #             )
        #             - 2
        #             <= BIG + multi_vertex_vars[multi_vertex]
        #         )

        # edge between V_C and real components
        shared_boundary_edges = [
            tuple(sorted([v, v + 1 % 6])) for v in position_of_real_components
        ]

        # bonus edges for each real component
        # |E_i_inner|
        # |V_i| >= 5  ->  6 - |E_i_inner|
        # |V_i| == 4  ->  3 - |E_i_inner|
        # |V_i| == 3  ->  0
        base_rc_bonus = {1: 0, 2: 3, 3: 6}
        neg_rc_bonus = [
            [-base_rc_bonus[case[i]]]
            + [chord_vars[inner_chord] for inner_chord in inner_rc_chords[i]]
            for i in range(number_of_real_components)
        ]
        neg_rc_bonus_flat = [item for rc_list in neg_rc_bonus for item in rc_list]

        # target function
        prob += (
            12
            - lpSum([chord_vars[c_cr] for c_cr in C_cr])
            # + lpSum(
            #     [
            #         encapsulated_vertices_vars[enc_vertex]
            #         for enc_vertex in encapsulated_vertices
            #     ]
            # )
            # + lpSum(
            #     [multi_vertex_vars[multi_vertex] for multi_vertex in multi_vertices]
            # )
            + number_of_real_components
            - lpSum(chord_vars[shared_edge] for shared_edge in shared_boundary_edges)
            - lpSum(neg_rc_bonus_flat)
        )

        prob.solve()

        # analyse solution
        chords_in_solution = [
            chord for chord in chords if chord_vars[chord].value() >= 1
        ]

        print(f"Case {case} solved with value {prob.objective.value()}")

        results[case] = (
            prob.objective.value(),
            ConnectingSolution(
                k, 6, chords_in_solution, connector_vertices=sum(V_i, [])
            ),
        )

    critical_case = sorted(results.values(), key=lambda item: item[0])[0]
    results["critical"] = critical_case

    print(f"solved all cases with critical case having value {critical_case[0]}")

    return results


def case_1rc():
    return analyse_real_components(k=6, position_of_real_components=(0,))


def case_2rc_connected():
    return analyse_real_components(k=6, position_of_real_components=(0, 1))


def case_2rc_adjacent():
    return analyse_real_components(k=6, position_of_real_components=(0, 2))


case_2rc_connected()
