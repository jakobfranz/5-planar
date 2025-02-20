from chord_solver import ChordSolver, ChordSolution, Edge

# homotopic edges of K_12 in the top and bottom Dodekagon
# in the construction of 6.2n lower bound for simple
# framed graphs

top_multiedges = [
    (0, 2),
    (1, 3),
    (5, 7),
    (6, 8),
    (8, 10),
    (9, 11),
    (0, 10),
    (8, 11),
    (0, 9),
    (0, 8),
]

bottom_multiedges = [
    (2, 4),
    (1, 3),
    (5, 7),
    (4, 6),
    (8, 10),
    (9, 11),
    (0, 10),
    (8, 11),
    (0, 9),
    (0, 8),
]


def optimal_dodecagons() -> list[ChordSolution]:
    # both optimal outer 5-planar graphs with n=12
    # only one will be used for simple graphs
    solver = ChordSolver(5, 12)
    solution1 = solver.solve()
    solver.remove_solution(solution1)
    solution2 = solver.solve()

    print(solution1.chords)
    print(solution2.chords)
    return [solution1, solution2]


def simple_dodecagon(multiedges: list[Edge]):
    solver = ChordSolver(5, 12)
    for edge in multiedges:
        solver.remove_edge(edge)
    solution = solver.solve()
    print(solution.chords)
    return solution


optimal_dodecagons()
simple_dodecagon(top_multiedges)
simple_dodecagon(bottom_multiedges)
