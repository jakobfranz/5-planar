from chord_solver import ChordSolver, ChordSolution, Edge

# homotopic edges of K_12 in the top and bottom Dodekagon
# in the construction of 6.2n lower bound for simple
# framed graphs

top_homotopic_edges = [
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

bottom_homotopic_edges = [
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
    solver = ChordSolver(5, 12)
    solution1 = solver.solve()
    solver.remove_solution(solution1)
    solution2 = solver.solve()

    print(solution1.chords)
    print(solution2.chords)
    return [solution1, solution2]


def simple_dodecagon(homotopic_edges: list[Edge]):
    solver = ChordSolver(5, 12)
    for edge in homotopic_edges:
        solver.remove_edge(edge)
    solution = solver.solve()
    print(solution.chords)
    return solution


optimal_dodecagons()
simple_dodecagon(top_homotopic_edges)
simple_dodecagon(bottom_homotopic_edges)
