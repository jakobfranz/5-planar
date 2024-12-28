from chord_solver import ChordSolver, ChordSolution, Edge, Vertex
from persistent_cache import persistent_cache
import os.path

FILE = os.path.basename(__file__)


@persistent_cache(0, FILE)
def solve(k: int, n: int) -> ChordSolution:
    """Get a solution of a plain polygon"""
    return ChordSolver(k, n).solve()


@persistent_cache((), FILE)
def number_of_chords(k: int, n: int) -> int:
    """Get maximal number of chords in a n-gon"""
    return solve(k, n).size


@persistent_cache(0, FILE)
def all_solutions(k: int, n: int, chords=-1) -> list[ChordSolution]:
    """Calculate all unique solutions with given number of chords.

    If ```number_of_chords``` is ```-1```, get all solutions with maximal number of chords.
    """

    def chord_num_solver(
        solution: ChordSolution, found_solutions: list[ChordSolution]
    ) -> tuple[bool, list[ChordSolution]]:
        if solution.size > chords:
            # shouldn't happen
            return False, found_solutions
        if solution.size == chords:
            found_solutions.append(solution)
            return False, found_solutions
        return True, found_solutions

    solver = ChordSolver(k, n)
    max_chords = number_of_chords(k, n)
    if chords == -1:
        # get maximal number of chords
        chords = max_chords
        solver.remove_solution(solve(k, n))  # solve(k, n) is in cache
    else:
        for chord_number in range(chords + 1, max_chords + 1):
            for better_solution in all_solutions(k, n, chord_number):
                solver.remove_solution(better_solution, True)

    _, solutions = solver.solve_until(chord_num_solver, [], True)
    return solutions + [solve(k, n)]


@persistent_cache(0, FILE)
def solve_special(
    k: int,
    n: int,
    removed_edges: tuple[Edge] = (),
    vertex_incidence: tuple[tuple[Vertex, int]] = (),
) -> ChordSolution:
    pass


# s18 = solve(5, 18)

# for i in range(4, 18):
#     sol = all_solutions(5, i)


# def plot_density(max: int):
