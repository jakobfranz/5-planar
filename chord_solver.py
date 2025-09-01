from __future__ import annotations
from typing import Callable, TypeVar, Dict
from pulp import LpProblem, LpMaximize, LpInteger, pulp
import networkx as nx
import math

Vertex = float
Edge = tuple[Vertex, Vertex]


def cross(a: Edge, b: Edge):
    """Checks, if two edges cross.

    Assumes both edges to be tuples of the adjacent
    vertices with the lower numbered vertex coming first.

    Vertices are assumed to lie on a convex face in the order
    of their numbering.

    Args:
        a (Edge): an edge
        b (Edge): an edge

    Returns:
        True, if a and b cross
    """
    if a[0] not in b[:2] and a[1] not in b[:2]:
        return (a[0] < b[0] and b[0] < a[1] and a[1] < b[1]) or (
            b[0] < a[0] and a[0] < b[1] and b[1] < a[1]
        )
    return False


def crossings(edge: Edge, other_edges: list[Edge]) -> list[Edge]:
    """List all edges that are crossing a specific edge.

    Crossings are determined by the :func:`~cross` function.

    Args:
        edge (Edge):
            An edge
        other_edges (list[Edge]):
            A list of edges. May or may not contain `edge`.

    Returns:
        All edges in `other_edges`, that cross `edge`."""
    return [other_edge for other_edge in other_edges if cross(edge, other_edge)]


class ChordSolution:
    def __init__(
        self, k: int, n: int, chords: list[Edge], draw_outer_edges: bool = True
    ) -> None:
        self.k = k
        self.n = n
        self.chords = chords
        self.draw_outer_edges = draw_outer_edges

        # crossing number of chords in solution
        self.crossings = {chord: self.crossing_number(chord) for chord in self.chords}

        # number of chords in solution
        self.size = len(self.chords)

        # number of chords at each vertex
        self.vertex_chord_degree = [
            len([chord for chord in self.chords if vertex in chord])
            for vertex in range(self.n)
        ]

        # saturated chords, i.e. chords with k crossings
        self.saturated_chords = [
            chord for chord in self.chords if self.crossings[chord] >= self.k
        ]

        # vertices, that have no adjacent chord with k crossings and no more than k-1 adjacent chords
        self.crossable_vertices: list[Vertex] = [vertex for vertex in range(self.n)]
        for saturated_chord in self.saturated_chords:
            for vertex in saturated_chord:
                if vertex in self.crossable_vertices:
                    self.crossable_vertices.remove(vertex)

        for vertex in range(self.n):
            if (
                self.vertex_chord_degree[vertex] > k - 1
                and vertex in self.crossable_vertices
            ):
                self.crossable_vertices.remove(vertex)

        # edge density of the best (?) (infinite) tiling of this polygon
        self.edge_density = (self.size + self.n / 2) / (self.n / 2 - 1)

        self.description = (
            f"{self.k}-planar drawing of {self.n}-gon with {self.size} chords"
        )
        if type(self) is ChordSolution:
            print(self.description)

    def crossing_number(self, edge: Edge) -> int:
        return len(crossings(edge, self.chords))

    def graph_crossing_number(self) -> int:
        return sum([self.crossing_number(edge) for edge in self.chords]) / 2

    def connect(self, removable_chords: int = 0) -> tuple[ChordSolution, Edge]:
        """Finds the optimal outer edge, on which another polygon may dock."""
        # disect solution on saturated edges
        sections = [list(range(self.n))]
        for a, b in self.saturated_chords:
            for section in sections:
                if len(section) == 1:
                    sections.remove(section)
                    continue
                if a in section and b in section:
                    sections.remove(section)
                    if len(section) == self.n:
                        # complete circle
                        sections.append(list(range(a, b + 1)))
                        sections.append(list(range(b, self.n)) + list(range(a + 1)))
                    else:
                        ind_a = section.index(a)
                        ind_b = section.index(b)
                        if ind_a > ind_b:
                            # in case a is further into the segment than b (if section contains n-1 - 0 edge)
                            temp = ind_a
                            ind_a = ind_b
                            ind_b = temp
                        sections.append(section[: ind_a + 1])
                        sections.append(section[ind_a : ind_b + 1])
                        sections.append(section[ind_b:])
                else:
                    for vertex in (vertex for vertex in (a, b) if vertex in section):
                        # split at vertex
                        sections.remove(section)
                        sections.append(section[: vertex + 1])
                        sections.append(section[vertex:])
        max_new_edges = 0
        for section in sections:
            if len(section <= max_new_edges - 2):
                continue
            for outer_edge in (
                sorted(section[i], section[i + 1]) for i in range(len(section) - 1)
            ):
                possible_vertices = section.copy()
                possible_vertices.remove(outer_edge[0])
                possible_vertices.remove(outer_edge[1])
                for vertex in possible_vertices:
                    if (
                        self.crossing_number(sorted(vertex, outer_edge[0] + 0.5))
                        > self.k
                    ):
                        pass
        raise NotImplementedError

    def to_nx_graph(self) -> nx.Graph:
        G = nx.Graph()
        G.add_edges_from(self.chords)
        return G

    def _point_positions(self) -> Dict[Vertex, tuple[float, float]]:
        n = self.n
        point_positions = {
            i: (math.sin((2 * math.pi / n) * i), math.cos((2 * math.pi / n) * i))
            for i in range(n)
        }
        return point_positions

    def visualize(self) -> None:
        """Visualizes a solution using networkx."""
        G = self.to_nx_graph()
        # G.add_edges_from([(i, (i + 1) % self.n) for i in range(self.n)])

        point_positions = self._point_positions()

        nx.draw(G, pos=point_positions, with_labels=True)
        nx.draw_networkx_edges(
            G,
            point_positions,
            edgelist=self.saturated_chords,
            width=8,
            alpha=0.5,
            edge_color="tab:red",
        )
        if self.draw_outer_edges:
            nx.draw_networkx_edges(
                G,
                point_positions,
                edgelist=[(i, (i + 1) % self.n) for i in range(self.n)],
                alpha=0.8,
                edge_color="tab:grey",
            )
        nx.draw_networkx_nodes(
            G,
            point_positions,
            self.crossable_vertices,
            node_size=800,
            node_color="tab:green",
        )

    def __eq__(self, other: ChordSolution) -> bool:
        """Checks, if two solutions are congruent.

        I.e. rotated or mirrored.

        Args:
            other: solution to be compared to this one.

        Returns:
            True, if the solutions are congruent
        """
        if self.n != other.n or self.size != other.size:
            return False

        # check for rotation
        for offset in range(self.n):
            if all(
                [
                    tuple(sorted(((vertex + offset) % self.n for vertex in chord)))
                    in other.chords
                    for chord in self.chords
                ]
            ):
                return True

        # check for mirrored rotation
        for offset in range(self.n):
            if all(
                [
                    tuple(
                        sorted(((-start + offset) & self.n, (-end + offset) % self.n))
                    )
                    in other.chords
                    for (start, end) in self.chords
                ]
            ):
                return True

        return False


class ConnectingSolution(ChordSolution):
    def split_edges(
        edges: list[Edge],
        connector_vertices: list[Vertex],
        is_neighbor: Callable[[Vertex, Vertex], bool] = None,
    ) -> tuple[list[Edge], ...]:
        if is_neighbor is None:
            max_vertex = math.floor(max([max(edge[:2]) for edge in edges]))

            def is_neighbor(a, b):
                return abs(a - b) == 1 or abs(a - b) == max_vertex

        real_chords = [
            chord
            for chord in edges
            if chord[0] not in connector_vertices
            and chord[1] not in connector_vertices
            and not is_neighbor(chord[0], chord[1])
        ]
        outer_edges = [
            edge
            for edge in edges
            if edge[0] not in connector_vertices
            and edge[1] not in connector_vertices
            and is_neighbor(edge[0], edge[1])
        ]
        connecting_chords = [
            chord
            for chord in edges
            if (chord[0] in connector_vertices) ^ (chord[1] in connector_vertices)
        ]
        through_chords = [
            chord
            for chord in edges
            if chord[0] in connector_vertices and chord[1] in connector_vertices
        ]
        return outer_edges, real_chords, connecting_chords, through_chords

    def __init__(
        self,
        k,
        n,
        edges,
        connector_vertices,
        draw_outer_edges=False,
    ):
        (
            outer_edges,
            real_chords,
            connecting_chords,
            through_chords,
        ) = ConnectingSolution.split_edges(edges, connector_vertices)
        print(connector_vertices)
        print(ConnectingSolution.split_edges(edges, connector_vertices))
        super().__init__(k, n, real_chords, draw_outer_edges)

        self.connector_vertices = connector_vertices

        self.outer_edges = outer_edges
        self.connecting_chords = connecting_chords
        self.through_chords = through_chords

        self.weighted_edge_number = (
            self.size + (len(self.outer_edges) + len(connecting_chords)) / 2
        )

        self.edge_density = self.weighted_edge_number / (n / 2 - 1)

        self.description = f"{k}-planar Connecting Solution in an {n}-gon with {len(connecting_chords)} connections and {self.size} chords. Weighted edge number: {self.weighted_edge_number}"
        print(self.description)

    def to_nx_graph(self) -> nx.Graph:
        G = super().to_nx_graph()
        G.add_edges_from(self.outer_edges)
        G.add_edges_from(self.connecting_chords)
        G.add_edges_from(self.through_chords)
        return G

    def _point_positions(self) -> Dict[Vertex, tuple[float, float]]:
        point_positions = super()._point_positions()

        for connector_vertex in self.connector_vertices:
            low = math.floor(connector_vertex)
            offset = connector_vertex - low
            high = math.ceil(connector_vertex) % self.n

            low = point_positions[low]
            high = point_positions[high]

            point_positions[connector_vertex] = tuple(
                (offset * high[dim] + (1 - offset) * low[dim] for dim in range(2))
            )

        return point_positions

    def visualize(self):
        super().visualize()
        G = self.to_nx_graph()
        point_positions = self._point_positions()

        nx.draw_networkx_edges(
            G,
            point_positions,
            edgelist=self.connecting_chords,
            width=8,
            alpha=0.5,
            edge_color="tab:blue",
        )

        nx.draw_networkx_edges(
            G,
            point_positions,
            edgelist=self.through_chords,
            width=8,
            alpha=0.3,
            edge_color="tab:blue",
        )

        nx.draw_networkx_edges(
            G,
            point_positions,
            edgelist=self.outer_edges,
            alpha=0.8,
            edge_color="tab:grey",
        )

        nx.draw_networkx_nodes(
            G, point_positions, self.connector_vertices, node_color="#666666"
        )


class ChordSolver:
    def crossable_vertex_predicate(
        max_iterations: int,
    ) -> tuple[Callable[[ChordSolution, int], tuple[bool, int]], int]:
        def _predicate(solution: ChordSolution, iteration: int) -> tuple[bool, int]:
            return (
                len(solution.crossable_vertices) > 0 or iteration > max_iterations,
                iteration + 1,
            )

        return (_predicate, 0)

    def __init__(self, k: int, n: int) -> None:
        """Initialise ChordSolver.

        Issues all chords of complete graph with n vertices.

        Args:
            k: maximum number of crossings per edge,
                as in k-planar
            n: size of polygon"""
        self.k = k
        self.n = n
        self.forbidden_chord_configurations: list[list[Edge]] = []

        # issue chords of complete graph K_n
        self.chords: list[Edge] = []
        for i in range(0, n - 2):
            for j in range(i + 2, n):
                self.chords.append((i, j))
        self.chords.remove((0, n - 1))

        self.max_crossings: dict[Edge, int] = {chord: k for chord in self.chords}
        self.max_chord_degree: dict[Vertex, int] = {}

    def remove_edge(self, edge: Edge) -> bool:
        """Remove a chord before solving.

        Args:
            edge: Edge to be removed

        Returns:
            True if edge existed and was removed
        """
        if edge in self.chords:
            self.chords.remove(edge)
            return True
        else:
            return False

    def set_crossable_vertex(self, vertex: Vertex, chord_degree: int = -1) -> None:
        """Require a vertex to be crossable.

        If called multiple times, each vertex will be crossable at least as
        many times as this function was called on it at the same time.

        Args:
            vertex: The vertex that should be crossable.
            chord_degree: Maximum number of adjacent chords.
                If chord_degree is -1, then no restriction will be made.
        """
        # Flag adjacent chords to be crossable
        for chord in self.chords:
            if vertex in chord:
                self.max_crossings[chord] -= 1

        # Flag vertex to have given chord degree
        if chord_degree >= 0:
            self.max_chord_degree[vertex] = chord_degree

    def solve(self) -> ChordSolution:
        """Solves this chord configuration via ILP.

        Returns:
            Solution as a dictionary
        """
        # Create ILP
        prob = LpProblem("Chords in Polygon", LpMaximize)
        n_squared = self.n * self.n
        x = pulp.LpVariable.dicts(
            "chords", self.chords, lowBound=0, upBound=1, cat=LpInteger
        )
        prob += pulp.lpSum([x[chord] for chord in self.chords])
        for chord in self.chords:
            prob += (
                pulp.lpSum(
                    [x[otherChord] for otherChord in crossings(chord, self.chords)]
                )
                + n_squared * x[chord]
                <= self.max_crossings[chord] + n_squared
            )

        # limit chord_degree, where required
        for vertex, max_chord_degree in self.max_chord_degree.items():
            prob += (
                pulp.lpSum([x[chord] for chord in self.chords if vertex in chord])
                <= max_chord_degree
            )

        # prohibit forbidden solutions
        for chord_configuration in self.forbidden_chord_configurations:
            if all((chord in self.chords for chord in chord_configuration)):
                # chords of forbidden configurations might have been removed
                prob += (
                    pulp.lpSum([x[chord] for chord in chord_configuration])
                    <= len(chord_configuration) - 1
                )

        # solve
        prob.solve()

        # analyse solution

        # chords in solution
        chords_in_solution = [chord for chord in self.chords if x[chord].value() == 1]

        solution = ChordSolution(self.k, self.n, chords_in_solution)
        return solution

    def remove_solution(
        self, solution: ChordSolution, remove_congruents: bool = True
    ) -> None:
        """Removes a optimal solution from the solution space of the ILP.

        Useful to get multiple solutions of the same problem.
        IMPORTANT: Assumes the solution to be removed to be optimal.
        Otherwise better solutions, that have the chords of the removed
        solution as their subset can't be reached as well.

        Args:
            solution: The solution dictionary of a solutions
                that should be excluded from the solution space
            remove_congruents:
                Congruent (i.e. rotated or mirrored) solutions
                will be removed as well
        """
        if not remove_congruents:
            self.forbidden_chord_configurations.append(solution.chords)
        else:
            # remove rotations
            print("Removing Congruents")
            for offset in range(solution.n):
                for mirrored in (1, -1):
                    self.forbidden_chord_configurations.append(
                        [
                            tuple(
                                sorted(
                                    (
                                        (mirrored * start + offset) % self.n,
                                        (mirrored * end + offset) % self.n,
                                    )
                                )
                            )
                            for (start, end) in solution.chords
                        ]
                    )

    T = TypeVar("T")

    def solve_until(
        self,
        predicate: Callable[[ChordSolution, T], tuple[bool, T]],
        aggregator: T = 0,
        remove_congruents: bool = False,
    ) -> tuple[ChordSolution, T]:
        """Solves this problem until a requirement is met.

        Repeatedly solves the ILP and removes the solution from the
        solution space, until the requirement given by ```predicate``` is met.

        Args:
            predicate: A function ```(solution, aggregator) -> (stop_solving, aggregator)```
                Iterating the solve-delete cycle stops, when ```stop_solving``` is true. ```predicate``` has
                the last calculated solution and an aggregator as arguments.
                The aggregator may be modified by ```predicate``` and
                will be used when ```predicate``` is called
                in the next iteration.
            aggregator: The starting value for the aggregator used
                by ```predicate```.
            remove_congruents: remove congruent solutions of the found solution
        Returns:
            Last calculated solution
            aggregator
        """
        stop_iterating = False
        while not stop_iterating:
            solution = self.solve()
            self.remove_solution(solution, remove_congruents)
            stop_iterating, aggregator = predicate(solution, aggregator)
        return solution, aggregator
