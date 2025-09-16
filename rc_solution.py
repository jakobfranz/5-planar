from __future__ import annotations
import math
from abc import ABC, abstractmethod

import shapely as shp
import shapely.plotting as spt
import matplotlib.pyplot as plt

from chord_solver import ChordSolution, Vertex, Edge


class Point:
    def __init__(self, x: float, y: float):
        self.x = x
        self.y = y
        self.point = shp.Point(self.x, self.y)

    def __add__(self, other: Point) -> Point:
        return Point(self.x + other.x, self.y + other.y)

    def __sub__(self, other: Point) -> Point:
        return Point(self.x - other.x, self.y - other.y)

    def __mul__(self, scalar: float) -> Point:
        return Point(self.x * scalar, self.y * scalar)

    def __truediv__(self, scalar: float) -> Point:
        return Point(self.x / scalar, self.y / scalar)

    def __eq__(self, other: Point) -> bool:
        accuracy = 1e-6
        dif = self - other
        return dif.x < accuracy and dif.y < accuracy

    def __repr__(self) -> tuple:
        return f"Point({self.x}, {self.y})"

    def normalize(self) -> Point:
        length = math.sqrt(self.x**2 + self.y**2)
        return self / length

    def plot(self, axes, color="green", markersize=20, label=""):
        spt.plot_points(self.point, ax=axes, color=color, markersize=markersize)
        axes.text(
            self.x,
            self.y,
            label,
            horizontalalignment="center",
            verticalalignment="center",
            weight="semibold",
        )


class Line(ABC):
    line_width = 1
    offset_length = 0.02
    multiedge_fan_size = 0.2

    @abstractmethod
    def intersections(self, other: Line) -> list[Point]:
        raise NotImplementedError

    @abstractmethod
    def plot(self, color, multiplicity: int):
        raise NotImplementedError

    def __and__(self, other: Line) -> list[Point]:
        """Alias to crosses."""
        return self.intersections(other)


class Straight(Line):
    def __init__(self, a: Point, b: Point):
        self.a = a
        self.b = b
        self.linestring = shp.LineString([a.point, b.point])

    def intersections(self, other: Line) -> list[Point]:
        if any([isinstance(other, typ) for typ in [Straight, PolyLine]]):
            intersection = self.linestring.intersection(other.linestring)
            if type(intersection) is shp.Point:
                return [Point(intersection.x, intersection.y)]
            else:
                return []
        else:
            raise NotImplementedError

    def plot(self, axes, color="black", multiplicity=1):
        PolyLine(self.a, self.b).plot(axes, color=color, multiplicity=multiplicity)


class PolyLine(Line):
    def __init__(self, *points):
        self.points = points
        self.linestring = shp.LineString([[p.x, p.y] for p in points])

    def intersections(self, other: Line) -> list[Point]:
        if any([isinstance(other, typ) for typ in [Straight, PolyLine]]):
            intersection = self.linestring.intersection(other.linestring)
            if type(intersection) is shp.Point:
                return [Point(intersection.x, intersection.y)]
            else:
                return []
        else:
            raise NotImplementedError

    def plot(self, axes, color="black", multiplicity=1):
        if multiplicity == 1:
            spt.plot_line(
                self.linestring,
                ax=axes,
                add_points=False,
                color=color,
                linewidth=Line.line_width,
            )
        else:
            diff = self.points[-1] - self.points[0]
            if diff.x == 0:
                offset_unit = Point(Line.offset_length, 0)
            else:
                offset_m = -1 / (diff.y / diff.x)
                offset_x = math.sqrt(Line.offset_length**2 / (offset_m**2 + 1))
                offset_unit = Point(offset_x, offset_x * offset_m)
            for i in range(multiplicity):
                offset = offset_unit * (i - (multiplicity - 1) / 2)
                control_points = [point + offset for point in self.points]
                control_points[0] += (
                    control_points[1] - control_points[0]
                ).normalize() * Line.multiedge_fan_size
                control_points[-1] += (
                    control_points[-2] - control_points[-1]
                ).normalize() * Line.multiedge_fan_size
                control_points = [self.points[0]] + control_points + [self.points[-1]]
                PolyLine(*control_points).plot(axes, color, 1)


class RCSolution(ChordSolution):
    layer_depth = 0.3
    max_bend = 0.8

    def __init__(
        self,
        k: int,
        vertices: list[Vertex],
        chords: list[Edge],
        deltas: dict[str | int, float],
    ):
        print(chords)
        self.k = k
        self.vertices = vertices
        self.chords = [tuple(chord) for chord in chords]
        self.deltas = deltas
        max_vertex = max(vertices)
        if max_vertex % 1 == 0:
            self.n = int(max_vertex) + 1
        else:
            self.n = math.ceil(max_vertex)
        print(self._vertex_points())

    def _circle_point_positions(self) -> dict[Vertex, Point]:
        return {
            vertex: Point(
                math.sin((2 * math.pi / self.n) * vertex),
                math.cos((2 * math.pi / self.n) * vertex),
            )
            for vertex in self.vertices
        }

    def _vertex_points(self) -> dict[Vertex, Point]:
        points = dict()
        order = 1
        for vertex in self.vertices:
            base = math.floor(vertex)
            partial = vertex % 1
            if partial == 0:
                order = 1
            elif partial > 0 and partial < 1 / 2:
                order = 2
                partial -= 0.1
            elif partial == 1 / 2:
                order += 0.5
            elif partial > 1 / 2 and partial < 1:
                order = 2
                partial += 0.1

            shifted_vertex = base + partial
            radial_factor = 1 + (order - 1) * RCSolution.layer_depth

            points[vertex] = (
                Point(
                    math.sin((2 * math.pi / self.n) * shifted_vertex),
                    math.cos((2 * math.pi / self.n) * shifted_vertex),
                )
                * radial_factor
            )
        return points

    def _control_points(self) -> dict[Edge, list[Point]]:
        circle_positions = self._circle_point_positions()
        vertex_positions = self._vertex_points()

        bendlines = []
        for i in range(self.n):
            a = circle_positions[i]
            c = circle_positions[(i + 1) % self.n]
            b = (a + c) * RCSolution.max_bend / 2
            bendlines.append(PolyLine(a, b, c))

        control_points = dict()
        for chord in self.chords:
            old_line = Straight(circle_positions[chord[0]], circle_positions[chord[1]])
            c_points = [vertex_positions[chord[0]]]
            c_points += [
                old_line.intersections(bendlines[math.floor(inc_vertex)])[0]
                for inc_vertex in chord[:2]
                # non-inner-rc-vertices dont get shifted -> no bend necessary for it
                if inc_vertex % 1 != 0
                # check for inner-rc-chords, they don't collide with a bendline
                and math.floor(chord[0]) != math.floor(chord[1])
                and math.ceil(chord[0]) != math.ceil(chord[1])
            ]
            c_points += [vertex_positions[chord[1]]]
            control_points[chord] = c_points

        return control_points

    def visualize(self):
        fig, ax = plt.subplots()
        ax.axis("off")
        ax.set_aspect("equal", adjustable="box")

        vertex_positions = self._vertex_points()
        chord_control_points = self._control_points()

        for vertex, point in vertex_positions.items():
            if vertex % 1 == 0:
                color = "green"
                label = vertex
            else:
                color = "lightblue"
                label = f"{vertex:2.1f}"
            point.plot(ax, color=color, label=label)

        for chord, control_points in chord_control_points.items():
            if len(chord) == 2:
                mult = 1
            else:
                mult = chord[2]
            PolyLine(*control_points).plot(ax, multiplicity=mult)

        delta_min = self.deltas[0]
        delta = self.deltas[1]
        delta_C = self.deltas[2]
        if delta >= delta_min:
            comp = "\\geq"
            color = "green"
        else:
            comp = "<"
            color = "red"

        real_components = list(
            {math.floor(vertex) for vertex in self.vertices if vertex % 1 != 0}
        )
        for i in range(len(real_components)):
            rc_vertices = [
                vertex
                for vertex in self.vertices
                if vertex > real_components[i] and vertex < real_components[i] + 1
            ]
            number_of_rc_vertices = len(rc_vertices)
            if number_of_rc_vertices % 2 == 0:
                delta_pos = (
                    vertex_positions[rc_vertices[0]] + vertex_positions[rc_vertices[1]]
                ) / 2
            else:
                delta_pos = vertex_positions[
                    rc_vertices[(number_of_rc_vertices - 1) // 2]
                ]
            delta_pos = delta_pos * 1.15
            ax.text(
                delta_pos.x,
                delta_pos.y,
                f"$\\Delta_{{{real_components[i]}}} \\geq {self.deltas[3 + real_components[i]]:.2g}$",
                va="center",
                ha="center",
            )

        # borders
        print(self.vertices)
        for i in range(len(self.vertices)):
            Straight(
                vertex_positions[self.vertices[i]],
                vertex_positions[self.vertices[(i + 1) % len(self.vertices)]],
            ).plot(ax, "grey", 1)

        ax.set_title(
            f"Outer-{self.k}-planar configuration\n$\\Delta = {delta:.2g} {comp} {delta_min:.2g} = \\Delta_{{min}}, \\Delta_C = {delta_C:.2g}$",
            color=color,
            pad=15,
        )

        plt.subplots_adjust(top=0.85)
        plt.show()


test_solution = RCSolution(
    6,
    [0, 0.5, 1, 1.5, 2, 3],
    [[0, 1], [0, 1.5, 2], [0.5, 1.5, 3], [0.5, 3], [1.5, 3], [1, 2]],
    None,
)
print(str(test_solution._vertex_points()))

#
# test_solution.visualize()
