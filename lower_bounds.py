from typing import Callable
from chords import number_of_chords
import os

FILE = os.path.basename(__file__)


def density_outerplanar(n, chords):
    """Calculate the density of outer-planar graphs made by tiling it with a smaller graph.

    Args:
        n (int):
            Size of the outer-planar graph tile
        chords (int):
            Number of chords in the outer-planar graph tile
    Returns
        Density of these tiled graphs with
        - bound_linear (float): Linear factor of density
        - bound_constant (float): Constant
        The density is ```bound_linear * n - bound_constant"""
    bound_linear = (chords + 1) / (n - 2) + 1
    bound_constant = 2 * bound_linear - 1
    return (bound_linear, bound_constant)


def density_h_framed(n, chords):
    """Calculate the density of k-planar h-framed graphs made by tiling it with a smaller graph.

    Args:
        n (int):
            Size of the outer-planar graph tile
        chords (int):
            Number of chords in the outer-planar graph tile
    Returns
        Density of these tiled graphs with
        - bound_linear (float): Linear factor of density
        - bound_constant (float): Constant
        The density is ```bound_linear * n - bound_constant"""
    bound_linear = (n + 2 * chords) / (n - 2)
    bound_constant = 2 * bound_linear
    return (bound_linear, bound_constant)


def lower_bound(
    k: int, max_n: int, density_function: Callable[[int, int], float]
) -> tuple[float, float]:
    chords = {n: number_of_chords(k, n) for n in range(3, max_n + 1)}
    density = {n: density_function(n, chords[n]) for n in chords.keys()}
    densest_n, highest_density = sorted(
        density.items(), key=lambda kv: kv[1], reverse=True
    )[0]
    print(f"Highest density at {highest_density}")
    print(
        f"Reached with optimal {k}-planar {densest_n}-gon with {chords[densest_n]} chords"
    )
    return highest_density
