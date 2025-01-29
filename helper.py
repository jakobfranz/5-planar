# usefull functions when dealing with graph drawings


def density_formula(t: float, highest_cell: int = 6) -> str:
    """Gives the density formula for a parameter t"""
    params = [-(t - 1) / 4 * size + t for size in range(highest_cell + 1)]
    sum = " ".join([f"{param:+.4g}·C{i}" for i, param in list(enumerate(params))[3:]])
    return f"|E| ≤ {t}·|V| {sum} - |X| - {2 * t}"


def density_formula_t(coefficient: float, cell_size: int) -> float:
    return (coefficient - 1 / 4) / (1 - cell_size / 4)
