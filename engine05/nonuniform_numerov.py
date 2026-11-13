"""Brunetti, Speciale and Rudan nonuniform-grid Numerov, 2021.

Journal of Computational Electronics 20, 1105–1113.
DOI 10.1007/s10825-021-01699-3.

This is a boundary-value system. It is not the 2005 variable-step
recurrence, and it is not a call to numerov_variable_step.
Equation solved: -u'' = c(x) u + s(x), with u known at the endpoints.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


class NonuniformNumerovError(RuntimeError):
    """Raised when the nonuniform-grid system cannot be solved."""


def _midpoint_weights(h: float, c_left: float, c_mid: float, c_right: float, s_left: float, s_mid: float, s_right: float) -> tuple[float, float, float]:
    """Uniform Numerov on the half-interval, solved for the midpoint value."""
    h2 = h * h
    denominator = 96.0 - 10.0 * h2 * c_mid
    if abs(denominator) < 1e-12:
        raise NonuniformNumerovError("Midpoint Numerov denominator is too small")
    a_left = (48.0 + h2 * c_left) / denominator
    a_right = (48.0 + h2 * c_right) / denominator
    b = (s_left + 10.0 * s_mid + s_right) / denominator
    return a_left, a_right, b


def assemble_nonuniform_system(nodes: np.ndarray, c_at, s_at, u_left: float, u_right: float) -> tuple[np.ndarray, np.ndarray]:
    """Tridiagonal system for interior nodes. The matrix is not forced to be symmetric."""
    x = np.asarray(nodes, dtype=float)
    if x.ndim != 1 or len(x) < 4 or np.any(np.diff(x) <= 0.0):
        raise ValueError("Need at least four strictly increasing nodes")
    if not all(math.isfinite(float(v)) for v in (u_left, u_right)):
        raise ValueError("Boundary values must be finite")
    n_internal = len(x) - 2
    matrix = np.zeros((n_internal, n_internal), dtype=float)
    rhs = np.zeros(n_internal, dtype=float)
    for row, i in enumerate(range(1, len(x) - 1)):
        h_i = float(x[i] - x[i - 1])
        h_next = float(x[i + 1] - x[i])
        left_mid = 0.5 * (x[i - 1] + x[i])
        right_mid = 0.5 * (x[i] + x[i + 1])
        c_im1, c_left_mid, c_i = float(c_at(x[i - 1])), float(c_at(left_mid)), float(c_at(x[i]))
        c_right_mid, c_ip1 = float(c_at(right_mid)), float(c_at(x[i + 1]))
        s_im1, s_left_mid, s_i = float(s_at(x[i - 1])), float(s_at(left_mid)), float(s_at(x[i]))
        s_right_mid, s_ip1 = float(s_at(right_mid)), float(s_at(x[i + 1]))
        a_prev, a_here_left, b_left = _midpoint_weights(h_i, c_im1, c_left_mid, c_i, s_im1, s_left_mid, s_i)
        a_here_right, a_next, b_right = _midpoint_weights(
            h_next, c_i, c_right_mid, c_ip1, s_i, s_right_mid, s_ip1
        )
        beta = 1.0 / h_i + a_prev * h_i * c_left_mid / 3.0
        gamma = 1.0 / h_next + a_next * h_next * c_right_mid / 3.0
        alpha = (
            1.0 / h_i
            + 1.0 / h_next
            - a_here_left * h_i * c_left_mid / 3.0
            - (h_i + h_next) * c_i / 6.0
            - a_here_right * h_next * c_right_mid / 3.0
        )
        sigma = h_i * s_left_mid / 3.0 + (h_i + h_next) * s_i / 6.0 + h_next * s_right_mid / 3.0
        eta = b_left * h_i**3 * c_left_mid / 3.0 + b_right * h_next**3 * c_right_mid / 3.0
        if row > 0:
            matrix[row, row - 1] = -beta
        else:
            rhs[row] += beta * u_left
        matrix[row, row] = alpha
        if row + 1 < n_internal:
            matrix[row, row + 1] = -gamma
        else:
            rhs[row] += gamma * u_right
        rhs[row] += sigma + eta
    if not np.all(np.isfinite(matrix)) or not np.all(np.isfinite(rhs)):
        raise NonuniformNumerovError("Non-finite matrix entries")
    return matrix, rhs


def solve_nonuniform_bvp(nodes: np.ndarray, c_at, s_at, u_left: float, u_right: float) -> np.ndarray:
    """Nodal solution of -u'' = c u + s on a strictly increasing mesh."""
    x = np.asarray(nodes, dtype=float)
    matrix, rhs = assemble_nonuniform_system(x, c_at, s_at, u_left, u_right)
    try:
        internal = np.linalg.solve(matrix, rhs)
    except np.linalg.LinAlgError as exc:
        raise NonuniformNumerovError("Nonuniform-grid matrix is singular") from exc
    if not np.all(np.isfinite(internal)):
        raise NonuniformNumerovError("Non-finite nodal solution")
    solution = np.empty(len(x), dtype=float)
    solution[0] = u_left
    solution[-1] = u_right
    solution[1:-1] = internal
    return solution


def brunetti_model(k: int = 2, p: int = 5):
    """Manufactured solution (20). c and s are the closed form of that solution.

    The extracted text of their (21) did not satisfy -u'' = c u + s, so those
    printed factors are not used.
    """
    if k <= 0 or p <= 0:
        raise ValueError("k and p must be positive")
    amplitude = k * math.pi * (1.0 + p)
    mu = p / amplitude

    def lam(x: float) -> float:
        return amplitude / (1.0 + p * x)

    def exact(x: float) -> float:
        value = lam(x)
        return value * math.sin(value)

    def coefficient(x: float) -> float:
        value = lam(x)
        return (mu * value) ** 2 * (value * value - 2.0)

    def source(x: float) -> float:
        value = lam(x)
        return -4.0 * (mu * value) ** 2 * value * value * math.cos(value)

    return exact, coefficient, source


def manufactured_identity_is_exact() -> bool:
    """Symbolic check that rebuilt c and s satisfy -u'' = c u + s for solution (20)."""
    import sympy as sp

    x = sp.symbols("x", real=True)
    amplitude = 2 * sp.pi * 6
    mu = 5 / amplitude
    lam = amplitude / (1 + 5 * x)
    u = lam * sp.sin(lam)
    coefficient = (mu * lam) ** 2 * (lam**2 - 2)
    source = -4 * (mu * lam) ** 2 * lam**2 * sp.cos(lam)
    residual = sp.simplify(-sp.diff(u, x, 2) - (coefficient * u + source))
    return residual == 0


def shifted_nodes(n_internal: int, p: int = 5) -> np.ndarray:
    """Left-shifted nodes from their transformation (24)."""
    if n_internal < 2:
        raise ValueError("Need at least two internal nodes")
    uniform = np.linspace(0.0, 1.0, n_internal + 2)
    inside = 1.0 + p * (p + 2.0) * (1.0 - uniform)
    shifted = (1.0 + p - np.sqrt(inside)) / p
    shifted[0] = 0.0
    shifted[-1] = 1.0
    if np.any(np.diff(shifted) <= 0.0):
        raise NonuniformNumerovError("Shifted nodes are not increasing")
    return shifted


@dataclass(frozen=True)
class ManufacturedRun:
    n_internal: int
    max_error: float
    finite_difference_max_error: float
    grid: str
    reduction_ratio: float | None
    formal_order: str = "NOT_ASSIGNED"


def finite_difference_model(nodes: np.ndarray, c_at, s_at) -> np.ndarray:
    """Equation (22): the paper's finite-difference comparison, not Numerov."""
    x = np.asarray(nodes, dtype=float)
    n_internal = len(x) - 2
    matrix = np.zeros((n_internal, n_internal), dtype=float)
    rhs = np.zeros(n_internal, dtype=float)
    for row, i in enumerate(range(1, len(x) - 1)):
        h_i = float(x[i] - x[i - 1])
        h_next = float(x[i + 1] - x[i])
        beta = 1.0 / h_i
        gamma = 1.0 / h_next
        alpha = beta + gamma - 0.5 * (h_i + h_next) * c_at(x[i])
        if row > 0:
            matrix[row, row - 1] = -beta
        matrix[row, row] = alpha
        if row + 1 < n_internal:
            matrix[row, row + 1] = -gamma
        rhs[row] = 0.5 * (h_i + h_next) * s_at(x[i])
    return np.linalg.solve(matrix, rhs)


def manufactured_convergence(n_values: tuple[int, ...] = (160, 320, 640)) -> list[ManufacturedRun]:
    exact, coefficient, source = brunetti_model()
    runs = []
    previous = None
    for n_internal in n_values:
        nodes = np.linspace(0.0, 1.0, n_internal + 2)
        solved = solve_nonuniform_bvp(nodes, coefficient, source, 0.0, 0.0)
        truth = np.array([exact(float(x)) for x in nodes])
        fd_internal = finite_difference_model(nodes, coefficient, source)
        fd_error = float(np.max(np.abs(fd_internal - truth[1:-1])))
        error = float(np.max(np.abs(solved - truth)))
        ratio = None if previous is None else previous / error
        runs.append(
            ManufacturedRun(
                n_internal=n_internal,
                max_error=error,
                finite_difference_max_error=fd_error,
                grid="uniform",
                reduction_ratio=ratio,
            )
        )
        previous = error
    return runs
