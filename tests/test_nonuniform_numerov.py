import numpy as np
import pytest

from engine05.nonuniform_numerov import (
    assemble_nonuniform_system,
    brunetti_model,
    manufactured_convergence,
    manufactured_identity_is_exact,
    shifted_nodes,
    solve_nonuniform_bvp,
)
from engine05.numerov_solver import numerov_forward
from engine05.ode_reference import integrate_second_order
from engine05.variable_numerov import numerov_variable_step


def test_rebuilt_source_satisfies_the_differential_equation():
    assert manufactured_identity_is_exact()


def test_varying_coefficient_keeps_the_matrix_asymmetric():
    nodes = np.linspace(0.0, 1.0, 10)
    varying, _rhs = assemble_nonuniform_system(nodes, lambda z: np.exp(8.0 * z), lambda _z: 0.0, 0.0, 0.0)
    constant, _rhs = assemble_nonuniform_system(nodes, lambda _z: 2.0, lambda _z: 0.0, 0.0, 0.0)
    assert np.max(np.abs(varying - varying.T)) > 1e-4
    assert np.max(np.abs(constant - constant.T)) < 1e-10


def test_error_reduction_is_recorded_without_an_order_claim():
    runs = manufactured_convergence((160, 320, 640))
    assert runs[0].reduction_ratio is None
    assert runs[1].reduction_ratio > 8.0
    assert runs[2].reduction_ratio > 8.0
    assert runs[2].max_error < runs[1].max_error
    assert all(run.formal_order == "NOT_ASSIGNED" for run in runs)
    assert runs[1].max_error < runs[1].finite_difference_max_error


def test_shifted_grid_matches_the_manufactured_solution():
    exact, coefficient, source = brunetti_model()
    nodes = shifted_nodes(320, p=5)
    solved = solve_nonuniform_bvp(nodes, coefficient, source, 0.0, 0.0)
    truth = np.array([exact(float(x)) for x in nodes])
    assert float(np.max(np.abs(solved - truth))) < 0.1


def test_dop853_and_fixed_numerov_agree_on_the_oscillator():
    h = 0.01
    x = np.arange(0.0, 1.0 + 0.5 * h, h)
    reference = integrate_second_order(lambda _z: 1.0, (0.0, 1.0), [0.0, 1.0], x)
    values = np.zeros_like(x)
    values[1] = np.sin(h)
    for n in range(1, len(x) - 1):
        factor = 1.0 + h * h / 12.0
        values[n + 1] = numerov_forward(factor, factor, factor, values[n - 1], values[n])
    assert np.max(np.abs(reference.y[0] - np.sin(x))) < 1e-10
    assert np.max(np.abs(values - np.sin(x))) < 1e-8


def test_variable_step_propagation_tracks_the_oscillator():
    h = 0.01
    c = 1.25
    xs = [0.0, h]
    while xs[-1] < 1.0:
        step = min(c * (xs[-1] - xs[-2]), 1.0 - xs[-1])
        if step <= 0.0:
            break
        xs.append(xs[-1] + step)
    values = [0.0, np.sin(xs[1])]
    for n in range(1, len(xs) - 1):
        values.append(
            numerov_variable_step(values[-2], values[-1], 1.0, 1.0, 1.0, xs[n] - xs[n - 1], xs[n + 1] - xs[n])
        )
    exact = np.sin(np.array(xs))
    assert np.max(np.abs(np.array(values) - exact)) < 1e-5
    with pytest.raises(ValueError):
        integrate_second_order(lambda _z: 1.0, (0.0, 1.0), [0.0, 1.0], np.array([0.2, 0.1]))
