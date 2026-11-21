"""Outer-boundary comparison. The matrix and the Numerov tail are not the same condition."""

from __future__ import annotations

from engine05.gaussian_matrix import ground_state
from engine05.numerov_solver import GaussianModel, solve_bound_state
from engine05.ode_matching import OdeMatchError, match_eigenvalue


def matrix_versus_numerov_domains(
    model: GaussianModel | None = None,
    *,
    h_over_a: float = 0.02,
    domains: tuple[float, ...] = (4.0, 6.0, 10.0),
    match_over_a: float = 0.17,
) -> dict:
    """Dirichlet matrix energies beside exponential-tail Numerov energies.

    The difference is not a pure boundary error. It still contains the
    finite-difference truncation.
    """
    model = model or GaussianModel()
    rows = []
    for domain in domains:
        matrix = ground_state(model, h_over_a=h_over_a, rmax_over_a=domain)
        numerov = solve_bound_state(model, h_over_a, domain, match_over_a)
        rows.append(
            {
                "rmax_over_a": domain,
                "h_over_a": h_over_a,
                "matrix_boundary": "u(0)=u(R)=0",
                "numerov_boundary": "u(0)=0 and a decaying exponential at R",
                "matrix_energy_eV": matrix.energy_eV,
                "matrix_nodes": matrix.nodes,
                "numerov_energy_eV": numerov.energy_eV,
                "numerov_nodes": numerov.nodes,
                "difference_eV": matrix.energy_eV - numerov.energy_eV,
                "like_for_like": False,
            }
        )
    matrix_values = [row["matrix_energy_eV"] for row in rows]
    numerov_values = [row["numerov_energy_eV"] for row in rows]
    return {
        "rows": rows,
        "matrix_domain_spread_eV": max(matrix_values) - min(matrix_values),
        "numerov_domain_spread_eV": max(numerov_values) - min(numerov_values),
        "cross_solver_difference_at_largest_domain_eV": rows[-1]["difference_eV"],
        "classification": "COMPUTED",
        "note": "Domain spread and cross-solver difference are separate. Neither is an experimental uncertainty.",
    }


def ode_domain_and_match_sensitivity(
    model: GaussianModel | None = None,
    *,
    domains: tuple[float, ...] = (4.0, 6.0, 8.0),
    matches: tuple[float, ...] = (0.17, 1.5),
) -> list[dict]:
    """DOP853 ground-state search. A nonzero-node root is a failure, not a result."""
    model = model or GaussianModel()
    rows = []
    for domain in domains:
        for match in matches:
            row = {
                "solver": "dop853-two-sided",
                "rmax_over_a": domain,
                "match_over_a": match,
                "boundary": "u(0)=0, u'(0)=1, inward u'(R)=-kappa u(R)",
            }
            if not match < domain:
                row.update({"status": "NOT_RUN", "reason": "MATCH_NOT_INSIDE_DOMAIN", "energy_eV": None, "residual": None})
                rows.append(row)
                continue
            try:
                matched = match_eigenvalue(
                    model,
                    rmax_over_a=domain,
                    match_over_a=match,
                    n_samples=16,
                    max_step=0.15,
                    residual_atol=1e-7,
                )
            except OdeMatchError as exc:
                row.update({"status": "FAILED", "reason": str(exc), "energy_eV": None, "residual": None, "nodes": None})
                rows.append(row)
                continue
            row.update(
                {
                    "status": "COMPUTED",
                    "reason": matched.termination_reason,
                    "energy_eV": matched.energy_eV,
                    "residual": matched.residual,
                    "nodes": 0,
                }
            )
            rows.append(row)
    return rows
