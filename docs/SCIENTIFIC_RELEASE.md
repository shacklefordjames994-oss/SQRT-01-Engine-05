# Scientific release

`SCIENTIFIC_RELEASE_STATUS = NOT_CERTIFIED`

This process does not sign an independent review. A passing test suite does not close that gate.

| Gate | Status | Evidence |
| --- | --- | --- |
| Clean checkout | NOT_RUN | The campaign was not tied to a commit at generation time |
| Environment | PASS | Python, NumPy, and SciPy imported and the campaign finished |
| Tests | PASS locally, not on GitHub CI | `66 passed in 6.77s` in `evidence/validation/pytest-followup.txt`. The repository reports no CI check runs. |
| Ground-state reproduction | PASS | Fixed-step `h/a = 0.0025` reproduced the repository archive |
| Solver agreement | PASS | DOP853 minus that Numerov energy is `3.0e-10 eV` at match `r/a = 0.17` |
| Wavefunction | PASS | Zero nodes, unit norm, matrix overlap `0.999998` |
| Uncertainty | PASS | Components are listed and not combined |
| Physical review | PASS | [PHYSICAL_MODEL_REVIEW.md](PHYSICAL_MODEL_REVIEW.md) |
| Reproducibility | IN_PROGRESS | Run JSON, manifest, and source hashes exist; checkout gate is still open |
| Independent review | NOT_RUN | No external reviewer |

The matrix solver uses `u(R) = 0`. Numerov and DOP853 use a decaying tail. At `h/a = 0.05`, changing the domain from `4a` to `8a` moves the matrix energy by about `2×10⁻¹⁵ eV` and does not move the Numerov energy. The cross-solver gap at `8a` is about `3.6×10⁻³ eV`. That gap is not a boundary-truncation effect, and it is not a like-for-like comparison.

DOP853 at match `r/a = 0.17` gives `-4.5741727208 eV` on domains `6a` and `8a`. The same search at match `r/a = 1.5` fails because it does not find a zero-node root. That failure is kept.

GitHub CI was not run for this branch. The local test log is not a substitute for a check run.

Historical energies in `evidence/SQRT-D02-GAUSS-001.json` were not overwritten.
