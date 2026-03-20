# SPGMR Reuse And Tuning Cloud Model Performance

Date: 2026-03-20

## Scope

- Branch purpose: reuse SPGMR preconditioner setup work when CVODE reports that the Jacobian can be reused and `gamma` has not materially changed
- `MXSTEP`: unchanged
- Benchmark model: single cloud model from `main_cloud.py` against `benchmark_outputs/stage1`
- Environment: `uclchem_cvode`

## Code Change

- Added preconditioner reuse state in `cvode_bridge.c`
- Reused the existing sparse Jacobian and preconditioner when `jok=.TRUE.` and `gamma` stayed within a configurable relative tolerance
- Added tuning env var: `UCLCHEM_PRECOND_GAMMA_REUSE_RTOL`

## Benchmark Results

### Stable Runs

- `spgmr + diag`:
  - command env: `UCLCHEM_CVODE_LINEAR_SOLVER=spgmr UCLCHEM_SPGMR_PRECONDITIONER=diag`
  - output: `/tmp/uclchem_spgmr_diag_reuse_single.pkl`
  - success flag: `0`
  - wall time: `168.24 s`

- `spgmr + klu` with default reuse tolerance:
  - command env: `UCLCHEM_CVODE_LINEAR_SOLVER=spgmr UCLCHEM_SPGMR_PRECONDITIONER=klu`
  - output: `/tmp/uclchem_spgmr_klu_reuse_single.pkl`
  - success flag: `0`
  - wall time: `187.36 s`

- `spgmr + klu` with `UCLCHEM_PRECOND_GAMMA_REUSE_RTOL=1e-8`:
  - command env: `UCLCHEM_CVODE_LINEAR_SOLVER=spgmr UCLCHEM_SPGMR_PRECONDITIONER=klu UCLCHEM_PRECOND_GAMMA_REUSE_RTOL=1e-8`
  - output: `/tmp/uclchem_spgmr_klu_reuse_rtol1e8_single.pkl`
  - success flag: `0`
  - wall time: `171.21 s`

- `spgmr + klu` with `UCLCHEM_PRECOND_GAMMA_REUSE_RTOL=1e-6`:
  - command env: `UCLCHEM_CVODE_LINEAR_SOLVER=spgmr UCLCHEM_SPGMR_PRECONDITIONER=klu UCLCHEM_PRECOND_GAMMA_REUSE_RTOL=1e-6`
  - output: `/tmp/uclchem_spgmr_klu_reuse_rtol1e6_single.pkl`
  - success flag: `0`
  - wall time: `175.91 s`

### Unsuccessful Or Rejected Runs

- `spgmr + diag + linear_solution_scaling=1`:
  - command env: `UCLCHEM_CVODE_LINEAR_SOLVER=spgmr UCLCHEM_SPGMR_PRECONDITIONER=diag UCLCHEM_CVODE_LINEAR_SOLUTION_SCALING=1`
  - output: `/tmp/uclchem_spgmr_diag_linscale_single.pkl`
  - success flag: `-6`
  - wall time: `75.45 s`
  - behavior: much faster, but it entered repeated `ISTATE -4` shortening and ended with `Ran out of timepoints in arrays`

- `dense`:
  - command env: `UCLCHEM_CVODE_LINEAR_SOLVER=dense`
  - output: no completed result file
  - wall time: aborted after roughly `165 s`
  - behavior: remained stuck in repeated `CV_TOO_MUCH_WORK` / `t + h = t` behavior late in the run, so it was not competitive enough to justify a hybrid startup implementation yet

## Interpretation

- The reuse change helps the KLU-preconditioner path materially:
  - previous branch without reuse: `222.07 s`
  - this branch with default reuse: `187.36 s`
  - best stable reuse setting tested: `171.21 s` with `UCLCHEM_PRECOND_GAMMA_REUSE_RTOL=1e-8`

- The current best stable configuration on this branch is still `spgmr + diag` at `168.24 s`.

- A looser `gamma` reuse threshold helps KLU until about `1e-8`; going looser to `1e-6` was slightly worse, but still stable.

- `linear_solution_scaling=1` is not acceptable in the current setup because it speeds the run up by failing early rather than by solving the chemistry correctly.

- Dense CVODE is not a promising next step for a hybrid startup at this stage because it already misses the time target badly before completing the same single-model run.
