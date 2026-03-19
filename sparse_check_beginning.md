# Sparse CVODE Beginning Check

Date: 2026-03-19

## Scope

This note records what was verified about the CVODE sparse backend during the early-stage convergence investigation for the cloud benchmark model run from:

- `/Users/binjia/Desktop/UCLCHEM_CVODE_test/main_cloud.py`
- benchmark input: `benchmark_outputs/stage1/s1cv/abundance/3.0/m*.csv`
- environment: `uclchem_cvode`

The original question was whether the sparse Jacobian structure or sparse Jacobian callback was wrong, and if not, what is happening at the very beginning of the CVODE run.

## What Was Checked

### 1. Sparse solver wiring in the C bridge

Checked in:

- `src/fortran_src/cvode_bridge.c`

Verified:

- sparse path uses `SUNSparseMatrix(..., CSC_MAT, ...)`
- pattern is loaded once through `uclchem_sparse_jacobian_pattern`
- sparse linear solver is `SUNLinSol_KLU`
- sparse Jacobian callback is registered through `CVodeSetJacFn`
- dense path also exists and uses the same chemistry/Jacobian source for comparison

Conclusion:

- the high-level sparse CVODE wiring is correct

### 2. Sparse Jacobian generation and structure

Checked in:

- `src/fortran_src/jacobian.f90`
- `src/fortran_src/chemistry.f90`

Verified:

- sparse CSC pattern arrays are internally consistent
- no duplicate row indices within any column
- row indices are strictly sorted within each column
- on this SUNDIALS build, `sunindextype` is 32-bit, so the `int *` handoff used by the pattern callback is not the source of corruption

Conclusion:

- the sparse structure itself does not appear to be wrong

### 3. Sparse Jacobian values vs dense Jacobian values

Checked using the existing dense-vs-sparse debug path in:

- `src/fortran_src/chemistry.f90`

Observed in `cvode_solver.log`:

- `Jacobian compare call: 1..4`
- `mismatch_count: 0`
- `max_abs_diff: 0`
- `max_rel_diff: 0`

Conclusion:

- the sparse Jacobian values match the dense Jacobian exactly for the first callbacks
- the early failure is not caused by a wrong sparse Jacobian formula

## Reproduced Behavior

### Sparse CVODE

The sparse CVODE run fails almost immediately at the beginning.

Typical observed messages:

- `CVODE ERROR ... At t = 2.09e-08, mxstep steps taken before reaching tout`
- then repeated `CV_TOO_MUCH_WORK` returns
- then the outer target time is repeatedly reduced, but the solver still only crawls forward

Observed behavior:

- after repeated retries, the solution had only advanced to about `2.26e-07` while still hitting `mxstep = 10000`
- this is already far too slow for the required 1 minute runtime budget

### Dense CVODE

Dense CVODE behaves better than sparse in the startup region, but it is also not healthy enough yet.

Observed behavior:

- dense does not get trapped as early as sparse
- however, dense also still hits `CV_TOO_MUCH_WORK` and shows strong step-count pressure
- this means the sparse problem is not a pure sparse-structure bug; there is also a broader startup stiffness/scaling issue in the CVODE path

## What the Solver Diagnostics Suggest

From `cvode_c_debug.log` for the sparse run:

- `last_lin_flag = CVLS_SUCCESS`
- KLU reports successful factorization status
- structural rank is full
- however, there are many nonlinear convergence failures (`nncf`) and very large nonlinear iteration counts (`nni`)

Interpretation:

- KLU is not reporting an outright linear-solver failure
- CVODE is instead spending huge effort in Newton iteration at the beginning
- the startup Jacobian is numerically difficult, even though its sparse pattern and entries are correct

Important:

- the KLU diagnostics do not prove the startup matrix is numerically well-conditioned
- they only show that the factorization path did not fail in a simple structural sense

## Main Conclusion

The evidence from this check is:

1. The sparse Jacobian pattern is not the problem.
2. The sparse Jacobian values are not the problem.
3. The early-stage failure happens because the startup Newton solve with sparse KLU behaves much worse than dense LU on this very stiff, badly scaled initial state.

In short:

- sparse is structurally correct
- sparse is numerically much less robust than dense at the very beginning of this model

## Likely Physical/Numerical Reason

At the beginning of the run:

- many species are near the abundance floor
- a small number of species dominate
- reaction timescales differ strongly
- the Jacobian is therefore highly stiff and badly scaled

Dense LU appears to tolerate this startup transient better.
KLU preserves sparsity and is potentially faster later, but at startup it appears to be less robust for this problem.

## Code Changes Tried During This Check

### Kept

In `src/fortran_src/cvode_bridge.c`:

- removed the hard-coded `CVodeSetMaxOrd(..., 3)`
- replaced it with configurable behavior:
  - `UCLCHEM_CVODE_MAX_ORDER`
  - default now set to `5`

Reason:

- hard-limiting BDF to order 3 was unnecessarily restrictive
- this is a safe cleanup even though it did not solve the sparse startup problem

### Tried and Reverted

- preserving solver history across immediate `CV_TOO_MUCH_WORK` retries

Reason for revert:

- it did not materially improve the sparse early-stage behavior

## Practical Outcome

At the end of this check:

- sparse CVODE still does not converge acceptably at the beginning
- the problem is not in sparse Jacobian generation
- the bottleneck is the early sparse Newton/KLU behavior

## Most Promising Next Step

The most plausible next engineering approach is:

- start with dense CVODE during the initial transient
- switch to sparse KLU only after the solution leaves the pathological startup regime

That directly targets the observed weakness:

- dense is more robust at the beginning
- sparse may still become useful later when the system is less numerically hostile

