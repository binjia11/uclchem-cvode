# CVODE Performance Improvement Options

This note summarizes practical ways to improve ODE solve performance in the current CVODE-based UCLCHEM fork.

## Current situation

The present CVODE bridge is minimal:

- It creates and destroys the CVODE solver, `SUNContext`, dense matrix, and dense linear solver on every integration call in `src/fortran_src/cvode_bridge.c`.
- It uses a dense linear solver: `SUNDenseMatrix` + `SUNLinSol_Dense`.
- It does not currently provide an analytic Jacobian to CVODE.

## Improvement options

### 1. Switch away from dense linear solves

For these chemical networks, the Jacobian is sparse. Dense `N x N` factorizations usually scale poorly for larger networks.

Best candidates:

- sparse direct solve with `SUNMatrix_Sparse` + `SUNLinSol_KLU`
- iterative solve with `SUNLinSol_SPGMR` plus a preconditioner

This is likely the highest-impact improvement for larger networks.

### 2. Reuse CVODE memory instead of rebuilding every solve

Each integration currently recreates:

- `SUNContext`
- `N_Vector`
- absolute-tolerance vector
- `CVode` memory
- linear solver
- matrix

This should be replaced with persistent solver state plus `CVodeReInit(...)` where possible.

### 3. Provide an analytic Jacobian or at least Jacobian sparsity

UCLCHEM already has Jacobian-related generation logic in Makerates and an old Jacobian stub in `src/fortran_src/chemistry.f90`.

Benefits:

- sparse direct solvers become viable
- Newton iterations get cheaper and more stable
- finite-difference Jacobian estimation cost is avoided

Even just supplying the sparsity structure would help.

### 4. Add a preconditioner for iterative linear solves

If the solver is changed to GMRES/SPGMR, performance will depend strongly on preconditioning.

Candidate approaches:

- diagonal preconditioner
- block-diagonal gas/surface/bulk preconditioner
- sparse ILU-style preconditioner if practical

### 5. Cache RHS work that does not need recalculation inside nonlinear iterations

Inside `F(...)`, some quantities are recomputed every RHS evaluation. Some are state-dependent and must stay, but others may be cached more aggressively per macro-step.

Potential targets:

- photorate-related bookkeeping
- surface/bulk bookkeeping where reuse is safe
- repeated scalar setup work

### 6. Expose and tune more CVODE controls

The current bridge only sets tolerances and `MaxNumSteps`.

Additional controls worth testing:

- initial step size
- min/max step size
- nonlinear iteration limits
- Jacobian/setup reuse behavior
- error test limits

### 7. Revisit tolerance strategy

Current defaults:

- `reltol = 1d-8`
- `abstol_factor = 1d-14`
- `abstol_min = 1d-25`

These may be stricter than necessary for some workloads. Careful relaxation can speed up solves, but this must be validated scientifically.

### 8. Reduce expensive retry behavior

Current failure handling often shrinks the target step by a factor of 10 and retries. That helps robustness, but it can also hide solver configuration issues and create many extra solves.

Once the linear algebra path is improved, this logic may need less intervention.

### 9. Exploit block structure in the network

The network naturally separates into gas, surface, and bulk species. Reordering species and designing sparse/preconditioner logic around that structure may help:

- reduce bandwidth
- reduce sparse fill-in
- improve block preconditioning

### 10. Profile before micro-optimizing chemistry kernels

Time is likely split between:

- RHS/rate evaluation
- Jacobian/linear solve
- solver setup overhead

For larger networks, dense linear algebra is likely the main bottleneck. For smaller networks, repeated solver setup may dominate.

## Recommended implementation order

1. Reuse CVODE objects across calls.
2. Move to sparse linear algebra, ideally KLU first.
3. Supply Jacobian values or at least sparsity structure.
4. If needed, try GMRES + preconditioner.
5. Then revisit tolerances and RHS micro-optimizations.
