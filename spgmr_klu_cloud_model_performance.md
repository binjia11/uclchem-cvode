# SPGMR KLU Cloud Model Performance

Date: 2026-03-20

## Configuration

- Solver backend: `UCLCHEM_CVODE_LINEAR_SOLVER=spgmr`
- Preconditioner: sparse KLU left preconditioner on `P = I - gamma J`
- Jacobian-times-vector: exact sparse `J*v`
- `MXSTEP`: unchanged from the project default used by the cloud model

## Benchmark Command

```bash
env UCLCHEM_CVODE_LINEAR_SOLVER=spgmr UCLCHEM_SPGMR_PRECONDITIONER=klu /usr/bin/time -p \
  conda run --no-capture-output -n uclchem_cvode \
  python /Users/binjia/Desktop/UCLCHEM_CVODE_test/main_cloud.py \
    --solver cvode \
    --benchmark-root /Users/binjia/Desktop/UCLCHEM_CVODE_test/benchmark_outputs/stage1 \
    --output /tmp/uclchem_spgmr_klu_single.pkl
```

## Result

- Model set size: 1 (`model_0`)
- Output file: `/tmp/uclchem_spgmr_klu_single.pkl`
- Success flag: `0`
- Wall time: `222.07 s`

## Interpretation

- The model completes successfully with the SPGMR backend and KLU preconditioner.
- This remains robust through the startup transient, so it does not reproduce the original sparse-KLU direct-solver failure mode.
- On this cloud model it is slower than the diagonal SPGMR preconditioner benchmark (`162.94 s`).
- Runtime is still above the 1 minute target because the run repeatedly hits `CV_TOO_MUCH_WORK` and later emits many `t + h = t` warnings.
