# SPGMR Cloud Model Performance

Date: 2026-03-20

## Configuration

- Solver backend: `UCLCHEM_CVODE_LINEAR_SOLVER=spgmr`
- Preconditioner: diagonal left preconditioner built from `I - gamma J`
- Jacobian-times-vector: exact sparse `J*v`
- `MXSTEP`: unchanged from the project default used by the cloud model

## Benchmark Command

```bash
env UCLCHEM_CVODE_LINEAR_SOLVER=spgmr /usr/bin/time -p \
  conda run --no-capture-output -n uclchem_cvode \
  python /Users/binjia/Desktop/UCLCHEM_CVODE_test/main_cloud.py \
    --solver cvode \
    --benchmark-root /Users/binjia/Desktop/UCLCHEM_CVODE_test/benchmark_outputs/stage1 \
    --output /tmp/uclchem_spgmr_single.pkl
```

## Result

- Model set size: 1 (`model_0`)
- Output file: `/tmp/uclchem_spgmr_single.pkl`
- Success flag: `0`
- Wall time: `162.94 s`

## Interpretation

- The model completes successfully with the SPGMR backend.
- This removes the earlier sparse-KLU behavior where CVODE stalled very near the startup transient.
- The current limitation is performance, not early convergence failure.
- Runtime is still above the 1 minute target because the cloud run repeatedly hits `CV_TOO_MUCH_WORK` and later emits many `t + h = t` warnings.
