# UCLCHEM CVODE Migration and Installation Guide

## Scope
This document explains:

1. What was changed to migrate UCLCHEM from DVODE to CVODE.
2. How to install and run the CVODE-enabled build in `conda` environment `uclchem_cvode`.
3. How to verify Python and Fortran CLI workflows.
4. Current limitations and troubleshooting.

Date: 2026-02-17

---

## 1. Repository Workflow (as implemented)

The current workflow in this repository has two main execution paths:

1. Python-first path (primary/maintained):
   - Build via `meson-python` (`pip install .`).
   - Runtime via `uclchem.model.*` Python API (calls Fortran extension `uclchemwrap`).

2. Fortran CLI path (legacy):
   - Build from `src/fortran_src/Makefile`.
   - Run binary as `./uclchem <MODEL> <input_file> [model args]`.

The solver integration point is centralized in:

- `src/fortran_src/chemistry.f90`

That module previously called DVODE directly. Migration targeted this integration point while minimizing unrelated code changes.

---

## 2. What Was Changed

### 2.1 Solver migration (DVODE -> CVODE)

#### Added
- `src/fortran_src/cvode_bridge.c`

This file provides a C bridge that:

1. Creates SUNDIALS context.
2. Wraps Fortran state arrays as CVODE serial `N_Vector`.
3. Configures CVODE with BDF and dense linear solver.
4. Calls `CVode(...)`.
5. Returns CVODE status code to Fortran.

#### Updated
- `src/fortran_src/chemistry.f90`

Key updates:

1. Removed direct dependency on `DVODE_F90_M`.
2. Added ISO C bindings and C interface to `uclchem_cvode_integrate`.
3. Added CVODE status code mapping constants.
4. Replaced DVODE call in `integrateODESystem` with `uclchem_cvode_integrate(...)`.
5. Added C-compatible RHS callback (`cvode_rhs_c`) forwarding to existing Fortran `F(...)`.
6. Preserved original error-recovery behavior pattern (reduce `targetTime`, relax tolerance factor, unrecoverable error paths).

### 2.2 Python extension build/link changes

#### Updated
- `meson.build`

Key updates:

1. Added C compiler usage (`cc`).
2. Removed DVODE static library build/link.
3. Added SUNDIALS library resolution from `CONDA_PREFIX/lib`:
   - `sundials_cvode`
   - `sundials_nvecserial`
   - `sundials_sunlinsoldense`
   - `sundials_sunmatrixdense`
   - `sundials_generic`
4. Added SUNDIALS include path (`CONDA_PREFIX/include`) to extension compile flags.
5. Added `src/fortran_src/cvode_bridge.c` to extension sources.

### 2.3 Legacy Fortran Makefile changes

#### Updated
- `src/fortran_src/Makefile`

Key updates:

1. Added C compiler rules for `cvode_bridge.c`.
2. Added SUNDIALS include/lib flags using `SUNDIALS_PREFIX` (defaults to `CONDA_PREFIX`).
3. Removed DVODE object from link list.
4. Added SUNDIALS libs to `main` and `python` targets.

### 2.4 Fortran CLI compatibility update

#### Updated
- `src/fortran_src/main.f90`

Key updates:

1. Updated CLI calls for `CLOUD`, `HOTCORE`, `CSHOCK` to use current wrapper signatures via keyword args.
2. Fixed `nspec` import to `f2py_constants`.
3. Added explicit behavior for `POSTPROCESS` CLI case:
   - prints unsupported message
   - exits with code `2`
   - rationale: postprocess path requires trajectory arrays provided by Python API.

---

## 3. Installation in `uclchem_cvode`

## 3.1 Prerequisites

1. `conda` installed.
2. Existing env named `uclchem_cvode` (or create it with Python 3.12+).
3. Compiler toolchain available in env (`gfortran`, clang toolchain).

## 3.2 Install CVODE/SUNDIALS and build tooling

```bash
conda install -n uclchem_cvode -c conda-forge -y sundials
conda install -n uclchem_cvode -c conda-forge -y meson-python meson ninja
```

## 3.3 Install UCLCHEM package from source

Use no-build-isolation so install works offline against already-installed env packages:

```bash
conda run -n uclchem_cvode pip install --no-build-isolation --no-deps .
```

Notes:

1. If your environment defaults to user-site install and fails on permissions, run install in a writable env context.
2. `--no-build-isolation` avoids downloading build dependencies from network.

---

## 4. Build and Run Fortran CLI

## 4.1 Build binary

```bash
PATH=/Users/binjia/anaconda3/envs/uclchem_cvode/bin:$PATH \
CONDA_PREFIX=/Users/binjia/anaconda3/envs/uclchem_cvode \
make -C src/fortran_src main
```

This produces binary:

- `./uclchem`

## 4.2 Run examples (minimal)

Example minimal input dictionary file:

```json
{"finaltime": 10.0, "points": 1, "outputfile": "/tmp/uclchem_full.dat"}
```

Run:

```bash
DYLD_LIBRARY_PATH=/Users/binjia/anaconda3/envs/uclchem_cvode/lib \
./uclchem CLOUD /tmp/uclchem_cli_minimal.inp
```

Hot core:

```bash
DYLD_LIBRARY_PATH=/Users/binjia/anaconda3/envs/uclchem_cvode/lib \
./uclchem HOTCORE /tmp/uclchem_cli_minimal.inp 3 300.0
```

C-shock:

```bash
DYLD_LIBRARY_PATH=/Users/binjia/anaconda3/envs/uclchem_cvode/lib \
./uclchem CSHOCK /tmp/uclchem_cli_minimal.inp 20.0 0.01 10.0
```

POSTPROCESS (CLI):

```bash
DYLD_LIBRARY_PATH=/Users/binjia/anaconda3/envs/uclchem_cvode/lib \
./uclchem POSTPROCESS /tmp/uclchem_cli_minimal.inp
```

Expected: explicit unsupported message and exit code 2.

---

## 5. Python API Verification

Smoke test command:

```bash
MPLCONFIGDIR=/tmp/mpl_uclchem_cvode \
/Users/binjia/anaconda3/envs/uclchem_cvode/bin/python - << 'PY'
from uclchem import model
physics, chem, rates, start, flag = model.cloud(
    param_dict={'finalTime': 1.0e1, 'points': 1},
    return_array=True,
    timepoints=50,
)
print("flag", flag)
print("physics_shape", physics.shape)
print("chem_shape", chem.shape)
PY
```

Expected:

1. `flag 0`
2. sensible non-empty output shapes.

---

## 6. Files Modified

1. `meson.build`
2. `src/fortran_src/Makefile`
3. `src/fortran_src/chemistry.f90`
4. `src/fortran_src/main.f90`
5. `src/fortran_src/cvode_bridge.c` (new)

---

## 7. Known Limitations

1. CLI `POSTPROCESS` remains unsupported in `main.f90` by design, because current `postprocess` wrapper requires trajectory arrays.
2. If macOS runtime cannot locate SUNDIALS shared libs for `./uclchem`, set:
   - `DYLD_LIBRARY_PATH=$CONDA_PREFIX/lib`
3. Legacy warning in `wrap.f90` about label placement may still appear during build; it does not block successful build.

---

## 8. Quick End-to-End Commands

```bash
# 1) Install deps
conda install -n uclchem_cvode -c conda-forge -y sundials meson-python meson ninja

# 2) Install package
conda run -n uclchem_cvode pip install --no-build-isolation --no-deps .

# 3) Build CLI binary
PATH=/Users/binjia/anaconda3/envs/uclchem_cvode/bin:$PATH \
CONDA_PREFIX=/Users/binjia/anaconda3/envs/uclchem_cvode \
make -C src/fortran_src main

# 4) Run a minimal model
cat > /tmp/uclchem_cli_minimal.inp << 'EOF'
{"finaltime": 10.0, "points": 1, "outputfile": "/tmp/uclchem_full.dat"}
EOF

DYLD_LIBRARY_PATH=/Users/binjia/anaconda3/envs/uclchem_cvode/lib \
./uclchem CLOUD /tmp/uclchem_cli_minimal.inp
```

