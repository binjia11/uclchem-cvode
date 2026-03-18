#!/usr/bin/env python3
"""Check that the generated sparse Jacobian pattern is self-consistent."""

from __future__ import annotations

import re
import sys
from pathlib import Path


def parse_integer_parameter(source: str, name: str) -> list[int]:
    match = re.search(
        rf"INTEGER,\s+PARAMETER\s+::\s+{name}\(\d+\)\s*=\s*\(/\s*(.*?)\s*/\)",
        source,
        re.DOTALL,
    )
    if match is None:
        raise RuntimeError(f"Could not find {name} in generated Jacobian module")
    values = match.group(1).replace("&", " ")
    return [int(value.strip()) for value in values.split(",") if value.strip()]


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]
    jacobian_path = repo_root / "src" / "fortran_src" / "jacobian.f90"
    source = jacobian_path.read_text()
    nnz_match = re.search(r"INTEGER,\s+PARAMETER\s+::\s+jac_nnz\s*=\s*(\d+)", source)
    if nnz_match is None:
        raise RuntimeError("Could not find jac_nnz in generated Jacobian module")
    jac_nnz = int(nnz_match.group(1))
    col_ptr = parse_integer_parameter(source, "jac_col_ptr")
    row_ind = parse_integer_parameter(source, "jac_row_ind")
    neq = len(col_ptr) - 1
    monotonic = all(col_ptr[idx] <= col_ptr[idx + 1] for idx in range(len(col_ptr) - 1))
    row_bounds_ok = all(0 <= row < neq for row in row_ind)
    length_ok = len(row_ind) == jac_nnz
    last_ptr_ok = col_ptr[-1] == jac_nnz

    print(f"neq:        {neq}")
    print(f"jac_nnz:    {jac_nnz}")
    print(f"row_ind:    {len(row_ind)}")
    print(f"last colptr:{col_ptr[-1]}")

    return 0 if monotonic and row_bounds_ok and length_ok and last_ptr_ok else 1


if __name__ == "__main__":
    sys.exit(main())
