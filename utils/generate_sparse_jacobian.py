#!/usr/bin/env python3
"""Generate a sparse Jacobian Fortran module from jacobian.f90 assignments."""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path


ASSIGNMENT_RE = re.compile(r"^J\((\d+),(\d+)\)=(.*)$")
SUM_RE = re.compile(r"^SUM\(J\((\w+),(\d+)\)\)$")


def split_terms(expr: str) -> list[str]:
    expr = expr.strip()
    if not expr:
        return []
    if expr[0] not in "+-":
        expr = "+" + expr

    terms: list[str] = []
    current = [expr[0]]
    for idx, ch in enumerate(expr[1:], start=1):
        prev = expr[idx - 1]
        if ch in "+-" and prev not in "dDeE":
            terms.append("".join(current).strip())
            current = [ch]
        else:
            current.append(ch)
    terms.append("".join(current).strip())
    return terms


def wrap_values(values: list[str], indent: str, per_line: int = 8) -> list[str]:
    lines = []
    for idx in range(0, len(values), per_line):
        chunk = values[idx : idx + per_line]
        suffix = ", &" if idx + per_line < len(values) else " &"
        lines.append(f"{indent}" + ", ".join(chunk) + suffix)
    return lines


def parse_list_from_fortran(module_text: str, list_name: str) -> list[int]:
    marker = f":: {list_name} "
    start = module_text.find(marker)
    if start == -1:
        raise RuntimeError(f"Could not find {list_name} in network.f90")

    list_start = module_text.find("(/", start)
    list_end = module_text.find("/)", list_start)
    if list_start == -1 or list_end == -1:
        raise RuntimeError(f"Could not parse {list_name} values in network.f90")

    raw_values = module_text[list_start + 2 : list_end].replace("&", " ")
    return [int(value.strip()) for value in raw_values.split(",") if value.strip()]


def render_sum_assignment(index: int, dependency_indices: list[int]) -> list[str]:
    if not dependency_indices:
        return [f"    data({index}) = 0.0_dp"]

    terms = [f"data({dep_index})" for dep_index in dependency_indices]
    lines = [f"    data({index}) = &"]
    for pos in range(0, len(terms), 8):
        chunk = terms[pos : pos + 8]
        suffix = " + &" if pos + 8 < len(terms) else ""
        lines.append(f"      {' + '.join(chunk)}{suffix}")
    return lines


def render_assignment(index: int, expr: str) -> list[str]:
    terms = split_terms(expr)
    if not terms:
        return [f"    data({index}) = 0.0_dp"]

    lines = [f"    data({index}) = &"]
    for pos, term in enumerate(terms):
        suffix = " &" if pos < len(terms) - 1 else ""
        lines.append(f"      {term}{suffix}")
    return lines


def main():
    repo_root = Path(__file__).resolve().parents[1]
    input_path = repo_root / "src" / "fortran_src" / "jacobian.f90"
    network_path = repo_root / "src" / "fortran_src" / "network.f90"
    output_path = repo_root / "src" / "fortran_src" / "jacobian_sparse.f90"

    entries: dict[tuple[int, int], str] = {}
    max_row = 0
    max_col = 0

    for line in input_path.read_text().splitlines():
        match = ASSIGNMENT_RE.match(line.strip())
        if not match:
            continue
        row = int(match.group(1))
        col = int(match.group(2))
        expr = match.group(3).strip()
        entries[(row, col)] = expr
        max_row = max(max_row, row)
        max_col = max(max_col, col)

    if max_row != max_col:
        raise RuntimeError(f"Jacobian is not square: rows={max_row}, cols={max_col}")

    network_text = network_path.read_text()
    list_values = {
        "surfaceList": parse_list_from_fortran(network_text, "surfaceList"),
        "bulkList": parse_list_from_fortran(network_text, "bulkList"),
    }

    by_col: dict[int, list[tuple[int, str]]] = defaultdict(list)
    for (row, col), expr in entries.items():
        by_col[col].append((row, expr))
    for col in by_col:
        by_col[col].sort()

    row_ind: list[int] = []
    col_ptr: list[int] = [0]
    ordered_exprs: list[str] = []
    ordered_keys: list[tuple[int, int]] = []
    for col in range(1, max_col + 1):
        column_entries = by_col.get(col, [])
        for row, expr in column_entries:
            row_ind.append(row - 1)
            ordered_exprs.append(expr)
            ordered_keys.append((row, col))
        col_ptr.append(len(row_ind))

    entry_index = {key: index for index, key in enumerate(ordered_keys, start=1)}
    aggregate_dependencies: dict[int, list[int]] = {}
    for index, ((row, col), expr) in enumerate(zip(ordered_keys, ordered_exprs), start=1):
        match = SUM_RE.fullmatch(expr)
        if not match:
            continue

        list_name = match.group(1)
        ref_col = int(match.group(2))
        if ref_col != col:
            raise RuntimeError(
                f"Unexpected dense Jacobian sum mismatch at row {row}, col {col}: {expr}"
            )
        if list_name not in list_values:
            raise RuntimeError(f"Unsupported dense Jacobian aggregate list {list_name}")

        dependency_indices = [
            entry_index[(dep_row, col)]
            for dep_row in list_values[list_name]
            if dep_row <= max_row and (dep_row, col) in entry_index
        ]
        aggregate_dependencies[index] = dependency_indices

    lines: list[str] = []
    lines.append("MODULE jacobian_sparse")
    lines.append("USE constants")
    lines.append("USE iso_c_binding, ONLY: c_int")
    lines.append("USE physicscore, ONLY: densdot")
    lines.append("IMPLICIT NONE")
    lines.append("")
    lines.append(f"INTEGER, PARAMETER :: jac_neq = {max_row}")
    lines.append(f"INTEGER, PARAMETER :: jac_nnz = {len(row_ind)}")
    lines.append(
        f"INTEGER(c_int), PARAMETER :: jac_col_ptr(jac_neq + 1) = (/ &"
    )
    lines.extend(wrap_values([str(v) for v in col_ptr], "    "))
    lines.append("    /)")
    lines.append(
        f"INTEGER(c_int), PARAMETER :: jac_row_ind(jac_nnz) = (/ &"
    )
    lines.extend(wrap_values([str(v) for v in row_ind], "    "))
    lines.append("    /)")
    lines.append("")
    lines.append("CONTAINS")
    lines.append("")
    lines.append(
        "SUBROUTINE fill_sparse_jacobian_values(RATE, Y, safeMantle, D, bulkLayersReciprocal, totalSwap, data)"
    )
    lines.append(
        "    REAL(dp), INTENT(IN) :: RATE(:), Y(:), safeMantle, D, bulkLayersReciprocal, totalSwap"
    )
    lines.append("    REAL(dp), INTENT(OUT) :: data(jac_nnz)")
    lines.append("")
    for idx, expr in enumerate(ordered_exprs, start=1):
        if idx in aggregate_dependencies:
            lines.extend(render_sum_assignment(idx, aggregate_dependencies[idx]))
        else:
            lines.extend(render_assignment(idx, expr))
    lines.append("END SUBROUTINE fill_sparse_jacobian_values")
    lines.append("")
    lines.append("REAL(dp) FUNCTION ddensdensdot(D)")
    lines.append("    REAL(dp), INTENT(IN) :: D")
    lines.append("    REAL(dp) :: safeD")
    lines.append("    safeD = MAX(D, 1.0d-30)")
    lines.append("    ddensdensdot = (densdot(1.0001d0 * safeD) - densdot(safeD)) / (0.0001d0 * safeD)")
    lines.append("END FUNCTION ddensdensdot")
    lines.append("")
    lines.append("END MODULE jacobian_sparse")
    lines.append("")

    output_path.write_text("\n".join(lines))
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
