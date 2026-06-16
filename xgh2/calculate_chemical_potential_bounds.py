#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Calculate chemical-potential bounds from competing phase total energies.

Assumptions:
  - chemical/<formula>/scf/ contains REPORT and atom.config.
  - directory names are chemical formulas, such as O2, Si, SiO2, H2O.
  - atom.config first line starts with total atom count.
  - REPORT contains E_tot(eV) or an Etot,Ep,Ek line.

The target phase is taken from calculate/bulk/scf by default and normalized by
the target formula. Competing phases are taken from chemical/.
"""

from __future__ import annotations

import argparse
import itertools
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import yaml

if sys.platform == "win32":
    import io

    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")


NUMBER_RE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
FORMULA_RE = re.compile(r"([A-Z][a-z]?)(\d*)")
ATOM_COUNT_RE = re.compile(r"^\s*(\d+)\s*(?:atoms?)?")
ETOT_PATTERNS = [
    re.compile(rf"E_tot\s*\(eV\)\s*[:=]\s*({NUMBER_RE})", re.IGNORECASE),
    re.compile(rf"Etot\s*,\s*Ep\s*,\s*Ek\s*=\s*({NUMBER_RE})", re.IGNORECASE),
    re.compile(rf"\bEtot\b\s*[:=]\s*({NUMBER_RE})", re.IGNORECASE),
    re.compile(rf"Total\s+energy.*?[:=]\s*({NUMBER_RE})", re.IGNORECASE),
]


@dataclass
class Phase:
    name: str
    formula: str
    counts: Dict[str, int]
    report_path: Path
    atom_config_path: Path
    atoms_total: int
    atoms_per_formula: int
    n_formula: int
    e_tot_eV: float
    e_per_formula_eV: float
    role: str = "competing"
    dHf_eV: Optional[float] = None


@dataclass
class Constraint:
    kind: str
    phase: str
    formula: str
    coeffs: Dict[str, float]
    op: str
    rhs: float
    note: str = ""


@dataclass
class SolveResult:
    feasible: bool
    vertices: List[Dict[str, float]] = field(default_factory=list)
    bounds: Dict[str, Dict[str, Optional[float]]] = field(default_factory=dict)
    message: str = ""


def find_project_root(start: Optional[Path] = None) -> Path:
    current = (start or Path(__file__).resolve().parent).resolve()
    fallback = None
    while True:
        if (current / "calculate").is_dir() or (current / "chemical").is_dir():
            return current
        if fallback is None and ((current / "config.yaml").exists() or (current / "input").exists()):
            fallback = current
        parent = current.parent
        if parent == current:
            return fallback or (start or Path(__file__).resolve().parent).resolve()
        current = parent


def parse_formula(formula: str) -> Dict[str, int]:
    pos = 0
    counts: Dict[str, int] = {}
    for match in FORMULA_RE.finditer(formula):
        if match.start() != pos:
            raise ValueError(f"无法解析化学式: {formula}")
        elem, num = match.groups()
        counts[elem] = counts.get(elem, 0) + (int(num) if num else 1)
        pos = match.end()
    if pos != len(formula) or not counts:
        raise ValueError(f"无法解析化学式: {formula}")
    return counts


def formula_to_text(counts: Dict[str, int]) -> str:
    parts = []
    for elem in sorted(counts):
        n = counts[elem]
        parts.append(elem if n == 1 else f"{elem}{n}")
    return "".join(parts)


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="ignore")


def read_atom_count(path: Path) -> int:
    if not path.exists():
        raise FileNotFoundError(f"atom.config 不存在: {path}")
    first_line = path.open("r", encoding="utf-8", errors="ignore").readline()
    match = ATOM_COUNT_RE.search(first_line)
    if not match:
        raise ValueError(f"无法从首行读取原子数: {path}")
    return int(match.group(1))


def read_energy_from_text(text: str, source: Path) -> float:
    for pattern in ETOT_PATTERNS:
        match = pattern.search(text)
        if match:
            return float(match.group(1))
    raise ValueError(f"无法读取总能: {source}")


def read_total_energy(report_path: Path, fallback_config: Optional[Path] = None) -> float:
    if report_path.exists():
        try:
            return read_energy_from_text(read_text(report_path), report_path)
        except ValueError:
            pass
    if fallback_config and fallback_config.exists():
        return read_energy_from_text(read_text(fallback_config), fallback_config)
    raise FileNotFoundError(f"无法从 REPORT 或 atom.config 读取总能: {report_path}")


def phase_from_paths(
    name: str,
    formula: str,
    report_path: Path,
    atom_config_path: Path,
    role: str,
) -> Phase:
    counts = parse_formula(formula)
    atoms_per_formula = sum(counts.values())
    atoms_total = read_atom_count(atom_config_path)
    if atoms_total % atoms_per_formula != 0:
        raise ValueError(
            f"{name}: 原子总数 {atoms_total} 不能被 formula 原子数 {atoms_per_formula} 整除"
        )
    n_formula = atoms_total // atoms_per_formula
    e_tot = read_total_energy(report_path, atom_config_path)
    return Phase(
        name=name,
        formula=formula,
        counts=counts,
        report_path=report_path,
        atom_config_path=atom_config_path,
        atoms_total=atoms_total,
        atoms_per_formula=atoms_per_formula,
        n_formula=n_formula,
        e_tot_eV=e_tot,
        e_per_formula_eV=e_tot / n_formula,
        role=role,
    )


def scan_chemical_phases(chemical_dir: Path) -> List[Phase]:
    phases: List[Phase] = []
    if not chemical_dir.is_dir():
        raise FileNotFoundError(f"chemical 目录不存在: {chemical_dir}")
    for item in sorted(chemical_dir.iterdir(), key=lambda p: p.name):
        if not item.is_dir():
            continue
        formula = item.name
        scf_dir = item / "scf"
        report = scf_dir / "REPORT"
        atom_config = scf_dir / "atom.config"
        phases.append(phase_from_paths(formula, formula, report, atom_config, "competing"))
    return phases


def parse_ref_options(values: Iterable[str]) -> Dict[str, str]:
    refs: Dict[str, str] = {}
    for value in values:
        for part in str(value).split(","):
            part = part.strip()
            if not part:
                continue
            if "=" not in part:
                raise ValueError(f"--ref 格式应为 Element=Formula，收到: {part}")
            elem, formula = part.split("=", 1)
            refs[elem.strip()] = formula.strip()
    return refs


def default_ref_formula(element: str) -> str:
    molecular = {"H": "H2", "N": "N2", "O": "O2", "F": "F2", "Cl": "Cl2"}
    return molecular.get(element, element)


def choose_reference_phases(
    elements: Sequence[str],
    phases_by_formula: Dict[str, Phase],
    refs: Dict[str, str],
) -> Dict[str, Phase]:
    ref_phases: Dict[str, Phase] = {}
    for elem in elements:
        formula = refs.get(elem, default_ref_formula(elem))
        phase = phases_by_formula.get(formula)
        if phase is None:
            available = ", ".join(sorted(phases_by_formula))
            raise ValueError(f"找不到元素 {elem} 的参考相 {formula}。可用 chemical 相: {available}")
        if elem not in phase.counts:
            raise ValueError(f"参考相 {formula} 不含元素 {elem}")
        ref_phases[elem] = phase
    return ref_phases


def reference_mu0(ref_phases: Dict[str, Phase]) -> Dict[str, float]:
    mu0: Dict[str, float] = {}
    for elem, phase in ref_phases.items():
        mu0[elem] = phase.e_per_formula_eV / phase.counts[elem]
    return mu0


def formation_energy(phase: Phase, mu0: Dict[str, float]) -> float:
    missing = [elem for elem in phase.counts if elem not in mu0]
    if missing:
        raise ValueError(f"{phase.formula}: 缺少元素参考化学势: {', '.join(missing)}")
    return phase.e_per_formula_eV - sum(phase.counts[elem] * mu0[elem] for elem in phase.counts)


def coeff_vector(coeffs: Dict[str, float], elements: Sequence[str]) -> List[float]:
    return [float(coeffs.get(elem, 0.0)) for elem in elements]


def dot(a: Sequence[float], b: Sequence[float]) -> float:
    return sum(x * y for x, y in zip(a, b))


def solve_linear_system(matrix: List[List[float]], rhs: List[float], tol: float = 1e-12) -> Optional[List[float]]:
    n = len(rhs)
    aug = [row[:] + [rhs_i] for row, rhs_i in zip(matrix, rhs)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(aug[r][col]))
        if abs(aug[pivot][col]) < tol:
            return None
        if pivot != col:
            aug[col], aug[pivot] = aug[pivot], aug[col]
        pivot_val = aug[col][col]
        aug[col] = [x / pivot_val for x in aug[col]]
        for row in range(n):
            if row == col:
                continue
            factor = aug[row][col]
            if abs(factor) < tol:
                continue
            aug[row] = [x - factor * y for x, y in zip(aug[row], aug[col])]
    return [aug[i][-1] for i in range(n)]


def solve_bounds(
    elements: Sequence[str],
    equality: Constraint,
    inequalities: List[Constraint],
    tol: float = 1e-8,
) -> SolveResult:
    n = len(elements)
    eq_vec = coeff_vector(equality.coeffs, elements)
    eq_rhs = equality.rhs
    ineq_vecs = [(constraint, coeff_vector(constraint.coeffs, elements), constraint.rhs) for constraint in inequalities]

    vertices: List[List[float]] = []
    if n == 1:
        if abs(eq_vec[0]) < tol:
            return SolveResult(False, message="目标相等式无法约束单元素体系")
        candidate = [eq_rhs / eq_vec[0]]
        if all(dot(vec, candidate) <= rhs + tol for _, vec, rhs in ineq_vecs):
            vertices.append(candidate)
    else:
        for active in itertools.combinations(ineq_vecs, n - 1):
            matrix = [eq_vec] + [vec for _, vec, _ in active]
            rhs = [eq_rhs] + [rhs_i for _, _, rhs_i in active]
            candidate = solve_linear_system(matrix, rhs)
            if candidate is None:
                continue
            if abs(dot(eq_vec, candidate) - eq_rhs) > 1e-6:
                continue
            if all(dot(vec, candidate) <= rhs_i + 1e-7 for _, vec, rhs_i in ineq_vecs):
                if not any(max(abs(a - b) for a, b in zip(candidate, old)) < 1e-7 for old in vertices):
                    vertices.append(candidate)

    if not vertices:
        return SolveResult(False, message="No feasible chemical potential region.")

    vertex_dicts = [
        {elem: float(value) for elem, value in zip(elements, vertex)}
        for vertex in sorted(vertices, key=lambda row: tuple(row))
    ]
    bounds: Dict[str, Dict[str, Optional[float]]] = {}
    for idx, elem in enumerate(elements):
        vals = [vertex[idx] for vertex in vertices]
        bounds[elem] = {"dmu_min_eV": min(vals), "dmu_max_eV": max(vals)}
    return SolveResult(True, vertices=vertex_dicts, bounds=bounds)


def constraint_text(constraint: Constraint, elements: Sequence[str]) -> str:
    terms = []
    for elem in elements:
        coeff = constraint.coeffs.get(elem, 0.0)
        if abs(coeff) < 1e-14:
            continue
        if abs(coeff - 1.0) < 1e-14:
            terms.append(f"dmu_{elem}")
        else:
            terms.append(f"{coeff:g}*dmu_{elem}")
    lhs = " + ".join(terms) if terms else "0"
    return f"{lhs} {constraint.op} {constraint.rhs:.8f}"


def fmt(value: Optional[float]) -> str:
    if value is None:
        return "nan"
    return f"{float(value):.8f}"


def write_phase_log(phases: List[Phase], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        handle.write(
            "phase formula atoms_total atoms_per_formula n_formula "
            "E_tot_eV E_per_formula_eV dHf_eV role report atom_config\n"
        )
        for phase in phases:
            handle.write(
                f"{phase.name} {phase.formula} {phase.atoms_total} {phase.atoms_per_formula} "
                f"{phase.n_formula} {phase.e_tot_eV:.8f} {phase.e_per_formula_eV:.8f} "
                f"{fmt(phase.dHf_eV)} {phase.role} {phase.report_path} {phase.atom_config_path}\n"
            )


def write_constraints_log(
    constraints: List[Constraint],
    elements: Sequence[str],
    output: Path,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        handle.write("kind phase formula op rhs_eV expression note\n")
        for constraint in constraints:
            handle.write(
                f"{constraint.kind} {constraint.phase} {constraint.formula} {constraint.op} "
                f"{constraint.rhs:.8f} \"{constraint_text(constraint, elements)}\" {constraint.note}\n"
            )


def write_bounds_log(
    elements: Sequence[str],
    mu0: Dict[str, float],
    bounds: Dict[str, Dict[str, Optional[float]]],
    output: Path,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        handle.write("element mu0_eV dmu_min_eV dmu_max_eV mu_min_eV mu_max_eV\n")
        for elem in elements:
            b = bounds.get(elem, {})
            dmin = b.get("dmu_min_eV")
            dmax = b.get("dmu_max_eV")
            mu_min = mu0[elem] + dmin if dmin is not None else None
            mu_max = mu0[elem] + dmax if dmax is not None else None
            handle.write(
                f"{elem} {mu0[elem]:.8f} {fmt(dmin)} {fmt(dmax)} {fmt(mu_min)} {fmt(mu_max)}\n"
            )


def write_yaml_output(
    target: Phase,
    ref_phases: Dict[str, Phase],
    mu0: Dict[str, float],
    solve: SolveResult,
    constraints: List[Constraint],
    output: Path,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "target": {
            "formula": target.formula,
            "E_per_formula_eV": target.e_per_formula_eV,
            "dHf_eV": target.dHf_eV,
        },
        "references": {
            elem: {
                "formula": phase.formula,
                "mu0_eV": mu0[elem],
                "E_per_formula_eV": phase.e_per_formula_eV,
            }
            for elem, phase in ref_phases.items()
        },
        "feasible": solve.feasible,
        "bounds": {
            elem: {
                **vals,
                "mu_min_eV": mu0[elem] + vals["dmu_min_eV"] if vals.get("dmu_min_eV") is not None else None,
                "mu_max_eV": mu0[elem] + vals["dmu_max_eV"] if vals.get("dmu_max_eV") is not None else None,
            }
            for elem, vals in solve.bounds.items()
        },
        "vertices": solve.vertices,
        "constraints": [
            {
                "kind": c.kind,
                "phase": c.phase,
                "formula": c.formula,
                "op": c.op,
                "rhs_eV": c.rhs,
                "coeffs": c.coeffs,
                "note": c.note,
            }
            for c in constraints
        ],
    }
    with output.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(data, handle, allow_unicode=True, sort_keys=False)


def parse_args() -> argparse.Namespace:
    root = find_project_root()
    parser = argparse.ArgumentParser(
        description="根据 chemical/ 杂相总能计算化学势限制范围。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  python standalone/calculate_chemical_potential_bounds.py -i . -o result/chemical_potential_bounds.yaml --target SiO2
  python standalone/calculate_chemical_potential_bounds.py -i chemical -o bounds.yaml --target SiO2 --ref Si=Si --ref O=O2
        """,
    )
    parser.add_argument(
        "-i",
        "--input",
        default=None,
        help="输入项目根目录或 chemical 目录。若为 chemical 目录，默认项目根目录取其父目录。",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="主输出 YAML 文件。未显式指定 sidecar log 时，会按该文件名派生 phase/constraints/bounds log。",
    )
    parser.add_argument("--target", required=True, help="目标主相化学式，例如 SiO2。")
    parser.add_argument("--project-root", default=str(root), help="项目根目录。")
    parser.add_argument("--chemical-dir", default=None, help="chemical 目录，默认 <project-root>/chemical。")
    parser.add_argument("--target-report", default=None, help="目标主相 REPORT，默认 calculate/bulk/scf/REPORT。")
    parser.add_argument("--target-config", default=None, help="目标主相 atom.config，默认 calculate/bulk/scf/atom.config。")
    parser.add_argument(
        "--ref",
        action="append",
        default=[],
        help="元素参考相，格式 Element=Formula，例如 --ref O=O2。可多次使用或逗号分隔。",
    )
    parser.add_argument("--phase-log", default=None, help="phase energy log 输出路径。")
    parser.add_argument("--constraints-log", default=None, help="constraints log 输出路径。")
    parser.add_argument("--bounds-log", default=None, help="bounds log 输出路径。")
    parser.add_argument("--yaml", default=None, help="YAML 输出路径。")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_path = Path(args.input).expanduser().resolve() if args.input else None
    if input_path and input_path.name == "chemical":
        project_root = input_path.parent
        chemical_dir = input_path
    elif input_path:
        project_root = input_path
        chemical_dir = Path(args.chemical_dir).expanduser().resolve() if args.chemical_dir else project_root / "chemical"
    else:
        project_root = Path(args.project_root).expanduser().resolve()
        chemical_dir = Path(args.chemical_dir).expanduser().resolve() if args.chemical_dir else project_root / "chemical"
    result_dir = project_root / "result"

    yaml_candidate = args.output or args.yaml
    yaml_path = Path(yaml_candidate).expanduser().resolve() if yaml_candidate else result_dir / "chemical_potential_bounds.yaml"
    output_stem = yaml_path.with_suffix("")

    phase_log = (
        Path(args.phase_log).expanduser().resolve()
        if args.phase_log
        else (output_stem.parent / f"{output_stem.name}_phase_energies.log" if args.output else result_dir / "chemical_phase_energies.log")
    )
    constraints_log = (
        Path(args.constraints_log).expanduser().resolve()
        if args.constraints_log
        else (output_stem.parent / f"{output_stem.name}_constraints.log" if args.output else result_dir / "chemical_potential_constraints.log")
    )
    bounds_log = (
        Path(args.bounds_log).expanduser().resolve()
        if args.bounds_log
        else (output_stem.parent / f"{output_stem.name}_bounds.log" if args.output else result_dir / "chemical_potential_bounds.log")
    )

    target_formula = args.target
    target_counts = parse_formula(target_formula)
    elements = sorted(target_counts)

    chemical_phases = scan_chemical_phases(chemical_dir)
    phases_by_formula = {phase.formula: phase for phase in chemical_phases}
    ref_options = parse_ref_options(args.ref)
    ref_phases = choose_reference_phases(elements, phases_by_formula, ref_options)
    mu0 = reference_mu0(ref_phases)

    target_report = Path(args.target_report).expanduser().resolve() if args.target_report else project_root / "calculate" / "bulk" / "scf" / "REPORT"
    target_config = Path(args.target_config).expanduser().resolve() if args.target_config else project_root / "calculate" / "bulk" / "scf" / "atom.config"
    target_phase = phase_from_paths("target", target_formula, target_report, target_config, "target")
    target_phase.dHf_eV = formation_energy(target_phase, mu0)

    for phase in chemical_phases:
        phase_elements = set(phase.counts)
        if phase in ref_phases.values():
            ref_for = [elem for elem, ref_phase in ref_phases.items() if ref_phase is phase]
            phase.role = "reference:" + ",".join(ref_for)
        elif phase_elements.issubset(set(elements)):
            phase.role = "competing"
            phase.dHf_eV = formation_energy(phase, mu0)
        else:
            phase.role = "ignored"

    constraints: List[Constraint] = []
    target_constraint = Constraint(
        kind="target",
        phase=target_phase.name,
        formula=target_phase.formula,
        coeffs={elem: float(coeff) for elem, coeff in target_phase.counts.items()},
        op="=",
        rhs=float(target_phase.dHf_eV),
        note="target stability equality",
    )
    constraints.append(target_constraint)

    inequalities: List[Constraint] = []
    for elem in elements:
        c = Constraint(
            kind="element",
            phase=ref_phases[elem].name,
            formula=ref_phases[elem].formula,
            coeffs={elem: 1.0},
            op="<=",
            rhs=0.0,
            note="element-rich upper bound",
        )
        inequalities.append(c)
        constraints.append(c)

    for phase in chemical_phases:
        if phase.role != "competing":
            continue
        c = Constraint(
            kind="competing",
            phase=phase.name,
            formula=phase.formula,
            coeffs={elem: float(coeff) for elem, coeff in phase.counts.items()},
            op="<=",
            rhs=float(phase.dHf_eV),
            note="avoid competing phase precipitation",
        )
        inequalities.append(c)
        constraints.append(c)

    solve = solve_bounds(elements, target_constraint, inequalities)
    all_phase_rows = [target_phase] + chemical_phases

    write_phase_log(all_phase_rows, phase_log)
    write_constraints_log(constraints, elements, constraints_log)
    if solve.feasible:
        write_bounds_log(elements, mu0, solve.bounds, bounds_log)
    else:
        bounds_log.parent.mkdir(parents=True, exist_ok=True)
        bounds_log.write_text(f"No feasible chemical potential region.\n{solve.message}\n", encoding="utf-8")
    write_yaml_output(target_phase, ref_phases, mu0, solve, constraints, yaml_path)

    print("=" * 70)
    print("化学势限制范围计算")
    print("=" * 70)
    print(f"target:       {target_formula}")
    print(f"chemical dir: {chemical_dir}")
    print(f"elements:     {', '.join(elements)}")
    print(f"feasible:     {solve.feasible}")
    print(f"phase log:    {phase_log}")
    print(f"constraints:  {constraints_log}")
    print(f"bounds log:   {bounds_log}")
    print(f"yaml:         {yaml_path}")
    if solve.feasible:
        print()
        for elem in elements:
            b = solve.bounds[elem]
            print(
                f"{elem}: dmu=[{b['dmu_min_eV']:.8f}, {b['dmu_max_eV']:.8f}] eV, "
                f"mu=[{mu0[elem] + b['dmu_min_eV']:.8f}, {mu0[elem] + b['dmu_max_eV']:.8f}] eV"
            )
        return 0

    print(f"错误: {solve.message}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
