#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Extract reorganization energies from explicitly specified input files.

Supported methods:
  structural_relaxation:
    S(0->q) = E_first - E_END from -i S1_RELAXSTEPS S2_RELAXSTEPS
    S(q->0) = E_first - E_END from -i S1_RELAXSTEPS S2_RELAXSTEPS

  fixed_charge_static:
    E0        = -i E0_REPORT E0_PRIME_REPORT EQ_REPORT EQ_PRIME_REPORT
    E0_prime  = -i E0_REPORT E0_PRIME_REPORT EQ_REPORT EQ_PRIME_REPORT
    Eq        = -i E0_REPORT E0_PRIME_REPORT EQ_REPORT EQ_PRIME_REPORT
    Eq_prime  = -i E0_REPORT E0_PRIME_REPORT EQ_REPORT EQ_PRIME_REPORT
    S(0->q) = E0_prime - Eq
    S(q->0) = Eq_prime - E0
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path
from typing import List, Optional, Tuple

if sys.platform == "win32":
    import io

    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")


NUMBER_RE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
ENERGY_RE = re.compile(rf"\bE\s*=\s*({NUMBER_RE})")
NUMBER_FIND_RE = re.compile(NUMBER_RE)
ETOT_EV_LABEL_RE = re.compile(r"E_tot\s*\(eV\)", re.IGNORECASE)
REPORT_ETOT_PATTERNS = [
    re.compile(rf"Etot\s*,\s*Ep\s*,\s*Ek\s*=\s*({NUMBER_RE})", re.IGNORECASE),
    re.compile(rf"\bEtot\b\s*[:=]\s*({NUMBER_RE})", re.IGNORECASE),
    re.compile(rf"Total\s+energy.*?[:=]\s*({NUMBER_RE})", re.IGNORECASE),
]


def charge_label(charge: int) -> str:
    return f"+{charge}" if charge > 0 else str(charge)


def defect_index(name: str) -> str:
    match = re.match(r"^(\d+)(?:_|$)", str(name))
    return match.group(1) if match else str(name)


def format_float(value: Optional[float]) -> str:
    if value is None or not math.isfinite(float(value)):
        return "nan"
    return f"{float(value):.8f}"


def tsv_value(value) -> str:
    if value is None:
        return "nan"
    return str(value).replace("\t", " ").replace("\r", " ").replace("\n", " ")


def read_relaxsteps_energy(path: Path, require_end: bool = True) -> dict:
    if not path.exists():
        return {
            "ok": False,
            "message": "missing",
            "first_eV": None,
            "end_eV": None,
            "lambda_eV": None,
            "end_found": False,
        }

    energies: List[Tuple[float, str]] = []
    end_energy: Optional[float] = None
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        match = ENERGY_RE.search(line)
        if not match:
            continue
        energy = float(match.group(1))
        energies.append((energy, line.strip()))
        if "*END" in line:
            end_energy = energy

    if not energies:
        return {
            "ok": False,
            "message": "no_energy",
            "first_eV": None,
            "end_eV": None,
            "lambda_eV": None,
            "end_found": False,
        }

    first_energy = energies[0][0]
    end_found = end_energy is not None
    if end_energy is None:
        if require_end:
            return {
                "ok": False,
                "message": "missing_END",
                "first_eV": first_energy,
                "end_eV": None,
                "lambda_eV": None,
                "end_found": False,
            }
        end_energy = energies[-1][0]

    return {
        "ok": True,
        "message": "ok" if end_found else "fallback_last_step",
        "first_eV": first_energy,
        "end_eV": end_energy,
        "lambda_eV": first_energy - end_energy,
        "end_found": end_found,
    }


def read_report_total_energy(path: Path) -> dict:
    if not path.exists():
        return {
            "ok": False,
            "message": "missing",
            "path": path,
            "E_eV": None,
            "line_number": None,
            "line_text": "",
            "pattern": "",
        }

    last_etot_ev: Optional[dict] = None
    fallback: Optional[dict] = None
    for line_number, line in enumerate(path.read_text(encoding="utf-8", errors="ignore").splitlines(), start=1):
        label_match = ETOT_EV_LABEL_RE.search(line)
        if label_match:
            tail = line[label_match.end():]
            number_match = NUMBER_FIND_RE.search(tail)
            if number_match:
                last_etot_ev = {
                    "ok": True,
                    "message": "ok",
                    "path": path,
                    "E_eV": float(number_match.group(0)),
                    "line_number": line_number,
                    "line_text": line.strip(),
                    "pattern": "E_tot(eV)",
                }

        for pattern in REPORT_ETOT_PATTERNS:
            match = pattern.search(line)
            if match:
                fallback = {
                    "ok": True,
                    "message": "ok:fallback",
                    "path": path,
                    "E_eV": float(match.group(1)),
                    "line_number": line_number,
                    "line_text": line.strip(),
                    "pattern": pattern.pattern,
                }
                break

    if last_etot_ev is not None:
        return last_etot_ev
    if fallback is not None:
        return fallback
    return {
        "ok": False,
        "message": "no_total_energy",
        "path": path,
        "E_eV": None,
        "line_number": None,
        "line_text": "",
        "pattern": "",
    }


def structural_relaxation_row(args: argparse.Namespace) -> dict:
    s1_path = Path(args.s1_relaxsteps).expanduser().resolve()
    s2_path = Path(args.s2_relaxsteps).expanduser().resolve()
    s1 = read_relaxsteps_energy(s1_path, require_end=not args.allow_unconverged)
    s2 = read_relaxsteps_energy(s2_path, require_end=not args.allow_unconverged)
    return {
        "defect": args.defect,
        "charge": charge_label(args.charge),
        "lambda_0_to_q_eV": s1["lambda_eV"],
        "lambda_q_to_0_eV": s2["lambda_eV"],
        "method": "structural_relaxation",
        "S1_source": s1_path,
        "S2_source": s2_path,
        "S1_first_eV": s1["first_eV"],
        "S1_end_eV": s1["end_eV"],
        "S2_first_eV": s2["first_eV"],
        "S2_end_eV": s2["end_eV"],
        "S1_status": s1["message"],
        "S2_status": s2["message"],
        "E0_eV": None,
        "E0_prime_eV": None,
        "Eq_eV": None,
        "Eq_prime_eV": None,
        "E0_report": "",
        "E0_prime_report": "",
        "Eq_report": "",
        "Eq_prime_report": "",
        "E0_line_number": "",
        "E0_prime_line_number": "",
        "Eq_line_number": "",
        "Eq_prime_line_number": "",
        "E0_status": "",
        "E0_prime_status": "",
        "Eq_status": "",
        "Eq_prime_status": "",
        "E0_line_text": "",
        "E0_prime_line_text": "",
        "Eq_line_text": "",
        "Eq_prime_line_text": "",
    }


def fixed_charge_static_row(args: argparse.Namespace) -> dict:
    e0_path = Path(args.E0).expanduser().resolve()
    e0_prime_path = Path(args.E0_prime).expanduser().resolve()
    eq_path = Path(args.Eq).expanduser().resolve()
    eq_prime_path = Path(args.Eq_prime).expanduser().resolve()

    e0 = read_report_total_energy(e0_path)
    e0_prime = read_report_total_energy(e0_prime_path)
    eq = read_report_total_energy(eq_path)
    eq_prime = read_report_total_energy(eq_prime_path)

    s_0_to_q = None
    s_q_to_0 = None
    s1_status = "ok"
    s2_status = "ok"

    if eq["E_eV"] is not None and e0_prime["E_eV"] is not None:
        s_0_to_q = float(e0_prime["E_eV"]) - float(eq["E_eV"])
    else:
        s1_status = f"{eq['message']}:{eq_path}" if eq["E_eV"] is None else f"{e0_prime['message']}:{e0_prime_path}"

    if e0["E_eV"] is not None and eq_prime["E_eV"] is not None:
        s_q_to_0 = float(eq_prime["E_eV"]) - float(e0["E_eV"])
    else:
        s2_status = f"{e0['message']}:{e0_path}" if e0["E_eV"] is None else f"{eq_prime['message']}:{eq_prime_path}"

    return {
        "defect": args.defect,
        "charge": charge_label(args.charge),
        "lambda_0_to_q_eV": s_0_to_q,
        "lambda_q_to_0_eV": s_q_to_0,
        "method": "fixed_charge_static",
        "S1_source": e0_prime_path,
        "S2_source": eq_prime_path,
        "S1_first_eV": None,
        "S1_end_eV": None,
        "S2_first_eV": None,
        "S2_end_eV": None,
        "S1_status": s1_status,
        "S2_status": s2_status,
        "E0_eV": e0["E_eV"],
        "E0_prime_eV": e0_prime["E_eV"],
        "Eq_eV": eq["E_eV"],
        "Eq_prime_eV": eq_prime["E_eV"],
        "E0_report": e0_path,
        "E0_prime_report": e0_prime_path,
        "Eq_report": eq_path,
        "Eq_prime_report": eq_prime_path,
        "E0_line_number": e0["line_number"],
        "E0_prime_line_number": e0_prime["line_number"],
        "Eq_line_number": eq["line_number"],
        "Eq_prime_line_number": eq_prime["line_number"],
        "E0_status": e0["message"],
        "E0_prime_status": e0_prime["message"],
        "Eq_status": eq["message"],
        "Eq_prime_status": eq_prime["message"],
        "E0_line_text": e0["line_text"],
        "E0_prime_line_text": e0_prime["line_text"],
        "Eq_line_text": eq["line_text"],
        "Eq_prime_line_text": eq_prime["line_text"],
    }


def write_summary_log(row: dict, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = ["idx", "S[0->+1]", "S[+1->0]", "S[0->-1]", "S[-1->0]"]
    charge = int(str(row.get("charge", "")).replace("+", ""))
    values = {
        "idx": defect_index(str(row.get("defect", ""))),
        "S[0->+1]": None,
        "S[+1->0]": None,
        "S[0->-1]": None,
        "S[-1->0]": None,
    }
    if charge == 1:
        values["S[0->+1]"] = row.get("lambda_0_to_q_eV")
        values["S[+1->0]"] = row.get("lambda_q_to_0_eV")
    elif charge == -1:
        values["S[0->-1]"] = row.get("lambda_0_to_q_eV")
        values["S[-1->0]"] = row.get("lambda_q_to_0_eV")

    with output.open("w", encoding="utf-8") as handle:
        handle.write(" ".join(fields) + "\n")
        handle.write(" ".join([str(values["idx"])] + [format_float(values[field]) for field in fields[1:]]) + "\n")


def write_detail_log(row: dict, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "defect",
        "charge",
        "method",
        "S0_to_q_eV",
        "Sq_to_0_eV",
        "E0_eV",
        "E0_prime_eV",
        "Eq_eV",
        "Eq_prime_eV",
        "E0_report",
        "E0_prime_report",
        "Eq_report",
        "Eq_prime_report",
        "S1_source",
        "S2_source",
        "S1_first_eV",
        "S1_end_eV",
        "S2_first_eV",
        "S2_end_eV",
        "E0_line_number",
        "E0_prime_line_number",
        "Eq_line_number",
        "Eq_prime_line_number",
        "S1_status",
        "S2_status",
        "E0_status",
        "E0_prime_status",
        "Eq_status",
        "Eq_prime_status",
        "E0_line_text",
        "E0_prime_line_text",
        "Eq_line_text",
        "Eq_prime_line_text",
    ]
    values = {
        **row,
        "S0_to_q_eV": format_float(row.get("lambda_0_to_q_eV")),
        "Sq_to_0_eV": format_float(row.get("lambda_q_to_0_eV")),
        "E0_eV": format_float(row.get("E0_eV")),
        "E0_prime_eV": format_float(row.get("E0_prime_eV")),
        "Eq_eV": format_float(row.get("Eq_eV")),
        "Eq_prime_eV": format_float(row.get("Eq_prime_eV")),
    }
    with output.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(fields) + "\n")
        handle.write("\t".join(tsv_value(values.get(field, "")) for field in fields) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="从显式指定的输入文件提取单个缺陷、单个电荷态的重组能。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # structural_relaxation：两个 RELAXSTEPS
  python standalone/extract_reorganization_energy.py --method structural_relaxation --defect 0_v_Si --charge 1 -i q1/relax/RELAXSTEPS q1/S21/RELAXSTEPS -o reorganization_energy.log

  # fixed_charge_static：四个 REPORT
  python standalone/extract_reorganization_energy.py --method fixed_charge_static --defect 0_v_Si --charge 1 -i q0/scf/REPORT q0/hse-meta-q1/REPORT q1/scf/REPORT q1/hse-meta-q0/REPORT -o reorganization_energy.log
        """,
    )
    parser.add_argument("--method", choices=["structural_relaxation", "fixed_charge_static"], required=True, help="重组能提取方法。")
    parser.add_argument("--defect", required=True, help="缺陷名称或编号，例如 0_v_Si 或 0。")
    parser.add_argument("--charge", required=True, type=int, help="目标电荷态，例如 1 或 -1。")
    parser.add_argument(
        "-i",
        "--input",
        nargs="+",
        default=None,
        help=(
            "显式输入文件。structural_relaxation 需要两个 RELAXSTEPS: S1 S2；"
            "fixed_charge_static 需要四个 REPORT: E0 E0_prime Eq Eq_prime。"
        ),
    )
    parser.add_argument("-o", "--output", default="reorganization_energy.log", help="summary log 输出路径，默认当前目录 reorganization_energy.log。")
    parser.add_argument("--detail-output", default=None, help="详细核对 log 输出路径，默认与 summary 同目录 detail_reorg.log。")
    parser.add_argument("--s1-relaxsteps", default=None, help="structural_relaxation 的 S(0->q) RELAXSTEPS；若同时给 -i，则以 -i 为准。")
    parser.add_argument("--s2-relaxsteps", default=None, help="structural_relaxation 的 S(q->0) RELAXSTEPS；若同时给 -i，则以 -i 为准。")
    parser.add_argument("--allow-unconverged", action="store_true", help="RELAXSTEPS 没有 *END 时使用最后一个能量步。")
    parser.add_argument("--E0", default=None, help="fixed_charge_static: q=0 在 Q0 的 REPORT；若同时给 -i，则以 -i 为准。")
    parser.add_argument("--E0-prime", dest="E0_prime", default=None, help="fixed_charge_static: q=0 电荷在 Qq 构型的 REPORT；若同时给 -i，则以 -i 为准。")
    parser.add_argument("--Eq", default=None, help="fixed_charge_static: q 电荷在 Qq 的 REPORT；若同时给 -i，则以 -i 为准。")
    parser.add_argument("--Eq-prime", dest="Eq_prime", default=None, help="fixed_charge_static: q 电荷在 Q0 构型的 REPORT；若同时给 -i，则以 -i 为准。")
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    if args.charge == 0:
        raise SystemExit("错误: --charge 不能为 0。")

    if args.input:
        if args.method == "structural_relaxation":
            if len(args.input) != 2:
                raise SystemExit("错误: structural_relaxation 的 -i/--input 需要两个 RELAXSTEPS: S1 S2。")
            args.s1_relaxsteps, args.s2_relaxsteps = args.input
        else:
            if len(args.input) != 4:
                raise SystemExit("错误: fixed_charge_static 的 -i/--input 需要四个 REPORT: E0 E0_prime Eq Eq_prime。")
            args.E0, args.E0_prime, args.Eq, args.Eq_prime = args.input

    if args.method == "structural_relaxation":
        missing = [name for name in ["s1_relaxsteps", "s2_relaxsteps"] if getattr(args, name) is None]
    else:
        missing = [name for name in ["E0", "E0_prime", "Eq", "Eq_prime"] if getattr(args, name) is None]
    if missing:
        readable = ", ".join("--" + name.replace("_", "-") for name in missing)
        raise SystemExit(f"错误: {args.method} 缺少必需参数: {readable}")


def main() -> int:
    args = parse_args()
    validate_args(args)

    if args.method == "structural_relaxation":
        row = structural_relaxation_row(args)
    else:
        row = fixed_charge_static_row(args)

    output = Path(args.output).expanduser().resolve()
    detail_output = Path(args.detail_output).expanduser().resolve() if args.detail_output else output.parent / "detail_reorg.log"
    write_summary_log(row, output)
    write_detail_log(row, detail_output)

    ok = row.get("lambda_0_to_q_eV") is not None and row.get("lambda_q_to_0_eV") is not None
    print("=" * 70)
    print("重组能提取")
    print("=" * 70)
    print(f"method:  {args.method}")
    print(f"defect:  {args.defect}")
    print(f"charge:  {charge_label(args.charge)}")
    print(f"S0->q:   {format_float(row.get('lambda_0_to_q_eV'))} eV")
    print(f"Sq->0:   {format_float(row.get('lambda_q_to_0_eV'))} eV")
    print(f"output:  {output}")
    print(f"detail:  {detail_output}")
    if not ok:
        print("警告: 至少一个重组能未成功提取，请检查 detail log。", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
