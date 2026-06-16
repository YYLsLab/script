#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot a two-state Configuration Coordinate Diagram (CCD)."""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import yaml

if sys.platform == "win32":
    import io

    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")


def find_project_root(start: Optional[Path] = None) -> Path:
    current = (start or Path(__file__).resolve().parent).resolve()
    fallback = None
    while True:
        if (current / "result").is_dir() or (current / "calculate").is_dir():
            return current
        if fallback is None and (current / "config.yaml").exists():
            fallback = current
        parent = current.parent
        if parent == current:
            return fallback or (start or Path(__file__).resolve().parent).resolve()
        current = parent


def defect_index(name: str) -> str:
    match = re.match(r"^(\d+)(?:_|$)", str(name))
    return match.group(1) if match else str(name)


def same_defect(name: str, query: str) -> bool:
    return str(name) == str(query) or defect_index(name) == str(query)


def row_matches_defect(row: dict, defect: str) -> bool:
    if "defect" in row and same_defect(row.get("defect", ""), defect):
        return True
    if "idx" in row:
        return str(row.get("idx", "")) in {str(defect), defect_index(defect)}
    return False


def parse_table(path: Path) -> List[dict]:
    if not path.exists():
        raise FileNotFoundError(f"输入文件不存在: {path}")
    lines = [
        line.strip()
        for line in path.read_text(encoding="utf-8", errors="ignore").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if not lines:
        return []
    headers = lines[0].split()
    rows = []
    for line in lines[1:]:
        parts = line.split()
        if len(parts) < len(headers):
            continue
        rows.append({key: value for key, value in zip(headers, parts)})
    return rows


def parse_charge(value) -> int:
    text = str(value).strip()
    if text.startswith("q="):
        text = text[2:]
    return int(text.replace("+", ""))


def charge_label(charge: int) -> str:
    return f"+{charge}" if charge > 0 else str(charge)


def read_deltaq(deltaq_log: Path, defect: str, charge: int) -> float:
    rows = parse_table(deltaq_log)
    if not rows:
        raise ValueError(f"ΔQ log 为空: {deltaq_log}")

    # New wide format:
    # defect deltaQ_sqrt_amu_A[q+1] deltaQ_sqrt_amu_A[q-1]
    wide_key = f"deltaQ_sqrt_amu_A[q{charge_label(charge)}]"
    for row in rows:
        if same_defect(row.get("defect", ""), defect) and wide_key in row:
            return float(row[wide_key])

    # Backward-compatible long format:
    # defect charge deltaQ_sqrt_amu_A ...
    for row in rows:
        if same_defect(row.get("defect", ""), defect) and "charge" in row:
            if parse_charge(row["charge"]) == charge:
                return float(row["deltaQ_sqrt_amu_A"])

    raise ValueError(f"无法在 {deltaq_log} 中找到 defect={defect}, charge={charge} 的 ΔQ")


def read_reorganization(reorg_log: Path, defect: str, charge: int) -> Tuple[float, float, str]:
    rows = parse_table(reorg_log)

    # New wide format:
    # idx S[0->+1] S[+1->0] S[0->-1] S[-1->0]
    if charge == 1:
        forward_key, backward_key = "S[0->+1]", "S[+1->0]"
    elif charge == -1:
        forward_key, backward_key = "S[0->-1]", "S[-1->0]"
    else:
        forward_key = backward_key = ""
    for row in rows:
        if not row_matches_defect(row, defect):
            continue
        if forward_key in row and backward_key in row:
            return float(row[forward_key]), float(row[backward_key]), "wide:S"

    # Backward-compatible long format:
    # defect charge lambda_0_to_q_eV lambda_q_to_0_eV ...
    for row in rows:
        if not row_matches_defect(row, defect):
            continue
        if "charge" not in row or parse_charge(row["charge"]) != charge:
            continue

        if "lambda_0_to_q_eV" in row and "lambda_q_to_0_eV" in row:
            return float(row["lambda_0_to_q_eV"]), float(row["lambda_q_to_0_eV"]), "asymmetric"
        if "S_0_to_q_eV" in row and "S_q_to_0_eV" in row:
            return float(row["S_0_to_q_eV"]), float(row["S_q_to_0_eV"]), "asymmetric"
        if "lambda_eV" in row:
            val = float(row["lambda_eV"])
            return val, val, "symmetric"
        if "S_eV" in row:
            val = float(row["S_eV"])
            return val, val, "symmetric"

    raise ValueError(f"无法在 {reorg_log} 中找到 defect={defect}, charge={charge} 的重组能")


def load_defect_results(path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(f"defect_results.yaml 不存在: {path}")
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def get_defect_entry(data: dict, defect: str) -> dict:
    for entry in data.get("defects", []) or []:
        if same_defect(entry.get("name", ""), defect):
            return entry
    raise ValueError(f"defect_results.yaml 中找不到缺陷: {defect}")


def charge_state_map(defect_entry: dict) -> Dict[int, dict]:
    states: Dict[int, dict] = {}
    for state in defect_entry.get("charge_states", []) or []:
        try:
            states[int(state.get("q"))] = state
        except (TypeError, ValueError):
            continue
    return states


def transition_keys(charge: int) -> List[str]:
    if charge == 1:
        return ["(0/+1)", "(0/+)", "(+1/0)"]
    if charge == -1:
        return ["(-1/0)", "(-/0)", "(0/-1)", "(0/-)"]
    return [f"(0/{charge_label(charge)})", f"({charge_label(charge)}/0)"]


def get_transition_level(defect_entry: dict, charge: int) -> Optional[float]:
    levels = defect_entry.get("transition_levels") or {}
    for key in transition_keys(charge):
        if key in levels:
            return float(levels[key])
    return None


def delta_e_from_transition(defect_entry: dict, charge: int, fermi_level: str) -> Tuple[float, Optional[float], str]:
    if fermi_level == "transition":
        eps = get_transition_level(defect_entry, charge)
        if eps is None:
            raise ValueError(f"缺少 q=0 与 q={charge_label(charge)} 的转变能级")
        return 0.0, eps, "transition_degenerate"

    ef = float(fermi_level)
    eps = get_transition_level(defect_entry, charge)
    if eps is not None:
        return charge * (ef - eps), eps, "transition_level"

    states = charge_state_map(defect_entry)
    if 0 not in states or charge not in states:
        raise ValueError(f"缺少 q=0 或 q={charge_label(charge)} 的 E_f0")
    delta_e = float(states[charge]["E_f0"]) + charge * ef - float(states[0]["E_f0"])
    return delta_e, None, "formation_energy"


def solve_crossing(
    delta_q: float,
    lambda_0_to_q: float,
    lambda_q_to_0: float,
    delta_e: float,
) -> dict:
    k0 = lambda_0_to_q / (delta_q ** 2)
    kq = lambda_q_to_0 / (delta_q ** 2)
    e0 = 0.0
    eq = delta_e

    # E0 + k0 Q^2 = Eq + kq (Q - D)^2
    a = k0 - kq
    b = 2.0 * kq * delta_q
    c = e0 - eq - kq * delta_q ** 2

    roots: List[float] = []
    if abs(a) < 1e-14:
        if abs(b) > 1e-14:
            roots = [-c / b]
    else:
        disc = b * b - 4.0 * a * c
        if disc >= -1e-12:
            disc = max(0.0, disc)
            sqrt_disc = math.sqrt(disc)
            roots = [(-b + sqrt_disc) / (2.0 * a), (-b - sqrt_disc) / (2.0 * a)]

    q_cross = None
    if roots:
        inside = [root for root in roots if -1e-10 <= root <= delta_q + 1e-10]
        if inside:
            q_cross = min(inside, key=lambda x: abs(x - 0.5 * delta_q))
        else:
            q_cross = min(roots, key=lambda x: abs(x - 0.5 * delta_q))

    if q_cross is None:
        return {
            "k0": k0,
            "kq": kq,
            "Q_cross": None,
            "E_cross": None,
            "barrier_0_to_q": None,
            "barrier_q_to_0": None,
            "crossing_inside": False,
        }

    e_cross = e0 + k0 * q_cross ** 2
    barrier_0_to_q = e_cross - e0
    barrier_q_to_0 = e_cross - eq
    return {
        "k0": k0,
        "kq": kq,
        "Q_cross": q_cross,
        "E_cross": e_cross,
        "barrier_0_to_q": barrier_0_to_q,
        "barrier_q_to_0": barrier_q_to_0,
        "crossing_inside": 0.0 <= q_cross <= delta_q,
    }


def format_float(value) -> str:
    if value is None:
        return "nan"
    return f"{float(value):.8f}"


def write_ccd_log(path: Path, row: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "defect",
        "charge",
        "mode",
        "deltaQ_sqrt_amu_A",
        "lambda_0_to_q_eV",
        "lambda_q_to_0_eV",
        "deltaE_eV",
        "transition_level_eV",
        "fermi_level",
        "k0_eV_per_Q2",
        "kq_eV_per_Q2",
        "Q_cross_sqrt_amu_A",
        "E_cross_eV",
        "barrier_0_to_q_eV",
        "barrier_q_to_0_eV",
        "crossing_inside",
        "source",
    ]
    with path.open("w", encoding="utf-8") as handle:
        handle.write(" ".join(fields) + "\n")
        values = []
        for field in fields:
            value = row.get(field)
            if isinstance(value, float) or value is None:
                values.append(format_float(value))
            else:
                values.append(str(value))
        handle.write(" ".join(values) + "\n")


def plot_ccd(
    output_base: Path,
    defect: str,
    charge: int,
    delta_q: float,
    lambda_0_to_q: float,
    lambda_q_to_0: float,
    delta_e: float,
    crossing: dict,
    formats: Iterable[str],
    dpi: int,
    show: bool,
) -> None:
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "错误: 缺少 matplotlib 或 numpy，无法绘图。\n"
            "请先安装: pip install matplotlib numpy\n"
            "如果在服务器 conda 环境中运行，也可用: conda install matplotlib numpy"
        ) from exc

    k0 = crossing["k0"]
    kq = crossing["kq"]
    e0 = 0.0
    eq = delta_e

    q_min = -0.25 * delta_q
    q_max = 1.25 * delta_q
    q_cross = crossing.get("Q_cross")
    if q_cross is not None:
        q_min = min(q_min, q_cross - 0.15 * delta_q)
        q_max = max(q_max, q_cross + 0.15 * delta_q)

    q_values = np.linspace(q_min, q_max, 500)
    v0 = e0 + k0 * q_values ** 2
    vq = eq + kq * (q_values - delta_q) ** 2

    fig, ax = plt.subplots(figsize=(6.2, 4.5))
    ax.plot(q_values, v0, linewidth=2.0, color="#1f77b4", label="D^0")
    ax.plot(q_values, vq, linewidth=2.0, color="#d62728", label=f"D^{charge_label(charge)}")
    ax.scatter([0.0, delta_q], [e0, eq], color=["#1f77b4", "#d62728"], zorder=4)
    ax.text(0.0, e0, "  Q0", va="bottom", fontsize=9)
    ax.text(delta_q, eq, f"  Q{charge_label(charge)}", va="bottom", fontsize=9)

    ax.annotate(
        "",
        xy=(delta_q, min(e0, eq) - 0.08 * max(lambda_0_to_q, lambda_q_to_0, 1e-6)),
        xytext=(0.0, min(e0, eq) - 0.08 * max(lambda_0_to_q, lambda_q_to_0, 1e-6)),
        arrowprops={"arrowstyle": "<->", "linewidth": 1.0, "color": "black"},
    )
    ax.text(0.5 * delta_q, min(e0, eq) - 0.12 * max(lambda_0_to_q, lambda_q_to_0, 1e-6), "ΔQ", ha="center", fontsize=9)

    if q_cross is not None:
        e_cross = crossing["E_cross"]
        ax.scatter([q_cross], [e_cross], color="black", marker="x", s=55, zorder=5)
        ax.text(q_cross, e_cross, "  crossing", fontsize=9, va="bottom")
        ax.axvline(q_cross, color="black", alpha=0.22, linewidth=0.9)

    ax.set_xlabel(r"Configuration coordinate Q ($\sqrt{\mathrm{amu}}\,\AA$)")
    ax.set_ylabel("Energy (eV)")
    ax.set_title(f"Configuration Coordinate Diagram: {defect}, q={charge_label(charge)}")
    ax.grid(True, alpha=0.25, linewidth=0.7)
    ax.legend(frameon=False)
    fig.tight_layout()

    output_base.parent.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        fmt = fmt.strip().lower().lstrip(".")
        if not fmt:
            continue
        out = output_base.with_suffix(f".{fmt}")
        fig.savefig(out, dpi=dpi if fmt in {"png", "jpg", "jpeg"} else None)
        print(f"图像已写入: {out}")

    if show:
        plt.show()
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    root = find_project_root()
    parser = argparse.ArgumentParser(
        description="绘制指定缺陷位点的二态 Configuration Coordinate Diagram (CCD)。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 只需要重组能log + 手动输入ΔQ；默认 ΔG=0，表示两态简并
  python standalone/plot_ccd.py -i result/reorganization_energy.log --deltaQ 2.4 --defect 0_v_Si --charge 1 -o result/ccd/0_v_Si.png

  # 非简并情况：额外手动指定 ΔG = E_q - E_0
  python standalone/plot_ccd.py -i result/reorganization_energy.log --deltaQ 2.4 --defect 0_v_Si --charge 1 --deltaG 0.25

  # 兼容旧方式：从deltaQ log读取ΔQ
  python standalone/plot_ccd.py -i result --deltaq-log result/deltaQ_scf_atom.log --defect 0_v_Si --charge 1
        """,
    )
    parser.add_argument(
        "-i",
        "--input",
        default=None,
        help="输入重组能log，或result目录/项目根目录。文件路径时按reorganization_energy.log解析。",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="输出图像文件或输出基名。带 .png/.svg/.pdf 后缀时只输出该格式；无后缀时使用 --formats。",
    )
    parser.add_argument("--defect", required=True, help="缺陷名或编号，例如 10_H_i 或 10。")
    parser.add_argument("--charge", required=True, type=int, help="目标电荷态，例如 1 或 -1。")
    parser.add_argument(
        "--mode",
        choices=["intrinsic", "transition"],
        default="intrinsic",
        help="能量差来源。默认 intrinsic，即直接使用手动输入的 --deltaG/--delta-e；transition 才读取 defect_results.yaml。",
    )
    parser.add_argument(
        "--deltaG",
        "--delta-g",
        "--delta-e",
        dest="delta_e",
        type=float,
        default=0.0,
        help="手动指定 CCD 中 E_q - E_0，即 ΔG，单位 eV；默认0，表示两态简并。",
    )
    parser.add_argument("--fermi-level", default="transition", help="transition 模式下的费米能级；可用 transition 或数值。")
    parser.add_argument("--deltaQ", "--delta-q", dest="delta_q", type=float, default=None, help="直接指定ΔQ，单位 sqrt(amu)*Å。")
    parser.add_argument("--S1", type=float, default=None, help="直接指定 λ_0→q，即 V0 在 Qq 处的重组能，单位 eV。")
    parser.add_argument("--S2", type=float, default=None, help="直接指定 λ_q→0，即 Vq 在 Q0 处的重组能，单位 eV。")
    parser.add_argument("--deltaq-log", default=None, help="ΔQ log 路径。只有未指定 --deltaQ 时才读取。")
    parser.add_argument("--reorg-log", default=None, help="重组能 log 路径。默认由 -i 推断。")
    parser.add_argument("--defect-results", default=None, help="defect_results.yaml 路径。仅 --mode transition 需要。")
    parser.add_argument("--output-dir", default=None, help="输出目录。默认 <input>/ccd。")
    parser.add_argument("--formats", default="png,svg", help="输出格式，默认 png,svg。")
    parser.add_argument("--dpi", type=int, default=300, help="PNG 分辨率，默认 300。")
    parser.add_argument("--show", action="store_true", help="绘图后显示窗口。")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    defect = args.defect
    charge = args.charge
    if charge == 0:
        print("错误: --charge 不能为 0。", file=sys.stderr)
        return 1

    root = find_project_root()
    input_path = Path(args.input).expanduser().resolve() if args.input else root / "result" / "reorganization_energy.log"
    if input_path.is_file():
        input_dir = input_path.parent
        inferred_reorg_log = input_path
    else:
        input_dir = input_path
        if (input_dir / "result").is_dir() and not (input_dir / "reorganization_energy.log").exists():
            input_dir = input_dir / "result"
        inferred_reorg_log = input_dir / "reorganization_energy.log"

    deltaq_log = Path(args.deltaq_log).expanduser().resolve() if args.deltaq_log else None
    reorg_log = Path(args.reorg_log).expanduser().resolve() if args.reorg_log else inferred_reorg_log
    defect_results = (
        Path(args.defect_results).expanduser().resolve()
        if args.defect_results
        else input_dir / "defect_results.yaml"
    )

    if args.delta_q is not None:
        delta_q = float(args.delta_q)
        deltaq_source = "cli:deltaQ"
    elif deltaq_log is not None:
        delta_q = read_deltaq(deltaq_log, defect, charge)
        deltaq_source = str(deltaq_log)
    else:
        print("错误: plot_ccd.py 需要 --deltaQ；或通过 --deltaq-log 指定可读取的ΔQ log。", file=sys.stderr)
        return 1
    if delta_q <= 0:
        print(f"错误: ΔQ 必须大于 0，当前为 {delta_q}", file=sys.stderr)
        return 1

    if args.S1 is not None or args.S2 is not None:
        if args.S1 is None or args.S2 is None:
            print("错误: --S1 和 --S2 必须同时指定。", file=sys.stderr)
            return 1
        lambda_0_to_q = float(args.S1)
        lambda_q_to_0 = float(args.S2)
        lambda_source = "cli:S1,S2"
    else:
        lambda_0_to_q, lambda_q_to_0, lambda_source = read_reorganization(reorg_log, defect, charge)

    transition_level = None
    energy_source = args.mode
    if args.mode == "intrinsic":
        delta_e = float(args.delta_e)
        energy_source = "deltaG"
    else:
        data = load_defect_results(defect_results)
        defect_entry = get_defect_entry(data, defect)
        delta_e, transition_level, energy_source = delta_e_from_transition(
            defect_entry, charge, str(args.fermi_level)
        )

    crossing = solve_crossing(delta_q, lambda_0_to_q, lambda_q_to_0, delta_e)
    label = f"{defect}_q0_q{charge_label(charge)}".replace("+", "p").replace("-", "m")
    formats = args.formats.split(",")
    if args.output:
        output_path = Path(args.output).expanduser().resolve()
        suffix = output_path.suffix.lower().lstrip(".")
        if suffix:
            output_base = output_path.with_suffix("")
            formats = [suffix]
        else:
            output_base = output_path
    else:
        output_dir = Path(args.output_dir).expanduser().resolve() if args.output_dir else input_dir / "ccd"
        output_base = output_dir / label
    plot_ccd(
        output_base,
        defect,
        charge,
        delta_q,
        lambda_0_to_q,
        lambda_q_to_0,
        delta_e,
        crossing,
        formats,
        args.dpi,
        args.show,
    )

    log_row = {
        "defect": defect,
        "charge": charge_label(charge),
        "mode": args.mode,
        "deltaQ_sqrt_amu_A": delta_q,
        "lambda_0_to_q_eV": lambda_0_to_q,
        "lambda_q_to_0_eV": lambda_q_to_0,
        "deltaE_eV": delta_e,
        "transition_level_eV": transition_level,
        "fermi_level": args.fermi_level,
        "k0_eV_per_Q2": crossing["k0"],
        "kq_eV_per_Q2": crossing["kq"],
        "Q_cross_sqrt_amu_A": crossing["Q_cross"],
        "E_cross_eV": crossing["E_cross"],
        "barrier_0_to_q_eV": crossing["barrier_0_to_q"],
        "barrier_q_to_0_eV": crossing["barrier_q_to_0"],
        "crossing_inside": crossing["crossing_inside"],
        "source": f"{energy_source};lambda={lambda_source};deltaQ={deltaq_source}",
    }
    log_path = output_base.with_suffix(".log")
    write_ccd_log(log_path, log_row)
    print(f"log 已写入: {log_path}")

    if crossing["Q_cross"] is None:
        print("警告: 两条抛物线没有实数交点。")
    elif not crossing["crossing_inside"]:
        print("警告: 交点位于 [Q0, Qq] 区间之外。")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
