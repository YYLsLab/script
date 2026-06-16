#!/usr/bin/env python3
"""
带电缺陷转变能级修正脚本 (ICIC + PA)。

功能:
  prepare — 生成 defect.input，输出 HPC 执行指南
  collect — 收集 ICIC/PA 修正结果，输出 result/correction_results.yaml

用法:
  python script/correction.py prepare
  python script/correction.py collect
  python script/correction.py all
"""

import os
import sys
import argparse
import numpy as np
import yaml

from defect_utils import (
    get_project_root, load_config, read_atom_config,
    find_defect_position, get_charge_dirs, format_charge,
    read_etot_from_report, read_etot_from_atom_config,
    read_ecoul_from_report, read_vatom_farthest,
    get_bulk_info, read_formation_out_energies,
    parse_atom_change, CleanDumper, FlowList, write_yaml,
)


# ============================================================
#  Prepare: 生成 defect.input
# ============================================================

def prepare(project_root, cfg, output_root=None):
    calc_dir = os.path.join(project_root, 'calculate')
    bulk_scf = os.path.join(calc_dir, 'bulk', 'scf')
    bulk_config = os.path.join(bulk_scf, 'atom.config')

    if not os.path.exists(bulk_config):
        print(f"错误: 找不到 bulk atom.config: {bulk_config}")
        return

    defect_dirs = []
    for item in sorted(os.listdir(calc_dir)):
        full = os.path.join(calc_dir, item)
        if os.path.isdir(full) and item != 'bulk':
            defect_dirs.append((item, full))

    if not defect_dirs:
        print("错误: calculate/ 下未找到缺陷目录")
        return

    print(f"{'='*60}")
    print(f"  [correction.py] 批量生成 defect.input")
    print(f"{'='*60}")
    print(f"检测到 {len(defect_dirs)} 个缺陷目录\n")

    total_generated = 0
    run_commands = []

    for defect_name, defect_dir in defect_dirs:
        print(f"--- {defect_name} ---")
        charge_dirs = get_charge_dirs(defect_dir)
        if not charge_dirs:
            print(f"  未找到 q_* 目录，跳过")
            continue

        charged_states = [q for q in charge_dirs if q != 0]
        if not charged_states:
            print(f"  仅有 q_0，无需修正")
            continue

        neutral_scf = os.path.join(defect_dir, 'q_0', 'scf')
        if not os.path.isdir(neutral_scf):
            print(f"  中性态目录不存在，跳过")
            continue

        defect_config = os.path.join(neutral_scf, 'atom.config')
        defect_pos = None
        if os.path.exists(defect_config):
            defect_pos = find_defect_position(bulk_config, defect_config)
        if defect_pos is None:
            centers = cfg.get('structure', {}).get('defect_center', [])
            if centers:
                defect_pos = tuple(centers[0])
        if defect_pos is None:
            print(f"  无法确定缺陷位置，跳过")
            continue

        bound_pos = tuple((c + 0.5) % 1.0 for c in defect_pos)
        print(f"  缺陷位置: {defect_pos[0]:.6f} {defect_pos[1]:.6f} {defect_pos[2]:.6f}")

        for q in sorted(charged_states):
            qdir = charge_dirs[q]
            scf_dir = os.path.join(qdir, 'scf')
            if not os.path.isdir(scf_dir):
                continue
            charge_str = format_charge(q)
            if output_root:
                rel_scf = os.path.relpath(scf_dir, project_root)
                output_path = os.path.join(output_root, rel_scf, 'defect.input')
                os.makedirs(os.path.dirname(output_path), exist_ok=True)
            else:
                output_path = os.path.join(scf_dir, 'defect.input')
            content = (
                f"defect {defect_pos[0]:.6f} {defect_pos[1]:.6f} {defect_pos[2]:.6f}\n"
                f"bound {bound_pos[0]:.6f} {bound_pos[1]:.6f} {bound_pos[2]:.6f}\n"
                f"charged state {charge_str}\n"
                f"bulk {os.path.abspath(bulk_scf)}\n"
                f"neutral {os.path.abspath(neutral_scf)}\n"
            )
            with open(output_path, 'w') as f:
                f.write(content)
            print(f"  生成: {output_path}")
            total_generated += 1
            run_commands.append((defect_name, charge_str, scf_dir))

    print(f"\n共生成 {total_generated} 个 defect.input")

    if run_commands:
        print(f"\n在 HPC 上执行:")
        script_root = os.path.abspath(os.path.dirname(__file__))
        for defect_name, charge_str, scf_dir in run_commands:
            neutral_dir = os.path.join(os.path.dirname(os.path.dirname(scf_dir)), 'q_0', 'scf')
            print(f"\n  # {defect_name} q={charge_str}")
            print(f"  cd {scf_dir}")
            print(f"  cp {neutral_dir}/OUT.OCC OUT.OCC0 && cp OUT.OCC OUT.OCC1")
            print(f"  bash {script_root}/1_get_rho.sh")
            print(f"  bash {script_root}/2_coulomb_integral.sh")
            print(f"  bash {script_root}/3_get_results.sh")
        print(f"\n修正完成后运行: python {script_root}/correction.py collect")


# ============================================================
#  Collect: 收集 ICIC + PA 修正
# ============================================================

def collect_correction(scf_dir, q, defect_pos, dielectric=None):
    """收集单个电荷态的修正数据"""
    sign = '+' if q > 0 else '-'
    q_abs = abs(q)
    corr_dir = os.path.join(scf_dir, f'image-corr_{sign}{q_abs}')

    E_icic = None
    delta_V = None

    report_inf = os.path.join(corr_dir, 'REPORT')
    report_0 = os.path.join(corr_dir, 'REPORT.0')
    if os.path.exists(report_inf) and os.path.exists(report_0):
        E_inf = read_ecoul_from_report(report_inf)
        E_P = read_ecoul_from_report(report_0)
        if E_inf is not None and E_P is not None and dielectric:
            E_icic = (E_inf - E_P) / dielectric

    defect_input_path = os.path.join(scf_dir, 'defect.input')
    if os.path.exists(defect_input_path):
        with open(defect_input_path, 'r') as f:
            lines = f.readlines()
        neutral_dir = bulk_dir = None
        for line in lines:
            if line.startswith('neutral'):
                neutral_dir = line.split(None, 1)[1].strip()
            if line.startswith('bulk'):
                bulk_dir = line.split(None, 1)[1].strip()
        if neutral_dir and bulk_dir:
            vatom_n = os.path.join(neutral_dir, 'OUT.VATOM')
            vatom_b = os.path.join(bulk_dir, 'OUT.VATOM')
            cfg_n = os.path.join(neutral_dir, 'final.config')
            cfg_b = os.path.join(bulk_dir, 'final.config')
            if not os.path.exists(cfg_n):
                cfg_n = os.path.join(neutral_dir, 'atom.config')
            if not os.path.exists(cfg_b):
                cfg_b = os.path.join(bulk_dir, 'atom.config')
            if all(os.path.exists(p) for p in [vatom_n, vatom_b, cfg_n, cfg_b]):
                V_n = read_vatom_farthest(vatom_n, defect_pos, cfg_n)
                V_b = read_vatom_farthest(vatom_b, defect_pos, cfg_b)
                if V_n is not None and V_b is not None:
                    delta_V = V_n - V_b

    return E_icic, delta_V


def compute_formation_energy(E_corrected, E_bulk, delta_n, q, VBM, chem_pot):
    """
    E_f(q, E_F) = E_corrected - E_bulk + Σ(-Δn_i)·μ_i + q·(VBM + E_F)
    返回 E_f(q, E_F=0)，即费米能级相对 VBM 为 0 时的形成能。
    """
    if E_corrected is None or E_bulk is None or VBM is None:
        return None
    chem_term = 0.0
    for sym, dn in delta_n.items():
        mu = chem_pot.get(sym)
        if mu is not None:
            chem_term += (-dn) * mu
        else:
            print(f"  警告: 元素 {sym} 的化学势未定义")
    return E_corrected - E_bulk + chem_term + q * VBM


def parse_charge_dir_options(values):
    charge_dirs = {}
    for value in values or []:
        if "=" not in value:
            raise ValueError(f"--charge-dir 格式应为 q=DIR，收到: {value}")
        q_text, path = value.split("=", 1)
        q_text = q_text.strip().lower().replace("q", "")
        if q_text in {"p1", "+1", "1"}:
            q = 1
        elif q_text in {"m1", "-1"}:
            q = -1
        else:
            q = int(q_text.replace("+", ""))
        if q != 0:
            charge_dirs[q] = os.path.abspath(path)
    return charge_dirs


def explicit_bulk_info(project_root, cfg, bulk_scf):
    info = get_bulk_info(project_root, cfg)
    system_cfg = cfg.get('system', {})
    if 'chemical_potential' not in info:
        info['chemical_potential'] = system_cfg.get('chemical_potential', {})
    if 'Gap' not in info and system_cfg.get('gap') is not None:
        info['Gap'] = system_cfg.get('gap')
    if 'dielectric' not in info and system_cfg.get('dielectric') is not None:
        info['dielectric'] = system_cfg.get('dielectric')
    bulk_report = os.path.join(bulk_scf, 'REPORT')
    bulk_config = os.path.join(bulk_scf, 'atom.config')
    if os.path.exists(bulk_report):
        e_bulk = read_etot_from_report(bulk_report)
        if e_bulk is not None:
            info['E_bulk'] = e_bulk
    if info.get('E_bulk') is None and os.path.exists(bulk_config):
        info['E_bulk'] = read_etot_from_atom_config(bulk_config)
    return info


def collect_explicit(project_root, cfg, bulk_scf, neutral_scf, charge_dirs, output_path=None, defect_name='defect'):
    """按显式 bulk/q0/qX 静态目录收集单个缺陷的修正结果。"""
    bulk_scf = os.path.abspath(bulk_scf)
    neutral_scf = os.path.abspath(neutral_scf)
    charge_dirs = {q: os.path.abspath(path) for q, path in charge_dirs.items() if q != 0}
    if not os.path.isdir(bulk_scf):
        print(f"错误: bulk 目录不存在: {bulk_scf}")
        return
    if not os.path.isdir(neutral_scf):
        print(f"错误: q=0 目录不存在: {neutral_scf}")
        return
    missing = [path for path in charge_dirs.values() if not os.path.isdir(path)]
    if missing:
        print("错误: 以下带电态目录不存在:")
        for path in missing:
            print(f"  {path}")
        return

    bulk_info = explicit_bulk_info(project_root, cfg, bulk_scf)
    dielectric = bulk_info.get('dielectric')
    bulk_config_path = os.path.join(bulk_scf, 'atom.config')
    defect_config = os.path.join(neutral_scf, 'atom.config')
    defect_pos = None
    if os.path.exists(defect_config) and os.path.exists(bulk_config_path):
        defect_pos = find_defect_position(bulk_config_path, defect_config)
    if defect_pos is None:
        centers = cfg.get('structure', {}).get('defect_center', [])
        if centers:
            defect_pos = tuple(centers[0])

    print(f"{'='*70}")
    print("  [correction.py] 显式目录收集 ICIC + PA 修正结果")
    print(f"{'='*70}\n")
    print(f"--- {defect_name} ---")

    charge_data = {}
    for q, scf_dir in sorted(charge_dirs.items()):
        E_raw = None
        report = os.path.join(scf_dir, 'REPORT')
        config = os.path.join(scf_dir, 'atom.config')
        if os.path.exists(report):
            E_raw = read_etot_from_report(report)
        if E_raw is None and os.path.exists(config):
            E_raw = read_etot_from_atom_config(config)

        E_icic = delta_V = None
        E_corr = 0.0
        if defect_pos is not None:
            E_icic, delta_V = collect_correction(scf_dir, q, defect_pos, dielectric=dielectric)
            if E_icic is not None:
                E_corr += E_icic
            if delta_V is not None:
                E_corr += q * delta_V

        charge_data[q] = {
            'E_raw': E_raw,
            'E_icic': E_icic,
            'delta_V': delta_V,
            'E_corr': E_corr,
            'E_corrected': (E_raw + E_corr) if E_raw is not None else None,
        }
        status = ""
        if E_icic is not None and delta_V is not None:
            status = "✓"
        elif E_icic is not None or delta_V is not None:
            status = "△ 部分"
        else:
            status = "✗ 无修正"
        e_s = f"{E_raw:14.5f}" if E_raw is not None else "N/A"
        print(f"  q={format_charge(q):>3s}  E_raw={e_s} eV  E_corr={E_corr:.5f} eV  {status}")

    # 加入中性态能量。
    E0_raw = None
    report0 = os.path.join(neutral_scf, 'REPORT')
    config0 = os.path.join(neutral_scf, 'atom.config')
    if os.path.exists(report0):
        E0_raw = read_etot_from_report(report0)
    if E0_raw is None and os.path.exists(config0):
        E0_raw = read_etot_from_atom_config(config0)
    charge_data[0] = {
        'E_raw': E0_raw,
        'E_icic': None,
        'delta_V': None,
        'E_corr': 0.0,
        'E_corrected': E0_raw,
    }

    delta_n = parse_atom_change(defect_name, bulk_config_path, defect_config)
    out = _build_correction_yaml(
        cfg,
        bulk_info,
        {defect_name: charge_data},
        calc_dir=None,
        defect_delta_n={defect_name: delta_n},
    )
    output_path = output_path or os.path.join(os.getcwd(), 'correction_results.yaml')
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    write_yaml(out, output_path, header_lines=[
        "============================================================",
        " ICIC + PA 修正、形成能与转变能级",
        " 由 correction.py collect 显式目录模式生成",
        "============================================================",
    ])
    print(f"\n结果已写入: {output_path}")
    for defect in out.get('defects', []):
        levels = defect.get('transition_levels') or {}
        if levels:
            text = ", ".join(f"ε{k}={v:.4f} eV" for k, v in levels.items())
            print(f"  {defect.get('name')}: {text}")


def collect(project_root, cfg, output_path=None):
    """收集修正结果，输出 correction_results.yaml"""
    calc_dir = os.path.join(project_root, 'calculate')
    bulk_info = get_bulk_info(project_root, cfg)
    formation_data = read_formation_out_energies(project_root)
    dielectric = bulk_info.get('dielectric')

    print(f"{'='*70}")
    print(f"  [correction.py] 收集 ICIC + PA 修正结果")
    print(f"{'='*70}\n")

    defect_dirs = []
    for item in sorted(os.listdir(calc_dir)):
        full = os.path.join(calc_dir, item)
        if os.path.isdir(full) and item != 'bulk':
            defect_dirs.append((item, full))

    bulk_config_path = os.path.join(calc_dir, 'bulk', 'scf', 'atom.config')
    all_results = {}

    for defect_name, defect_dir in defect_dirs:
        print(f"--- {defect_name} ---")
        charge_dirs = get_charge_dirs(defect_dir)
        if not charge_dirs:
            continue

        neutral_scf = os.path.join(defect_dir, 'q_0', 'scf')
        defect_config = os.path.join(neutral_scf, 'atom.config')
        defect_pos = None
        if os.path.exists(defect_config) and os.path.exists(bulk_config_path):
            defect_pos = find_defect_position(bulk_config_path, defect_config)
        if defect_pos is None:
            centers = cfg.get('structure', {}).get('defect_center', [])
            if centers:
                defect_pos = tuple(centers[0])

        charge_data = {}
        for q in sorted(charge_dirs.keys()):
            scf_dir = os.path.join(charge_dirs[q], 'scf')
            if not os.path.isdir(scf_dir):
                continue

            # 读取 E_tot
            E_raw = None
            report = os.path.join(scf_dir, 'REPORT')
            config = os.path.join(scf_dir, 'atom.config')
            if os.path.exists(report):
                E_raw = read_etot_from_report(report)
            if E_raw is None and os.path.exists(config):
                E_raw = read_etot_from_atom_config(config)
            if E_raw is None and defect_name in formation_data:
                fd = formation_data[defect_name]['charges'].get(q)
                if fd:
                    E_raw = fd['E_raw']

            E_icic = delta_V = None
            E_corr = 0.0
            if q != 0 and defect_pos is not None:
                E_icic, delta_V = collect_correction(scf_dir, q, defect_pos, dielectric=dielectric)
                if E_icic is not None:
                    E_corr += E_icic
                if delta_V is not None:
                    E_corr += q * delta_V

            charge_data[q] = {
                'E_raw': E_raw,
                'E_icic': E_icic,
                'delta_V': delta_V,
                'E_corr': E_corr,
                'E_corrected': (E_raw + E_corr) if E_raw is not None else None,
            }

            status = ""
            if q != 0:
                if E_icic is not None and delta_V is not None:
                    status = "✓"
                elif E_icic is not None or delta_V is not None:
                    status = "△ 部分"
                else:
                    status = "✗ 无修正"
            e_s = f"{E_raw:14.5f}" if E_raw is not None else "N/A"
            print(f"  q={format_charge(q):>3s}  E_raw={e_s} eV  E_corr={E_corr:.5f} eV  {status}")

        all_results[defect_name] = charge_data

    # ── 输出 YAML ──
    out = _build_correction_yaml(cfg, bulk_info, all_results, calc_dir)
    output_path = output_path or os.path.join(project_root, 'result', 'correction_results.yaml')
    result_dir = os.path.dirname(os.path.abspath(output_path))
    os.makedirs(result_dir, exist_ok=True)
    write_yaml(out, output_path, header_lines=[
        "============================================================",
        " ICIC + PA 修正结果",
        " 由 correction.py collect 生成",
        " 作为 formation_energy.py 的输入",
        "============================================================",
    ])
    print(f"\n结果已写入: {output_path}")
    for defect in out.get('defects', []):
        levels = defect.get('transition_levels') or {}
        if levels:
            text = ", ".join(f"ε{k}={v:.4f} eV" for k, v in levels.items())
            print(f"  {defect.get('name')}: {text}")


def _build_correction_yaml(cfg, bulk_info, all_results, calc_dir=None, defect_delta_n=None):
    """构建修正结果 YAML"""
    defect_delta_n = defect_delta_n or {}
    bulk_config_path = os.path.join(calc_dir, 'bulk', 'scf', 'atom.config') if calc_dir else None
    E_bulk = bulk_info.get('E_bulk')
    VBM = bulk_info.get('VBM')
    Gap = bulk_info.get('Gap')
    chem_pot = bulk_info.get('chemical_potential', {})
    out = {
        'system': {
            'name': cfg.get('system', {}).get('name', ''),
            'E_bulk': E_bulk,
            'VBM': VBM,
            'E_VBM': VBM,
            'E_CBM': (VBM + Gap) if (VBM is not None and Gap is not None) else None,
            'gap': Gap,
            'dielectric': cfg.get('system', {}).get('dielectric'),
            'chemical_potential': chem_pot,
        },
        'band_alignment': cfg.get('band_alignment', {}),
        'defects': [],
    }

    for defect_name, charge_data in all_results.items():
        # 获取 atom_change
        delta_n = defect_delta_n.get(defect_name)
        if delta_n is None:
            neutral_scf = os.path.join(calc_dir, defect_name, 'q_0', 'scf')
            dc_path = os.path.join(neutral_scf, 'atom.config')
            if not os.path.exists(dc_path):
                dc_path = os.path.join(neutral_scf, 'final.config')
            delta_n = parse_atom_change(defect_name, bulk_config_path, dc_path)

        sorted_qs = sorted(charge_data.keys())
        cs_list = []
        formation_energies = {}
        for q in sorted_qs:
            d = charge_data[q]
            E_f0 = compute_formation_energy(d['E_corrected'], E_bulk, delta_n, q, VBM, chem_pot)
            formation_energies[q] = E_f0
            cs_list.append({
                'q': int(q),
                'label': f"q={format_charge(q)}",
                'E_f0': float(E_f0) if E_f0 is not None else None,
                'E_raw': float(d['E_raw']) if d['E_raw'] is not None else None,
                'E_corrected': float(d['E_corrected']) if d['E_corrected'] is not None else None,
                'corrections': {
                    'E_icic': float(d['E_icic']) if d['E_icic'] is not None else None,
                    'delta_V_PA': float(d['delta_V']) if d['delta_V'] is not None else None,
                    'E_corr_total': float(d['E_corr']),
                },
            })

        tl = {}
        for i in range(len(sorted_qs)):
            for j in range(i + 1, len(sorted_qs)):
                q1, q2 = sorted_qs[i], sorted_qs[j]
                Ef1, Ef2 = formation_energies.get(q1), formation_energies.get(q2)
                if Ef1 is not None and Ef2 is not None and (q2 - q1) != 0:
                    eps = (Ef1 - Ef2) / (q2 - q1)
                    label = f"({format_charge(q1)}/{format_charge(q2)})"
                    tl[label] = float(eps)

        fe_table = None
        if Gap is not None and any(v is not None for v in formation_energies.values()):
            ef_vals = np.linspace(0, Gap, 101)
            ef_min_list = []
            stable_q_list = []
            for ef in ef_vals:
                best_ef, best_q = None, None
                for q in sorted_qs:
                    Ef0 = formation_energies.get(q)
                    if Ef0 is None:
                        continue
                    val = Ef0 + q * ef
                    if best_ef is None or val < best_ef:
                        best_ef = val
                        best_q = q
                ef_min_list.append(float(best_ef) if best_ef is not None else None)
                stable_q_list.append(int(best_q) if best_q is not None else None)
            fe_table = {
                'E_F_values': FlowList([float(round(x, 5)) for x in ef_vals]),
                'E_f_min': FlowList([float(round(x, 5)) if x is not None else None for x in ef_min_list]),
                'stable_q': FlowList(stable_q_list),
            }

        defect_entry = {
            'name': defect_name,
            'atom_change': delta_n if delta_n else None,
            'charge_states': cs_list,
            'transition_levels': tl if tl else None,
        }
        if fe_table:
            defect_entry['formation_energy_table'] = fe_table
        out['defects'].append(defect_entry)

    return out


# ============================================================
#  Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="生成 defect.input 或收集 ICIC/PA 修正结果。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""示例:
  python standalone/correction.py prepare -i /path/to/project -o /tmp/defect_inputs
  python standalone/correction.py collect -i /path/to/project -o /path/to/correction_results.yaml
  python standalone/correction.py collect -ibulk /path/to/bulk/scf -iq0 /path/to/q0/scf -iq1 /path/to/q1/scf -iqm1 /path/to/qm1/scf -o /path/to/correction_results.yaml
  python standalone/correction.py all -i /path/to/project -o /path/to/correction_results.yaml
""",
    )
    parser.add_argument("mode", choices=["prepare", "collect", "all"], help="运行模式。")
    parser.add_argument("-i", "--input", default=None, help="兼容模式输入项目根目录；该目录下应包含 calculate/。")
    parser.add_argument("-ibulk", "--input-bulk", default=None, help="显式模式 bulk 静态计算目录。")
    parser.add_argument("-iq0", "--input-q0", default=None, help="显式模式 q=0 静态计算目录。")
    parser.add_argument("-iq1", "--input-q1", default=None, help="显式模式 q=+1 静态计算目录。")
    parser.add_argument("-iqm1", "--input-qm1", default=None, help="显式模式 q=-1 静态计算目录。")
    parser.add_argument(
        "--charge-dir",
        action="append",
        default=[],
        help="显式模式额外带电态目录，格式 q=DIR，例如 --charge-dir 2=/path/to/q2/scf。",
    )
    parser.add_argument("--defect-name", default="defect", help="显式模式输出 YAML 中的缺陷名称。")
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="prepare 模式下为 defect.input 镜像输出目录；collect/all 模式下为 correction_results.yaml。",
    )
    parser.add_argument("--config", default=None, help="显式指定 config.yaml。")
    args = parser.parse_args()

    explicit_mode = bool(args.input_bulk or args.input_q0 or args.input_q1 or args.input_qm1 or args.charge_dir)

    if args.input:
        input_path = os.path.abspath(args.input)
        project_root = os.path.dirname(input_path) if os.path.basename(input_path) == "calculate" else input_path
    elif explicit_mode:
        project_root = os.getcwd()
    else:
        project_root = get_project_root(__file__)
    cfg = load_config(project_root)
    if args.config:
        with open(args.config, 'r', encoding='utf-8') as f:
            cfg = yaml.safe_load(f) or {}
    mode = args.mode.lower()

    if mode == 'prepare':
        if explicit_mode:
            print("错误: prepare 显式目录模式请使用 generate_defect_input.py。")
            sys.exit(1)
        prepare(project_root, cfg, output_root=os.path.abspath(args.output) if args.output else None)
    elif mode == 'collect':
        if explicit_mode:
            if not args.input_bulk or not args.input_q0:
                print("错误: collect 显式目录模式需要 -ibulk 和 -iq0。")
                sys.exit(1)
            charge_dirs = parse_charge_dir_options(args.charge_dir)
            if args.input_q1:
                charge_dirs[1] = os.path.abspath(args.input_q1)
            if args.input_qm1:
                charge_dirs[-1] = os.path.abspath(args.input_qm1)
            if not charge_dirs:
                print("错误: collect 显式目录模式至少需要一个非零电荷态目录。")
                sys.exit(1)
            collect_explicit(
                project_root,
                cfg,
                args.input_bulk,
                args.input_q0,
                charge_dirs,
                output_path=os.path.abspath(args.output) if args.output else None,
                defect_name=args.defect_name,
            )
        else:
            collect(project_root, cfg, output_path=os.path.abspath(args.output) if args.output else None)
    elif mode == 'all':
        if explicit_mode:
            print("错误: all 模式不支持显式目录输入；请分别运行 generate_defect_input.py 和 correction.py collect。")
            sys.exit(1)
        prepare(project_root, cfg)
        print("\n")
        collect(project_root, cfg, output_path=os.path.abspath(args.output) if args.output else None)
    else:
        print(f"未知模式: {mode}")
        sys.exit(1)


if __name__ == '__main__':
    main()
