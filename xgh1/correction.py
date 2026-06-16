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

def prepare(project_root, cfg):
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


def collect(project_root, cfg):
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
    result_dir = os.path.join(project_root, 'result')
    os.makedirs(result_dir, exist_ok=True)
    output_path = os.path.join(result_dir, 'correction_results.yaml')
    write_yaml(out, output_path, header_lines=[
        "============================================================",
        " ICIC + PA 修正结果",
        " 由 correction.py collect 生成",
        " 作为 formation_energy.py 的输入",
        "============================================================",
    ])
    print(f"\n结果已写入: {output_path}")
    script_root = os.path.abspath(os.path.dirname(__file__))
    print(f"接下来运行: python {script_root}/formation_energy.py")


def _build_correction_yaml(cfg, bulk_info, all_results, calc_dir):
    """构建修正结果 YAML"""
    bulk_config_path = os.path.join(calc_dir, 'bulk', 'scf', 'atom.config')
    out = {
        'system': {
            'name': cfg.get('system', {}).get('name', ''),
            'E_bulk': bulk_info.get('E_bulk'),
            'VBM': bulk_info.get('VBM'),
            'gap': bulk_info.get('Gap'),
            'dielectric': cfg.get('system', {}).get('dielectric'),
            'chemical_potential': bulk_info.get('chemical_potential', {}),
        },
        'band_alignment': cfg.get('band_alignment', {}),
        'defects': [],
    }

    for defect_name, charge_data in all_results.items():
        # 获取 atom_change
        neutral_scf = os.path.join(calc_dir, defect_name, 'q_0', 'scf')
        dc_path = os.path.join(neutral_scf, 'atom.config')
        if not os.path.exists(dc_path):
            dc_path = os.path.join(neutral_scf, 'final.config')
        delta_n = parse_atom_change(defect_name, bulk_config_path, dc_path)

        sorted_qs = sorted(charge_data.keys())
        cs_list = []
        for q in sorted_qs:
            d = charge_data[q]
            cs_list.append({
                'q': int(q),
                'E_raw': float(d['E_raw']) if d['E_raw'] is not None else None,
                'E_corrected': float(d['E_corrected']) if d['E_corrected'] is not None else None,
                'corrections': {
                    'E_icic': float(d['E_icic']) if d['E_icic'] is not None else None,
                    'delta_V_PA': float(d['delta_V']) if d['delta_V'] is not None else None,
                    'E_corr_total': float(d['E_corr']),
                },
            })

        out['defects'].append({
            'name': defect_name,
            'atom_change': delta_n if delta_n else None,
            'charge_states': cs_list,
        })

    return out


# ============================================================
#  Main
# ============================================================

def main():
    project_root = get_project_root(__file__)
    if len(sys.argv) < 2:
        print("用法:")
        print("  python script/correction.py prepare   # 生成 defect.input")
        print("  python script/correction.py collect   # 收集修正 → correction_results.yaml")
        print("  python script/correction.py all")
        sys.exit(1)

    cfg = load_config(project_root)
    mode = sys.argv[1].lower()

    if mode == 'prepare':
        prepare(project_root, cfg)
    elif mode == 'collect':
        collect(project_root, cfg)
    elif mode == 'all':
        prepare(project_root, cfg)
        print("\n")
        collect(project_root, cfg)
    else:
        print(f"未知模式: {mode}")
        sys.exit(1)


if __name__ == '__main__':
    main()
