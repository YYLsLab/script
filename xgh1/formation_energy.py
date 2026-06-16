#!/usr/bin/env python3
"""
缺陷形成能与转变能级计算脚本。

输入:
  result/correction_results.yaml — 由 script/correction.py collect 生成
  （或直接从 calculate/ 目录读取，如果修正结果不存在）

输出:
  result/defect_results.yaml — 形成能、转变能级汇总（兼容 CCD 脚本）
  result/E_forms/E_formation_corrected_<defect>.txt — E_f vs E_F 表格

用法:
  python script/formation_energy.py
  python script/formation_energy.py --input result/correction_results.yaml
"""

import os
import sys
import numpy as np
import yaml

from defect_utils import (
    get_project_root, load_config, format_charge,
    get_bulk_info, parse_atom_change, read_atom_config,
    CleanDumper, FlowList, write_yaml,
)


# ============================================================
#  形成能计算
# ============================================================

def compute_formation_energy(E_corrected, E_bulk, delta_n, q, VBM, chem_pot):
    """
    E_f(q, E_F) = E_corrected - E_bulk + Σ(-Δn_i)·μ_i + q·(VBM + E_F)
    返回 E_f(q, E_F=0)
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


def load_correction_results(project_root):
    """加载 correction_results.yaml"""
    path = os.path.join(project_root, 'result', 'correction_results.yaml')
    if not os.path.exists(path):
        return None
    with open(path, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)


def run(project_root, cfg):
    """主计算流程"""
    # 尝试加载修正结果
    corr_data = load_correction_results(project_root)

    if corr_data is None:
        print("未找到 result/correction_results.yaml")
        print("请先运行: python script/correction.py collect")
        print("或确保 calculate/ 目录下有完整的计算结果。")
        sys.exit(1)

    # 从修正结果中提取 bulk 信息
    sys_info = corr_data.get('system', {})
    E_bulk = sys_info.get('E_bulk')
    VBM = sys_info.get('VBM')
    Gap = sys_info.get('gap')
    chem_pot = sys_info.get('chemical_potential', {})

    # 如果修正结果中缺少信息，从 config / formation.out 补充
    bulk_info = get_bulk_info(project_root, cfg)
    if E_bulk is None:
        E_bulk = bulk_info.get('E_bulk')
    if VBM is None:
        VBM = bulk_info.get('VBM')
    if Gap is None:
        Gap = bulk_info.get('Gap')
    if not chem_pot:
        chem_pot = bulk_info.get('chemical_potential', {})

    print(f"{'='*70}")
    print(f"  [formation_energy.py] 计算形成能与转变能级")
    print(f"{'='*70}")
    print(f"E_bulk = {E_bulk} eV")
    print(f"VBM = {VBM} eV")
    print(f"Gap = {Gap} eV")
    print(f"化学势: {chem_pot}\n")

    # ── 逐缺陷计算 ──
    defects_in = corr_data.get('defects', [])
    defects_out = []

    for dentry in defects_in:
        defect_name = dentry['name']
        delta_n = dentry.get('atom_change') or {}
        charge_states = dentry.get('charge_states', [])

        print(f"--- {defect_name} ---")
        if delta_n:
            print(f"  原子数变化: {delta_n}")

        sorted_cs = sorted(charge_states, key=lambda x: x['q'])
        formation_energies = {}
        cs_out = []

        for cs in sorted_cs:
            q = cs['q']
            E_corrected = cs.get('E_corrected')
            E_raw = cs.get('E_raw')

            E_f0 = compute_formation_energy(E_corrected, E_bulk, delta_n, q, VBM, chem_pot)
            formation_energies[q] = E_f0

            cs_out_entry = {
                'q': int(q),
                'label': f"q={format_charge(q)}",
                'E_f0': float(E_f0) if E_f0 is not None else None,
                'E_raw': float(E_raw) if E_raw is not None else None,
                'E_corrected': float(E_corrected) if E_corrected is not None else None,
                'corrections': cs.get('corrections', {}),
            }
            cs_out.append(cs_out_entry)

            # 打印表达式
            if E_f0 is not None:
                q_str = format_charge(q)
                if q == 0:
                    print(f"  E_f(q={q_str}, E_F) = {E_f0:.5f} eV")
                else:
                    sign = "+" if q > 0 else "-"
                    coeff = "" if abs(q) == 1 else f"{abs(q)}·"
                    print(f"  E_f(q={q_str}, E_F) = {E_f0:.5f} {sign} {coeff}E_F  (eV)")

        # ── 转变能级 ──
        sorted_qs = sorted(formation_energies.keys())
        tl = {}
        for i in range(len(sorted_qs)):
            for j in range(i + 1, len(sorted_qs)):
                q1, q2 = sorted_qs[i], sorted_qs[j]
                Ef1, Ef2 = formation_energies.get(q1), formation_energies.get(q2)
                if Ef1 is not None and Ef2 is not None and (q2 - q1) != 0:
                    eps = (Ef1 - Ef2) / (q2 - q1)
                    label = f"({format_charge(q1)}/{format_charge(q2)})"
                    tl[label] = float(eps)
                    in_gap = ""
                    if Gap is not None:
                        in_gap = " [带隙内]" if 0 <= eps <= Gap else " [带隙外]"
                    print(f"  ε{label} = {eps:.4f} eV{in_gap}")

        # ── 包络线表格 ──
        fe_table = None
        if Gap is not None and any(v is not None for v in formation_energies.values()):
            n_pts = 101
            ef_vals = np.linspace(0, Gap, n_pts)
            ef_min_list, stable_q_list = [], []
            for ef in ef_vals:
                best_ef, best_q = None, None
                for q in sorted_qs:
                    Ef0 = formation_energies.get(q)
                    if Ef0 is not None:
                        val = Ef0 + q * ef
                        if best_ef is None or val < best_ef:
                            best_ef = val
                            best_q = q
                ef_min_list.append(float(best_ef))
                stable_q_list.append(int(best_q))
            fe_table = {
                'E_F_values': FlowList([float(round(x, 5)) for x in ef_vals]),
                'E_f_min': FlowList([float(round(x, 5)) for x in ef_min_list]),
                'stable_q': FlowList(stable_q_list),
            }

        defect_out = {
            'name': defect_name,
            'atom_change': delta_n if delta_n else None,
            'charge_states': cs_out,
            'transition_levels': tl if tl else None,
        }
        if fe_table:
            defect_out['formation_energy_table'] = fe_table
        defects_out.append(defect_out)

    # ── 构建输出 YAML ──
    out = {
        'system': {
            'name': cfg.get('system', {}).get('name', ''),
            'E_VBM': float(VBM) if VBM is not None else None,
            'E_CBM': float(VBM + Gap) if (VBM is not None and Gap is not None) else None,
            'gap': float(Gap) if Gap is not None else None,
            'E_bulk': float(E_bulk) if E_bulk is not None else None,
            'chemical_potential': {k: float(v) for k, v in chem_pot.items()},
            'dielectric': cfg.get('system', {}).get('dielectric'),
        },
        'band_alignment': corr_data.get('band_alignment', cfg.get('band_alignment', {})),
        'defects': defects_out,
    }

    # 写入主 YAML
    result_dir = os.path.join(project_root, 'result')
    os.makedirs(result_dir, exist_ok=True)
    yaml_path = os.path.join(result_dir, 'defect_results.yaml')
    write_yaml(out, yaml_path, header_lines=[
        "============================================================",
        " 缺陷形成能与转变能级计算结果",
        " 由 formation_energy.py 生成",
        " 格式兼容 CCD 脚本输入",
        "============================================================",
    ])
    print(f"\n  YAML 结果: {yaml_path}")

    # 写入 E_f vs E_F 表格
    ef_dir = os.path.join(result_dir, 'E_forms')
    os.makedirs(ef_dir, exist_ok=True)
    for dentry in defects_out:
        ft = dentry.get('formation_energy_table')
        if not ft:
            continue
        dname = dentry['name']
        sorted_qs = sorted([cs['q'] for cs in dentry['charge_states']])
        form_en = {cs['q']: cs['E_f0'] for cs in dentry['charge_states']}
        table_path = os.path.join(ef_dir, f'E_formation_corrected_{dname}.txt')
        with open(table_path, 'w', encoding='utf-8') as f:
            q_headers = [f"q={format_charge(q)}(eV)" for q in sorted_qs]
            f.write(f"{'E_F(eV)':>12s}")
            for qh in q_headers:
                f.write(f"  {qh:>14s}")
            f.write(f"  {'E_f_min(eV)':>14s}  {'stable_q':>10s}\n")
            for idx, ef_val in enumerate(ft['E_F_values']):
                f.write(f"{ef_val:12.5f}")
                for q in sorted_qs:
                    Ef0 = form_en.get(q)
                    if Ef0 is not None:
                        f.write(f"  {Ef0 + q * ef_val:14.5f}")
                    else:
                        f.write(f"  {'N/A':>14s}")
                f.write(f"  {ft['E_f_min'][idx]:14.5f}  {format_charge(ft['stable_q'][idx]):>10s}\n")
        print(f"  E_f 表格: {table_path}")

    print(f"\n{'='*70}")
    print(f"  完成")
    print(f"{'='*70}")


# ============================================================
#  Main
# ============================================================

def main():
    project_root = get_project_root(__file__)
    if len(sys.argv) > 1 and sys.argv[1] in ('-h', '--help'):
        print("用法: python script/formation_energy.py")
        print("  读取 result/correction_results.yaml，计算形成能和转变能级")
        print("  输出 result/defect_results.yaml")
        sys.exit(0)

    cfg = load_config(project_root)
    run(project_root, cfg)


if __name__ == '__main__':
    main()
