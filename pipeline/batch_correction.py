#!/usr/bin/env python3
"""
批量缺陷能级修正与转变能级计算脚本。

功能:
  1. 自动扫描 calculate/ 下所有缺陷目录（排除 bulk）
  2. 对每个缺陷的每个带电态:
     - 自动生成 defect.input
     - 生成 ICIC+PA 修正所需的脚本
  3. 收集已完成的修正结果 (E_∞, E_P, ΔV_PA)
  4. 计算修正后的形成能和转变能级
  5. 输出 YAML 格式结果到 result/defect_results.yaml（兼容 CCD 脚本）

输入:
  config.yaml — 统一配置文件（替代旧的 input / input.json）

用法:
  python script/batch_correction.py prepare   # 生成 defect.input 和修正脚本
  python script/batch_correction.py collect   # 收集结果并计算转变能级
  python script/batch_correction.py all       # 先 prepare 再 collect
"""

import os
import sys
import re
import json
import numpy as np
import yaml
from itertools import combinations


# ============================================================
#  配置读取
# ============================================================

def get_project_root():
    """Find the project root from the script location.

    The script may be placed either in the project root or in a script/
    subdirectory on HPC.  Use the nearest parent containing project inputs.
    """
    start = os.path.abspath(os.path.dirname(__file__))
    current = start
    while True:
        if (
            os.path.isdir(os.path.join(current, 'calculate')) or
            os.path.exists(os.path.join(current, 'config.yaml')) or
            os.path.exists(os.path.join(current, 'input'))
        ):
            return current
        parent = os.path.dirname(current)
        if parent == current:
            return start
        current = parent


def load_config(project_root):
    """
    加载配置，优先读取 config.yaml，回退到旧格式 (input + prepare/input.json + formation.out)。
    返回统一的 dict。
    """
    yaml_path = os.path.join(project_root, 'config.yaml')
    if os.path.exists(yaml_path):
        with open(yaml_path, 'r', encoding='utf-8') as f:
            cfg = yaml.safe_load(f)
        print(f"已加载配置: {yaml_path}")
        return cfg

    # ---- 回退: 从旧格式构建等效 config ----
    print("未找到 config.yaml，尝试从旧格式 (input / input.json / formation.out) 加载...")
    cfg = _load_legacy_config(project_root)
    return cfg


def _load_legacy_config(project_root):
    """从旧的 input + input.json + formation.out 构建等效配置 dict"""
    cfg = {}

    # ---- input.json ----
    json_path = os.path.join(project_root, 'prepare', 'input.json')
    if os.path.exists(json_path):
        with open(json_path, 'r') as f:
            jdata = json.load(f)
    else:
        jdata = {}

    # job
    js = jdata.get('job_set', {})
    cfg['job'] = {
        'queue': js.get('queue', 'slurm'),
        'partition': js.get('partition', ''),
        'nodes': js.get('nodes', 1),
        'ntasks_per_node': js.get('ntasks-per-node', 4),
        'ngpu': js.get('ngpu', 4),
    }

    # calculate
    cs = jdata.get('calculate_set', {})
    cfg['calculate'] = {
        'psp_path': cs.get('psp_path', ''),
        'psp_name': cs.get('psp_name', []),
        'e_cut': cs.get('e_cut', 60),
        'spin': cs.get('spin', 1),
        'level': cs.get('level', 1),
    }

    # system / bulk
    bs = jdata.get('bulk_set', {})
    si = jdata.get('structure_information', {})
    elements = si.get('defect_elements', [])
    chem_list = si.get('chemical_potential', [])
    chem_pot = {}
    for i, elem in enumerate(elements):
        if i < len(chem_list):
            chem_pot[elem] = float(chem_list[i])

    cfg['system'] = {
        'name': '',
        'num_atoms': [int(x) for x in bs.get('num_atoms', ['64', '320'])],
        'dielectric': float(bs.get('diec', 13.0)),
        'refractive': float(bs.get('refractive', 3.9766)),
        'effective_mass': {
            'electron': float(bs.get('effectivemass', ['1.0', '1.0'])[0]),
            'hole': float(bs.get('effectivemass', ['1.0', '1.0'])[1]),
        },
        'gap': float(bs.get('gap', 0.0)),
        'elements': elements,
        'chemical_potential': chem_pot,
    }

    # structure
    cfg['structure'] = {
        'defect_center': si.get('defect_center', []),
        'bulk_elements': si.get('bulk_elements', []),
        'defect_elements': elements,
    }

    # formation_energy
    fe = jdata.get('formation_energy_calculation', {})
    cfg['formation_energy'] = {
        'charge': [int(x) for x in fe.get('charge', [])],
        'e_corr': fe.get('e_corr', 'T') in ('T', 'True', True),
        't_range': [float(x) for x in fe.get('t_range', ['300', '1000'])],
        'fermi_range': [float(x) for x in fe.get('ferimi_range', [])],
        'mu_range': [float(x) for x in fe.get('mu_range', [])],
    }

    # band_alignment (placeholder)
    cfg['band_alignment'] = {'E_VBM': 0.0}

    # ---- formation.out 补充 VBM / Gap / E_bulk ----
    formation_path = os.path.join(project_root, 'result', 'formation.out')
    if os.path.exists(formation_path):
        vbm, gap = _read_vbm_gap_from_formation(formation_path)
        if vbm is not None:
            cfg.setdefault('_bulk', {})['VBM'] = vbm
        if gap is not None:
            cfg['system']['gap'] = gap
        with open(formation_path, 'r') as f:
            for line in f:
                if 'bulk Energy' in line:
                    m = re.search(r'([-\d.]+)', line.split(':')[-1])
                    if m:
                        cfg.setdefault('_bulk', {})['E_bulk'] = float(m.group(1))

    return cfg


# ============================================================
#  工具函数
# ============================================================

def read_atom_config(filepath):
    """读取 atom.config / final.config，返回 (原子数, 晶格矢量, [(元素序号, x, y, z), ...])"""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    num_atoms = int(lines[0].strip().split()[0])
    lattice = []
    lat_start = None
    for i, line in enumerate(lines):
        if 'lattice' in line.lower():
            lat_start = i + 1
            break
    for i in range(lat_start, lat_start + 3):
        parts = lines[i].split()
        lattice.append([float(parts[0]), float(parts[1]), float(parts[2])])
    lattice = np.array(lattice)
    atoms = []
    pos_start = None
    for i, line in enumerate(lines):
        if 'position' in line.lower():
            pos_start = i + 1
            break
    for i in range(pos_start, pos_start + num_atoms):
        parts = lines[i].split()
        elem = int(parts[0])
        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        atoms.append((elem, x, y, z))
    return num_atoms, lattice, atoms


def find_vacancy_position(bulk_atoms, defect_atoms, lattice):
    G = np.dot(lattice, lattice.T)
    def frac_distance(f1, f2):
        delta = np.array(f2) - np.array(f1)
        delta -= np.round(delta)
        return np.sqrt(np.dot(delta, np.dot(G, delta)))
    defect_coords = [(a[1], a[2], a[3]) for a in defect_atoms]
    vacancies = []
    for batom in bulk_atoms:
        bcoord = (batom[1], batom[2], batom[3])
        matched = any(frac_distance(bcoord, dc) < 0.1 for dc in defect_coords)
        if not matched:
            vacancies.append(bcoord)
    return vacancies


def find_defect_position(bulk_config_path, defect_config_path):
    n_bulk, lat_bulk, atoms_bulk = read_atom_config(bulk_config_path)
    n_def, lat_def, atoms_def = read_atom_config(defect_config_path)
    G = np.dot(lat_bulk, lat_bulk.T)
    def frac_distance(f1, f2):
        delta = np.array(f2) - np.array(f1)
        delta -= np.round(delta)
        return np.sqrt(np.dot(delta, np.dot(G, delta)))
    def count_by_element(atoms):
        counts = {}
        for atom in atoms:
            counts[atom[0]] = counts.get(atom[0], 0) + 1
        return counts
    if n_bulk > n_def:
        vacs = find_vacancy_position(atoms_bulk, atoms_def, lat_bulk)
        if vacs:
            return vacs[0]
    elif n_bulk < n_def:
        extras = find_vacancy_position(atoms_def, atoms_bulk, lat_def)
        if extras:
            return extras[0]
    else:
        # Substitution with preserved atom order: one line changes element.
        changed_sites = [
            (ba, da) for ba, da in zip(atoms_bulk, atoms_def)
            if ba[0] != da[0]
        ]
        if len(changed_sites) == 1:
            ba, _ = changed_sites[0]
            return (ba[1], ba[2], ba[3])

        # Substitution after atom reordering: infer the removed site from
        # element counts, then find the bulk atom of the removed species that
        # no longer has a nearby same-element counterpart in the defect cell.
        bulk_counts = count_by_element(atoms_bulk)
        defect_counts = count_by_element(atoms_def)
        all_elements = set(bulk_counts) | set(defect_counts)
        removed = [
            z for z in all_elements
            if defect_counts.get(z, 0) - bulk_counts.get(z, 0) == -1
        ]
        added = [
            z for z in all_elements
            if defect_counts.get(z, 0) - bulk_counts.get(z, 0) == 1
        ]
        if len(removed) == 1 and len(added) == 1:
            removed_z = removed[0]
            defect_same = [
                (a[1], a[2], a[3]) for a in atoms_def if a[0] == removed_z
            ]
            candidates = []
            for ba in atoms_bulk:
                if ba[0] != removed_z:
                    continue
                bc = (ba[1], ba[2], ba[3])
                nearest = min(
                    (frac_distance(bc, dc) for dc in defect_same),
                    default=float('inf'))
                candidates.append((nearest, bc))
            if candidates:
                nearest, coord = max(candidates, key=lambda item: item[0])
                if nearest > 0.5:
                    return coord

        for ba in atoms_bulk:
            bc = np.array([ba[1], ba[2], ba[3]])
            for da in atoms_def:
                dc = np.array([da[1], da[2], da[3]])
                if frac_distance(bc, dc) < 0.1 and ba[0] != da[0]:
                    return (ba[1], ba[2], ba[3])
    return None


def get_charge_dirs(defect_dir):
    charge_dirs = {}
    for item in sorted(os.listdir(defect_dir)):
        m = re.match(r'^q_([-+]?\d+)$', item)
        if m:
            q = int(m.group(1))
            charge_dirs[q] = os.path.join(defect_dir, item)
    return charge_dirs


def format_charge(q):
    return f"+{q}" if q > 0 else str(q)


def format_float_list(values, ndigits=6):
    return '[' + ', '.join(f'{float(v):.{ndigits}f}' for v in values) + ']'


def charge_input_name(q):
    return f"defect_{format_charge(q)}.input"


def parse_defect_input(defect_input_path):
    """Parse paths and defect coordinate from an ICIC defect.input file."""
    info = {}
    if not defect_input_path or not os.path.exists(defect_input_path):
        return info
    with open(defect_input_path, 'r') as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            key = parts[0]
            if key == 'defect' and len(parts) >= 4:
                try:
                    info['defect_pos'] = tuple(float(x) for x in parts[1:4])
                except ValueError:
                    pass
            elif key in ('bulk', 'neutral') and len(parts) >= 2:
                info[key] = parts[1]
            elif key == 'charged' and len(parts) >= 3 and parts[1] == 'state':
                info['charge_state'] = parts[2]
    return info


def defect_index(defect_name):
    m = re.match(r'^(\d+)(?:_|$)', str(defect_name))
    return m.group(1) if m else ''


def parse_defect_filter(args):
    filters = []
    i = 0
    while i < len(args):
        arg = args[i]
        if arg in ('--idx', '--defect', '--only'):
            if i + 1 >= len(args):
                raise ValueError(f"{arg} 需要一个参数")
            values = args[i + 1]
            filters.extend(v.strip() for v in values.split(',') if v.strip())
            i += 2
        elif arg.startswith('--idx=') or arg.startswith('--defect=') or arg.startswith('--only='):
            values = arg.split('=', 1)[1]
            filters.extend(v.strip() for v in values.split(',') if v.strip())
            i += 1
        else:
            filters.extend(v.strip() for v in arg.split(',') if v.strip())
            i += 1
    return set(filters)


def normalize_correction_mode(mode):
    """Normalize user/config correction mode names.

    Returned values:
      both: use ICIC + PA
      icic: use image-charge correction only
      pa: use potential alignment only
      none: use raw energies without correction
    """
    if mode is None:
        return None
    key = str(mode).strip().lower().replace('_', '-')
    aliases = {
        'both': 'both',
        'all': 'both',
        'full': 'both',
        'icic+pa': 'both',
        'pa+icic': 'both',
        'icic-pa': 'both',
        'pa-icic': 'both',
        'true': 'both',
        't': 'both',
        'yes': 'both',
        'y': 'both',
        '1': 'both',
        'icic': 'icic',
        'image': 'icic',
        'image-charge': 'icic',
        'image-correction': 'icic',
        'mirror': 'icic',
        'pa': 'pa',
        'potential': 'pa',
        'potential-alignment': 'pa',
        'alignment': 'pa',
        'delta-v': 'pa',
        'dv': 'pa',
        'none': 'none',
        'raw': 'none',
        'off': 'none',
        'false': 'none',
        'f': 'none',
        'no': 'none',
        'n': 'none',
        '0': 'none',
    }
    if key not in aliases:
        raise ValueError(
            f"未知修正模式: {mode}  (可用: both / icic / pa / none)"
        )
    return aliases[key]


def get_correction_mode(cfg, override=None):
    """Return correction mode from CLI override or config.yaml."""
    if override is not None:
        return normalize_correction_mode(override)

    fe = cfg.get('formation_energy', {}) or {}
    e_corr = fe.get('e_corr', True)
    if normalize_correction_mode(e_corr) == 'none':
        return 'none'

    mode = (
        fe.get('correction_mode',
        fe.get('correction',
        fe.get('correction_type',
        fe.get('correction_method', 'both'))))
    )
    return normalize_correction_mode(mode)


def parse_cli_options(args):
    """Parse common options for prepare/collect/all.

    Unknown positional tokens are kept as defect filters for backward
    compatibility, e.g. ``collect 108``.
    """
    filters = []
    correction_mode = None
    i = 0
    while i < len(args):
        arg = args[i]
        if arg in ('--idx', '--defect', '--only'):
            if i + 1 >= len(args):
                raise ValueError(f"{arg} 需要一个参数")
            values = args[i + 1]
            filters.extend(v.strip() for v in values.split(',') if v.strip())
            i += 2
        elif arg.startswith('--idx=') or arg.startswith('--defect=') or arg.startswith('--only='):
            values = arg.split('=', 1)[1]
            filters.extend(v.strip() for v in values.split(',') if v.strip())
            i += 1
        elif arg in ('--corr', '--correction', '--correction-mode'):
            if i + 1 >= len(args):
                raise ValueError(f"{arg} 需要一个参数")
            correction_mode = normalize_correction_mode(args[i + 1])
            i += 2
        elif arg.startswith('--corr=') or arg.startswith('--correction=') or arg.startswith('--correction-mode='):
            correction_mode = normalize_correction_mode(arg.split('=', 1)[1])
            i += 1
        else:
            filters.extend(v.strip() for v in arg.split(',') if v.strip())
            i += 1

    return set(filters), correction_mode


def match_defect_filter(defect_name, filters):
    if not filters:
        return True
    return defect_name in filters or defect_index(defect_name) in filters


def read_etot_from_report(report_path):
    etot = None
    with open(report_path, 'r') as f:
        for line in f:
            if 'E_tot(eV)' in line and 'Ryd' not in line:
                for p in line.split():
                    try:
                        val = float(p)
                        if abs(val) > 100:
                            etot = val
                    except ValueError:
                        continue
    return etot


def read_etot_from_atom_config(config_path):
    with open(config_path, 'r') as f:
        line = f.readline()
    m = re.search(r'Etot,Ep,Ek\s*=\s*([-\d.E+]+)', line)
    if m:
        return float(m.group(1))
    return None


def read_ecoul_from_report(report_path):
    ecoul = None
    with open(report_path, 'r') as f:
        for line in f:
            if 'E_Coul(eV)' in line:
                for p in line.split():
                    try:
                        ecoul = float(p)
                    except ValueError:
                        continue
    return ecoul


def periodic_distance(frac1, frac2, lattice):
    G = np.dot(np.array(lattice), np.array(lattice).T)
    delta = np.array(frac1) - np.array(frac2)
    delta -= np.round(delta)
    return float(np.sqrt(np.dot(delta, np.dot(G, delta))))


def read_vatom_records(vatom_path, config_path):
    """Read OUT.VATOM rows with lattice-aware metadata."""
    _, lattice, _ = read_atom_config(config_path)
    with open(vatom_path, 'r') as f:
        lines = f.readlines()
    records = []
    for fallback_idx, line in enumerate(lines[1:], start=1):
        parts = line.split()
        if len(parts) >= 5:
            records.append({
                'potential': float(parts[4]),
                'line_index': fallback_idx,
                'element_z': int(parts[0]) if parts[0].lstrip('+-').isdigit() else None,
                'frac_coord': [float(parts[1]), float(parts[2]), float(parts[3])],
            })
    return records, np.array(lattice)


def annotate_pa_record(record, bound, lattice, match_distance=None):
    result = dict(record)
    result['distance_to_bound'] = periodic_distance(record['frac_coord'], bound, lattice)
    if match_distance is not None:
        result['match_distance'] = float(match_distance)
    result['bound'] = np.array(bound).tolist()
    return result


def read_vatom_far_atom(vatom_path, defect_center, config_path):
    """Read the atom nearest to the PA opposite point, without element matching."""
    records, lattice = read_vatom_records(vatom_path, config_path)
    if not records:
        return None
    bound = (np.array(defect_center) + 0.5) % 1.0
    far_record = min(records, key=lambda rec: periodic_distance(rec['frac_coord'], bound, lattice))
    return annotate_pa_record(far_record, bound, lattice)


def match_pa_atoms_by_element(neutral_vatom, neutral_config, bulk_vatom, bulk_config,
                              defect_center, max_match_distance=2.0):
    """Choose PA atoms with the same element in neutral and bulk cells.

    First choose neutral candidates by distance to the opposite point
    ``bound = defect + 0.5``.  For each neutral candidate, find the nearest bulk
    atom with the same element.  The first candidate with a reasonable match is
    used, preventing O-Hf local-potential comparisons.
    """
    neutral_records, neutral_lattice = read_vatom_records(neutral_vatom, neutral_config)
    bulk_records, bulk_lattice = read_vatom_records(bulk_vatom, bulk_config)
    if not neutral_records or not bulk_records:
        return None, None

    bound = (np.array(defect_center) + 0.5) % 1.0
    neutral_candidates = sorted(
        neutral_records,
        key=lambda rec: periodic_distance(rec['frac_coord'], bound, neutral_lattice),
    )
    for neutral_rec in neutral_candidates:
        z = neutral_rec.get('element_z')
        if z is None:
            continue
        same_element_bulk = [rec for rec in bulk_records if rec.get('element_z') == z]
        if not same_element_bulk:
            continue
        bulk_rec = min(
            same_element_bulk,
            key=lambda rec: periodic_distance(rec['frac_coord'], neutral_rec['frac_coord'], bulk_lattice),
        )
        match_distance = periodic_distance(
            bulk_rec['frac_coord'], neutral_rec['frac_coord'], bulk_lattice)
        if max_match_distance is not None and match_distance > max_match_distance:
            continue
        neutral_info = annotate_pa_record(neutral_rec, bound, neutral_lattice)
        bulk_info = annotate_pa_record(bulk_rec, bound, bulk_lattice, match_distance)
        return neutral_info, bulk_info

    return None, None


def read_vatom_farthest(vatom_path, defect_center, config_path):
    info = read_vatom_far_atom(vatom_path, defect_center, config_path)
    return info['potential'] if info else None


def _read_vbm_gap_from_formation(formation_path):
    vbm = None
    gap = None
    with open(formation_path, 'r') as f:
        for line in f:
            if 'VBM' in line and 'value' in line:
                m = re.search(r'([-\d.]+)\s*$', line.strip())
                if m:
                    vbm = float(m.group(1))
            if 'Gap' in line and 'value' in line:
                m = re.search(r'([-\d.]+)\s*$', line.strip())
                if m:
                    gap = float(m.group(1))
    return vbm, gap


# ── 元素序号 → 符号 ──
Z_TO_SYMBOL = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O',
    9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',
    16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti',
    23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu',
    30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr',
    37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 44: 'Ru',
    45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn', 51: 'Sb',
    52: 'Te', 53: 'I', 55: 'Cs', 56: 'Ba', 57: 'La', 72: 'Hf', 73: 'Ta',
    74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg',
    81: 'Tl', 82: 'Pb', 83: 'Bi',
}


def get_bulk_info(project_root, cfg):
    """
    获取 bulk 信息: E_bulk, VBM, Gap, chemical_potential。
    优先从 config.yaml 的 _bulk 字段读取（旧格式回退时填充），
    否则从 calculate/bulk/scf 读取。
    """
    info = {}
    # 从 cfg 中获取（旧格式回退时 _load_legacy_config 已填充）
    _b = cfg.get('_bulk', {})
    if 'E_bulk' in _b:
        info['E_bulk'] = _b['E_bulk']
    if 'VBM' in _b:
        info['VBM'] = _b['VBM']

    system_cfg = cfg.get('system', {})
    band_cfg = cfg.get('band_alignment', {})

    vbm = system_cfg.get('VBM', system_cfg.get('E_VBM'))
    if vbm is not None:
        info['VBM'] = float(vbm)
    elif 'E_VBM' in band_cfg:
        info['VBM'] = float(band_cfg['E_VBM'])

    gap = system_cfg.get('gap', 0.0)
    if gap and gap > 0:
        info['Gap'] = float(gap)

    # 也尝试从 formation.out 读取（兼容）
    formation_path = os.path.join(project_root, 'result', 'formation.out')
    if os.path.exists(formation_path):
        vbm, fgap = _read_vbm_gap_from_formation(formation_path)
        if 'VBM' not in info and vbm is not None:
            info['VBM'] = vbm
        if 'Gap' not in info and fgap is not None:
            info['Gap'] = fgap
        with open(formation_path, 'r') as f:
            for line in f:
                if 'bulk Energy' in line:
                    m = re.search(r'([-\d.]+)', line.split(':')[-1])
                    if m and 'E_bulk' not in info:
                        info['E_bulk'] = float(m.group(1))

    # 从 REPORT / atom.config 读取
    bulk_report = os.path.join(project_root, 'calculate', 'bulk', 'scf', 'REPORT')
    if 'E_bulk' not in info and os.path.exists(bulk_report):
        info['E_bulk'] = read_etot_from_report(bulk_report)
    bulk_config = os.path.join(project_root, 'calculate', 'bulk', 'scf', 'atom.config')
    if 'E_bulk' not in info and os.path.exists(bulk_config):
        info['E_bulk'] = read_etot_from_atom_config(bulk_config)

    # chemical_potential
    info['chemical_potential'] = cfg.get('system', {}).get('chemical_potential', {})
    dielectric = cfg.get('system', {}).get('dielectric')
    if dielectric is not None:
        try:
            info['dielectric'] = float(dielectric)
        except (TypeError, ValueError):
            info['dielectric'] = None

    return info


def parse_atom_change(defect_name, bulk_config_path, defect_config_path):
    """比较 bulk/defect atom.config，返回 {元素符号: Δn}"""
    if os.path.exists(bulk_config_path) and os.path.exists(defect_config_path):
        _, _, bulk_atoms = read_atom_config(bulk_config_path)
        _, _, defect_atoms = read_atom_config(defect_config_path)
        bulk_count = {}
        for a in bulk_atoms:
            bulk_count[a[0]] = bulk_count.get(a[0], 0) + 1
        defect_count = {}
        for a in defect_atoms:
            defect_count[a[0]] = defect_count.get(a[0], 0) + 1
        all_z = set(list(bulk_count.keys()) + list(defect_count.keys()))
        delta = {}
        for z in all_z:
            diff = defect_count.get(z, 0) - bulk_count.get(z, 0)
            if diff != 0:
                sym = Z_TO_SYMBOL.get(z, str(z))
                delta[sym] = diff
        return delta

    # 从名称推断
    delta = {}
    parts = defect_name.split('_')
    if len(parts) >= 3:
        if parts[1] == 'v':
            delta[parts[2]] = -1
        elif parts[-1] == 'i':
            delta[parts[1]] = 1
        else:
            delta[parts[1]] = 1
            delta[parts[2]] = -1
    return delta


def compute_formation_energy(E_corrected, E_bulk, delta_n, q, VBM, chem_pot):
    """
    E_f(q, E_F) = E_corrected - E_bulk + Σ(-Δn_i)·μ_i + q·(VBM + E_F)
    返回 E_f(q, E_F=0)，即 E_F 相对 VBM 为零时的形成能。
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


def read_formation_out_energies(project_root):
    """从旧 formation.out 读取已有能量（兼容回退）"""
    formation_path = os.path.join(project_root, 'result', 'formation.out')
    defects = {}
    if not os.path.exists(formation_path):
        return defects
    current_defect = None
    with open(formation_path, 'r') as f:
        for line in f:
            if line.startswith('defect_type:'):
                current_defect = line.split(':')[1].strip()
                defects[current_defect] = {'charges': {}}
            if current_defect and re.match(
                    r'^[-+]?\d+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+\s+[-\d.]+', line.strip()):
                parts = line.split()
                if len(parts) >= 5:
                    q = int(parts[0].replace('+', ''))
                    defects[current_defect]['charges'][q] = {
                        'E_raw': float(parts[1]),
                    }
    return defects


# ============================================================
#  YAML 输出
# ============================================================

def _float_representer(dumper, value):
    """控制 YAML 中浮点数的格式"""
    if value != value:  # NaN
        return dumper.represent_scalar('tag:yaml.org,2002:float', '.nan')
    if value == float('inf'):
        return dumper.represent_scalar('tag:yaml.org,2002:float', '.inf')
    if value == float('-inf'):
        return dumper.represent_scalar('tag:yaml.org,2002:float', '-.inf')
    s = f"{value:.5f}"
    return dumper.represent_scalar('tag:yaml.org,2002:float', s)


class CleanDumper(yaml.SafeDumper):
    """自定义 Dumper: 浮点数保留 5 位小数，不输出 Python 对象标签"""
    pass

CleanDumper.add_representer(float, _float_representer)
# numpy float → Python float
CleanDumper.add_representer(np.float64,
    lambda d, v: _float_representer(d, float(v)))
CleanDumper.add_representer(np.int64,
    lambda d, v: d.represent_int(int(v)))


def build_output_yaml(cfg, bulk_info, all_results, defect_delta_n, correction_mode=None):
    """
    构建与 CCD 脚本兼容的输出 YAML dict。

    结构:
      system:
        name, E_VBM, E_CBM, E_bulk, chemical_potential
      defects:
        - name, atom_change, charge_states[], transition_levels{}
          charge_states:
            - q, label, E_f0, E_raw, E_corrected, corrections{}
          transition_levels:
            "(q1/q2)": value
          formation_energy_table:
            E_F_values: [...]
            E_f_min: [...]
            stable_q: [...]
      band_alignment: (从 config 透传)
    """
    E_bulk = bulk_info.get('E_bulk')
    VBM = bulk_info.get('VBM')
    Gap = bulk_info.get('Gap')
    chem_pot = bulk_info.get('chemical_potential', {})
    dielectric = bulk_info.get('dielectric')

    out = {}

    # ── system ──
    sys_name = cfg.get('system', {}).get('name', '')
    out['system'] = {
        'name': sys_name,
        'E_VBM': float(VBM) if VBM is not None else None,
        'E_CBM': float(VBM + Gap) if (VBM is not None and Gap is not None) else None,
        'gap': float(Gap) if Gap is not None else None,
        'E_bulk': float(E_bulk) if E_bulk is not None else None,
        'chemical_potential': {k: float(v) for k, v in chem_pot.items()},
        'dielectric': cfg.get('system', {}).get('dielectric'),
    }
    if correction_mode is not None:
        out['system']['correction_mode'] = correction_mode

    # ── band_alignment (透传) ──
    ba = cfg.get('band_alignment', {})
    out['band_alignment'] = {k: float(v) if isinstance(v, (int, float)) else v
                             for k, v in ba.items()}

    # ── defects ──
    defects_list = []
    for defect_name, charge_data in all_results.items():
        delta_n = defect_delta_n.get(defect_name, {})
        sorted_qs = sorted(charge_data.keys())

        # charge_states
        cs_list = []
        formation_energies = {}
        for q in sorted_qs:
            d = charge_data[q]
            E_f0 = compute_formation_energy(
                d['E_corrected'], E_bulk, delta_n, q, VBM, chem_pot)
            formation_energies[q] = E_f0

            cs_entry = {
                'q': int(q),
                'label': f"q={format_charge(q)}",
                'E_f0': float(E_f0) if E_f0 is not None else None,
                'E_raw': float(d['E_raw']) if d['E_raw'] is not None else None,
                'E_corrected': float(d['E_corrected']) if d['E_corrected'] is not None else None,
                'corrections': {
                    'correction_mode': d.get('correction_mode', correction_mode),
                    'E_icic': float(d['E_icic']) if d['E_icic'] is not None else None,
                    'E_icic_raw': float(d['E_icic_raw']) if d.get('E_icic_raw') is not None else None,
                    'E_icic_dielectric': float(d['E_icic_dielectric']) if d.get('E_icic_dielectric') is not None else None,
                    'delta_V_PA': float(d['delta_V']) if d['delta_V'] is not None else None,
                    'E_corr_total': float(d['E_corr']),
                },
            }
            pa_info = d.get('pa_info')
            if pa_info:
                cs_entry['corrections']['PA_site'] = pa_info
            cs_list.append(cs_entry)

        # transition_levels
        tl = {}
        for i in range(len(sorted_qs)):
            for j in range(i + 1, len(sorted_qs)):
                q1, q2 = sorted_qs[i], sorted_qs[j]
                Ef1, Ef2 = formation_energies.get(q1), formation_energies.get(q2)
                if Ef1 is not None and Ef2 is not None and (q2 - q1) != 0:
                    eps = (Ef1 - Ef2) / (q2 - q1)
                    label = f"({format_charge(q1)}/{format_charge(q2)})"
                    tl[label] = float(eps)

        # formation_energy_table (包络线)
        fe_table = None
        if Gap is not None and any(v is not None for v in formation_energies.values()):
            n_pts = 101
            ef_vals = np.linspace(0, Gap, n_pts)
            ef_min_list = []
            stable_q_list = []
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
                'E_F_values': [float(round(x, 5)) for x in ef_vals],
                'E_f_min': [float(round(x, 5)) for x in ef_min_list],
                'stable_q': stable_q_list,
            }

        defect_entry = {
            'name': defect_name,
            'atom_change': delta_n if delta_n else None,
            'charge_states': cs_list,
            'transition_levels': tl if tl else None,
        }
        if fe_table:
            defect_entry['formation_energy_table'] = fe_table

        defects_list.append(defect_entry)

    out['defects'] = defects_list
    return out


class FlowList(list):
    """标记为 flow style 输出的 list"""
    pass

def _flow_list_representer(dumper, data):
    return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)

CleanDumper.add_representer(FlowList, _flow_list_representer)


def write_output_yaml(out_dict, output_path):
    """写入 YAML 输出文件，formation_energy_table 用 flow style 紧凑输出"""
    # 将 formation_energy_table 中的长列表转为 FlowList
    for dentry in out_dict.get('defects', []):
        ft = dentry.get('formation_energy_table')
        if ft:
            ft['E_F_values'] = FlowList(ft['E_F_values'])
            ft['E_f_min'] = FlowList(ft['E_f_min'])
            ft['stable_q'] = FlowList(ft['stable_q'])

    with open(output_path, 'w', encoding='utf-8') as f:
        f.write("# ============================================================\n")
        f.write("#  缺陷形成能与转变能级计算结果\n")
        f.write("#  由 batch_correction.py collect 自动生成\n")
        f.write("#  格式兼容 CCD 脚本输入\n")
        f.write("# ============================================================\n\n")
        yaml.dump(out_dict, f, Dumper=CleanDumper, default_flow_style=False,
                  allow_unicode=True, sort_keys=False, width=200)
    print(f"  YAML 结果已写入: {output_path}")


# ============================================================
#  Phase 1: Prepare
# ============================================================

def prepare(project_root, cfg, defect_filter=None):
    """生成 defect.input 并输出修正计算指南"""
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
    if defect_filter:
        defect_dirs = [
            (name, path) for name, path in defect_dirs
            if match_defect_filter(name, defect_filter)
        ]

    if not defect_dirs:
        print("错误: calculate/ 下未找到匹配的缺陷目录")
        return

    print(f"{'='*60}")
    print(f"  批量生成 defect.input")
    print(f"{'='*60}")
    print(f"检测到 {len(defect_dirs)} 个缺陷目录")
    print()

    total_generated = 0
    run_commands = []

    for defect_name, defect_dir in defect_dirs:
        print(f"--- 缺陷: {defect_name} ---")
        charge_dirs = get_charge_dirs(defect_dir)
        if not charge_dirs:
            print(f"  未找到 q_* 目录，跳过")
            continue

        charged_states = [q for q in charge_dirs if q != 0]
        if not charged_states:
            print(f"  仅有 q_0，无需修正，跳过")
            continue

        neutral_scf = os.path.join(defect_dir, 'q_0', 'scf')
        if not os.path.isdir(neutral_scf):
            print(f"  警告: 中性态目录不存在: {neutral_scf}，跳过")
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
            print(f"  错误: 无法确定缺陷位置，跳过")
            continue

        bound_pos = tuple((c + 0.5) % 1.0 for c in defect_pos)
        print(f"  缺陷位置: {defect_pos[0]:.6f} {defect_pos[1]:.6f} {defect_pos[2]:.6f}")
        print(f"  带电态: {sorted(charged_states)}")

        for q in sorted(charged_states):
            qdir = charge_dirs[q]
            scf_dir = os.path.join(qdir, 'scf')
            if not os.path.isdir(scf_dir):
                continue
            charge_str = format_charge(q)
            content = (
                f"defect {defect_pos[0]:.6f} {defect_pos[1]:.6f} {defect_pos[2]:.6f}\n"
                f"bound {bound_pos[0]:.6f} {bound_pos[1]:.6f} {bound_pos[2]:.6f}\n"
                f"charged state {charge_str}\n"
                f"bulk {os.path.abspath(bulk_scf)}\n"
                f"neutral {os.path.abspath(neutral_scf)}\n"
            )
            output_paths = [
                os.path.join(scf_dir, 'defect.input'),
                os.path.join(neutral_scf, charge_input_name(q)),
            ]
            output_path = output_paths[0]
            if q < 0:
                output_paths.append(os.path.join(neutral_scf, f"defect_m{abs(q)}.input"))
            elif q > 0:
                output_paths.append(os.path.join(neutral_scf, f"defect_p{q}.input"))
            seen_paths = []
            for path in output_paths:
                if path not in seen_paths:
                    seen_paths.append(path)
            for path in seen_paths:
                with open(path, 'w') as f:
                    f.write(content)
            print(f"  已生成: {', '.join(seen_paths)}")
            total_generated += 1
            run_commands.append((defect_name, charge_str, scf_dir))

    print(f"\n共生成 {total_generated} 个 defect.input 文件")

    if run_commands:
        print(f"\n接下来请在 HPC 上依次执行:")
        script_root = os.path.abspath(os.path.dirname(__file__))
        for defect_name, charge_str, scf_dir in run_commands:
            print(f"\n  # {defect_name} q={charge_str}")
            print(f"  cd {scf_dir}")
            neutral_dir = os.path.join(os.path.dirname(os.path.dirname(scf_dir)), 'q_0', 'scf')
            print(f"  cp {neutral_dir}/OUT.OCC OUT.OCC0")
            print(f"  cp OUT.OCC OUT.OCC1")
            print(f"  bash {script_root}/1_get_rho.sh")
            print(f"  bash {script_root}/2_coulomb_integral.sh")
            print(f"  bash {script_root}/3_get_results.sh")
        print(f"\n全部修正完成后，运行: python {script_root}/batch_correction.py collect")


# ============================================================
#  Phase 2: Collect
# ============================================================

def collect_correction(scf_dir, q, defect_pos, neutral_scf=None, bulk_scf=None, dielectric=None):
    """收集 ICIC + PA 修正。

    返回 (E_icic, delta_V, E_icic_raw, E_icic_dielectric)，严格使用图中公式：
    E_icic_raw = E_inf - E_P，E_icic = E_icic_raw / dielectric。
    如果 dielectric 缺失或无效，不使用未屏蔽的 E_icic_raw 作为最终修正。

    新流程会在 neutral q_0/scf 目录下生成 image-corr_*；旧流程可能在
    charged q_*/scf 目录下生成。这里同时兼容两种目录布局。
    """
    sign = '+' if q > 0 else '-'
    q_abs = abs(q)
    corr_name = f'image-corr_{sign}{q_abs}'

    E_icic = None
    E_icic_raw = None
    delta_V = None
    if dielectric is not None and dielectric <= 0:
        dielectric = None

    defect_input_path = os.path.join(scf_dir, 'defect.input')
    neutral_charge_input = (
        os.path.join(neutral_scf, charge_input_name(q)) if neutral_scf else None
    )
    neutral_defect_input = (
        os.path.join(neutral_scf, 'defect.input') if neutral_scf else None
    )
    neutral_dir = bulk_dir = None
    input_for_paths = None
    for candidate in [neutral_charge_input, defect_input_path, neutral_defect_input]:
        if candidate and os.path.exists(candidate):
            input_for_paths = candidate
            break
    if input_for_paths:
        parsed_input = parse_defect_input(input_for_paths)
        neutral_dir = parsed_input.get('neutral')
        bulk_dir = parsed_input.get('bulk')
        defect_pos = parsed_input.get('defect_pos', defect_pos)
    if neutral_dir is None:
        neutral_dir = neutral_scf
    if bulk_dir is None:
        bulk_dir = bulk_scf

    corr_candidates = [os.path.join(scf_dir, corr_name)]
    if neutral_dir:
        corr_candidates.insert(0, os.path.join(neutral_dir, corr_name))
    seen = set()
    for corr_dir in corr_candidates:
        if corr_dir in seen:
            continue
        seen.add(corr_dir)
        report_inf = os.path.join(corr_dir, 'REPORT')
        report_0 = os.path.join(corr_dir, 'REPORT.0')
        if os.path.exists(report_inf) and os.path.exists(report_0):
            E_inf = read_ecoul_from_report(report_inf)
            E_P = read_ecoul_from_report(report_0)
            if E_inf is not None and E_P is not None:
                E_icic_raw = E_inf - E_P
                if dielectric:
                    E_icic = E_icic_raw / dielectric
                else:
                    E_icic = None
                break

    pa_info = None
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
            neutral_pa, bulk_pa = match_pa_atoms_by_element(
                vatom_n, cfg_n, vatom_b, cfg_b, defect_pos)
            pa_method = 'same_element_nearest'
            if neutral_pa and bulk_pa:
                delta_V = neutral_pa['potential'] - bulk_pa['potential']
                pa_info = {
                    'method': pa_method,
                    'defect_pos': list(defect_pos),
                    'bound': neutral_pa['bound'],
                    'neutral': neutral_pa,
                    'bulk': bulk_pa,
                }

    return E_icic, delta_V, E_icic_raw, dielectric, pa_info


def collect(project_root, cfg, defect_filter=None, correction_mode=None):
    """收集修正结果，计算形成能和转变能级，输出 YAML"""
    calc_dir = os.path.join(project_root, 'calculate')
    correction_mode = get_correction_mode(cfg, correction_mode)
    use_icic = correction_mode in ('both', 'icic')
    use_pa = correction_mode in ('both', 'pa')
    bulk_info = get_bulk_info(project_root, cfg)
    formation_data = read_formation_out_energies(project_root)

    E_bulk = bulk_info.get('E_bulk')
    VBM = bulk_info.get('VBM')
    Gap = bulk_info.get('Gap')
    chem_pot = bulk_info.get('chemical_potential', {})
    dielectric = bulk_info.get('dielectric')

    print(f"{'='*70}")
    print(f"  收集修正结果 & 计算形成能 / 转变能级")
    print(f"{'='*70}")
    print(f"E_bulk = {E_bulk} eV")
    print(f"VBM = {VBM} eV")
    print(f"Gap = {Gap} eV")
    print(f"化学势: {chem_pot}")
    print(f"介电常数: {dielectric}")
    print(f"修正模式: {correction_mode}  ({'ICIC' if use_icic else 'no ICIC'} / {'PA' if use_pa else 'no PA'})")
    print()

    defect_dirs = []
    for item in sorted(os.listdir(calc_dir)):
        full = os.path.join(calc_dir, item)
        if os.path.isdir(full) and item != 'bulk':
            defect_dirs.append((item, full))
    if defect_filter:
        defect_dirs = [
            (name, path) for name, path in defect_dirs
            if match_defect_filter(name, defect_filter)
        ]

    all_results = {}
    defect_delta_n = {}
    bulk_config_path = os.path.join(calc_dir, 'bulk', 'scf', 'atom.config')

    for defect_name, defect_dir in defect_dirs:
        print(f"\n{'─'*50}")
        print(f"缺陷: {defect_name}")
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

        # 原子数变化
        dc_path = defect_config if os.path.exists(defect_config) else \
            os.path.join(neutral_scf, 'final.config')
        delta_n = parse_atom_change(defect_name, bulk_config_path, dc_path)
        defect_delta_n[defect_name] = delta_n
        if delta_n:
            print(f"  原子数变化: {delta_n}")

        charge_data = {}
        for q in sorted(charge_dirs.keys()):
            qdir = charge_dirs[q]
            scf_dir = os.path.join(qdir, 'scf')
            if not os.path.isdir(scf_dir):
                continue

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

            E_icic = delta_V = E_icic_raw = E_icic_dielectric = pa_info = None
            E_corr = 0.0
            if q != 0 and defect_pos is not None:
                E_icic, delta_V, E_icic_raw, E_icic_dielectric, pa_info = collect_correction(
                    scf_dir,
                    q,
                    defect_pos,
                    neutral_scf=neutral_scf,
                    bulk_scf=os.path.join(calc_dir, 'bulk', 'scf'),
                    dielectric=dielectric,
                )
                if use_icic and E_icic is not None:
                    E_corr += E_icic
                if use_pa and delta_V is not None:
                    E_corr += q * delta_V

            charge_data[q] = {
                'E_raw': E_raw,
                'E_icic': E_icic,
                'E_icic_raw': E_icic_raw,
                'E_icic_dielectric': E_icic_dielectric,
                'delta_V': delta_V,
                'pa_info': pa_info,
                'E_corr': E_corr,
                'E_corrected': E_raw + E_corr if E_raw is not None else None,
                'correction_mode': correction_mode,
            }

            # 打印状态
            corr_status = ""
            if q != 0:
                if correction_mode == 'none':
                    corr_status = "无修正(已禁用)"
                elif correction_mode == 'icic':
                    corr_status = "✓ ICIC" if E_icic is not None else "✗ ICIC缺失"
                elif correction_mode == 'pa':
                    corr_status = "✓ PA" if delta_V is not None else "✗ PA缺失"
                else:
                    if E_icic is not None and delta_V is not None:
                        corr_status = "✓ ICIC+PA"
                    elif E_icic is not None:
                        corr_status = "△ 仅ICIC"
                    elif delta_V is not None:
                        corr_status = "△ 仅PA"
                    else:
                        corr_status = "✗ 无修正"
                if delta_V is not None and abs(delta_V) > 1.0:
                    corr_status += " ⚠ PA异常"
            e_s = f"{E_raw:14.5f}" if E_raw is not None else f"{'N/A':>14s}"
            print(f"  q={format_charge(q):>3s}  E_raw={e_s} eV  E_corr={E_corr:10.5f} eV  {corr_status}")
            if q != 0 and delta_V is not None and abs(delta_V) > 1.0 and pa_info:
                n_info = pa_info.get('neutral', {})
                b_info = pa_info.get('bulk', {})
                print(f"       PA警告: ΔV={delta_V:.5f} eV, bound={format_float_list(pa_info.get('bound', []), 5)}")
                print(f"       neutral line {n_info.get('line_index')} Z={n_info.get('element_z')} "
                      f"coord={format_float_list(n_info.get('frac_coord', []), 5)} "
                      f"V={n_info.get('potential'):.6f}")
                print(f"       bulk    line {b_info.get('line_index')} Z={b_info.get('element_z')} "
                      f"coord={format_float_list(b_info.get('frac_coord', []), 5)} "
                      f"V={b_info.get('potential'):.6f}")

        all_results[defect_name] = charge_data

    # ── 计算形成能并打印 ──
    print(f"\n{'='*70}")
    print(f"  形成能 & 转变能级")
    print(f"{'='*70}")

    for defect_name, charge_data in all_results.items():
        delta_n = defect_delta_n.get(defect_name, {})
        sorted_qs = sorted(charge_data.keys())
        print(f"\n--- {defect_name} ---")

        for q in sorted_qs:
            d = charge_data[q]
            E_f0 = compute_formation_energy(
                d['E_corrected'], E_bulk, delta_n, q, VBM, chem_pot)
            q_str = format_charge(q)
            if E_f0 is not None:
                if q == 0:
                    print(f"  E_f(q={q_str}, E_F) = {E_f0:.5f} eV")
                else:
                    sign = "+" if q > 0 else "-"
                    q_abs = abs(q)
                    coeff = "" if q_abs == 1 else f"{q_abs}·"
                    print(f"  E_f(q={q_str}, E_F) = {E_f0:.5f} {sign} {coeff}E_F  (eV)")

        # 转变能级
        formation_energies = {}
        for q in sorted_qs:
            d = charge_data[q]
            formation_energies[q] = compute_formation_energy(
                d['E_corrected'], E_bulk, delta_n, q, VBM, chem_pot)
        for i in range(len(sorted_qs)):
            for j in range(i + 1, len(sorted_qs)):
                q1, q2 = sorted_qs[i], sorted_qs[j]
                Ef1, Ef2 = formation_energies.get(q1), formation_energies.get(q2)
                if Ef1 is not None and Ef2 is not None and (q2 - q1) != 0:
                    eps = (Ef1 - Ef2) / (q2 - q1)
                    in_gap = ""
                    if Gap is not None:
                        in_gap = " [带隙内]" if 0 <= eps <= Gap else " [带隙外]"
                    print(f"  ε({format_charge(q1)}/{format_charge(q2)}) = {eps:.4f} eV{in_gap}")

    # ── 输出 YAML ──
    out_dict = build_output_yaml(
        cfg, bulk_info, all_results, defect_delta_n, correction_mode
    )

    result_dir = os.path.join(project_root, 'result')
    os.makedirs(result_dir, exist_ok=True)
    yaml_path = os.path.join(result_dir, 'defect_results.yaml')
    write_output_yaml(out_dict, yaml_path)

    # 同时输出 E_f vs E_F 表格 (txt，方便绘图)
    ef_dir = os.path.join(result_dir, 'E_forms')
    os.makedirs(ef_dir, exist_ok=True)
    for dentry in out_dict.get('defects', []):
        ft = dentry.get('formation_energy_table')
        if not ft:
            continue
        dname = dentry['name']
        table_path = os.path.join(ef_dir, f'E_formation_corrected_{dname}.txt')
        sorted_qs = sorted([cs['q'] for cs in dentry['charge_states']])
        form_en = {cs['q']: cs['E_f0'] for cs in dentry['charge_states']}
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
        print(f"  E_f 表格已写入: {table_path}")

    print(f"\n{'='*70}")
    print(f"  完成")
    print(f"{'='*70}")


# ============================================================
#  Main
# ============================================================

def main():
    project_root = get_project_root()

    if len(sys.argv) < 2:
        print("用法:")
        print("  python script/batch_correction.py prepare [--idx 108|--defect 108_v_O]")
        print("  python script/batch_correction.py collect [--idx 108|--defect 108_v_O] [--corr both|icic|pa|none]")
        print("  python script/batch_correction.py all [--idx 108|--defect 108_v_O] [--corr both|icic|pa|none]")
        sys.exit(1)

    cfg = load_config(project_root)
    mode = sys.argv[1].lower()
    try:
        defect_filter, correction_mode = parse_cli_options(sys.argv[2:])
    except ValueError as exc:
        print(f"参数错误: {exc}")
        sys.exit(1)
    if defect_filter:
        print(f"仅处理缺陷: {', '.join(sorted(defect_filter))}")
    if correction_mode is not None:
        print(f"命令行指定修正模式: {correction_mode}")

    if mode == 'prepare':
        prepare(project_root, cfg, defect_filter)
    elif mode == 'collect':
        collect(project_root, cfg, defect_filter, correction_mode)
    elif mode == 'all':
        prepare(project_root, cfg, defect_filter)
        print("\n\n")
        collect(project_root, cfg, defect_filter, correction_mode)
    else:
        print(f"未知模式: {mode}  (可用: prepare / collect / all)")
        sys.exit(1)


if __name__ == '__main__':
    main()
