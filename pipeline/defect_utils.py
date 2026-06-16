#!/usr/bin/env python3
"""
缺陷计算共享工具模块。

提供配置加载、原子结构读取、能量解析等公共函数，
供 correction.py 和 formation_energy.py 共同使用。
"""

import os
import re
import json
import numpy as np
import yaml


# ============================================================
#  配置读取
# ============================================================

def get_project_root(ref_file=None):
    """Find the project root from a script/module path.

    Scripts are often copied into a script/ subdirectory on HPC, while the
    project data remain one level above.  Walk upward until project inputs are
    found instead of assuming dirname(__file__) is always the root.
    """
    start = os.path.abspath(os.path.dirname(ref_file or __file__))
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
    """加载 config.yaml，回退到旧格式。"""
    yaml_path = os.path.join(project_root, 'config.yaml')
    if os.path.exists(yaml_path):
        with open(yaml_path, 'r', encoding='utf-8') as f:
            cfg = yaml.safe_load(f)
        print(f"已加载配置: {yaml_path}")
        return cfg
    print("未找到 config.yaml，尝试从旧格式加载...")
    return _load_legacy_config(project_root)


def _load_legacy_config(project_root):
    cfg = {}
    json_path = os.path.join(project_root, 'prepare', 'input.json')
    jdata = {}
    if os.path.exists(json_path):
        with open(json_path, 'r') as f:
            jdata = json.load(f)

    js = jdata.get('job_set', {})
    cfg['job'] = {
        'queue': js.get('queue', 'slurm'),
        'partition': js.get('partition', ''),
        'nodes': js.get('nodes', 1),
        'ntasks_per_node': js.get('ntasks-per-node', 4),
        'ngpu': js.get('ngpu', 4),
    }
    cs = jdata.get('calculate_set', {})
    cfg['calculate'] = {
        'psp_path': cs.get('psp_path', ''),
        'psp_name': cs.get('psp_name', []),
        'e_cut': cs.get('e_cut', 60),
        'spin': cs.get('spin', 1),
        'level': cs.get('level', 1),
    }
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
    cfg['structure'] = {
        'defect_center': si.get('defect_center', []),
        'bulk_elements': si.get('bulk_elements', []),
        'defect_elements': elements,
    }
    fe = jdata.get('formation_energy_calculation', {})
    cfg['formation_energy'] = {
        'charge': [int(x) for x in fe.get('charge', [])],
        'e_corr': fe.get('e_corr', 'T') in ('T', 'True', True),
        't_range': [float(x) for x in fe.get('t_range', ['300', '1000'])],
        'fermi_range': [float(x) for x in fe.get('ferimi_range', [])],
        'mu_range': [float(x) for x in fe.get('mu_range', [])],
    }
    cfg['band_alignment'] = {'E_VBM': 0.0}

    formation_path = os.path.join(project_root, 'result', 'formation.out')
    if os.path.exists(formation_path):
        vbm, gap = read_vbm_gap_from_formation(formation_path)
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
#  原子结构读取
# ============================================================

def read_atom_config(filepath):
    """读取 atom.config / final.config"""
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
        if not any(frac_distance(bcoord, dc) < 0.1 for dc in defect_coords):
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
        return vacs[0] if vacs else None
    elif n_bulk < n_def:
        extras = find_vacancy_position(atoms_def, atoms_bulk, lat_def)
        return extras[0] if extras else None
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


# ============================================================
#  目录扫描与能量读取
# ============================================================

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
    return float(m.group(1)) if m else None


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


def read_vatom_farthest(vatom_path, defect_center, config_path):
    """Read potential at atom nearest to (defect_center + 0.5) mod 1."""
    _, lattice, _ = read_atom_config(config_path)
    G = np.dot(np.array(lattice), np.array(lattice).T)
    with open(vatom_path, 'r') as f:
        lines = f.readlines()
    coords, vatom_vals = [], []
    for line in lines[1:]:
        parts = line.split()
        if len(parts) >= 5:
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
            vatom_vals.append(float(parts[4]))
    min_dist, max_idx = float('inf'), 0
    bound = (np.array(defect_center) + 0.5) % 1.0
    for i, fc in enumerate(coords):
        delta = np.array(fc) - bound
        delta -= np.round(delta)
        dist = np.sqrt(np.dot(delta, np.dot(G, delta)))
        if dist < min_dist:
            min_dist = dist
            max_idx = i
    return vatom_vals[max_idx]


def read_vbm_gap_from_formation(formation_path):
    vbm, gap = None, None
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


def read_formation_out_energies(project_root):
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
                    defects[current_defect]['charges'][q] = {'E_raw': float(parts[1])}
    return defects


# ============================================================
#  Bulk 信息
# ============================================================

def get_bulk_info(project_root, cfg):
    info = {}
    _b = cfg.get('_bulk', {})
    if 'E_bulk' in _b:
        info['E_bulk'] = _b['E_bulk']
    if 'VBM' in _b:
        info['VBM'] = _b['VBM']
    system_cfg = cfg.get('system', {})
    dielectric = system_cfg.get('dielectric')
    if dielectric is not None:
        try:
            info['dielectric'] = float(dielectric)
        except (TypeError, ValueError):
            info['dielectric'] = None
    gap = system_cfg.get('gap', 0.0)
    if gap and gap > 0:
        info['Gap'] = gap

    formation_path = os.path.join(project_root, 'result', 'formation.out')
    if os.path.exists(formation_path):
        vbm, fgap = read_vbm_gap_from_formation(formation_path)
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

    bulk_report = os.path.join(project_root, 'calculate', 'bulk', 'scf', 'REPORT')
    if 'E_bulk' not in info and os.path.exists(bulk_report):
        info['E_bulk'] = read_etot_from_report(bulk_report)
    bulk_config = os.path.join(project_root, 'calculate', 'bulk', 'scf', 'atom.config')
    if 'E_bulk' not in info and os.path.exists(bulk_config):
        info['E_bulk'] = read_etot_from_atom_config(bulk_config)

    info['chemical_potential'] = system_cfg.get('chemical_potential', {})
    return info


# ============================================================
#  元素与原子数变化
# ============================================================

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
                delta[Z_TO_SYMBOL.get(z, str(z))] = diff
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


# ============================================================
#  YAML 输出工具
# ============================================================

def _float_representer(dumper, value):
    if value != value:
        return dumper.represent_scalar('tag:yaml.org,2002:float', '.nan')
    if value == float('inf'):
        return dumper.represent_scalar('tag:yaml.org,2002:float', '.inf')
    if value == float('-inf'):
        return dumper.represent_scalar('tag:yaml.org,2002:float', '-.inf')
    return dumper.represent_scalar('tag:yaml.org,2002:float', f"{value:.5f}")


class CleanDumper(yaml.SafeDumper):
    pass

CleanDumper.add_representer(float, _float_representer)
CleanDumper.add_representer(np.float64, lambda d, v: _float_representer(d, float(v)))
CleanDumper.add_representer(np.int64, lambda d, v: d.represent_int(int(v)))


class FlowList(list):
    pass

def _flow_list_representer(dumper, data):
    return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)

CleanDumper.add_representer(FlowList, _flow_list_representer)


def write_yaml(out_dict, output_path, header_lines=None):
    """写入 YAML 文件"""
    with open(output_path, 'w', encoding='utf-8') as f:
        if header_lines:
            for h in header_lines:
                f.write(f"# {h}\n")
            f.write("\n")
        yaml.dump(out_dict, f, Dumper=CleanDumper, default_flow_style=False,
                  allow_unicode=True, sort_keys=False, width=200)
