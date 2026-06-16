#!/usr/bin/env python3
"""
自动生成 defect.input 文件的脚本。

用法:
    python script/generate_defect_input.py <defect_dir>

示例:
    python script/generate_defect_input.py calculate/0_v_Si

脚本会:
  1. 通过比较 bulk 和缺陷超胞的 atom.config，自动找到缺陷位置（空位坐标）
  2. 计算 bound 坐标（缺陷坐标 + 0.5，取模到 [0,1)）
  3. 扫描 q_* 目录，找到所有带电态
  4. 在每个非零带电态的 scf/ 目录下生成 defect.input
"""

import os
import sys
import glob
import re
import numpy as np


def get_project_root():
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


def read_atom_config(filepath):
    """读取 atom.config 或 final.config，返回 (原子数, 晶格矢量, [(元素序号, x, y, z), ...])"""
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # 第一行：原子数
    first_line = lines[0].strip()
    num_atoms = int(first_line.split()[0])

    # 读取晶格矢量（跳过 "Lattice vector" 行）
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

    # 读取原子坐标（跳过 "Position" 行）
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
    """
    比较 bulk 和 defect 的原子坐标，找到空位位置。
    对于每个 bulk 原子，检查是否能在 defect 中找到匹配（考虑周期性），
    找不到匹配的就是空位所在位置。
    """
    threshold = 0.1  # 分数坐标匹配阈值（Å）

    G = np.dot(lattice, lattice.T)

    def frac_distance(f1, f2):
        delta = np.array(f2) - np.array(f1)
        delta -= np.round(delta)
        return np.sqrt(np.dot(delta, np.dot(G, delta)))

    defect_coords = [(a[1], a[2], a[3]) for a in defect_atoms]
    vacancy_positions = []

    for batom in bulk_atoms:
        bcoord = (batom[1], batom[2], batom[3])
        matched = False
        for dcoord in defect_coords:
            dist = frac_distance(bcoord, dcoord)
            if dist < threshold:
                matched = True
                break
        if not matched:
            vacancy_positions.append(bcoord)

    return vacancy_positions


def find_defect_position_from_json(project_root):
    """从 input.json 中读取 defect_center"""
    json_path = os.path.join(project_root, 'prepare', 'input.json')
    if not os.path.exists(json_path):
        return None

    import json
    with open(json_path, 'r') as f:
        data = json.load(f)

    centers = data.get('structure_information', {}).get('defect_center', [])
    if centers:
        return centers  # [[x, y, z], ...]
    return None


def get_charge_dirs(defect_dir):
    """扫描 defect_dir 下的 q_* 目录，返回 {charge_int: dir_path}"""
    charge_dirs = {}
    for item in sorted(os.listdir(defect_dir)):
        m = re.match(r'^q_([-+]?\d+)$', item)
        if m:
            q = int(m.group(1))
            charge_dirs[q] = os.path.join(defect_dir, item)
    return charge_dirs


def format_charge_state(q):
    """将整数电荷转为 defect.input 格式，如 +1, -1, +2"""
    if q > 0:
        return f"+{q}"
    else:
        return str(q)


def main():
    if len(sys.argv) < 2:
        print("用法: python script/generate_defect_input.py <defect_dir>")
        print("示例: python script/generate_defect_input.py calculate/0_v_Si")
        sys.exit(1)

    defect_dir = os.path.abspath(sys.argv[1])
    project_root = get_project_root()

    if not os.path.isdir(defect_dir):
        print(f"错误: 目录不存在: {defect_dir}")
        sys.exit(1)

    # ---- 1. 确定 bulk scf 路径 ----
    bulk_scf = os.path.join(project_root, 'calculate', 'bulk', 'scf')
    if not os.path.isdir(bulk_scf):
        print(f"错误: 找不到 bulk scf 目录: {bulk_scf}")
        sys.exit(1)

    # ---- 2. 确定 neutral (q_0) scf 路径 ----
    neutral_scf = os.path.join(defect_dir, 'q_0', 'scf')
    if not os.path.isdir(neutral_scf):
        print(f"错误: 找不到中性态目录: {neutral_scf}")
        sys.exit(1)

    # ---- 3. 找到缺陷位置 ----
    # 优先从 atom.config 比较得到，如果失败则尝试 input.json
    bulk_config = os.path.join(bulk_scf, 'atom.config')
    defect_config = os.path.join(neutral_scf, 'atom.config')

    defect_pos = None

    if os.path.exists(bulk_config) and os.path.exists(defect_config):
        print(f"读取 bulk atom.config: {bulk_config}")
        print(f"读取 defect atom.config: {defect_config}")

        n_bulk, lattice_bulk, atoms_bulk = read_atom_config(bulk_config)
        n_defect, lattice_defect, atoms_defect = read_atom_config(defect_config)

        if n_bulk > n_defect:
            # vacancy: bulk 有更多原子
            vacancies = find_vacancy_position(atoms_bulk, atoms_defect, lattice_bulk)
            if vacancies:
                defect_pos = vacancies[0]
                print(f"检测到空位位置 (分数坐标): {defect_pos[0]:.6f} {defect_pos[1]:.6f} {defect_pos[2]:.6f}")
                if len(vacancies) > 1:
                    print(f"  警告: 检测到 {len(vacancies)} 个空位，使用第一个")
            else:
                print("警告: 未通过原子比较找到空位位置")
        elif n_bulk < n_defect:
            # interstitial: defect 有更多原子，找多出来的原子
            extras = find_vacancy_position(atoms_defect, atoms_bulk, lattice_defect)
            if extras:
                defect_pos = extras[0]
                print(f"检测到间隙原子位置 (分数坐标): {defect_pos[0]:.6f} {defect_pos[1]:.6f} {defect_pos[2]:.6f}")
            else:
                print("警告: 未通过原子比较找到间隙位置")
        else:
            # antisite: 原子数相同，找元素不同的原子
            print("原子数相同，可能是反位缺陷，尝试查找元素变化的原子...")
            G = np.dot(lattice_bulk, lattice_bulk.T)
            for batom in atoms_bulk:
                bcoord = np.array([batom[1], batom[2], batom[3]])
                for datom in atoms_defect:
                    dcoord = np.array([datom[1], datom[2], datom[3]])
                    delta = dcoord - bcoord
                    delta -= np.round(delta)
                    dist = np.sqrt(np.dot(delta, np.dot(G, delta)))
                    if dist < 0.1 and batom[0] != datom[0]:
                        defect_pos = (batom[1], batom[2], batom[3])
                        print(f"检测到反位位置 (分数坐标): {defect_pos[0]:.6f} {defect_pos[1]:.6f} {defect_pos[2]:.6f}")
                        break
                if defect_pos:
                    break

    # 如果原子比较失败，尝试从 input.json 读取
    if defect_pos is None:
        print("尝试从 input.json 读取缺陷位置...")
        centers = find_defect_position_from_json(project_root)
        if centers:
            defect_pos = tuple(centers[0])
            print(f"从 input.json 读取到缺陷位置: {defect_pos[0]:.6f} {defect_pos[1]:.6f} {defect_pos[2]:.6f}")
        else:
            print("错误: 无法确定缺陷位置，请手动指定")
            sys.exit(1)

    # ---- 4. 计算 bound 坐标 ----
    bound_pos = tuple((c + 0.5) % 1.0 for c in defect_pos)
    print(f"bound 坐标: {bound_pos[0]:.6f} {bound_pos[1]:.6f} {bound_pos[2]:.6f}")

    # ---- 5. 扫描带电态目录并生成 defect.input ----
    charge_dirs = get_charge_dirs(defect_dir)
    if not charge_dirs:
        print(f"错误: 在 {defect_dir} 下未找到 q_* 目录")
        sys.exit(1)

    print(f"\n检测到电荷态: {sorted(charge_dirs.keys())}")

    generated = 0
    for q, qdir in sorted(charge_dirs.items()):
        if q == 0:
            continue  # 中性态不需要修正

        scf_dir = os.path.join(qdir, 'scf')
        if not os.path.isdir(scf_dir):
            print(f"  跳过 q={q}: scf 目录不存在 ({scf_dir})")
            continue

        charge_str = format_charge_state(q)
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

        print(f"  已生成: {output_path}")
        print(f"    charge state: {charge_str}")
        generated += 1

    print(f"\n完成! 共生成 {generated} 个 defect.input 文件")


if __name__ == '__main__':
    main()
