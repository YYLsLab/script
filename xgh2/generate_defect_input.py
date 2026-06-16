#!/usr/bin/env python3
"""
自动生成 defect.input 文件的脚本。

用法:
    python standalone/generate_defect_input.py -ibulk <bulk_dir> -iq0 <q0_dir> -iq1 <q1_dir>

示例:
    python standalone/generate_defect_input.py -ibulk bulk/scf -iq0 neutral/scf -iq1 positive/scf

脚本会:
  1. 通过比较 bulk 和缺陷超胞的 atom.config，自动找到缺陷位置（空位坐标）
  2. 计算 bound 坐标（缺陷坐标 + 0.5，取模到 [0,1)）
  3. 根据显式给定的带电态静态计算目录生成 defect.input
"""

import os
import sys
import glob
import re
import argparse
import numpy as np


def get_project_root(start_file=None):
    start_path = os.path.abspath(start_file or __file__)
    start = os.path.dirname(start_path) if os.path.isfile(start_path) else start_path
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


def find_defect_position_from_json(project_root, json_path=None):
    """从 input.json 中读取 defect_center"""
    json_path = json_path or os.path.join(project_root, 'prepare', 'input.json')
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


def parse_charge_dir(values):
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
        charge_dirs[q] = os.path.abspath(path)
    return charge_dirs


def main():
    parser = argparse.ArgumentParser(
        description="根据显式给定的 bulk、中性态和带电态静态计算目录生成 defect.input。",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""示例:
  python standalone/generate_defect_input.py -ibulk /path/to/bulk/scf -iq0 /path/to/q0/scf -iq1 /path/to/q1/scf -o /tmp/defect_inputs
  python standalone/generate_defect_input.py -ibulk /path/to/bulk/scf -iq0 /path/to/q0/scf -iqm1 /path/to/qm1/scf -o /path/to/defect.input
  python standalone/generate_defect_input.py -ibulk /path/to/bulk/scf -iq0 /path/to/q0/scf --charge-dir 2=/path/to/q2/scf

说明:
  - 不假设任何 q_*/scf 目录架构。
  - 多个非零电荷态时，-o 应为目录。
  - 只有一个非零电荷态且 -o 带文件名时，直接写入该文件。
""",
    )
    parser.add_argument("-ibulk", "--input-bulk", required=True, help="bulk 静态计算目录，目录内应有 atom.config。")
    parser.add_argument("-iq0", "--input-q0", required=True, help="中性态 q=0 静态计算目录，目录内应有 atom.config。")
    parser.add_argument("-iq1", "--input-q1", default=None, help="q=+1 静态计算目录。")
    parser.add_argument("-iqm1", "--input-qm1", default=None, help="q=-1 静态计算目录。")
    parser.add_argument(
        "--charge-dir",
        action="append",
        default=[],
        help="额外带电态目录，格式 q=DIR，例如 --charge-dir 2=/path/to/q2/scf 或 --charge-dir -2=/path/to/qm2/scf。",
    )
    parser.add_argument("-o", "--output", default=None, help="输出目录或单个 defect.input 文件。默认写回各带电态输入目录。")
    parser.add_argument("--input-json", default=None, help="显式指定包含 defect_center 的 input.json。")
    parser.add_argument("--defect-center", nargs=3, type=float, metavar=("X", "Y", "Z"), help="手动指定缺陷分数坐标。")
    args = parser.parse_args()

    bulk_scf = os.path.abspath(args.input_bulk)
    neutral_scf = os.path.abspath(args.input_q0)
    if not os.path.isdir(bulk_scf):
        print(f"错误: 找不到 bulk scf 目录: {bulk_scf}")
        sys.exit(1)
    if not os.path.isdir(neutral_scf):
        print(f"错误: 找不到中性态目录: {neutral_scf}")
        sys.exit(1)

    charge_dirs = parse_charge_dir(args.charge_dir)
    if args.input_q1:
        charge_dirs[1] = os.path.abspath(args.input_q1)
    if args.input_qm1:
        charge_dirs[-1] = os.path.abspath(args.input_qm1)
    charge_dirs = {q: path for q, path in charge_dirs.items() if q != 0}
    if not charge_dirs:
        print("错误: 至少需要通过 -iq1、-iqm1 或 --charge-dir 指定一个非零电荷态目录")
        sys.exit(1)
    missing_dirs = [path for path in charge_dirs.values() if not os.path.isdir(path)]
    if missing_dirs:
        print("错误: 以下带电态目录不存在:")
        for path in missing_dirs:
            print(f"  {path}")
        sys.exit(1)

    bulk_config = os.path.join(bulk_scf, 'atom.config')
    defect_config = os.path.join(neutral_scf, 'atom.config')

    defect_pos = None
    if args.defect_center:
        defect_pos = tuple(args.defect_center)

    if defect_pos is None and os.path.exists(bulk_config) and os.path.exists(defect_config):
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
        centers = find_defect_position_from_json(
            os.getcwd(),
            json_path=os.path.abspath(args.input_json) if args.input_json else None,
        )
        if centers:
            defect_pos = tuple(centers[0])
            print(f"从 input.json 读取到缺陷位置: {defect_pos[0]:.6f} {defect_pos[1]:.6f} {defect_pos[2]:.6f}")
        else:
            print("错误: 无法确定缺陷位置，请用 --defect-center X Y Z 手动指定")
            sys.exit(1)

    # ---- 4. 计算 bound 坐标 ----
    bound_pos = tuple((c + 0.5) % 1.0 for c in defect_pos)
    print(f"bound 坐标: {bound_pos[0]:.6f} {bound_pos[1]:.6f} {bound_pos[2]:.6f}")

    # ---- 5. 生成 defect.input ----
    print(f"\n输入电荷态: {sorted(charge_dirs.keys())}")
    charged_items = sorted(charge_dirs.items())
    output_target = os.path.abspath(args.output) if args.output else None
    output_is_file = bool(output_target and os.path.splitext(output_target)[1])
    if output_is_file and len(charged_items) != 1:
        print("错误: 检测到多个非零电荷态，-o 不能是单个文件；请指定输出目录")
        sys.exit(1)

    generated = 0
    for q, scf_dir in charged_items:
        charge_str = format_charge_state(q)
        if output_target:
            if output_is_file:
                output_path = output_target
            else:
                q_label = f"q_{q}" if q < 0 else f"q_+{q}"
                output_path = os.path.join(output_target, q_label, 'defect.input')
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

        print(f"  已生成: {output_path}")
        print(f"    charge state: {charge_str}")
        generated += 1

    print(f"\n完成! 共生成 {generated} 个 defect.input 文件")


if __name__ == '__main__':
    main()
