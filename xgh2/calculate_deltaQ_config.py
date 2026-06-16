#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
计算deltaQ - atom.config版本
基于正确的计算方法，完全按照dQ_calculator.py的公式
"""

import numpy as np
import sys
import math
import argparse
import json
from pathlib import Path
from typing import List, Tuple

# 修复Windows控制台编码问题
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')


# ================== 原子质量表（完整118元素） ==================
def get_mass_table() -> dict:
    """返回完整的原子质量表（单位：amu）"""
    return {
        1: 1.008, 2: 4.002602, 3: 6.94, 4: 9.012182, 5: 10.81,
        6: 12.017, 7: 14.006, 8: 15.999, 9: 18.9984032, 10: 20.1797,
        11: 22.98976928, 12: 24.305, 13: 26.9815385, 14: 28.085, 15: 30.973762,
        16: 32.06, 17: 35.45, 18: 39.948, 19: 39.0983, 20: 40.078,
        21: 44.955908, 22: 47.867, 23: 50.9415, 24: 51.9961, 25: 54.938044,
        26: 55.845, 27: 58.933194, 28: 58.6934, 29: 63.546, 30: 65.38,
        31: 69.723, 32: 72.63, 33: 74.921595, 34: 78.971, 35: 79.904,
        36: 83.798, 37: 85.4678, 38: 87.62, 39: 88.90584, 40: 91.224,
        41: 92.90638, 42: 95.95, 43: 98.0, 44: 101.07, 45: 102.9055,
        46: 106.42, 47: 107.8682, 48: 112.411, 49: 114.818, 50: 118.71,
        51: 121.76, 52: 127.6, 53: 126.90447, 54: 131.293, 55: 132.90545196,
        56: 137.327, 57: 138.90547, 58: 140.116, 59: 140.90766, 60: 144.242,
        61: 145.0, 62: 150.36, 63: 151.964, 64: 157.25, 65: 158.92535,
        66: 162.5, 67: 164.93032, 68: 167.259, 69: 168.93421, 70: 173.054,
        71: 174.9668, 72: 178.49, 73: 180.94788, 74: 183.84, 75: 186.207,
        76: 190.23, 77: 192.217, 78: 195.084, 79: 196.966569, 80: 200.592,
        81: 204.38, 82: 207.2, 83: 208.9804, 84: 208.982, 85: 208.982,
        86: 222.0, 87: 223.0, 88: 226.0, 89: 227.0, 90: 231.03806,
        91: 238.02891, 92: 238.02891, 93: 242.0, 94: 247.0, 95: 247.0,
        96: 247.0, 97: 247.0, 98: 251.0, 99: 252.0, 100: 257.0,
        101: 258.0, 102: 259.0, 103: 262.0, 104: 265.0, 105: 268.0,
        106: 271.0, 107: 270.0, 108: 277.0, 109: 276.0, 110: 285.0,
        111: 280.0, 112: 285.0, 113: 284.0, 114: 289.0, 115: 288.0,
        116: 293.0, 117: 294.0, 118: 294.0
    }


# ================== 读取atom.config文件 ==================
def read_atom_config(filename: str) -> Tuple[int, List[int], np.ndarray, np.ndarray, np.ndarray]:
    """
    读取atom.config格式的结构文件

    参数:
        filename: 文件名

    返回:
        natoms: 原子数
        atomic_numbers: 原子序数列表
        positions: 分数坐标 (n_atoms, 3)
        lattice: 晶格向量 (3, 3)
        fix_tags: 固定标签 (n_atoms, 3)
    """
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    # 读取原子数
    natoms = int(lines[0].split()[0])

    # 读取晶格向量
    lattice = np.array([
        [float(x) for x in lines[2].split()],
        [float(x) for x in lines[3].split()],
        [float(x) for x in lines[4].split()]
    ])

    # 读取原子位置
    atomic_numbers, positions, fix_tags = [], [], []
    for line in lines[6:6 + natoms]:
        parts = line.split()
        atomic_numbers.append(int(parts[0]))
        coords = [float(x) for x in parts[1:4]]
        # 处理周期性边界条件（原始dQ_calculator.py中的处理）
        # coords = [x - 1.0 if x > 0.99 else x for x in coords]  # 注释掉，使用我们的方法
        positions.append(coords)
        fix = [int(x) for x in parts[4:7]]
        fix_tags.append(fix)

    return natoms, atomic_numbers, np.array(positions), lattice, np.array(fix_tags)


# ================== 写入atom.config文件 ==================
def write_atom_config(filename: str, natoms: int, atomic_numbers: List[int],
                      positions: np.ndarray, lattice: np.ndarray, fix_tags: np.ndarray):
    """
    写入atom.config格式的结构文件

    参数:
        filename: 文件名
        natoms: 原子数
        atomic_numbers: 原子序数列表
        positions: 分数坐标 (n_atoms, 3)
        lattice: 晶格向量 (3, 3)
        fix_tags: 固定标签 (n_atoms, 3)
    """
    with open(filename, 'w') as f:
        f.write(f"{natoms}\n")
        f.write("LATTICE\n")
        for vec in lattice:
            f.write(f"{vec[0]:15.9f} {vec[1]:15.9f} {vec[2]:15.9f}\n")
        f.write("POSITION\n")
        for z, pos, fix in zip(atomic_numbers, positions, fix_tags):
            f.write(f"{z:4d} {pos[0]:15.9f} {pos[1]:15.9f} {pos[2]:15.9f} "
                  f"{fix[0]:d} {fix[1]:d} {fix[2]:d}\n")


# ================== deltaQ计算核心（完全按照dQ_calculator.py） ==================
def calculate_dQ(ref_numbers: List[int], ref_pos_frac: np.ndarray, ref_lattice: np.ndarray,
                 target_numbers: List[int], target_pos_frac: np.ndarray, target_lattice: np.ndarray,
                 method: str = 'normalized') -> dict:
    """
    计算两个结构之间的deltaQ（完全按照dQ_calculator.py的公式）

    参数:
        ref_numbers: 参考结构的原子序数列表
        ref_pos_frac: 参考结构的分数坐标 (n_atoms, 3)
        ref_lattice: 参考结构的晶格向量 (3, 3)
        target_numbers: 目标结构的原子序数列表
        target_pos_frac: 目标结构的分数坐标 (n_atoms, 3)
        target_lattice: 目标结构的晶格向量 (3, 3)
        method: 计算方法 ('normalized' 或 'simple')

    返回:
        results: 包含详细结果的字典
    """
    # 验证原子一致性
    if ref_numbers != target_numbers:
        raise ValueError(f"原子序数不匹配: {len(ref_numbers)} vs {len(target_numbers)}")

    natoms = len(ref_numbers)

    # 获取原子质量表
    mass_table = get_mass_table()

    # ================== 考虑周期性边界条件 ==================
    # 计算分数坐标位移（应用最小镜像约定）
    delta_frac = target_pos_frac - ref_pos_frac
    delta_frac = delta_frac - np.round(delta_frac)  # 最小镜像约定

    # 转换为笛卡尔坐标位移
    displacements_cart = np.dot(delta_frac, ref_lattice.T)  # (n_atoms, 3)

    # 计算参考和目标的笛卡尔坐标
    ref_pos_cart = np.dot(ref_pos_frac, ref_lattice.T)
    target_pos_cart_adjusted = ref_pos_cart + displacements_cart

    # ================== 计算deltaQ ==================
    if method == 'simple':
        # 方法1：简单方法（不做归一化）
        dQ0 = 0.0
        individual_displacements = []

        for z, r_ref, r_target in zip(ref_numbers, ref_pos_cart, target_pos_cart_adjusted):
            displacement = r_ref - r_target
            di_sq = np.sum(displacement ** 2)
            individual_displacements.append(np.sqrt(di_sq))
            mass = mass_table[z]
            dQ0 += mass * di_sq

        dQ = np.sqrt(dQ0)  # 单位: √amu * Å
        dQ_au = None

    elif method == 'normalized':
        # 方法2：归一化方法（完全按照dQ_calculator.py）
        # 公式: dQ = [Σ √(mi * me^(-1) * me/mp) * di^2] / √(Σ di^2)

        dQ0 = 0.0
        norm_sum = 0.0
        individual_displacements = []

        # 原子单位转换：Å → a.u.
        au_conversion = 1.0 / 0.5291772083

        for z, r_ref, r_target in zip(ref_numbers, ref_pos_cart, target_pos_cart_adjusted):
            # 计算位移并转换为原子单位
            displacement = (r_ref - r_target) * au_conversion  # di in a.u.
            di_sq = np.sum(displacement ** 2)  # di^2 in a.u.^2

            # 记录实际位移（Å）
            individual_displacements.append(np.sqrt(np.sum((r_ref - r_target) ** 2)))

            # 归一化因子
            norm_sum += di_sq

            # 质量加权（完全按照dQ_calculator.py，使用精确值）
            # mass * (1 amu / me) = mass * (1.660539e-27 / 9.109383e-31)
            mass = mass_table[z]
            mass_weighted = np.sqrt(mass * 1.660539e-27 / 9.109383e-31)
            dQ0 += mass_weighted * di_sq

        # 归一化
        dQ_au = dQ0 / np.sqrt(norm_sum)  # in a.u.

        # 转换为常用单位
        # dQ (√amu * Å) = dQ_au * √(1/1822.888...) * 0.5291772083
        dQ = dQ_au * math.sqrt(1.0 / 1822.888486209) * 0.5291772083

    else:
        raise ValueError(f"未知的method: {method}")

    # ================== 准备结果 ==================
    results = {
        'n_atoms': natoms,
        'method': method,
        'dQ': dQ,
        'dQ_au': dQ_au if method == 'normalized' else None,
        'individual_displacements': individual_displacements,
        'mean_displacement': float(np.mean(individual_displacements)),
        'max_displacement': float(np.max(individual_displacements)),
        'min_displacement': float(np.min(individual_displacements)),
        'rmsd': float(np.sqrt(np.mean(np.array(individual_displacements) ** 2)))
    }

    return results


# ================== 主程序 ==================
def write_results(output_path: str, results: dict) -> None:
    """Write results as JSON for .json paths, otherwise as a compact text log."""
    path = Path(output_path).expanduser()
    if path.parent != Path("."):
        path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix.lower() == ".json":
        with path.open("w", encoding="utf-8") as f:
            json.dump(results, f, indent=2)
        return

    with path.open("w", encoding="utf-8") as f:
        f.write("deltaQ atom.config calculation\n")
        f.write(f"n_atoms {results['n_atoms']}\n")
        f.write(f"method {results['method']}\n")
        f.write(f"deltaQ_sqrt_amu_A {results['dQ']:.10g}\n")
        if results.get("dQ_au") is not None:
            f.write(f"deltaQ_au {results['dQ_au']:.10g}\n")
        f.write(f"mean_displacement_A {results['mean_displacement']:.10g}\n")
        f.write(f"max_displacement_A {results['max_displacement']:.10g}\n")
        f.write(f"min_displacement_A {results['min_displacement']:.10g}\n")
        f.write(f"rmsd_A {results['rmsd']:.10g}\n")


def main():
    parser = argparse.ArgumentParser(
        description='计算deltaQ - atom.config版本（基于正确的计算方法）',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:

  # 归一化方法（推荐）
  python calculate_deltaQ_config.py 0atom.config 1atom.config

  # 简单方法
  python calculate_deltaQ_config.py 0atom.config 1atom.config --method simple

  # 默认保存到当前目录 dQ.log
  python calculate_deltaQ_config.py 0atom.config 1atom.config

  # 保存JSON结果
  python calculate_deltaQ_config.py 0atom.config 1atom.config --output results.json

  # 显示帮助
  python calculate_deltaQ_config.py --help

注意:
  - 自动考虑周期性边界条件（最小镜像约定）
  - normalized方法: 归一化deltaQ（推荐，与dQ_calculator.py一致）
  - simple方法: 不做归一化
  - 使用精确的质量转换因子（完全按照dQ_calculator.py）
        """
    )

    parser.add_argument('files', nargs='*', help='兼容旧用法: 参考结构文件 目标结构文件')
    parser.add_argument('-i', '--input', nargs=2, metavar=('REF_CONFIG', 'TARGET_CONFIG'),
                       help='输入结构文件: 参考 atom.config 和目标 atom.config')
    parser.add_argument('--method', type=str, default='normalized',
                       choices=['normalized', 'simple'],
                       help='计算方法 (默认: normalized)')
    parser.add_argument('-o', '--output', type=str, default='dQ.log',
                       help='输出文件，默认当前目录 dQ.log；后缀为 .json 时输出JSON')
    parser.add_argument('--verbose', action='store_true',
                       help='显示详细信息')

    args = parser.parse_args()

    if args.input:
        args.file1, args.file2 = args.input
    elif len(args.files) == 2:
        args.file1, args.file2 = args.files
    else:
        parser.error("需要通过 -i REF_CONFIG TARGET_CONFIG 或旧位置参数指定两个输入结构文件")

    print("=" * 70)
    print("deltaQ 计算 - atom.config版本")
    print("基于正确的计算方法（完全按照dQ_calculator.py）")
    print("=" * 70)
    print(f"\n输入文件:")
    print(f"  参考结构: {args.file1}")
    print(f"  目标结构: {args.file2}")
    print(f"  计算方法: {args.method}")

    try:
        # ================== 读取结构文件 ==================
        natoms1, ref_numbers, ref_pos_frac, ref_lattice, ref_fix = read_atom_config(args.file1)
        natoms2, target_numbers, target_pos_frac, target_lattice, target_fix = read_atom_config(args.file2)

        if args.verbose:
            print(f"\n读取结构:")
            print(f"  原子数: {natoms1}")
            print(f"  元素组成: {set(ref_numbers)}")

            print(f"\n晶格参数:")
            print(f"  参考: a={ref_lattice[0,0]:.3f}, b={ref_lattice[1,1]:.3f}, c={ref_lattice[2,2]:.3f} Å")
            print(f"  目标: a={target_lattice[0,0]:.3f}, b={target_lattice[1,1]:.3f}, c={target_lattice[2,2]:.3f} Å")

        # ================== 计算deltaQ ==================
        results = calculate_dQ(
            ref_numbers, ref_pos_frac, ref_lattice,
            target_numbers, target_pos_frac, target_lattice,
            method=args.method
        )

        # ================== 显示结果 ==================
        print("\n" + "=" * 70)
        print("计算结果")
        print("=" * 70)

        if args.method == 'normalized':
            print(f"\ndeltaQ (归一化):")
            print(f"  {results['dQ_au']:15.6f} a.u. (原子单位)")
            print(f"  {results['dQ']:15.6f} √amu * Å")
        elif args.method == 'simple':
            print(f"\ndeltaQ (简单方法):")
            print(f"  {results['dQ']:15.6f} √amu * Å")

        print(f"\n原子位移统计:")
        print(f"  平均位移: {results['mean_displacement']:15.6f} Å")
        print(f"  最大位移: {results['max_displacement']:15.6f} Å")
        print(f"  最小位移: {results['min_displacement']:15.6f} Å")
        print(f"  RMSD: {results['rmsd']:15.6f} Å")

        # 显示位移最大的几个原子
        if args.verbose:
            disp_array = np.array(results['individual_displacements'])
            max_indices = np.argsort(disp_array)[-5:][::-1]
            print(f"\n位移最大的5个原子:")
            for i, idx in enumerate(max_indices, 1):
                z = ref_numbers[idx]
                symbol = {1: 'H', 8: 'O', 72: 'Hf'}.get(z, f'Z{z}')
                print(f"  {i}. 原子 {idx:3d} ({symbol}): {disp_array[idx]:.6f} Å")

        # 结果解释
        print(f"\n结果解释:")
        if args.method == 'normalized':
            if abs(results['dQ']) < 1.0:
                print(f"  deltaQ较小，结构变化微弱")
            elif abs(results['dQ']) < 3.0:
                print(f"  deltaQ中等，结构变化中等")
            else:
                print(f"  deltaQ较大，结构变化显著")

            print(f"\n物理意义:")
            print(f"  归一化deltaQ考虑了:")
            print(f"  - 原子质量加权（√(mi * amu/me)）")
            print(f"  - 周期性边界条件（最小镜像约定）")
            print(f"  - 位移归一化")
            print(f"  - 使用精确的物理常数")

        # 保存结果
        if args.output:
            write_results(args.output, results)
            print(f"\n结果已保存到: {args.output}")

        print("\n" + "=" * 70)
        print("计算完成！")
        print("=" * 70)

    except FileNotFoundError as e:
        print(f"\n错误: 文件未找到 - {e}")
        print("请检查文件名是否正确")
        sys.exit(1)
    except Exception as e:
        print(f"\n错误: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
