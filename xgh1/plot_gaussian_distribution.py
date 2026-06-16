#!/usr/bin/env python3
"""Plot a histogram and Gaussian fit for a list of discrete values.

Examples:
  python script/plot_gaussian_distribution.py values.txt
  python script/plot_gaussian_distribution.py values.csv -o dist.png --bins 10
  python script/plot_gaussian_distribution.py --values 4.84 4.95 5.68 5.42
"""

import argparse
import math
import re
from pathlib import Path

import numpy as np
try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError as exc:
    raise SystemExit(
        "错误: 缺少 matplotlib，无法绘图。\n"
        "请先安装: pip install matplotlib numpy\n"
        "如果在服务器 conda 环境中运行，也可用: conda install matplotlib numpy"
    ) from exc


NUMBER_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?")


def read_values(path):
    """Extract all numeric tokens from a text/csv/log file."""
    text = Path(path).read_text(encoding="utf-8", errors="ignore")
    return [float(m.group(0)) for m in NUMBER_RE.finditer(text)]


def gaussian_pdf(x, mean, sigma):
    if sigma <= 0:
        return np.zeros_like(x)
    coef = 1.0 / (sigma * math.sqrt(2.0 * math.pi))
    return coef * np.exp(-0.5 * ((x - mean) / sigma) ** 2)


def default_bin_count(values):
    """Freedman-Diaconis bin count with a small-data fallback."""
    n = len(values)
    if n < 2:
        return 1
    q25, q75 = np.percentile(values, [25, 75])
    iqr = q75 - q25
    data_range = max(values) - min(values)
    if iqr <= 0 or data_range <= 0:
        return max(3, int(round(math.sqrt(n))))
    bin_width = 2.0 * iqr / (n ** (1.0 / 3.0))
    if bin_width <= 0:
        return max(3, int(round(math.sqrt(n))))
    return max(3, int(math.ceil(data_range / bin_width)))


def plot_distribution(values, output, bins=None, title=None, xlabel="Value",
                      ylabel="Density", dpi=300, show=False):
    values = np.array(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        raise ValueError("没有可绘制的有效数字。")

    mean = float(np.mean(values))
    sigma = float(np.std(values, ddof=1)) if values.size > 1 else 0.0
    median = float(np.median(values))
    bins = int(bins) if bins else default_bin_count(values)

    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.hist(
        values,
        bins=bins,
        density=True,
        color="#4C78A8",
        edgecolor="white",
        linewidth=0.9,
        alpha=0.82,
        label=f"Histogram (n={values.size})",
    )

    if sigma > 0:
        x_pad = 0.12 * (values.max() - values.min() or 1.0)
        x = np.linspace(values.min() - x_pad, values.max() + x_pad, 500)
        y = gaussian_pdf(x, mean, sigma)
        ax.plot(
            x,
            y,
            color="#D95F02",
            linewidth=2.2,
            label=f"Gaussian fit: mu={mean:.3f}, sigma={sigma:.3f}",
        )
        ax.axvline(mean, color="#D95F02", linestyle="--", linewidth=1.2, alpha=0.85)

    ax.axvline(median, color="#2F855A", linestyle=":", linewidth=1.4,
               label=f"Median={median:.3f}")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title or "Distribution with Gaussian Fit")
    ax.grid(axis="y", alpha=0.25, linewidth=0.8)
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=dpi)
    if show:
        plt.show()
    plt.close(fig)

    return {
        "n": int(values.size),
        "mean": mean,
        "median": median,
        "sigma_sample": sigma,
        "min": float(values.min()),
        "max": float(values.max()),
        "bins": bins,
        "output": str(output),
    }


def parse_args():
    parser = argparse.ArgumentParser(
        description="对一组离散数值绘制统计分布柱状图和高斯拟合曲线。"
    )
    parser.add_argument(
        "input",
        nargs="?",
        help="输入文件，支持 txt/csv/log 等文本文件；脚本会提取其中所有数字。",
    )
    parser.add_argument(
        "--values",
        nargs="+",
        type=float,
        help="直接在命令行输入数字；若指定该项，则不需要 input 文件。",
    )
    parser.add_argument("-o", "--output", default="result/gaussian_distribution.png",
                        help="输出图片路径，默认 result/gaussian_distribution.png。")
    parser.add_argument("--bins", type=int, default=None,
                        help="直方图 bin 数；默认用 Freedman-Diaconis 规则自动估计。")
    parser.add_argument("--title", default=None, help="图标题。")
    parser.add_argument("--xlabel", default="Value", help="x 轴标题。")
    parser.add_argument("--ylabel", default="Density", help="y 轴标题。")
    parser.add_argument("--dpi", type=int, default=300, help="输出图片分辨率。")
    parser.add_argument("--show", action="store_true", help="保存后弹出显示窗口。")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.values:
        values = args.values
    elif args.input:
        values = read_values(args.input)
    else:
        raise SystemExit("错误: 请提供 input 文件或使用 --values 输入数字。")

    stats = plot_distribution(
        values,
        output=args.output,
        bins=args.bins,
        title=args.title,
        xlabel=args.xlabel,
        ylabel=args.ylabel,
        dpi=args.dpi,
        show=args.show,
    )

    print("完成: 已输出", stats["output"])
    print(
        "统计: "
        f"n={stats['n']}, mean={stats['mean']:.6f}, "
        f"median={stats['median']:.6f}, sigma={stats['sigma_sample']:.6f}, "
        f"min={stats['min']:.6f}, max={stats['max']:.6f}, bins={stats['bins']}"
    )


if __name__ == "__main__":
    main()
