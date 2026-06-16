#!/usr/bin/env python3
"""Export defect transition levels from defect_results.yaml to CSV."""

import argparse
import os
import re
import sys

import yaml


def find_project_root():
    start = os.path.abspath(os.path.dirname(__file__))
    current = start
    while True:
        if (
            os.path.isdir(os.path.join(current, "result"))
            or os.path.isdir(os.path.join(current, "calculate"))
            or os.path.exists(os.path.join(current, "config.yaml"))
        ):
            return current
        parent = os.path.dirname(current)
        if parent == current:
            return start
        current = parent


def defect_index(defect_name):
    """Return the leading numeric defect index in names such as 149_H_O."""
    match = re.match(r"^(\d+)(?:_|$)", str(defect_name))
    return match.group(1) if match else ""


def load_yaml(path):
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def get_level(levels, *labels):
    for label in labels:
        if label in levels:
            value = levels[label]
            if isinstance(value, (int, float)):
                return f"{float(value):.8f}"
            return value
    return "nan"


def transition_rows(data, include_empty=False):
    for defect in data.get("defects", []) or []:
        name = defect.get("name", "")
        idx = defect_index(name)
        levels = defect.get("transition_levels") or {}

        if not levels:
            if include_empty:
                yield {
                    "idx": idx,
                    "tl(-/0)": "nan",
                    "tl(0/+)": "nan",
                    "tl(-/+)": "nan",
                }
            continue

        yield {
            "idx": idx,
            "tl(-/0)": get_level(levels, "(-1/0)", "(-/0)"),
            "tl(0/+)": get_level(levels, "(0/+1)", "(0/+)"),
            "tl(-/+)": get_level(levels, "(-1/+1)", "(-/+)"),
        }


def main():
    root = find_project_root()
    parser = argparse.ArgumentParser(
        description="Export transition levels from result/defect_results.yaml to CSV."
    )
    parser.add_argument(
        "-i",
        "--input",
        default=os.path.join(root, "result", "defect_results.yaml"),
        help="Input defect_results.yaml path.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=os.path.join(root, "result", "transition_levels.log"),
        help="Output log/text path.",
    )
    parser.add_argument(
        "--include-empty",
        action="store_true",
        help="Also write defects with no transition levels as empty rows.",
    )
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"ERROR: input YAML not found: {args.input}", file=sys.stderr)
        return 1

    data = load_yaml(args.input) or {}
    rows = list(transition_rows(data, include_empty=args.include_empty))

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    def sort_key(row):
        idx = row.get("idx", "")
        return (0, int(idx)) if str(idx).isdigit() else (1, str(idx))

    rows.sort(key=sort_key)
    fieldnames = ["idx", "tl(-/0)", "tl(0/+)", "tl(-/+)"]
    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write(" ".join(fieldnames) + "\n")
        for row in rows:
            handle.write(" ".join(str(row.get(name, "")) for name in fieldnames) + "\n")

    print(f"Wrote {len(rows)} rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
