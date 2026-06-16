# batch_correction / batch_submit_corrections 使用文档

本文档汇总：

- `script/batch_correction.py`
- `script/batch_submit_corrections.sh`

按三类流程说明：**只 PA**、**只 ICIC**、**both = ICIC + PA**。

## 1. 脚本分工

| 脚本 | 主要用途 |
|---|---|
| `batch_correction.py` | 生成 `defect.input`；收集修正；计算形成能和转变能级。 |
| `batch_submit_corrections.sh` | 批量执行 ICIC 流程：生成差分电荷、提交库仑积分、收集 ICIC 结果。 |

**只做 PA 时一般不用运行 `batch_submit_corrections.sh`。**

## 2. 通用目录要求

项目根目录建议保持如下树结构。**`config.yaml` 必须放在项目根目录，和 `calculate/` 同一级**，否则脚本可能无法正确读取体系参数、化学势、介电常数、带隙和电荷态设置。

```text
project_root/
├── config.yaml
├── calculate/
│   ├── bulk/
│   │   └── scf/
│   │       ├── REPORT
│   │       ├── atom.config 或 final.config
│   │       └── OUT.VATOM
│   └── 1_H_i/
│       ├── q_0/
│       │   └── scf/
│       │       ├── REPORT
│       │       ├── atom.config 或 final.config
│       │       ├── OUT.OCC
│       │       └── OUT.VATOM
│       ├── q_-1/
│       │   └── scf/
│       │       ├── REPORT
│       │       ├── atom.config
│       │       └── OUT.OCC
│       └── q_1/
│           └── scf/
│               ├── REPORT
│               ├── atom.config
│               └── OUT.OCC
└── script/
    ├── batch_correction.py
    ├── batch_submit_corrections.sh
    ├── 1_get_rho.sh
    ├── 2_coulomb_integral.sh
    └── 3_get_results.sh
```

当前仓库根目录下已经提供了一个模板：

```text
config.yaml
```

新项目可先复制这个模板到目标项目根目录，再按体系修改其中的：

```text
system.name
system.dielectric
system.gap
system.elements
system.chemical_potential
structure.defect_center
formation_energy.charge
band_alignment.E_VBM
```

## 3. 通用第一步：生成 defect.input

所有模式都建议先运行：

```bash
python script/batch_correction.py prepare
```

只处理单个缺陷：

```bash
python script/batch_correction.py prepare --defect 10_H_i
```

只处理某个编号：

```bash
python script/batch_correction.py prepare --idx 10
```

该步骤会生成：

```text
calculate/<defect>/q_<charge>/scf/defect.input
calculate/<defect>/q_0/scf/defect_<charge>.input
```

## 4. 只 PA

### 需要什么

PA 需要 bulk 和中性缺陷的 `OUT.VATOM`：

```text
calculate/bulk/scf/OUT.VATOM
calculate/bulk/scf/atom.config 或 final.config
calculate/<defect>/q_0/scf/OUT.VATOM
calculate/<defect>/q_0/scf/atom.config 或 final.config
```

带电态还需要能读到总能：

```text
calculate/<defect>/q_<charge>/scf/REPORT
```

或：

```text
calculate/<defect>/q_<charge>/scf/atom.config
```

### 执行命令

```bash
python script/batch_correction.py prepare
python script/batch_correction.py collect --corr pa
```

只处理单个缺陷：

```bash
python script/batch_correction.py collect --defect 10_H_i --corr pa
```

### 重要说明

`--corr pa` **不需要** ICIC 的 `image-corr_*`、`REPORT`、`REPORT.0`。

如果出现：

```text
✗ PA缺失
```

优先检查：

```bash
ls calculate/bulk/scf/OUT.VATOM
ls calculate/<defect>/q_0/scf/OUT.VATOM
```

## 5. 只 ICIC

### 需要什么

ICIC 需要中性态和带电态的占据文件：

```text
calculate/<defect>/q_0/scf/OUT.OCC
calculate/<defect>/q_<charge>/scf/OUT.OCC
```

运行后需要生成：

```text
calculate/<defect>/q_0/scf/image-corr_<charge>/REPORT
calculate/<defect>/q_0/scf/image-corr_<charge>/REPORT.0
```

例如：

```text
calculate/10_H_i/q_0/scf/image-corr_-1/REPORT
calculate/10_H_i/q_0/scf/image-corr_-1/REPORT.0
```

### 执行命令

1. 生成输入：

```bash
python script/batch_correction.py prepare
```

2. 提交 ICIC 相关计算：

```bash
bash script/batch_submit_corrections.sh submit
```

3. 查看状态：

```bash
bash script/batch_submit_corrections.sh status
```

4. Slurm 任务完成后收集 ICIC：

```bash
bash script/batch_submit_corrections.sh collect
```

5. 只使用 ICIC 计算形成能 / 转变能级：

```bash
python script/batch_correction.py collect --corr icic
```

只处理单个缺陷：

```bash
bash script/batch_submit_corrections.sh submit --defect 10_H_i
bash script/batch_submit_corrections.sh status --defect 10_H_i
bash script/batch_submit_corrections.sh collect --defect 10_H_i
python script/batch_correction.py collect --defect 10_H_i --corr icic
```

## 6. both：ICIC + PA

### 需要什么

both 同时需要：

- PA 的 `OUT.VATOM`
- ICIC 的 `OUT.OCC`
- ICIC 生成的 `image-corr_<charge>/REPORT` 和 `REPORT.0`

### 执行命令

```bash
python script/batch_correction.py prepare
bash script/batch_submit_corrections.sh submit
bash script/batch_submit_corrections.sh status
bash script/batch_submit_corrections.sh collect
python script/batch_correction.py collect --corr both
```

只处理单个缺陷：

```bash
python script/batch_correction.py prepare --defect 10_H_i
bash script/batch_submit_corrections.sh submit --defect 10_H_i
bash script/batch_submit_corrections.sh status --defect 10_H_i
bash script/batch_submit_corrections.sh collect --defect 10_H_i
python script/batch_correction.py collect --defect 10_H_i --corr both
```

## 7. 常用选项

### batch_correction.py

```bash
python script/batch_correction.py prepare
python script/batch_correction.py collect --corr pa
python script/batch_correction.py collect --corr icic
python script/batch_correction.py collect --corr both
python script/batch_correction.py collect --corr none
```

筛选缺陷：

```bash
--idx 10
--defect 10_H_i
--only 10_H_i
```

### batch_submit_corrections.sh

```bash
bash script/batch_submit_corrections.sh submit
bash script/batch_submit_corrections.sh status
bash script/batch_submit_corrections.sh collect
```

强制重新提交或重新收集：

```bash
bash script/batch_submit_corrections.sh submit --defect 10_H_i --force
bash script/batch_submit_corrections.sh collect --defect 10_H_i --force
```

## 8. 输出文件

最终结果由 `batch_correction.py collect` 写出：

```text
result/defect_results.yaml
result/E_forms/E_formation_corrected_<defect>.txt
```

终端会打印每个缺陷、电荷态的：

```text
E_raw
E_corr
修正状态
形成能
转变能级
```

## 9. 最小流程速查

| 目标 | 命令 |
|---|---|
| 只 PA | `python script/batch_correction.py prepare` 后运行 `python script/batch_correction.py collect --corr pa` |
| 只 ICIC | `prepare` → `bash ... submit` → `bash ... status` → `bash ... collect` → `python ... collect --corr icic` |
| both | `prepare` → `bash ... submit` → `bash ... status` → `bash ... collect` → `python ... collect --corr both` |
| 不加修正 | `python script/batch_correction.py collect --corr none` |
