# PWmatScript

PWmat 第一性原理缺陷计算后处理工具集。

---

## Git 基础操作

### 仓库操作

#### 克隆仓库（Clone）
```bash
# HTTPS 方式（推荐，通用性最好）
git clone https://github.com/your-org/PWmatScript.git

# SSH 方式（需要先配置 SSH Key）
git clone git@github.com:your-org/PWmatScript.git
```

#### 初始化仓库（Init）
```bash
# 在已有目录中初始化 git
git init

# 初始化并指定默认分支名
git init -b main
```

#### 管理远程仓库（Remote）
```bash
# 查看当前远程仓库
git remote -v

# 添加远程仓库
git remote add origin https://github.com/your-org/PWmatScript.git

# 修改远程仓库地址
git remote set-url origin https://github.com/your-org/PWmatScript.git
```

---

### 查看状态

```bash
# 查看工作区状态（修改了哪些文件、暂存了哪些）
git status

# 简洁模式
git status -s
```

#### 查看提交历史（Log）
```bash
# 查看提交日志
git log

# 单行显示，更紧凑
git log --oneline

# 显示最近 N 条
git log --oneline -10

# 图形化显示分支历史
git log --oneline --graph --all

# 查看某个文件的提交历史
git log --oneline -- path/to/file.py
```

#### 查看差异（Diff）
```bash
# 查看工作区未暂存的修改
git diff

# 查看已暂存但未提交的修改
git diff --staged

# 查看某次提交与当前的区别
git diff HEAD~1

# 查看两个分支之间的区别
git diff main..feature-branch
```

---

### 分支操作

```bash
# 查看所有分支（带 * 的是当前分支）
git branch

# 查看所有分支（含远程）
git branch -a

# 创建新分支
git branch feature/new-calculation

# 切换分支
git switch feature/new-calculation

# 创建并切换到新分支（一步到位）
git switch -c feature/new-calculation
# 或旧语法
git checkout -b feature/new-calculation

# 切换到已有分支
git switch main

# 删除本地分支（已合并的分支）
git branch -d feature/old-branch

# 强制删除本地分支（即使未合并）
git branch -D feature/old-branch
```

#### 合并分支（Merge）
```bash
# 将 feature 分支合并到当前分支
git merge feature/new-calculation

# 如果发生冲突，手动解决后：
git add .
git commit -m "resolve merge conflicts"

# 放弃合并，回到合并前状态
git merge --abort
```

---

### 暂存与提交

```bash
# 暂存单个文件
git add README.md

# 暂存某个目录下所有文件
git add pipeline/

# 暂存所有修改过的文件（不包括新文件）
git add -u

# 暂存当前目录所有修改
git add .

# 提交已暂存的文件
git commit -m "docs: 更新 README"

# 暂存所有修改过的文件并直接提交（跳过 git add）
git commit -a -m "fix: 修复 calculate_deltaQ 公式错误"

# 修改最近一次提交（补充遗漏的文件或修改提交信息）
git add forgotten_file.py
git commit --amend

# 只修改最近一次提交信息（不改文件）
git commit --amend -m "docs: 更完善的提交信息"
```

#### 提交信息规范建议
```
feat:     新功能
fix:      修复 bug
docs:     文档修改
refactor: 代码重构（不改变功能）
style:    代码格式调整
test:     测试相关
chore:    构建或辅助工具变动
```

---

### 同步远程

```bash
# 拉取远程更新并自动合并到当前分支
git pull

# 等价于 git pull 的完整写法
git pull origin main

# 仅拉取远程更新，不合并（更安全）
git fetch origin

# 查看远程有什么更新（fetch 后查看）
git log HEAD..origin/main --oneline

# 推送本地提交到远程
git push origin main

# 首次推送新分支并设置上游
git push -u origin feature/new-calculation

# 之后可以直接推送
git push

# 强制推送（慎用！会覆盖远程历史）
git push --force
```

---

### 撤销操作

```bash
# 撤销工作区中某文件的修改（恢复到最近一次提交的状态）
git restore README.md

# 撤销所有工作区修改
git restore .

# 将已暂存的文件移出暂存区（保留工作区修改）
git restore --staged README.md

# 回退到之前的提交（保留工作区修改）
git reset HEAD~1

# 回退并丢弃所有修改（危险！不可恢复）
git reset --hard HEAD~1

# 将当前工作区修改暂存起来，清理工作区
git stash

# 查看 stash 列表
git stash list

# 恢复最近一次 stash 的内容
git stash pop

# 恢复指定 stash
git stash pop stash@{1}
```

---

### .gitignore 示例

在项目根目录创建 `.gitignore` 文件，忽略不需要版本控制的文件：

```gitignore
# Python
__pycache__/
*.py[cod]
*.egg-info/
dist/
build/
.eggs/

# IDE / 编辑器
.vscode/
.idea/
*.swp
*.swo
*~

# 计算结果（通常较大，不需要入库）
result/
calculate/
*.log

# 操作系统
.DS_Store
Thumbs.db
```

---

## 项目简介

本项目为 **PWmat 第一性原理缺陷计算** 提供完整的后处理工具链，涵盖：

- 缺陷电荷修正（ICIC + PA 方法）
- 形成能计算与费米能级绘图
- 热力学缺陷转变能级计算
- 化学势允许范围分析
- 重组能提取与组态坐标图（CCD）绘制

---

## 目录结构

| 目录 | 说明 |
|------|------|
| `hpc_workflow/` | HPC 集群自动化缺陷计算流程（Slurm 作业脚本） |
| `formation_energy/` | 独立形成能绘制脚本（包含使用说明 PDF） |
| `standalone/` | **可独立使用的通用脚本（推荐）**，通过 `-i`/`-o` 显式指定输入输出路径，适配任意目录结构 |
| `pipeline/` | 项目内工作流脚本集，依赖固定项目目录结构 |

> 详细使用说明见 [`standalone/README.md`](standalone/README.md)

---

## 环境要求

- **操作系统**：Linux（集群环境）或 macOS / Windows（本地分析）
- **Python 3.6+**
- **PWmat**：第一性原理计算软件
- **Slurm**（可选）：HPC 作业调度系统

### Python 依赖

```bash
pip install numpy pyyaml matplotlib pandas
```

| 包 | 用途 |
|----|------|
| `numpy` | 数值计算、矩阵运算 |
| `PyYAML` | 读写 YAML 配置文件 |
| `matplotlib` | 绘制形成能图、CCD 图 |
| `pandas` | 数据处理与导出 |

---

## 快速开始

### 1. 准备工作

确保已完成 PWmat 缺陷计算，目录结构如下：

```
project/
├── bulk/scf/          # 完美超胞 SCF 计算
├── q0/scf/            # 中性缺陷 SCF 计算
├── q1/scf/            # +1 带电缺陷 SCF 计算
├── qm1/scf/           # -1 带电缺陷 SCF 计算
├── config.yaml        # 项目配置文件
└── chemical/          # 竞争相总能数据（可选）
```

### 2. 生成缺陷修正输入文件

```bash
# 生成 defect.input（缺陷位置文件）
python standalone/generate_defect_input.py \
    -ibulk bulk/scf \
    -iq0 q0/scf \
    -iq1 q1/scf \
    -o defect_inputs/

# 生成 occ.input（占据数分析）
python standalone/generate_occ.py \
    -i q0/scf \
    -o occ.input
```

### 3. 计算 ICIC + PA 电荷修正

在 HPC 上依次运行 `1_get_rho.sh` → `2_coulomb_integral.sh` → `3_get_results.sh`，或使用独立版脚本：

```bash
python standalone/correction.py collect \
    -ibulk bulk/scf \
    -iq0 q0/scf \
    -iq1 q1/scf \
    -iqm1 qm1/scf \
    -o result/correction_results.yaml
```

### 4. 计算化学势允许范围（可选）

```bash
python standalone/calculate_chemical_potential_bounds.py \
    -i chemical \
    -o result/chemical_potential_bounds.yaml \
    --target SiO2
```

### 5. 计算形成能与转变能级

```bash
python standalone/formation_energy.py \
    -i result/correction_results.yaml \
    -o result/defect_results.yaml
```

### 6. 绘制组态坐标图（CCD）

```bash
# 提取重组能
python standalone/extract_reorganization_energy.py \
    --method structural_relaxation \
    --defect v_O --charge 1 \
    -i q1/relax/RELAXSTEPS q1/S21/RELAXSTEPS \
    -o result/reorganization_energy.log

# 计算 ΔQ
python standalone/calculate_deltaQ_config.py \
    -c0 atom_config_0 \
    -cq atom_config_q \
    -o result/deltaQ.log

# 绘制 CCD
python standalone/plot_ccd.py \
    -i result/reorganization_energy.log \
    --deltaQ 2.4 \
    --defect v_O --charge 1 \
    -o result/ccd/v_O.png
```

---

## 工作流程图

```
PWmat DFT 计算
      │
      ▼
generate_defect_input.py  ──→  defect.input
generate_occ.py           ──→  occ.input / IN.OCC
      │
      ▼
1_get_rho.sh → 2_coulomb_integral.sh → 3_get_results.sh
  （或 pipeline/ 中的 batch_correction.py 批量处理）
      │
      ▼
correction.py collect  ──→  correction_results.yaml
      │
      ▼
calculate_chemical_potential_bounds.py（可选）
      │
      ▼
formation_energy.py  ──→  defect_results.yaml  +  形成能图
      │
      ▼
extract_reorganization_energy.py  →  calculate_deltaQ_config.py  →  plot_ccd.py
```

---

## 核心脚本一览

| 脚本 | 功能 |
|------|------|
| `generate_defect_input.py` | 自动生成缺陷位置文件 `defect.input` |
| `generate_occ.py` | 分析占据数，生成 `occ.input` 和 `IN.OCC` |
| `correction.py` | ICIC + PA 电荷修正（prepare / collect / all） |
| `formation_energy.py` | 形成能计算与费米能级绘图 |
| `calculate_chemical_potential_bounds.py` | 化学势允许范围计算 |
| `extract_reorganization_energy.py` | 重组能提取（两种方法） |
| `calculate_deltaQ_config.py` | 组态坐标 ΔQ 计算 |
| `plot_ccd.py` | 二维组态坐标图绘制 |
| `batch_correction.py` | 批量缺陷修正处理 |
| `export_transition_levels.py` | 转变能级导出为 CSV |

---

## 使用提示

- **推荐优先使用 `standalone/` 中的脚本**，通过 `-i`/`-o` 显式指定输入输出，适配不同项目结构
- 每个脚本都支持 `--help` 查看参数说明
- 详细使用说明见 [`standalone/README.md`](standalone/README.md)
- `pipeline/` 中的脚本适用于已有固定项目目录结构的情况
