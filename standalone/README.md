# standalone 使用说明

`standalone/` 存放的是可迁移版本的后处理和辅助脚本。核心目标是通过 `-i/--input` 和 `-o/--output` 显式指定输入、输出，从而适配不同项目目录或服务器目录结构。

## 基本原则

- **按任务选择工作流**：不同目标需要不同前置产物。
- **显式路径优先**：优先使用 `-i` 指定输入文件或目录，使用 `-o` 指定输出文件、输出目录或输出基名。
- **输出含义按脚本区分**：有些脚本输出单个文件，有些脚本输出多个文件；因此 `-o` 在不同脚本中可能表示文件、目录或图像基名。
- **可继续使用 `--help`**：每个脚本都支持 `--help` 查看当前参数说明。

建议使用带有 `numpy`、`PyYAML`、`matplotlib` 的 Python 环境。

- `numpy`：结构、占据数和数值计算需要。
- `PyYAML`：`correction.py`、`formation_energy.py`、`plot_ccd.py`、`calculate_chemical_potential_bounds.py` 需要。
- `matplotlib`：`plot_ccd.py` 绘图需要。

## 任务一：形成能和缺陷转变能级

目标：

- 得到修正后的形成能。
- 得到热力学缺陷转变能级，例如 `(-1/0)`、`(0/+1)`。
- 输出含修正项、形成能和转变能级的 YAML；需要另存标准形成能结果时，可再运行 `formation_energy.py` 输出 `defect_results.yaml` 和 `E_forms/E_formation_corrected_<defect>.txt`。

推荐工作流：

```text
1. generate_defect_input.py / generate_occ.py
   -> 生成缺陷修正所需输入文件

2. correction.py prepare / collect
   -> 生成或收集 ICIC/PA 修正
   -> 在化学势、E_bulk、VBM 可用时，同时计算 中性形成能（E_f0） 和转变能级
   -> correction_results.yaml

3. calculate_chemical_potential_bounds.py
   -> 计算化学势允许范围
   -> 用结果决定 chemical_potential
   -> 在 collect 或 formation_energy 前写入 config.yaml，或写入 correction_results.yaml 的 system.chemical_potential

4. formation_energy.py
   -> 可选；从 correction_results.yaml 重新生成/另存 formation-energy 输出
   -> defect_results.yaml
   -> E_forms/E_formation_corrected_<defect>.txt
```

注意：热力学转变能级不是从输出文件中“读取”的独立量，而是由各电荷态形成能交点计算得到。`correction.py collect` 现在会在收集修正项后直接计算 `E_f0` 和 `transition_levels`；`formation_energy.py` 是独立后处理，用于从已有 `correction_results.yaml` 重新生成或另存形成能结果。

具体操作：

### 1. 生成缺陷修正所需的缺陷中心坐标

使用 `generate_defect_input.py` 时，显式指定 bulk、中性态和带电态静态计算目录：

```powershell
python standalone\generate_defect_input.py `
  -ibulk \path\to\bulk_static `
  -iq0 \path\to\neutral_static `
  -iq1 \path\to\positive_static `
  -iqm1 \path\to\negative_static `
  -o \path\to\out\defect_inputs
```

输入目录要求：

- `-ibulk` 指定 bulk 静态计算目录，目录内需要 `atom.config`。
- `-iq0` 指定中性态 q=0 静态计算目录，目录内需要 `atom.config`。
- `-iq1` 指定 q=+1 静态计算目录。
- `-iqm1` 指定 q=-1 静态计算目录。
- 其他电荷态用 `--charge-dir q=DIR` 指定，例如 `--charge-dir 2=\path\to\q2_static`。

如果缺陷中心不能通过 bulk 和中性态 `atom.config` 自动识别，可手动指定：

```powershell
python standalone\generate_defect_input.py `
  -ibulk \path\to\bulk_static `
  -iq0 \path\to\neutral_static `
  -iq1 \path\to\positive_static `
  --defect-center 0.25 0.50 0.75 `
  -o \path\to\out\defect_inputs
```

如果存在多个非零电荷态，`-o` 必须是目录；如果只有一个非零电荷态，`-o` 可以是单个 `defect.input` 文件。不指定 `-o` 时，脚本会把 `defect.input` 写回各带电态输入目录。

如果需要为 ICIC 流程生成占据输入，使用 `generate_occ.py`。输入目录必须包含：

```text
OUT.KPT
OUT.OCC0
OUT.OCC1
```

示例：

```powershell
python standalone\generate_occ.py `
  -i \path\to\q_1\scf `
  -o \path\to\out\occ
```

无自旋时通常输出：

```text
occ.input
IN.OCC
```

有自旋时可能输出：

```text
occ_spin1.input
occ_spin2.input
IN.OCC_1
IN.OCC_2
```

### 2. 收集修正结果并计算转变能级

推荐使用显式目录模式，不假设项目目录架构：

```powershell
python standalone\correction.py collect `
  -ibulk \path\to\bulk_static `
  -iq0 \path\to\neutral_static `
  -iq1 \path\to\positive_static `
  -iqm1 \path\to\negative_static `
  --defect-name 0_v_Si `
  --config \path\to\config.yaml `
  -o \path\to\out\correction_results.yaml
```

参数含义：

- `-ibulk`：bulk 静态计算目录，用于读取 bulk 总能和 PA 参考。
- `-iq0`：中性态 q=0 静态计算目录。
- `-iq1`：q=+1 静态计算目录。
- `-iqm1`：q=-1 静态计算目录。
- `--charge-dir q=DIR`：额外电荷态目录，例如 `--charge-dir 2=\path\to\q2_static`。
- `--config`：包含 `system.chemical_potential`、`system.gap`、`system.dielectric` 等信息的配置文件。

旧的项目扫描模式仍保留兼容：

```powershell
python standalone\correction.py collect `
  -i \path\to\project `
  -o \path\to\out\correction_results.yaml
```

但该模式会扫描 `<project>/calculate` 下的缺陷目录，不属于完全普适的显式输入模式。

输出：

```text
correction_results.yaml
```

`correction_results.yaml` 会包含：

- `system.E_bulk`
- `system.VBM` / `system.E_VBM`
- `system.gap`
- `system.chemical_potential`
- 每个电荷态的 `E_raw`、`E_corrected`、`corrections`
- 每个电荷态的 `E_f0`
- 每个缺陷的 `transition_levels`
- 若带隙可用，还会包含 `formation_energy_table`

如果 `system.chemical_potential` 为空或缺少某个元素，相关缺陷的形成能化学势项会缺失，脚本会打印警告。此时 `transition_levels` 可能仍被计算，但物理含义取决于缺失项是否会在比较的电荷态之间抵消；不能默认视为可靠。

### 3. 计算化学势范围

默认目录假设：

```text
<project-root>/chemical/<formula>/scf/REPORT
<project-root>/chemical/<formula>/scf/atom.config
<project-root>/calculate/bulk/scf/REPORT
<project-root>/calculate/bulk/scf/atom.config
```

从项目根目录运行：

```powershell
python standalone\calculate_chemical_potential_bounds.py `
  -i \path\to\project `
  -o \path\to\out\chemical_potential_bounds.yaml `
  --target SiO2 `
  --ref Si=Si `
  --ref O=O2
```

也可以直接指定 `chemical` 目录：

```powershell
python standalone\calculate_chemical_potential_bounds.py `
  -i \path\to\project\chemical `
  -o \path\to\out\chemical_potential_bounds.yaml `
  --target SiO2 `
  --ref Si=Si `
  --ref O=O2
```

如果 bulk 主相文件不在默认位置：

```powershell
python standalone\calculate_chemical_potential_bounds.py `
  -i \path\to\project\chemical `
  -o \path\to\out\chemical_potential_bounds.yaml `
  --target SiO2 `
  --target-report \path\to\bulk\scf\REPORT `
  --target-config \path\to\bulk\scf\atom.config `
  --ref Si=Si `
  --ref O=O2
```

输出：

- 主 YAML：`-o` 指定的文件。
- phase energy log：默认由主 YAML 文件名派生。
- constraints log：默认由主 YAML 文件名派生。
- bounds log：默认由主 YAML 文件名派生。

例如 `-o chemical_potential_bounds.yaml` 时，相关日志会写为：

```text
chemical_potential_bounds_phase_energies.log
chemical_potential_bounds_constraints.log
chemical_potential_bounds_bounds.log
```

`calculate_chemical_potential_bounds.py` 只计算化学势范围，不会自动改写 `config.yaml`。运行 `correction.py collect` 或 `formation_energy.py` 前，需要确认化学势已经写入输入数据。

当前脚本读取化学势的规则如下：

| 脚本 | 化学势来源 | 需要的文件 |
|---|---|---|
| `correction.py collect` | `config.yaml` 中的 `system.chemical_potential`，经 `get_bulk_info()` 写入输出 YAML 的 `system.chemical_potential` | 显式模式推荐用 `--config` 指定 YAML；兼容项目扫描模式可从输入项目根目录读取 `config.yaml` |
| `formation_energy.py` | 优先读取输入 `correction_results.yaml` 的 `system.chemical_potential`；如果为空，再从 `config.yaml` 的 `system.chemical_potential` 回退读取 | `-i correction_results.yaml`；若其中无化学势，还需要 `--project-root` 指向含 `config.yaml` 的项目根目录，或用 `--config` 指定配置 |

`config.yaml` 中推荐写法：

```yaml
system:
  chemical_potential:
    Si: -107.17
    O: -9.86
```

如果 `calculate_chemical_potential_bounds.py` 输出了 `chemical_potential_bounds.yaml`，它只是给出允许范围和参考信息；你仍需要人工选择一组具体元素化学势，并写入 `config.yaml` 或 `correction_results.yaml`。

### 4. 可选：重新生成形成能和热力学转变能级文件

```powershell
python standalone\formation_energy.py `
  -i \path\to\out\correction_results.yaml `
  -o \path\to\out\defect_results.yaml `
  --project-root \path\to\project
```

输入：

- `correction_results.yaml`
- 化学势：优先来自 `correction_results.yaml` 的 `system.chemical_potential`
- 如果输入 YAML 中缺少化学势，则需要 `--project-root` 指向包含 `config.yaml` 的项目根目录，或用 `--config` 显式指定配置文件
- 项目根目录中的 `result/formation.out` 或 `calculate/bulk/scf/REPORT`、`atom.config`，用于在输入 YAML 缺少 bulk 信息时补充 `E_bulk`、`VBM`、`gap`

输出：

- `defect_results.yaml`
- `E_forms/E_formation_corrected_<defect>.txt`

`defect_results.yaml` 包含：

- 各电荷态 `E_f0`
- 修正后的能量
- 形成能表
- 热力学转变能级 `transition_levels`

如果 `correction_results.yaml` 已经包含 `system.chemical_potential`、`system.E_bulk`、`system.VBM` 和 `system.gap`，`formation_energy.py` 基本可以只靠 `-i` 和 `-o` 运行。若这些字段不全，建议显式指定 `--project-root` 或 `--config`，避免从错误目录推断 bulk 信息或化学势。

## 任务二：绘制 CCD

目标：

- 得到两态 CCD 图像。
- 得到同名 CCD log，记录 `DeltaQ`、重组能、`DeltaG` 或 `deltaE`、交点和势垒。

推荐工作流：

```text
1. 提取重组能
   -> reorganization_energy.log

2. calculate_deltaQ_config.py
   -> dQ.log 或 dQ.json

3. 手动确定 DeltaG
   -> plot_ccd.py 用 --deltaG 手动输入
   -> 若只画两态简并 CCD，使用默认 --deltaG 0

4. plot_ccd.py
   -> CCD 图像
   -> CCD 同名 log
```

### 1. 提取重组能

使用 `extract_reorganization_energy.py` 从显式指定的输入文件提取单个缺陷、单个电荷态的重组能。该脚本不扫描项目目录，不假设 `calculate/<defect>/q_*/...` 结构。

支持两种提取方法：

| 方法 | 输入 | 公式 |
|---|---|---|
| `structural_relaxation` | 两个 `RELAXSTEPS` | `S(0->q)=E_first(S1)-E_END(S1)`；`S(q->0)=E_first(S2)-E_END(S2)` |
| `fixed_charge_static` | 四个 `REPORT` | `S(0->q)=E0_prime-Eq`；`S(q->0)=Eq_prime-E0` |

`structural_relaxation` 示例：

```powershell
python standalone\extract_reorganization_energy.py `
  --method structural_relaxation `
  --defect 0_v_Si `
  --charge 1 `
  -i \path\to\S1\RELAXSTEPS \path\to\S2\RELAXSTEPS `
  -o \path\to\out\reorganization_energy.log
```

其中 `-i` 后两个文件的顺序必须是：

```text
S(0->q) 的 RELAXSTEPS
S(q->0) 的 RELAXSTEPS
```

`fixed_charge_static` 示例：

```powershell
python standalone\extract_reorganization_energy.py `
  --method fixed_charge_static `
  --defect 0_v_Si `
  --charge 1 `
  -i \path\to\E0\REPORT \path\to\E0_prime\REPORT \path\to\Eq\REPORT \path\to\Eq_prime\REPORT `
  -o \path\to\out\reorganization_energy.log
```

其中 `-i` 后四个文件的顺序必须是：

```text
E0        : q=0 电荷在 Q0 构型的 REPORT
E0_prime  : q=0 电荷在 Qq 构型的 REPORT
Eq        : q 电荷在 Qq 构型的 REPORT
Eq_prime  : q 电荷在 Q0 构型的 REPORT
```

默认还会在 summary log 同目录生成：

```text
detail_reorg.log
```

如果需要指定详细日志路径：

```powershell
python standalone\extract_reorganization_energy.py `
  --method fixed_charge_static `
  --defect 0_v_Si `
  --charge -1 `
  -i E0_REPORT E0_prime_REPORT Eq_REPORT Eq_prime_REPORT `
  -o reorganization_energy.log `
  --detail-output detail_reorg_qm1.log
```

`extract_reorganization_energy.py` 输出的 summary log 可被 `plot_ccd.py` 直接读取。当前推荐格式为宽表：

```text
idx S[0->+1] S[+1->0] S[0->-1] S[-1->0]
0 0.82 0.76 0.91 0.88
```

读取规则：

- `--charge 1` 读取 `S[0->+1]` 和 `S[+1->0]`
- `--charge -1` 读取 `S[0->-1]` 和 `S[-1->0]`

也兼容旧长表字段，例如：

```text
defect charge lambda_0_to_q_eV lambda_q_to_0_eV
```

或：

```text
defect charge S_0_to_q_eV S_q_to_0_eV
```

### 2. 计算 DeltaQ

默认输出为当前运行目录的 `dQ.log`：

```powershell
python standalone\calculate_deltaQ_config.py `
  -i \path\to\q_0\scf\atom.config \path\to\q_1\scf\atom.config
```

输出：

```text
.\dQ.log
```

`dQ.log` 为文本格式，包含：

```text
deltaQ_sqrt_amu_A
deltaQ_au
mean_displacement_A
max_displacement_A
min_displacement_A
rmsd_A
```

如果需要 JSON：

```powershell
python standalone\calculate_deltaQ_config.py `
  -i \path\to\q_0\scf\atom.config \path\to\q_1\scf\atom.config `
  -o \path\to\out\dQ.json
```

可选方法：

```powershell
python standalone\calculate_deltaQ_config.py `
  -i ref.atom.config target.atom.config `
  --method simple `
  -o dQ_simple.log
```

`normalized` 是默认方法，与原 `dQ_calculator.py` 公式保持一致；`simple` 是未归一化质量加权位移。

注意：如果两个结构完全相同，归一化方法的位移范数为 0，原公式会产生 `nan`。这属于零位移边界，不代表路径参数失败。

### 3. 手动指定 DeltaG

`plot_ccd.py` 使用 `--deltaG` 手动指定 CCD 中两态的垂直能量差：

```text
--deltaG
```

它表示 CCD 中 `E_q - E_0`，单位 eV。脚本也保留 `--delta-e` 作为兼容别名，但推荐文档和新命令统一使用 `--deltaG`。

只用“重组能 + DeltaQ”可以确定两条抛物线的曲率和横向距离，但不能唯一确定两条曲线的相对垂直能量差。因此本通用流程不自动读取本征能级或自动推断 `DeltaG`，而是由用户手动赋值。

| 场景 | 推荐输入 |
|---|---|
| 两态简并 CCD | 不指定 `--deltaG`，默认 `--deltaG 0` |
| 指定非简并能量差 | 手动给 `--deltaG <数值>` |
| “本征能级 + 重组能” | 自行提取本征能级相关能量差后，手动给 `--deltaG <数值>` |
| 热力学转变能级相关 CCD | 可参考任务一得到的 `defect_results.yaml`，但推荐仍是确认数值后手动给 `--deltaG` |

### 4. 绘制 CCD

最小输入是：

- `reorganization_energy.log`
- 命令行给出的 `--deltaQ`
- 缺陷编号/名称
- 电荷态

两态简并 CCD：

```powershell
python standalone\plot_ccd.py `
  -i \path\to\result\reorganization_energy.log `
  --deltaQ 2.4 `
  --defect 0_v_Si `
  --charge 1 `
  -o \path\to\out\ccd\0_v_Si_qp1.png
```

上面命令默认：

```text
--mode intrinsic
--deltaG 0
```

非简并 CCD：

```powershell
python standalone\plot_ccd.py `
  -i \path\to\result\reorganization_energy.log `
  --deltaQ 2.4 `
  --defect 0_v_Si `
  --charge 1 `
  --deltaG 0.25 `
  -o \path\to\out\ccd\0_v_Si_qp1_deltaG025.png
```

兼容旧的 `--delta-e` 写法：

```powershell
python standalone\plot_ccd.py `
  -i \path\to\result\reorganization_energy.log `
  --deltaQ 2.4 `
  --defect 0_v_Si `
  --charge 1 `
  --delta-e 0.25 `
  -o \path\to\out\ccd\0_v_Si_qp1_deltaE025.png
```

兼容旧的 DeltaQ log 读取方式：

```powershell
python standalone\plot_ccd.py `
  -i \path\to\result `
  --deltaq-log \path\to\result\deltaQ_scf_atom.log `
  --defect 0_v_Si `
  --charge 1 `
  -o \path\to\out\ccd\0_v_Si_qp1.png
```

输出：

- 图像：由 `-o` 的后缀决定，例如 `.png`、`.svg`、`.pdf`
- 同名 `.log`：记录 DeltaQ、重组能、`deltaE_eV`、交点和势垒。这里的 `deltaE_eV` 对应手动输入的 `--deltaG`。

## 当前脚本是否能读取本征缺陷能级

**不能。当前 `standalone/` 中没有脚本会自动读取“本征缺陷能级”的能量值。**

需要区分三类量：

| 名称 | 当前脚本支持情况 | 说明 |
|---|---|---|
| 热力学转变能级 | 支持 | `formation_energy.py` 从形成能交点计算，写入 `defect_results.yaml` 的 `transition_levels` |
| 缺陷带编号/占据变化 | 支持 | `generate_occ.py` 从 `OUT.OCC0/OUT.OCC1` 判断缺陷带范围和占据差 |
| 本征缺陷能级能量值 | 不支持 | 当前没有脚本从 eigenvalue 输出中读取具体缺陷能级能量 |

`generate_defect_input.py` 的功能是：

- 比较 bulk 和 defect 的 `atom.config`
- 确定缺陷位置和 bound 位置
- 为显式指定的带电态静态计算目录生成 `defect.input`

它不读取本征能级。

`generate_occ.py` 的功能是：

- 读取 `OUT.KPT`
- 读取 `OUT.OCC0` 和 `OUT.OCC1`
- 根据占据数差异判断缺陷带编号
- 生成 `occ.input` 和 `IN.OCC`

它只处理占据数，不读取 eigenvalue 能量。因此它可以告诉你“哪些 band 可能是缺陷能级”，但不能告诉你“这些缺陷能级的能量是多少”。

如果 CCD 采用“本征能级 + 重组能”路线，目前需要你从对应计算输出中自行提取本征能级能量差，然后作为 `DeltaG` 传给：

```powershell
python standalone\plot_ccd.py `
  -i \path\to\result\reorganization_energy.log `
  --deltaQ 2.4 `
  --defect 0_v_Si `
  --charge 1 `
  --deltaG <手动指定的DeltaG> `
  -o \path\to\out\ccd\0_v_Si_intrinsic.png
```

若后续需要自动化，应新增一个专门的本征能级提取脚本，明确读取哪个文件、band index 如何绑定、能量参考是 VBM、真空能级还是绝对本征值。

## 脚本总览

| 脚本 | 主要功能 | `-i` 含义 | `-o` 含义 |
|---|---|---|---|
| `generate_defect_input.py` | 根据显式给定的 bulk/q0/q± 目录生成 `defect.input` | `-ibulk`、`-iq0`、`-iq1`、`-iqm1` 或 `--charge-dir` | 输出目录，或单个 `defect.input` 文件 |
| `generate_occ.py` | 根据 `OUT.OCC0/1` 生成 `occ.input` 和 `IN.OCC` | 包含 `OUT.KPT`、`OUT.OCC0`、`OUT.OCC1` 的目录 | 输出目录，或无自旋主 `occ.input` 文件 |
| `correction.py` | 收集 ICIC/PA 修正并计算 `E_f0`/`transition_levels` | 推荐显式 `-ibulk/-iq0/-iq1/-iqm1/--charge-dir`；兼容项目根目录扫描 | `collect` 输出 `correction_results.yaml` |
| `calculate_chemical_potential_bounds.py` | 由 `chemical/` 杂相总能计算化学势范围 | 项目根目录或 `chemical` 目录 | 主 YAML 输出；相关 log 自动派生 |
| `formation_energy.py` | 根据修正结果计算形成能和热力学转变能级 | `correction_results.yaml` | `defect_results.yaml` |
| `extract_reorganization_energy.py` | 从显式指定的 `RELAXSTEPS` 或 `REPORT` 提取重组能 | `structural_relaxation` 为两个 `RELAXSTEPS`；`fixed_charge_static` 为四个 `REPORT` | `reorganization_energy.log`；默认同时写 `detail_reorg.log` |
| `calculate_deltaQ_config.py` | 由两个 `atom.config` 计算 DeltaQ | 两个结构文件：参考结构、目标结构 | 输出文件，默认当前目录 `dQ.log`；`.json` 后缀时输出 JSON |
| `plot_ccd.py` | 根据重组能、DeltaQ 和 `deltaE` 绘制 CCD | 重组能 log，或 `result` 目录/项目根目录 | 图像文件或输出基名 |

## 输入输出关系

形成能和转变能级：

```text
generate_defect_input.py / generate_occ.py
  -> defect.input / occ.input / IN.OCC
  -> correction.py collect
      -> correction_results.yaml
         (含 corrections, E_f0, transition_levels, formation_energy_table)
      -> formation_energy.py  [可选重算/另存]
          -> defect_results.yaml
          -> E_forms/E_formation_corrected_<defect>.txt
```

化学势：

```text
calculate_chemical_potential_bounds.py
  -> chemical_potential_bounds.yaml
  -> phase/constraints/bounds logs
  -> 手动选择化学势并更新 config.yaml
```

CCD：

```text
extract_reorganization_energy.py -> reorganization_energy.log
calculate_deltaQ_config.py -> dQ.log 或 dQ.json
手动指定 DeltaG
  -> plot_ccd.py
      -> CCD 图像
      -> CCD 同名 log
```

## 路径适配建议

如果脚本运行在与数据不同的目录，建议始终使用绝对路径。例如：

```powershell
python \get_aligned_transition_level\standalone\formation_energy.py `
  -i \data\case01\result\correction_results.yaml `
  -o \data\case01\result\defect_results.yaml `
  --project-root \data\case01
```

如果输入文件来自其他目录，但还需要项目配置或 bulk 信息，显式给 `--project-root` 或相关文件参数，避免脚本从当前 `standalone/` 附近误推断。

## 快速检查命令

检查所有脚本能否通过语法编译：

```powershell
python -m py_compile `
  standalone\calculate_chemical_potential_bounds.py `
  standalone\calculate_deltaQ_config.py `
  standalone\correction.py `
  standalone\extract_reorganization_energy.py `
  standalone\formation_energy.py `
  standalone\generate_defect_input.py `
  standalone\generate_occ.py `
  standalone\plot_ccd.py
```

查看单个脚本帮助：

```powershell
python standalone\extract_reorganization_energy.py --help
python standalone\plot_ccd.py --help
python standalone\calculate_deltaQ_config.py --help
python standalone\formation_energy.py --help
python standalone\generate_occ.py --help
```

## 注意事项

- `correction.py collect` 会在收集修正项后计算 `E_f0` 和 `transition_levels`；前提是 `E_bulk`、`VBM` 和所需元素化学势已经可读取。
- `calculate_chemical_potential_bounds.py` 只给出化学势允许范围，不会自动修改 `config.yaml`。
- 如果先运行 `correction.py collect` 后才修改化学势，需要重新运行 `correction.py collect`，或用 `formation_energy.py` 基于更新后的化学势重新生成形成能和转变能级。
- `plot_ccd.py` 默认 `--deltaG 0`，表示两态简并 CCD；若要表达实际非简并能量差，必须手动指定 `--deltaG`。
- `calculate_deltaQ_config.py` 的 `normalized` 方法在两个结构完全相同时会因位移范数为 0 得到 `nan`；这是原公式的零位移边界。
- `correction.py prepare` 是兼容旧项目结构的批量模式；完全普适的 `defect.input` 生成推荐使用 `generate_defect_input.py` 的显式目录参数。
- `generate_defect_input.py` 处理多个非零电荷态时，`-o` 必须是目录。
- `formation_energy.py` 若输入 YAML 不在项目 `result/` 下，建议显式指定 `--project-root`。
