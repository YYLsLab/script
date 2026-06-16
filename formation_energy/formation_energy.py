import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# ==============================================================
# 缺陷形成能公式：
# ΔH(α, q) = E(α, q) - E(bulk) + Σn_i(E_i + μ_i) + q[E_VBM(bulk) + E_f + V_q_R-V_b_R]                       # 无ICIC修正
# ΔH(α, q) = E(α, q) - E(bulk) + Σn_i(E_i + μ_i) + q[E_VBM(bulk) + E_f + V_0_R-V_b_R] + （E_infi-E_p）/ε     # ICIC修正
#
# 其中：
#   ΔH(α, q)  : 缺陷形成能 （y轴）
#   E(α, q)   : 缺陷在带电态q下的总能
#   E(bulk)   : 完美超胞总能
#   n_i       : 原子变化数  (原子被放入：ni<0, 原子被拿出ni>0)
#   E_i       : 第i种原子的总能
#   μ_i       : 第i种原子的化学势
#   q         : 缺陷带电量
#   E_VBM(bulk): bulk价带顶能量
#   E_f       : 费米能级（相对于VBM）（自变量x轴，不是从体系中读取的费米能级）
#   V_q_R        : 带电缺陷体系中距离缺陷最远的原子势
#   V_0_R        : 电中性缺陷体系中距离缺陷最远的原子势
#   V_b_R        : 无缺陷bulk中距离缺陷最远的原子势
#
# ICIC修正相关：
#   E_infi-E_p   : 带电缺陷无限大超胞与周期性超胞电势能之差
#   ε            : 材料相对介电常数 static dielectric constant

# ==============================================================

# ========================1. 基础参数设置==============================
out_fig = "O_i.png"  # 输出图片名称
icic = 1 # ICIC修正开关  1开0关  参考官网ICIC教程
epsilon = 22 # ICIC修正需要参数 材料相对静态介电常数

# bulk参数
E_host = -.96706686534013E+05  # 无缺陷态体系总能
E_VBM_host = 3.906262 # 无缺陷态体系VBM数值
Eg = 4.410540 # 带隙大小
Ef_list = np.linspace(0, Eg+1, 500)  # 费米能级取点数目，一般不用修改

# poor-rich 条件控制
# 所考虑的缺陷中，每多出现一个元素，就多写一行，先计算好每种元素的化学势再往里面填数字
atom_ref = {
    # O-rich
    "Ta": {"E": -1594.77485, "mu": -10.52 },    # "mu": -10.52
    "O": {"E":-432.40019408622, "mu":0 },       # "mu": -4.21   mu=0代表相应元素rich
    # "Se": {"E": -771.24761805382/3, "mu": -3.13},
}

# 画图时使用的颜色， 每种缺陷一种色系，不同电荷态颜色用深浅区分，确保颜色数量大于等于缺陷数量即可
colors = ["#1f77b4", '#d62728', "#ff7f0e", "#2ca02c", "#d62728"]

# ====================== 缺陷参数 ======================

defect_list = [
    # 缺陷1
    {
        "name": "O_i17" ,  # 缺陷名称，自定义
        "n_i": {"O": -1 ,
                # "Ta": 0,
                },  # 缺陷包含的元素，从体系中拿出为+，放进体系中为-
        "V_0_R": -63.86553,  # ICIC修正用V_0_R,不用V_q_R，不使用时填0即可
        # 电荷态参数 
        "charges": [
            #  -0.186807970      0.795776283      0.011742475
            {"q": 0,  "E_defect": -.97138159203612E+05 , "V_q_R": -63.86553, "V_b_R": -63.79623, "E_infi": 0.0, "E_p": 0.0},
            # -0.186291440      0.796176834      0.012076734
            {"q": -1, "E_defect": -.97129803707830E+05 , "V_q_R":-44.90083, "V_b_R":-63.79623, "E_infi":1.43012, "E_p": 0.47952},
            # -0.171839877      0.793339823      0.017945685
            {"q": -2, "E_defect": -.97122165150792E+05, "V_q_R":-63.71493, "V_b_R": -63.79623, "E_infi": 5.67292, "E_p": 1.96052},
        ]
    },

 {
        "name": "O_i39" ,  # 缺陷名称，自定义
        "n_i": {"O": -1 ,
                # "Ta": 0,
                },  # 缺陷包含的元素，从体系中拿出为+，放进体系中为-
        "V_0_R": -63.18687,  # ICIC修正用V_0_R,不用V_q_R，不使用时填0即可
        # 电荷态参数 
        "charges": [
            #  -0.186807970      0.795776283      0.011742475
            {"q": 0,  "E_defect": -97138.05224 , "V_q_R": -63.18687, "V_b_R": -63.22735, "E_infi": 0.0, "E_p": 0.0},
            # -0.186291440      0.796176834      0.012076734
            {"q": -1, "E_defect": -97129.29359 , "V_q_R":-63.21549, "V_b_R":-63.22735, "E_infi":1.47666, "E_p": 0.35922},
            # -0.171839877      0.793339823      0.017945685
            {"q": -2, "E_defect": -97121.75, "V_q_R":-63.24683, "V_b_R": -63.22735, "E_infi": 5.88159, "E_p": 1.39775},
        ]
    },

 {
        "name": "O_i6" ,  # 缺陷名称，自定义
        "n_i": {"O": -1 ,
                # "Ta": 0,
                },  # 缺陷包含的元素，从体系中拿出为+，放进体系中为-
        "V_0_R":-63.78504,  # ICIC修正用V_0_R,不用V_q_R，不使用时填0即可
        # 电荷态参数 
        "charges": [
            #  -0.186807970      0.795776283      0.011742475
            {"q": 0,  "E_defect":-97138.13302 , "V_q_R": -63.78504, "V_b_R":-63.81282, "E_infi": 0.0, "E_p": 0.0},
            # -0.186291440      0.796176834      0.012076734
            {"q": -1, "E_defect":-97131.42783 , "V_q_R":0, "V_b_R":-63.81282, "E_infi":1.44333, "E_p": 0.3476},
            # -0.171839877      0.793339823      0.017945685
            {"q": -2, "E_defect":-97127.43561, "V_q_R":0, "V_b_R": -63.81282, "E_infi": 5.86143, "E_p": 1.38333},
        ]
    },

    # # 缺陷2
    # {
    #     "name": "V_Cl",
    #     "n_i": {"Se": 0,
    #             "Bi":+1
    #             }, 
    #     "V_0_R": 0.0,
    #     "charges": [
    #         {"q": 0,  "E_defect":-.61568874337909E+05, "V_q_R": -38.60185, "V_b_R": -38.51117, "E_infi": 0.0, "E_p": 0.0},
    #         {"q": +1, "E_defect":-.61570119129369E+05, "V_q_R":-39.17429, "V_b_R": -38.51117, "E_infi": 0.0, "E_p": 0.0},
    #         {"q": +2, "E_defect": -.61566253037408E+05 , "V_q_R": -39.89512, "V_b_R": -38.51117, "E_infi": 0.0, "E_p": 0.0},
    #     ]
    # },
]

# --------------------------
# 形成能计算函数
# --------------------------
def calc_formation_energy(defect_data, Ef):
    E_defect = defect_data["E_defect"]
    q = defect_data["q"]
    n_i = defect_data["n_i"]
    
    V_0_R = defect_data.get("V_0_R", 0.0)
    V_q_R = defect_data.get("V_q_R", 0.0)
    V_b_R = defect_data.get("V_b_R", 0.0)
    
    if icic == 0:
        deltaV = V_q_R - V_b_R
        charge_term = q * (E_VBM_host + Ef + deltaV)
    else:
        deltaV = V_0_R - V_b_R
        E_infi = defect_data.get("E_infi", 0.0)
        E_p = defect_data.get("E_p", 0.0)
        E_corr = (E_infi - E_p) / epsilon
        charge_term = q * (E_VBM_host + Ef + deltaV) + E_corr

    chem_pot_term = sum(n * (atom_ref[atom]["E"] + atom_ref[atom]["mu"]) for atom, n in n_i.items())
    return E_defect - E_host + chem_pot_term + charge_term

# --------------------------
# 数据整合和计算
# --------------------------
unique_defects = [d["name"] for d in defect_list] # 缺陷种类
defect_arrays = {}  
# 键为缺陷的名字，值为这个缺陷不同带电态下的形成能矩阵：第一列为费米能级，第二列开始为形成能，最后一列为形成能最小值

for defect_info in defect_list:
    defect_name = defect_info["name"]
    n_i = defect_info["n_i"]
    V_0_R = defect_info["V_0_R"]
    charge_states = defect_info["charges"]

    dH_cols = []
    dH_Ef0_cols = []
    q_values = []
    # 遍历该缺陷所有电荷态
    for chg in charge_states:
        q = chg["q"]
        q_values.append(q)
        # 合并公共参数+当前电荷态参数，传给计算函数
        full_data = {
            **chg,  # v_b_R 参数内容在这里
            "n_i": n_i,
            "V_0_R": V_0_R
        }
        dH = np.array([calc_formation_energy(full_data, Ef) for Ef in Ef_list])
        dH_0 = np.array([calc_formation_energy(full_data, Ef=0)])
        dH_Ef0_cols.append(dH_0.reshape(-1, 1))
        dH_cols.append(dH.reshape(-1, 1))

    # 计算各费米能级下的最小形成能
    dH_matrix = np.hstack(dH_cols)
    min_vals = np.min(dH_matrix, axis=1)
    min_col = min_vals.reshape(-1, 1)
    
    # 拼接形成能数据
    ef_col = Ef_list.reshape(-1, 1)
    data_array = np.hstack([ef_col] + dH_cols + [min_col])
    header = ["Ef"] + [f"q={q}" for q in q_values] + ["min_formation_energy"]
    final_array = np.vstack([header, data_array]).astype(object)

    defect_arrays[defect_name] = final_array

    # 使用Ef=0处的形成能计算转变能级（Eq-Eq'）/(q'-q）
    dH_Ef0_matrix = np.hstack(dH_Ef0_cols)
    Etrans_matrix = np.zeros((len(q_values), len(q_values)))
    for i, qi in enumerate(q_values):
        for j,qj in enumerate(q_values):
            if i != j:
                Etrans_matrix[i][j]= ( dH_Ef0_matrix[0,i]-dH_Ef0_matrix[0,j] ) /(qj - qi)
    Etrans_matrix = np.round(Etrans_matrix, 3) # 保留小数 ,相对于VBM的转变能级
    Etrans_CBM = np.round(Etrans_matrix - Eg,3) # 得到相对于CBM的矩阵
    np.fill_diagonal(Etrans_CBM, 0)
    
    df_VBM = pd.DataFrame(Etrans_matrix, index=q_values, columns=q_values)
    vbm_header = pd.DataFrame( [["VBM"] + [""]*(len(df_VBM.columns)-1)],  index=[""], columns=df_VBM.columns)
    df_VBM = pd.concat([vbm_header, df_VBM], axis=0)
    
    # 空行
    empty_row = pd.DataFrame([[""] * df_VBM.shape[1]], index=[""], columns=df_VBM.columns)

    df_CBM = pd.DataFrame(Etrans_CBM,index=q_values, columns=q_values)
    cbm_header = pd.DataFrame( [["CBM"] + [""]*(len(df_VBM.columns)-1)],  index=[""], columns=df_VBM.columns)
    df_CBM = pd.concat([cbm_header, df_CBM], axis=0)

    df_total = pd.concat([df_VBM, empty_row, df_CBM], axis=0)
    # df = pd.DataFrame(Etrans_matrix, index=q_values,columns=q_values)
    df_total.to_csv(f"Etrans_{defect_name}.csv", encoding='utf-8-sig')



# --------------------------
# 绘图，导出CSV
# --------------------------
plt.figure(figsize=(8, 6))

# 左子图：缺陷所有电荷态曲线
plt.subplot(1, 2, 1)
for i, (name, arr) in enumerate(defect_arrays.items()):
    df = pd.DataFrame(arr[1:], columns=arr[0])
    df.to_csv(f"Eform_{name}.csv", index=False, encoding='utf-8-sig')

    Ef = arr[1:, 0].astype(float)
    q_cols = arr[1:, 1:-1].astype(float)
    q_values = [int(q_str.replace("q=", "")) for q_str in arr[0, 1:-1]]
    num_charges = q_cols.shape[1]

    for j in range(q_cols.shape[1]):
        q = q_values[j]
        label = f"{name} q={q}" if q == 0 else f"{name} q={q:+d}"

        base_color = colors[i]
        brightness = 1 - (j / num_charges) * 0.7
        light_color = mcolors.rgb_to_hsv(mcolors.to_rgb(base_color))
        light_color[2] = brightness
        light_color = mcolors.hsv_to_rgb(light_color)

        plt.plot(Ef, q_cols[:, j],
                 color=light_color,
                 linestyle='-',
                 linewidth=1.5,
                 label=label)

plt.title("Charge States", fontsize=14)
plt.xlabel("Fermi Level (eV)", fontsize=12)
plt.ylabel("Formation Energy (eV)", fontsize=12)
plt.grid(alpha=0.3)
plt.legend(fontsize=9)
plt.xlim(min(Ef_list), max(Ef_list))
# ==== 加直线（左图）====
plt.axvline(x=0, color='gray', linestyle='--', linewidth=1)          # 价带顶
plt.text(0, plt.ylim()[1]*0.9, 'VBM', fontsize=10, ha='right', va='top')
plt.axvline(x=Eg, color='gray', linestyle='--', linewidth=1)        # 导带底
plt.text(Eg, plt.ylim()[1]*0.9, 'CBM', fontsize=10, ha='left', va='top')

# 右子图：最小形成能曲线
plt.subplot(1, 2, 2)
for i, (name, arr) in enumerate(defect_arrays.items()):
    Ef = arr[1:, 0].astype(float)
    min_line = arr[1:, -1].astype(float)
    plt.plot(Ef, min_line,
             color=colors[i],
             linewidth=1.5,
             label=f"{name}")

plt.title("Min Formation energy", fontsize=14)
plt.xlabel("Fermi Level (eV)", fontsize=12)
plt.ylabel("Formation Energy (eV)", fontsize=12)
plt.grid(alpha=0.3)
plt.legend(fontsize=10)
plt.xlim(0,Eg)

plt.tight_layout()
plt.savefig(out_fig, dpi=300)
plt.show()