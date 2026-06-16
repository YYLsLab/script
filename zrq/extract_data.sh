#!/bin/bash
# ======================================================
# 自动提取 PWmat 能级、重组能及转变能级 (Transition Level)
# ======================================================

# 1. 创建输出数据文件夹
OUTPUT_DIR="Data_Summary"
mkdir -p "$OUTPUT_DIR"
OUT_FILE="${OUTPUT_DIR}/Results_Summary.txt"

> "$OUT_FILE"

echo "======================================================" >> "$OUT_FILE"
echo "        PWmat 缺陷计算数据提取报告" >> "$OUT_FILE"
echo "        生成时间: $(date '+%Y-%m-%d %H:%M:%S')" >> "$OUT_FILE"
echo "======================================================" >> "$OUT_FILE"
echo "" >> "$OUT_FILE"

# ==========================================
# 工具函数：提取 REPORT 中的总能 (E_tot)
# ==========================================
get_total_energy() {
    local report_file=$1
    if [ -f "$report_file" ]; then
        # 精准查找 E_tot(eV) 行，以等号分割后，只截取第一个数值 (总能)，过滤掉第二个 dE 误差数值
        grep "E_tot(eV)" "$report_file" | tail -1 | awk -F'=' '{print $2}' | awk '{print $1}'
    else
        echo ""
    fi
}

# ==========================================
# 工具函数：提取能带带边 (OUT.OCC)
# ==========================================
extract_band_edge() {
    local occ_file=$1
    local title=$2
    echo "------------------------------------------------------" >> "$OUT_FILE"
    echo "📊 能级提取: $title" >> "$OUT_FILE"
    echo "------------------------------------------------------" >> "$OUT_FILE"
    if [ ! -f "$occ_file" ]; then
        echo "❌ 未找到文件: $occ_file" >> "$OUT_FILE"
        return
    fi
    awk '/==========  SPIN/ {print "\n" $0; found=0; unocc=0; o1=""; o2=""; o3=""}
    NF >= 3 && $1 ~ /^[0-9]+$/ {
        if ($3 > 0.5) { o1=o2; o2=o3; o3=$0 }
        else if ($3 <= 0.5 && found == 0) {
            if(o1!="")print o1; if(o2!="")print o2; if(o3!="")print o3;
            print "----------------------------------------- (费米能级)";
            print $0; found=1; unocc=1
        }
        else if (found == 1 && unocc < 3) { print $0; unocc++ }
    }' "$occ_file" >> "$OUT_FILE"
    echo "" >> "$OUT_FILE"
}

# ==========================================
# 工具函数：计算重组能 (RELAXSTEPS)
# ==========================================
calculate_reorg_energy() {
    local relax_file=$1
    local title=$2
    echo "------------------------------------------------------" >> "$OUT_FILE"
    echo "⚡ 重组能计算: $title" >> "$OUT_FILE"
    echo "------------------------------------------------------" >> "$OUT_FILE"
    if [ ! -f "$relax_file" ]; then
        echo "❌ 未找到文件: $relax_file" >> "$OUT_FILE"
        return
    fi
    awk '/E=/ {for(i=1;i<=NF;i++){if($i=="E="){ce=$(i+1); if(fe=="")fe=ce; le=ce}}}
    END {if(fe!=""){printf "初始能量: %s eV\n最终能量: %s eV\n==> 重组能 (λ): %.6f eV\n", fe, le, fe-le}
    else {print "❌ 读取能量失败"}}' "$relax_file" >> "$OUT_FILE"
    echo "" >> "$OUT_FILE"
}

# ==========================================
# 2. 开始执行提取任务
# ==========================================
echo "正在提取中性态基准数据..."

# 提取主目录中性态总能
NEUTRAL_E=$(get_total_energy "./static/REPORT")
echo "主目录中性态总能: $NEUTRAL_E eV" >> "$OUT_FILE"
extract_band_edge "./static/OUT.OCC" "主目录 (中性态)"

# 自动遍历 lambda 开头的文件夹
for dir in lambda*/; do
    if [ ! -d "$dir" ]; then continue; fi
    dname=${dir%/}
    echo "正在处理 $dname ..."

    echo "======================================================" >> "$OUT_FILE"
    echo "        📂 文件夹: $dname 结果汇总" >> "$OUT_FILE"
    echo "======================================================" >> "$OUT_FILE"

    # 1. 提取带电态总能并计算转变能级
    CHARGED_E=$(get_total_energy "./${dname}/static/REPORT")
    if [ -n "$NEUTRAL_E" ] && [ -n "$CHARGED_E" ]; then
        DIFF=$(awk -v n="$NEUTRAL_E" -v c="$CHARGED_E" 'BEGIN {printf "%.6f", c-n}')
        echo "💡 转变能级计算 ($dname - 主目录):" >> "$OUT_FILE"
        echo "   中性态总能: $NEUTRAL_E eV" >> "$OUT_FILE"
        echo "   带电态总能: $CHARGED_E eV" >> "$OUT_FILE"
        echo "   ==> 转变能级 (Etot_diff): $DIFF eV" >> "$OUT_FILE"
        echo "" >> "$OUT_FILE"
    fi

    # 2. 提取能带数据
    extract_band_edge "./${dname}/static/OUT.OCC" "${dname} 状态"

    # 3. 提取重组能数据
    calculate_reorg_energy "./${dname}/RELAXSTEPS" "${dname} 构型"
done

echo "✅ 所有数据处理完成！结果保存在: ${OUT_FILE}"
