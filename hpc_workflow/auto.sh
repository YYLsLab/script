#!/bin/sh
#SBATCH --partition=3090
#SBATCH --job-name=HfO2_defect
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --gpus-per-task=1

module load mkl mpi
module load cuda/11.6
module load pwmat/20231124
#module load mkl mpi
#module load cuda/12.1
#module load pwmat/20241205_4090
# ==========================================
# 错误追踪配置
# ==========================================
ERROR_LOG="${SLURM_SUBMIT_DIR}/job_error.log"
echo "=== Job Started at $(date) ===" > "$ERROR_LOG"

check_error() {
    local stage_name=$1
    local output_file=$2
    if [ $? -ne 0 ]; then
        echo "[$(date '+%H:%M:%S')] CRITICAL: $stage_name 执行失败 (Exit Code != 0)" >> "$ERROR_LOG"
    fi
    if [ -f "$output_file" ]; then
        grep -iE "error|fault|divergence|fail|not converged" "$output_file" >> "$ERROR_LOG" 2>/dev/null
    fi
}

# ==========================================
# 参数设置区
# ==========================================
GLOBAL_PBE="false"

generate_input() {
    local JOB_TYPE=$1
    local IN_WG=$2
    local NUM_E=$3
    local XC="HSE"
    
    if [ "$JOB_TYPE" = "RELAX" ] || [ "$GLOBAL_PBE" = "true" ]; then
        XC="PBE"
    fi

    cat << INNER_EOF > etot.input
$SLURM_NPROCS  1
IN.ATOM = atom.config
JOB = $JOB_TYPE
$( [ "$JOB_TYPE" = "RELAX" ] && echo "RELAX_DETAIL = 1 100 0.05 0 0.01" )
IN.PSP1 = Hf.SG15.PBE.UPF
IN.PSP2 = O.SG15.PBE.UPF
IN.PSP3 = H.SG15.PBE.UPF
$( [ "$XC" = "HSE" ] && echo "HSEMASK_PSP1 = 0.059 3.50#Hf" )
$( [ "$XC" = "HSE" ] && echo "HSEMASK_PSP2 = 0.059 3.50#O" )
MP_N123 = 1 1 1 0 0 0
XCFUNCTIONAL = $XC
ECUT = 50
ECUT2 = 200
SPIN = 2
IN.WG = $IN_WG
E_ERROR   =   0     
RHO_ERROR =   5.0E-5
INNER_EOF
    [ -n "$NUM_E" ] && echo "NUM_ELECTRON = $NUM_E" >> etot.input
}

# ==========================================
# 1. 主目录：PBE 结构优化
# ==========================================
generate_input "RELAX" "F"
mpirun -np $SLURM_NPROCS PWmat | tee output 2>> "$ERROR_LOG"
check_error "Main_RELAX" "output"

# ==========================================
# 2. 主目录：Static & DOS
# ==========================================
mkdir -p static                                   
cp final.config static/atom.config 2>> "$ERROR_LOG"
cp *.UPF static/ 2>> "$ERROR_LOG"
cd static/
generate_input "SCF" "F"
mpirun -np $SLURM_NPROCS PWmat | tee output 2>> "$ERROR_LOG"
check_error "Main_Static" "output"

mkdir -p dos                                      
cp atom.config dos/ && cp *.UPF dos/ 2>> "$ERROR_LOG"
cd dos/
generate_input "DOS" "T"
[ -f "../OUT.WG" ] && ln -s ../OUT.WG IN.WG
mpirun -np $SLURM_NPROCS PWmat | tee output 2>> "$ERROR_LOG"
check_error "Main_DOS" "output"

cd ../../

# ==========================================
# 提取电子数
# ==========================================
if [ -f "./static/REPORT" ]; then
    BASE_ELEC=$(grep "NUM_ELECTRON" ./static/REPORT | tail -1 | awk -F'=' '{print $2}' | tr -d ' ')
    ELEC_M1=$(awk "BEGIN {printf \"%.6f\", $BASE_ELEC - 1}")
    ELEC_P1=$(awk "BEGIN {printf \"%.6f\", $BASE_ELEC + 1}")
else
    echo "[$(date '+%H:%M:%S')] LOGIC ERROR: 找不到 static/REPORT，无法进行后续 lambda 计算" >> "$ERROR_LOG"
    exit 1  
fi

# ==========================================
# 3. lambda-1 流程
# ==========================================
mkdir -p lambda-1                                 
cp final.config lambda-1/atom.config 2>> "$ERROR_LOG"
cp *.UPF lambda-1/ 2>> "$ERROR_LOG"
cd lambda-1/
generate_input "RELAX" "F" "$ELEC_M1"
mpirun -np $SLURM_NPROCS PWmat | tee output 2>> "$ERROR_LOG"
check_error "Lambda-1_RELAX" "output"

mkdir -p static                                   
cp final.config static/atom.config && cp *.UPF static/
cd static/
generate_input "SCF" "F" "$ELEC_M1"
mpirun -np $SLURM_NPROCS PWmat | tee output 2>> "$ERROR_LOG"
check_error "Lambda-1_Static" "output"

mkdir -p dos                                      
cp atom.config dos/ && cp *.UPF dos/
cd dos/
generate_input "DOS" "T" "$ELEC_M1"
[ -f "../OUT.WG" ] && ln -s ../OUT.WG IN.WG
mpirun -np $SLURM_NPROCS PWmat | tee output 2>> "$ERROR_LOG"
check_error "Lambda-1_DOS" "output"

cd ../../../

# ==========================================
# 4. lambda+1 流程
# ==========================================
mkdir -p lambda+1                                 
cp final.config lambda+1/atom.config 2>> "$ERROR_LOG"
cp *.UPF lambda+1/ 2>> "$ERROR_LOG"
cd lambda+1/
generate_input "RELAX" "F" "$ELEC_P1"
mpirun -np $SLURM_NPROCS PWmat | tee output 2>> "$ERROR_LOG"
check_error "Lambda+1_RELAX" "output"

mkdir -p static                                   
cp final.config static/atom.config && cp *.UPF static/
cd static/
generate_input "SCF" "F" "$ELEC_P1"
mpirun -np $SLURM_NPROCS PWmat | tee output 2>> "$ERROR_LOG"
check_error "Lambda+1_Static" "output"

mkdir -p dos                                      
cp atom.config dos/ && cp *.UPF dos/
cd dos/
generate_input "DOS" "T" "$ELEC_P1"
[ -f "../OUT.WG" ] && ln -s ../OUT.WG IN.WG
mpirun -np $SLURM_NPROCS PWmat | tee output 2>> "$ERROR_LOG"
check_error "Lambda+1_DOS" "output"

cd ../../../

echo "=== Job Finished at $(date) ===" >> "$ERROR_LOG"

if [ -f "./extract_data.sh" ]; then
    sh extract_data.sh
fi
