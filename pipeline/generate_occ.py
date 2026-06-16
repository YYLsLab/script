import numpy as np
import os

class OccupationAnalyzer:
    def __init__(self):
        self.nkpt = 0
        self.spin_type = 1  # 1: 无自旋, 2: 有自旋
        
    def read_nkpt(self):
        """读取k点数量"""
        if os.path.exists('OUT.KPT'):
            with open('OUT.KPT', 'r') as f:
                self.nkpt = int(f.readline().strip())
        else:
            raise FileNotFoundError("OUT.KPT文件不存在")
        print(f"k-point: {self.nkpt}")
    
    def count_lines(self, filename):
        """计算文件行数"""
        with open(filename, 'r') as f:
            return len(f.readlines())
    
    def find_band_number(self, filename):
        """读取能带数量"""
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        band_numbers = []
        reading_bands = False
        
        for line in lines:
            if 'NO.' in line:
                reading_bands = True
                continue
                
            if reading_bands:
                data = line.split()
                if len(data) >= 3 and data[0].isdigit():
                    band_numbers.append(int(data[0]))
                else:
                    break
        
        return max(band_numbers) if band_numbers else 0
    
    def check_spin_type(self, filename):
        """检查自旋类型"""
        with open(filename, 'r') as f:
            first_line = f.readline()
            return 2 if "SPIN" in first_line else 1
    
    def parse_occupation_file(self, filename, nkpt, nb, spin=1):
        """解析占据数文件"""
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        # 确定起始行
        if spin == 1:
            start_line = 0
            for i, line in enumerate(lines):
                if 'NO. ' in line:
                    start_line = i + 1
                    break

        if spin == 2:
            start_line = 0
            for i, line in enumerate(lines):
                if 'SPIN 2' in line:
                    start_line = i + 3
                    break
        occupations = []
        line_idx = start_line
        
        for kpt in range(nkpt):
            kpt_occ = []
            if 'KPOINTS' in lines[line_idx]:
                line_idx += 1
            if 'NO. 'in lines[line_idx]:
                line_idx += 1
            for band in range(nb):
                if line_idx >= len(lines):
                    break
                line = lines[line_idx].strip()
                if line:
                    parts = line.split()
                    if len(parts) >= 3:
                        try:
                            occ_value = float(parts[2])
                            kpt_occ.append(occ_value)
                        except ValueError:
                            print(f"Error parsing line {lines[line_idx]}")
                            continue
                line_idx += 1
            occupations.append(kpt_occ)      
        return occupations
    
    def analyze_defect_levels(self, occ1, occ2, nb):
        """分析缺陷能级
        Args:
            occ1: 第一个状态的占据数列表，形状为 [nkpt, nb],具有更多电子占据
            occ2: 第二个状态的占据数列表，形状为 [nkpt, nb]
            nb: 能级数量
        Returns:
            start_band: 缺陷能级起始位置(1-based)
            end_band: 缺陷能级结束位置(1-based) 
            n_defect: 缺陷能级数量
            occ_diff: 每个k点缺陷能级的占据数差异
        """
        # 计算每个k点的占据数差异
        occ_diff = []
        for kpt in range(len(occ1)):
            kpt_diff = [occ1[kpt][i] - occ2[kpt][i] for i in range(nb)]
            occ_diff.append(kpt_diff)
        
        # 1. 判断每个能级的占据情况（通过不同k点叠加）
        # 计算每个能级在所有k点上的平均占据数
        avg_occupation1 = [0] * nb  # 状态1的平均占据
        avg_occupation2 = [0] * nb  # 状态2的平均占据
        
        for band in range(nb):
            total_occ1 = 0
            total_occ2 = 0
            for kpt in range(len(occ1)):
                total_occ1 += occ1[kpt][band]
                total_occ2 += occ2[kpt][band]
            avg_occupation1[band] = total_occ1 
            avg_occupation2[band] = total_occ2 
        
        # 2. 根据能级是否为分数占据，确定能级是否为缺陷能级
        defect_bands = []  # 存储缺陷能级的索引(0-based)
        threshold = 0.01  # 分数占据阈值
        
        for band in range(nb):
            # 判断是否为分数占据：既不完全占据(1)也不完全空置(0)
            is_fractional1 = (avg_occupation1[band] > threshold and 
                            avg_occupation1[band] < 1 - threshold)
            is_fractional2 = (avg_occupation2[band] > threshold and 
                            avg_occupation2[band] < 1 - threshold)
            
            # 如果任一状态在该能级是分数占据，则认为是缺陷能级
            if is_fractional1 or is_fractional2:
                defect_bands.append(band)
            if  abs(avg_occupation1[band] - avg_occupation2[band]) > threshold:
                defect_bands.append(band)
        #去重
        defect_bands = list(set(defect_bands))
        
        # 3. 对缺陷能级进行归一化处理
        if defect_bands:
            # 找到连续的缺陷能级范围
            defect_bands.sort()
            start_band = defect_bands[0] + 1  # 转为1-based索引
            end_band = defect_bands[-1] + 1
            n_defect = len(defect_bands)
            
            print(f"Defect bands identified: {defect_bands} (0-based)")
            print(f"Defect level range: {start_band}-{end_band}, {n_defect} levels")
            
            # 对每个k点的缺陷能级占据数差异进行归一化
            normalized_occ_diff = []
            for kpt in range(len(occ_diff)):
                kpt_normalized = []
                for band in defect_bands:
                    # 归一化到[-1, 1]范围
                    diff_value = occ_diff[kpt][band]
                    # 可以根据需要调整归一化方法
                    normalized_value = diff_value  # 这里保持原值，可根据需要修改
                    kpt_normalized.append(normalized_value)
                normalized_occ_diff.append(kpt_normalized)
            #打印每个缺陷能级的贡献
            for band in defect_bands:
                occ=0
                for kpt in range(len(occ_diff)):
                    occ+=occ_diff[kpt][band]
                print(f"Defect band occ {band+1}: {occ:.4f}")
            return start_band, end_band, n_defect, normalized_occ_diff
        else:
            print("No defect levels found")
            return None, None, 0, occ_diff
    
    def calculate_total_electrons(self, occupations):
        """计算总电子数"""
        total = 0
        for kpt_occ in occupations:
            kpt_occ = np.array(kpt_occ)
            total += sum(kpt_occ)
        
        return round(total, 4)
    
    def generate_occ_input(self, nkpt, defect_info, occ_diff, q_diff, spin_suffix=""):
        """生成occ.input文件"""
        if spin_suffix:
            filename = f'occ{spin_suffix}.input'
        else:
            filename = 'occ.input'
            
        start_band, end_band, n_defect, _ = defect_info
        
        with open(filename, 'w') as f:
            f.write('nkpt  ndefect_levels  defect_level_1  q\n')
            f.write(f'{nkpt}          {n_defect}              {start_band}        {q_diff}\n')
            
            for kpt in range(nkpt):
                defect_occ = occ_diff[kpt]
                occ_str = ' '.join(f'{val:.3f}' for val in defect_occ)
                f.write(occ_str + '\n')
    
    def generate_in_occ(self, occupations, defect_info, nb, spin_suffix=""):
        """生成IN.OCC文件"""
        if spin_suffix:
            filename = f'IN.OCC{spin_suffix}'
        else:
            filename = 'IN.OCC'
            
        start_band, end_band, n_defect, _ = defect_info
        
        with open(filename, 'w') as f:
            for kpt in range(len(occupations)):
                # 计算相对占据数
                ref_occ = occupations[kpt][0]  # 参考占据数（第一个能级）
                relative_occ = []
                for band in range(start_band-1, end_band):
                    if ref_occ > 0:
                        rel_occ = occupations[kpt][band] / ref_occ
                    else:
                        rel_occ = 0.0
                    relative_occ.append(f'{rel_occ:.6f}')
                
                relative_str = ' '.join(relative_occ)
                f.write(f'{start_band-1}*1.000000 {relative_str} {nb-end_band}*0.000000\n')
    
    def analyze(self):
        """主分析函数"""
        print('generate occ.input')
        print('we need OUT.OCC0 (for 0 state) and OUT.OCC1 (for q state)')
        
        # 读取k点数量
        self.read_nkpt()
        
        # 检查文件是否存在
        files = ['OUT.OCC0', 'OUT.OCC1']
        for f in files:
            if not os.path.exists(f):
                raise FileNotFoundError(f"{f}文件不存在")
        
        # 确定自旋类型（以第一个文件为准）
        self.spin_type = self.check_spin_type('OUT.OCC0')
        print(f"Spin type: {'spin=2' if self.spin_type == 2 else 'spin=1'}")
        
        # 读取能带数量
        nb0 = self.find_band_number('OUT.OCC0')
        nb1 = self.find_band_number('OUT.OCC1')

        print(f"Number of bands in OUT.OCC0 and OUT.OCC1: {nb0}, {nb1}")
        
        if self.spin_type == 1:
            self.analyze_non_spin(nb0,nb1)
        else:
            self.analyze_spin_polarized(nb0,nb1)

    def check_electron_conservation(self, occ_diff, q_diff, spin_suffix=""):
        """检查电子数守恒"""
        total_occ_diff = 0
        for kpt_occ in occ_diff:
            total_occ_diff += sum(kpt_occ)
        
        total_occ_diff = round(total_occ_diff, 4)
        
        print(f"Electron conservation check{spin_suffix}:")
        print(f"  Sum of occupation differences: {total_occ_diff}")
        print(f"  Expected q_diff: {q_diff}")
        
        if abs(total_occ_diff - q_diff) > 0.01:  # 允许0.01的误差
            print(f"  WARNING: Electron conservation violated! Difference: {abs(total_occ_diff - q_diff):.4f}")
            return False
        else:
            return True
    def analyze_non_spin(self, nb0,nb1):
        """分析无自旋情况"""
        # 读取占据数
        occ0 = self.parse_occupation_file('OUT.OCC0', self.nkpt, nb0)
        occ1 = self.parse_occupation_file('OUT.OCC1', self.nkpt, nb1)
        
        # 计算总电子数
        q0 = self.calculate_total_electrons(occ0)
        q1 = self.calculate_total_electrons(occ1)
        print(f"Total electrons - State 0: {q0}, State 1: {q1}")
        
        # 确定基态和激发态
        if q0 > q1:
            ground, excited = occ0, occ1
            q_diff = q0 - q1
        else:
            ground, excited = occ1, occ0
            q_diff = q1 - q0
        nb=min(nb0,nb1)
        # 分析缺陷能级 ground有更多电子占据
        defect_info = self.analyze_defect_levels(ground, excited, nb)
        start_band, end_band, n_defect, occ_diff = defect_info
        print(start_band,end_band,n_defect)
        if n_defect > 0:
            print(f"Defect levels: {start_band}-{end_band}, {n_defect} levels")
            
            # 生成输出文件
            self.generate_occ_input(self.nkpt, defect_info, occ_diff, q_diff)
            self.generate_in_occ(excited, defect_info, nb)
        else:
            print("No defect levels found")
    
    def analyze_spin_polarized(self, nb0, nb1):
        """分析自旋极化情况"""
        print("Reading OUT.OCC0:")
        occ0_up = self.parse_occupation_file('OUT.OCC0', self.nkpt, nb0,spin=1)
        occ0_down = self.parse_occupation_file('OUT.OCC0', self.nkpt, nb0,spin=2)  
        print("Reading OUT.OCC1:")
        occ1_up = self.parse_occupation_file('OUT.OCC1', self.nkpt, nb1, spin=1)
        occ1_down = self.parse_occupation_file('OUT.OCC1', self.nkpt, nb1,spin=2)  # 简化
        q0_up = self.calculate_total_electrons(occ0_up)
        q0_down = self.calculate_total_electrons(occ0_down)
        q0_total = q0_up + q0_down
        
        q1_up = self.calculate_total_electrons(occ1_up)
        q1_down = self.calculate_total_electrons(occ1_down)
        q1_total = q1_up + q1_down

        print(f"Total electrons - State 0: {q0_total}, State 1: {q1_total}")
        
        # 确定基态和激发态
        if q0_total > q1_total:
            ground_up, excited_up = occ0_up, occ1_up
            ground_down, excited_down = occ0_down, occ1_down
            q_diff = q0_total - q1_total
        else:
            ground_up, excited_up = occ1_up, occ0_up
            ground_down, excited_down = occ1_down, occ0_down
            q_diff = q1_total - q0_total
        
        # 分析自旋向上的缺陷能级
        nb=min(nb0,nb1)
        print(f"analyzing Defect levels in spin-up:")
        defect_up = self.analyze_defect_levels(ground_up, excited_up, nb)
        
        start_up, end_up, n_defect_up, occ_diff_up = defect_up

        if n_defect_up > 0:
            print(f"Defect levels in spin-up: {start_up}-{end_up}, {n_defect_up} levels")
            self.generate_occ_input(self.nkpt, defect_up, occ_diff_up, q_diff, "_spin1")
            self.generate_in_occ(excited_up, defect_up, nb, "_1")
        
        # 分析自旋向下的缺陷能级 
        print(f"analyzing Defect levels in spin-down: ")
        defect_down = self.analyze_defect_levels(ground_down, excited_down, nb)
        start_down, end_down, n_defect_down, occ_diff_down = defect_down
        
        if n_defect_down > 0:
            print(f"Defect levels in spin-down: {start_down}-{end_down}, {n_defect_down} levels")
            self.generate_occ_input(self.nkpt, defect_down, occ_diff_down, q_diff, "_spin2")
            self.generate_in_occ(excited_down, defect_down, nb, "_2")

        #检查电子数守恒
        total_occ_diff_up = sum([sum(occ) for occ in occ_diff_up])
        total_occ_diff_down = sum([sum(occ) for occ in occ_diff_down])
        total_occ_diff = round(total_occ_diff_up + total_occ_diff_down, 4)
        
        print(f"Total electron conservation check:")
        print(f"  Sum of spin-up occupation differences: {total_occ_diff_up:.4f}")
        print(f"  Sum of spin-down occupation differences: {total_occ_diff_down:.4f}")
        print(f"  Total sum: {total_occ_diff}")
        print(f"  Expected q_diff: {q_diff}")
        if n_defect_up > 0 and n_defect_down > 0:
            if abs(total_occ_diff - q_diff) > 0.1:
                print(f"WARNING: Total electron conservation violated! Difference: {abs(total_occ_diff - q_diff):.4f}")
        if n_defect_up>0 and n_defect_down==0:
            if abs(total_occ_diff_up-q_diff)>0.1:
                print(f"WARNING: Total electron conservation violated! Difference: {abs(total_occ_diff_up - q_diff):.4f}")
        if n_defect_up==0 and n_defect_down>0:
            if abs(total_occ_diff_down-q_diff)>0.1:
                print(f"WARNING: Total electron conservation violated! Difference: {abs(total_occ_diff_down - q_diff):.4f}")
                
def main():
    
    analyzer = OccupationAnalyzer()
    analyzer.analyze()
    print("Analysis completed successfully!")
    

if __name__ == "__main__":
    main()
